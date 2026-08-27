#include "webporpid/core.hpp"
#include "webporpid/binary.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

void check(bool condition, const char* message) {
  if (!condition) { std::cerr << "FAIL: " << message << '\n'; std::exit(1); }
}

std::string degap(std::string value) {
  value.erase(std::remove(value.begin(), value.end(), '-'), value.end()); return value;
}

webporpid::Config synthetic_config() {
  webporpid::Config config; config.dataset = "synthetic";
  config.parameters.min_length = 20; config.parameters.max_length = 300;
  config.parameters.primer_window = 150; config.parameters.primer_tolerance = 0;
  config.parameters.spool_partitions = 8; config.parameters.family_size_threshold = 1;
  webporpid::Sample sample; sample.name = "sample_1";
  sample.cdna_primer = "CCGCTacgtaaNNNNNNNNGTCA";
  sample.second_strand_primer = "TAGG"; sample.panel_name = "panel.fa";
  config.samples.push_back(std::move(sample)); return config;
}

void test_utilities() {
  check(webporpid::reverse_complement("ACGTRYMK") == "MKRYACGT", "IUPAC reverse complement");
  check(webporpid::stable_hash(7, "read", "UMI") == webporpid::stable_hash(7, "read", "UMI"), "stable hash");
  auto alignment = webporpid::seeded_global_align(
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
    "ACGTACGTACGTACGTTACGTACGTACGTACGTACGTACGT");
  check(alignment.reference.size() == alignment.query.size(), "rectangular seeded alignment");
  check(degap(alignment.reference) == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "seeded reference reconstruction");
  check(degap(alignment.query) == "ACGTACGTACGTACGTTACGTACGTACGTACGTACGTACGT", "seeded query reconstruction");
}

void test_spool_roundtrip() {
  webporpid::SpoolRecord input{2, "ACGTACGT", "seq12", "ACGTN", "IIIII", 1234567};
  std::vector<std::uint8_t> encoded; webporpid::encode_spool_record(input, encoded); std::string error;
  const auto decoded = webporpid::decode_spool(encoded, error);
  check(error.empty() && decoded.size() == 1, "spool record decodes");
  check(decoded[0].sample == input.sample && decoded[0].umi == input.umi && decoded[0].sequence == input.sequence
    && decoded[0].quality == input.quality && decoded[0].sampling_hash == input.sampling_hash, "spool round trip");
}

void test_preprocess() {
  auto config = synthetic_config(); webporpid::Stats stats;
  const std::string tag = "ACGTAA" "AACCGGTT" "GTCA", payload(65, 'A');
  const std::string sequence = "TAGG" + webporpid::reverse_complement(tag + payload) + webporpid::reverse_complement("CCGCT");
  const std::string fastq = "@synthetic\n" + sequence + "\n+\n" + std::string(sequence.size(), 'I') + "\n";
  const auto routed = webporpid::preprocess_batch(config, stats, fastq, 10);
  check(stats.total_reads == 1 && stats.quality_reads == 1 && stats.demultiplexed_reads == 1, "synthetic read demultiplexes");
  check(stats.bpb_rejects == 0 && routed.size() > 5, "synthetic tag passes BPB");
  std::string error; const auto records = webporpid::decode_spool(std::span(routed).subspan(5), error);
  check(error.empty() && records.size() == 1, "routed spool payload decodes");
  check(records[0].umi == "AACCGGTT", "BPB DP extracts the eight-base UMI");
  check(records[0].sequence.starts_with(tag), "oriented insert starts at sample ID");
}

void test_umi_and_consensus() {
  auto config = synthetic_config();
  std::vector<std::uint8_t> counts; webporpid::binary::magic(counts, "WPN1"); webporpid::binary::number(counts, std::uint32_t{2});
  webporpid::binary::number(counts, std::uint16_t{0}); webporpid::binary::string(counts, "AACCGGTT"); webporpid::binary::number(counts, std::uint32_t{100000});
  webporpid::binary::number(counts, std::uint16_t{0}); webporpid::binary::string(counts, "AACCGGTA"); webporpid::binary::number(counts, std::uint32_t{1});
  std::string error; const auto model = webporpid::build_family_model(counts, config, error);
  check(error.empty() && model.size() == 2, "sparse UMI likelihood model builds");
  const auto child = std::find_if(model.begin(), model.end(), [](const auto& value) { return value.umi == "AACCGGTA"; });
  check(child != model.end() && child->parent == "AACCGGTT", "low-count one-edit tag is assigned to its parent");

  const std::string expected = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  std::vector<std::string> reads{expected, expected, expected, expected.substr(0, 22) + "T" + expected.substr(22)};
  double minimum = 0.0; std::vector<webporpid::LowAgreementSite> low;
  const auto consensus = webporpid::family_consensus(reads, minimum, low);
  check(consensus == expected, "three-pass indel-aware family consensus preserves the majority sequence");
  check(minimum >= 0.75 && minimum <= 1.0, "minimum agreement is counted on final alignments");

  const std::string split_expected = "ACGTACGTACGTACGTACGT";
  const std::string split_mutant = "ACGTACGTATGTACGTACGT";
  const std::vector<std::string> split_family{
    split_expected, split_expected, split_expected, split_mutant, split_mutant,
  };
  const auto split_consensus = webporpid::family_consensus(split_family, minimum, low);
  check(split_consensus == split_expected, "small-edit split family preserves the modal base");
  check(std::abs(minimum - 0.6) < 1e-9, "small-edit split family reports exact minimum agreement");

  const std::vector<std::string> noisy_family{
    "ACGTAAAACCGGTTGTCAATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA",
    "ACGTAAAACCGGTTGTCAATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA",
    "ACGTAAAACCGGTTGTCAATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA",
    "ACGTAAAACCGGTTGTCAATGCCTTGGGCCATCGGACCATATGTTTAACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA",
    "ACGTAAAACCGGTTGTCAATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACGTCAATTCGTCTCAGAGAAGTAA",
  };
  const auto noisy = webporpid::family_consensus(noisy_family, minimum, low);
  check(noisy.starts_with("ACGTAAAACCGGTTGTCA"), "five-read noisy family consensus completes");
}

} // namespace

int main() {
  test_utilities(); test_spool_roundtrip(); test_preprocess(); test_umi_and_consensus();
  std::cout << "webPORPID native core tests passed\n"; return 0;
}
