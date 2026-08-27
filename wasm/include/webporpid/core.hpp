#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace webporpid {

struct Parameters {
  double error_rate = 0.05;
  std::uint32_t min_length = 2100;
  std::uint32_t max_length = 4300;
  std::uint32_t primer_tolerance = 1;
  std::uint32_t primer_window = 200;
  std::uint32_t primer_chop = 0;
  std::uint32_t max_reads_per_sample = 100000;
  std::uint32_t family_size_threshold = 1;
  double lda_threshold = 0.995;
  double contamination_cluster_threshold = 0.015;
  double contamination_proportion_threshold = 0.2;
  double contamination_distance_threshold = 0.015;
  bool contamination_filter = true;
  double agreement_threshold = 0.6;
  double artefact_fraction = 0.25;
  double outlier_quantile = 0.99;
  double panel_threshold = 50.0;
  double functional_match_threshold = 0.7;
  std::uint32_t spool_partitions = 64;
  std::uint64_t deterministic_seed = 0x504f52504944ULL;
};

struct Sample {
  std::string name;
  std::string cdna_primer;
  std::string second_strand_primer;
  std::string panel_name;
  std::string functional_reference_name;
  std::vector<std::pair<std::string, std::string>> panel_sequences;
  std::string functional_reference;
  std::int32_t family_size_override = -1;
  double artefact_fraction_override = std::numeric_limits<double>::quiet_NaN();
  double outlier_quantile_override = std::numeric_limits<double>::quiet_NaN();
  double agreement_override = std::numeric_limits<double>::quiet_NaN();
  double functional_match_override = std::numeric_limits<double>::quiet_NaN();
};

struct Config {
  std::string dataset;
  Parameters parameters;
  std::vector<Sample> samples;
};

struct Stats {
  std::uint64_t total_reads = 0;
  std::uint64_t quality_reads = 0;
  std::uint64_t bad_reads = 0;
  std::uint64_t short_reads = 0;
  std::uint64_t long_reads = 0;
  std::uint64_t primer_rejects = 0;
  std::uint64_t id_rejects = 0;
  std::uint64_t demultiplexed_reads = 0;
  std::uint64_t bpb_rejects = 0;
  std::uint64_t malformed_records = 0;
  std::vector<std::uint64_t> per_sample;
};

struct FastqRecord { std::string name, sequence, quality; };
struct SpoolRecord {
  std::uint16_t sample = 0;
  std::string umi, name, sequence, quality;
  std::uint64_t sampling_hash = 0;
};

enum class FamilyDisposition : std::uint8_t {
  likely_real = 0, bpb_reject = 1, heteroduplex = 2, lda_reject = 3,
  umi_length_reject = 4, family_size_reject = 5,
};

struct FamilyDecision {
  std::uint16_t sample = 0;
  std::string umi, parent;
  std::uint32_t family_size = 0;
  double probability = 0.0;
  FamilyDisposition disposition = FamilyDisposition::lda_reject;
};

struct LowAgreementSite {
  std::uint32_t position = 0;
  float agreement = 0.0f;
  char modal_base = '-';
  std::uint32_t run_length = 0;
};

struct ConsensusRecord {
  std::uint16_t sample = 0;
  std::string id, umi, sequence;
  std::uint32_t family_size = 0;
  double minimum_agreement = 0.0;
  std::vector<LowAgreementSite> low_sites;
};

struct Alignment { std::string reference, query; };

bool decode_config(std::span<const std::uint8_t> bytes, Config& output, std::string& error);
std::vector<FastqRecord> parse_fastq(std::string_view text, Stats& stats);
std::vector<std::uint8_t> preprocess_batch(const Config& config, Stats& stats,
                                           std::string_view fastq, std::uint64_t first_ordinal);
std::vector<SpoolRecord> decode_spool(std::span<const std::uint8_t> bytes, std::string& error);
void encode_spool_record(const SpoolRecord& record, std::vector<std::uint8_t>& output);

std::vector<std::uint8_t> partition_sample_hashes(std::span<const std::uint8_t> bytes,
                                                  const Config& config, std::string& error);
std::vector<std::uint8_t> count_families(std::span<const std::uint8_t> bytes,
                                         std::span<const std::uint8_t> cutoffs,
                                         const Config& config, std::string& error);
std::vector<FamilyDecision> build_family_model(std::span<const std::uint8_t> merged_counts,
                                               const Config& config, std::string& error);
std::vector<std::uint8_t> encode_family_model(const std::vector<FamilyDecision>& model);
std::vector<FamilyDecision> decode_family_model(std::span<const std::uint8_t> bytes, std::string& error);

Alignment needleman_wunsch(std::string_view reference, std::string_view query);
Alignment seeded_global_align(std::string_view reference, std::string_view query);
std::string family_consensus(const std::vector<std::string>& reads,
                             double& minimum_agreement,
                             std::vector<LowAgreementSite>& low_sites);
std::vector<std::uint8_t> process_consensus_partition(std::span<const std::uint8_t> bytes,
                                                      std::span<const std::uint8_t> cutoffs,
                                                      const std::vector<FamilyDecision>& model,
                                                      const Config& config, std::string& error);

char complement(char value);
std::string reverse_complement(std::string_view sequence);
std::uint64_t stable_hash(std::uint64_t seed, std::string_view first, std::string_view second = {});
std::string json_escape(std::string_view value);

} // namespace webporpid
