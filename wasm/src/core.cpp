#include "webporpid/core.hpp"
#include "webporpid/binary.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace webporpid {

namespace {

bool read_parameters(binary::Reader& reader, Parameters& p) {
  std::uint8_t contamination = 0;
  return reader.number(p.error_rate)
      && reader.number(p.min_length) && reader.number(p.max_length)
      && reader.number(p.primer_tolerance) && reader.number(p.primer_window)
      && reader.number(p.primer_chop) && reader.number(p.max_reads_per_sample)
      && reader.number(p.family_size_threshold) && reader.number(p.lda_threshold)
      && reader.number(p.contamination_cluster_threshold)
      && reader.number(p.contamination_proportion_threshold)
      && reader.number(p.contamination_distance_threshold)
      && reader.number(contamination)
      && (p.contamination_filter = contamination != 0, true)
      && reader.number(p.agreement_threshold) && reader.number(p.artefact_fraction)
      && reader.number(p.outlier_quantile) && reader.number(p.panel_threshold)
      && reader.number(p.functional_match_threshold) && reader.number(p.spool_partitions)
      && reader.number(p.deterministic_seed);
}

bool valid_probability(double value) { return std::isfinite(value) && value >= 0.0 && value <= 1.0; }

std::string lowercase_run(std::string_view primer) {
  std::string best, current;
  for (char value : primer) {
    if (std::islower(static_cast<unsigned char>(value))) current.push_back(static_cast<char>(std::toupper(value)));
    else { if (current.size() > best.size()) best = current; current.clear(); }
  }
  if (current.size() > best.size()) best = current;
  return best;
}

} // namespace

bool decode_config(std::span<const std::uint8_t> bytes, Config& output, std::string& error) {
  binary::Reader reader(bytes);
  std::uint32_t version = 0, sample_count = 0;
  if (!reader.magic("WPC1") || !reader.number(version) || version != 1 || !reader.string(output.dataset)
      || !read_parameters(reader, output.parameters) || !reader.number(sample_count)) {
    error = "The compiled configuration is truncated or has an unsupported version.";
    return false;
  }
  if (output.dataset.empty() || sample_count == 0 || sample_count > 65535) {
    error = "The configuration must contain a dataset and at least one sample.";
    return false;
  }
  output.samples.clear(); output.samples.reserve(sample_count);
  for (std::uint32_t index = 0; index < sample_count; ++index) {
    Sample sample;
    std::uint32_t panel_count = 0;
    if (!reader.string(sample.name) || !reader.string(sample.cdna_primer)
        || !reader.string(sample.second_strand_primer) || !reader.string(sample.panel_name)
        || !reader.string(sample.functional_reference_name)
        || !reader.number(sample.family_size_override)
        || !reader.number(sample.artefact_fraction_override)
        || !reader.number(sample.outlier_quantile_override)
        || !reader.number(sample.agreement_override)
        || !reader.number(sample.functional_match_override)
        || !reader.number(panel_count) || panel_count > 100000) {
      error = "A sample entry in the compiled configuration is invalid.";
      return false;
    }
    sample.panel_sequences.reserve(panel_count);
    for (std::uint32_t panel = 0; panel < panel_count; ++panel) {
      std::string name, sequence;
      if (!reader.string(name) || !reader.string(sequence)) {
        error = "A panel FASTA entry is truncated."; return false;
      }
      sample.panel_sequences.emplace_back(std::move(name), std::move(sequence));
    }
    if (!reader.string(sample.functional_reference)) {
      error = "A functional-reference entry is truncated."; return false;
    }
    if (sample.name.empty() || lowercase_run(sample.cdna_primer).empty()
        || sample.second_strand_primer.empty()) {
      error = "Every sample needs a name, a lower-case sample ID, and a second-strand primer.";
      return false;
    }
    output.samples.push_back(std::move(sample));
  }
  const auto& p = output.parameters;
  if (!reader.done() || !valid_probability(p.error_rate) || p.min_length >= p.max_length
      || p.primer_window < 16 || p.spool_partitions == 0 || p.spool_partitions > 256
      || (p.spool_partitions & (p.spool_partitions - 1)) != 0
      || !valid_probability(p.lda_threshold) || !valid_probability(p.agreement_threshold)
      || !valid_probability(p.artefact_fraction) || !valid_probability(p.outlier_quantile)
      || !valid_probability(p.functional_match_threshold)) {
    error = "The compiled configuration contains invalid parameters or trailing bytes.";
    return false;
  }
  return true;
}

std::vector<FastqRecord> parse_fastq(std::string_view text, Stats& stats) {
  std::vector<std::string_view> lines;
  std::size_t start = 0;
  for (std::size_t index = 0; index <= text.size(); ++index) {
    if (index != text.size() && text[index] != '\n') continue;
    auto line = text.substr(start, index - start);
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    if (!line.empty() || index != text.size()) lines.push_back(line);
    start = index + 1;
  }
  std::vector<FastqRecord> output;
  output.reserve(lines.size() / 4);
  for (std::size_t index = 0; index + 3 < lines.size(); index += 4) {
    if (lines[index].empty() || lines[index][0] != '@' || lines[index + 2].empty()
        || lines[index + 2][0] != '+' || lines[index + 1].size() != lines[index + 3].size()) {
      ++stats.malformed_records; continue;
    }
    FastqRecord record{std::string(lines[index].substr(1)), std::string(lines[index + 1]), std::string(lines[index + 3])};
    std::transform(record.sequence.begin(), record.sequence.end(), record.sequence.begin(),
                   [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
    output.push_back(std::move(record));
  }
  if (lines.size() % 4 != 0) ++stats.malformed_records;
  return output;
}

char complement(char value) {
  switch (std::toupper(static_cast<unsigned char>(value))) {
    case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': case 'U': return 'A';
    case 'R': return 'Y'; case 'Y': return 'R'; case 'S': return 'S'; case 'W': return 'W';
    case 'K': return 'M'; case 'M': return 'K'; case 'B': return 'V'; case 'D': return 'H';
    case 'H': return 'D'; case 'V': return 'B'; case '-': return '-'; default: return 'N';
  }
}

std::string reverse_complement(std::string_view sequence) {
  std::string output; output.reserve(sequence.size());
  for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) output.push_back(complement(*iterator));
  return output;
}

std::uint64_t stable_hash(std::uint64_t seed, std::string_view first, std::string_view second) {
  std::uint64_t value = 1469598103934665603ULL ^ seed;
  const auto consume = [&value](std::string_view text) {
    for (unsigned char byte : text) { value ^= byte; value *= 1099511628211ULL; }
    value ^= 0xff; value *= 1099511628211ULL;
  };
  consume(first); consume(second);
  value ^= value >> 33; value *= 0xff51afd7ed558ccdULL;
  value ^= value >> 33; value *= 0xc4ceb9fe1a85ec53ULL;
  return value ^ (value >> 33);
}

std::string json_escape(std::string_view value) {
  std::ostringstream output;
  for (unsigned char character : value) {
    switch (character) {
      case '"': output << "\\\""; break; case '\\': output << "\\\\"; break;
      case '\b': output << "\\b"; break; case '\f': output << "\\f"; break;
      case '\n': output << "\\n"; break; case '\r': output << "\\r"; break;
      case '\t': output << "\\t"; break;
      default:
        if (character < 0x20) output << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(character) << std::dec;
        else output << static_cast<char>(character);
    }
  }
  return output.str();
}

void encode_spool_record(const SpoolRecord& record, std::vector<std::uint8_t>& output) {
  const auto body = static_cast<std::uint32_t>(2 + 2 + 4 + 4 + 8 + record.umi.size() + record.name.size()
      + record.sequence.size() + record.quality.size());
  binary::number(output, body); binary::number(output, record.sample);
  binary::number(output, static_cast<std::uint16_t>(record.umi.size()));
  binary::number(output, static_cast<std::uint32_t>(record.name.size()));
  binary::number(output, static_cast<std::uint32_t>(record.sequence.size()));
  binary::number(output, record.sampling_hash);
  output.insert(output.end(), record.umi.begin(), record.umi.end());
  output.insert(output.end(), record.name.begin(), record.name.end());
  output.insert(output.end(), record.sequence.begin(), record.sequence.end());
  output.insert(output.end(), record.quality.begin(), record.quality.end());
}

std::vector<SpoolRecord> decode_spool(std::span<const std::uint8_t> bytes, std::string& error) {
  binary::Reader reader(bytes);
  std::vector<SpoolRecord> records;
  while (!reader.done()) {
    std::uint32_t body = 0, name_length = 0, sequence_length = 0;
    std::uint16_t umi_length = 0;
    SpoolRecord record;
    if (!reader.number(body) || body < 20 || body > reader.remaining()
        || !reader.number(record.sample) || !reader.number(umi_length)
        || !reader.number(name_length) || !reader.number(sequence_length)
        || !reader.number(record.sampling_hash)) {
      error = "A spool record has an invalid header."; return {};
    }
    const std::size_t fields = static_cast<std::size_t>(umi_length) + name_length + sequence_length * 2;
    if (body != 20 + fields) { error = "A spool record has inconsistent field lengths."; return {}; }
    auto umi = reader.bytes(umi_length); auto name = reader.bytes(name_length);
    auto sequence = reader.bytes(sequence_length); auto quality = reader.bytes(sequence_length);
    if (umi.size() != umi_length || name.size() != name_length || sequence.size() != sequence_length || quality.size() != sequence_length) {
      error = "A spool record is truncated."; return {};
    }
    record.umi.assign(reinterpret_cast<const char*>(umi.data()), umi.size());
    record.name.assign(reinterpret_cast<const char*>(name.data()), name.size());
    record.sequence.assign(reinterpret_cast<const char*>(sequence.data()), sequence.size());
    record.quality.assign(reinterpret_cast<const char*>(quality.data()), quality.size());
    records.push_back(std::move(record));
  }
  return records;
}

} // namespace webporpid
