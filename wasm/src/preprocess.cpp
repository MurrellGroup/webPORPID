#include "webporpid/core.hpp"
#include "webporpid/binary.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <optional>
#include <sstream>

namespace webporpid {
namespace {

std::uint8_t iupac_mask(char value) {
  switch (std::toupper(static_cast<unsigned char>(value))) {
    case 'A': return 1; case 'C': return 2; case 'G': return 4; case 'T': case 'U': return 8;
    case 'R': return 5; case 'Y': return 10; case 'S': return 6; case 'W': return 9;
    case 'K': return 12; case 'M': return 3; case 'B': return 14; case 'D': return 13;
    case 'H': return 11; case 'V': return 7; case 'N': return 15; default: return 0;
  }
}

bool compatible(char left, char right) { return (iupac_mask(left) & iupac_mask(right)) != 0; }

std::string upper(std::string_view value) {
  std::string output(value);
  std::transform(output.begin(), output.end(), output.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return output;
}

std::string sample_id(std::string_view primer) {
  std::string output;
  bool started = false;
  for (char value : primer) {
    if (std::islower(static_cast<unsigned char>(value))) {
      started = true; output.push_back(static_cast<char>(std::toupper(value)));
    } else if (started) break;
  }
  return output;
}

std::string conserved_adapter(const Config& config) {
  std::size_t length = config.samples.front().cdna_primer.size();
  for (const auto& sample : config.samples) length = std::min(length, sample.cdna_primer.size());
  std::size_t end = 0;
  for (std::size_t position = 0; position < length; ++position) {
    const char value = config.samples.front().cdna_primer[position];
    if (std::islower(static_cast<unsigned char>(value))) break;
    bool same = true;
    for (const auto& sample : config.samples) same = same && sample.cdna_primer[position] == value;
    if (!same) break;
    end = position + 1;
  }
  return upper(std::string_view(config.samples.front().cdna_primer).substr(0, end));
}

struct Hit { std::size_t start = 0, end = 0, errors = 0; };

// Semi-global Levenshtein search with a free text prefix. BioSequences'
// findbestlast repeatedly takes the first hit, then resumes immediately after
// that non-overlapping hit. This distinction matters for repeated primers.
std::optional<Hit> approximate_find_first(std::string_view text, std::string_view pattern, std::size_t tolerance) {
  if (pattern.empty()) return Hit{};
  std::vector<std::size_t> previous_cost(pattern.size() + 1), previous_start(pattern.size() + 1, 0);
  for (std::size_t row = 0; row <= pattern.size(); ++row) previous_cost[row] = row;
  for (std::size_t column = 0; column < text.size(); ++column) {
    std::vector<std::size_t> cost(pattern.size() + 1, 0), starts(pattern.size() + 1, column + 1);
    for (std::size_t row = 1; row <= pattern.size(); ++row) {
      const std::array<std::pair<std::size_t, std::size_t>, 3> choices{{
        {previous_cost[row - 1] + !compatible(pattern[row - 1], text[column]), previous_start[row - 1]},
        {previous_cost[row] + 1, previous_start[row]}, {cost[row - 1] + 1, starts[row - 1]},
      }};
      auto choice = std::min_element(choices.begin(), choices.end(), [](const auto& left, const auto& right) {
        return left.first != right.first ? left.first < right.first : left.second > right.second;
      });
      cost[row] = choice->first; starts[row] = choice->second;
    }
    if (cost.back() <= tolerance) {
      return Hit{starts.back(), column + 1, cost.back()};
    }
    previous_cost.swap(cost); previous_start.swap(starts);
  }
  return std::nullopt;
}

std::optional<Hit> approximate_find_last(std::string_view text, std::string_view pattern, std::size_t tolerance) {
  std::optional<Hit> best;
  std::size_t cursor = 0;
  while (cursor < text.size()) {
    const auto next = approximate_find_first(text.substr(cursor), pattern, tolerance);
    if (!next) break;
    best = Hit{cursor + next->start, cursor + next->end, next->errors};
    if (best->end <= cursor) break;
    cursor = best->end;
  }
  return best;
}

std::optional<std::size_t> exact_find(std::string_view text, std::string_view pattern) {
  if (pattern.empty() || pattern.size() > text.size()) return std::nullopt;
  for (std::size_t start = 0; start + pattern.size() <= text.size(); ++start) {
    bool matches = true;
    for (std::size_t index = 0; index < pattern.size(); ++index) matches = matches && compatible(text[start + index], pattern[index]);
    if (matches) return start;
  }
  return std::nullopt;
}

struct OrientedRead { std::string sequence, quality; };

struct PrimerGroup { std::string forward; std::vector<std::size_t> samples; };

std::vector<PrimerGroup> primer_groups(const Config& config) {
  std::vector<std::string> unique;
  for (const auto& sample : config.samples) {
    const auto primer = upper(sample.second_strand_primer);
    if (std::find(unique.begin(), unique.end(), primer) == unique.end()) unique.push_back(primer);
  }
  std::vector<std::string> maximal;
  for (const auto& primer : unique) {
    bool substring = false;
    for (const auto& other : unique) if (primer != other && other.find(primer) != std::string::npos) substring = true;
    if (!substring) maximal.push_back(primer);
  }
  std::vector<PrimerGroup> groups(maximal.size());
  for (std::size_t index = 0; index < config.samples.size(); ++index) {
    const auto primer = upper(config.samples[index].second_strand_primer);
    std::size_t selected = 0; bool found = false;
    for (std::size_t group = 0; group < maximal.size(); ++group) {
      if (maximal[group].find(primer) == std::string::npos) continue;
      if (!found || maximal[group].size() < maximal[selected].size()) { selected = group; found = true; }
    }
    if (!found) { maximal.push_back(primer); groups.emplace_back(); selected = groups.size() - 1; }
    groups[selected].samples.push_back(index);
  }
  for (auto& group : groups) if (!group.samples.empty())
    group.forward = upper(config.samples[group.samples.front()].second_strand_primer);
  groups.erase(std::remove_if(groups.begin(), groups.end(), [](const auto& group) { return group.samples.empty(); }), groups.end());
  return groups;
}

std::optional<OrientedRead> primer_pair_trim(const FastqRecord& record, std::string_view forward,
                                             std::string_view reverse, const Parameters& p) {
  if (p.primer_chop >= forward.size() || p.primer_chop >= reverse.size()) return std::nullopt;
  for (bool reverse_orientation : {false, true}) {
    const std::string sequence = reverse_orientation ? reverse_complement(record.sequence) : record.sequence;
    std::string quality = record.quality;
    if (reverse_orientation) std::reverse(quality.begin(), quality.end());
    const auto fwd = approximate_find_last(std::string_view(sequence).substr(0, std::min<std::size_t>(p.primer_window, sequence.size())),
                                           forward.substr(p.primer_chop), p.primer_tolerance);
    if (!fwd) continue;
    const std::size_t tail_start = sequence.size() > p.primer_window ? sequence.size() - p.primer_window : 0;
    const auto tail_rc = reverse_complement(std::string_view(sequence).substr(tail_start));
    const auto rev = approximate_find_last(tail_rc, reverse.substr(p.primer_chop), p.primer_tolerance);
    if (!rev) continue;
    const std::size_t left = fwd->end;
    const std::size_t right = sequence.size() - rev->end;
    if (left >= right) continue;
    return OrientedRead{sequence.substr(left, right - left), quality.substr(left, right - left)};
  }
  return std::nullopt;
}

enum class TemplateState : std::uint8_t { fixed, barcode, repeat_any };
struct TemplateSymbol { TemplateState state; char base = 'N'; };

std::vector<TemplateSymbol> tag_template(const Sample& sample) {
  std::vector<TemplateSymbol> output;
  const auto id = sample_id(sample.cdna_primer);
  for (char value : id) output.push_back({TemplateState::fixed, value});
  const auto start = sample.cdna_primer.find_first_of("Nn");
  if (start == std::string::npos) return {};
  for (char value : std::string_view(sample.cdna_primer).substr(start)) {
    if (value == 'N' || value == 'n') output.push_back({TemplateState::barcode, 'N'});
    else output.push_back({TemplateState::fixed, static_cast<char>(std::toupper(static_cast<unsigned char>(value)))});
  }
  output.push_back({TemplateState::repeat_any, '*'});
  return output;
}

struct TagResult { std::string tag; std::size_t errors = 0; double score = 0.0; };
enum class Trace : std::uint8_t { diagonal, deletion, insertion };

TagResult extract_tag(std::string_view sequence, const std::vector<TemplateSymbol>& pattern) {
  constexpr double insertion = -4.605170185988091;
  constexpr double deletion = -4.605170185988091;
  constexpr double barcode = -1.3862943611198906;
  constexpr double match = -0.0010005003335835344;
  constexpr double mismatch = -8.006367567650246;
  constexpr double epsilon = 1e-10;
  const std::size_t rows = pattern.size() + 1, columns = sequence.size() + 1;
  const auto at = [columns](std::size_t row, std::size_t column) { return row * columns + column; };
  std::vector<double> scores(rows * columns, 0.0);
  std::vector<Trace> trace(rows * columns, Trace::diagonal);
  for (std::size_t row = 1; row < rows; ++row) {
    scores[at(row, 0)] = scores[at(row - 1, 0)] + (pattern[row - 1].state == TemplateState::repeat_any ? 0.0 : deletion);
    trace[at(row, 0)] = Trace::deletion;
  }
  for (std::size_t column = 1; column < columns; ++column) {
    scores[at(0, column)] = scores[at(0, column - 1)] + insertion;
    trace[at(0, column)] = Trace::insertion;
  }
  for (std::size_t row = 1; row < rows; ++row) for (std::size_t column = 1; column < columns; ++column) {
    const auto symbol = pattern[row - 1];
    if (symbol.state == TemplateState::repeat_any) {
      const double ins = scores[at(row, column - 1)], del = scores[at(row - 1, column)];
      if (del > ins - epsilon) { scores[at(row, column)] = del; trace[at(row, column)] = Trace::deletion; }
      else { scores[at(row, column)] = ins; trace[at(row, column)] = Trace::insertion; }
      continue;
    }
    const double emission = symbol.state == TemplateState::barcode ? barcode
      : (compatible(symbol.base, sequence[column - 1]) ? match : mismatch);
    const double diagonal = scores[at(row - 1, column - 1)] + emission;
    const double ins = scores[at(row, column - 1)] + insertion;
    const double del = scores[at(row - 1, column)] + deletion;
    if (diagonal > del - epsilon && diagonal > ins - epsilon) { scores[at(row, column)] = diagonal; trace[at(row, column)] = Trace::diagonal; }
    else if (ins > del - epsilon) { scores[at(row, column)] = ins; trace[at(row, column)] = Trace::insertion; }
    else { scores[at(row, column)] = del; trace[at(row, column)] = Trace::deletion; }
  }
  TagResult result; result.score = scores.back();
  std::size_t row = rows - 1, column = columns - 1;
  while (row || column) {
    switch (trace[at(row, column)]) {
      case Trace::diagonal:
        if (pattern[row - 1].state == TemplateState::barcode) result.tag.push_back(sequence[column - 1]);
        if (pattern[row - 1].state == TemplateState::fixed && !compatible(pattern[row - 1].base, sequence[column - 1])) ++result.errors;
        --row; --column; break;
      case Trace::insertion: {
        const bool current = row > 0 && pattern[row - 1].state == TemplateState::barcode;
        const bool next = row < pattern.size() && pattern[row].state == TemplateState::barcode;
        if (current || next) result.tag.push_back(sequence[column - 1]);
        if (!(row > 0 && pattern[row - 1].state == TemplateState::repeat_any)) ++result.errors;
        --column; break;
      }
      case Trace::deletion:
        if (pattern[row - 1].state != TemplateState::repeat_any) ++result.errors;
        --row; break;
    }
  }
  std::reverse(result.tag.begin(), result.tag.end());
  result.tag = upper(result.tag);
  return result;
}

double mean_error(std::string_view quality) {
  if (quality.empty()) return 1.0;
  double total = 0.0;
  for (unsigned char value : quality) total += std::pow(10.0, -static_cast<double>(value >= 33 ? value - 33 : 0) / 10.0);
  return total / static_cast<double>(quality.size());
}

} // namespace

std::vector<std::uint8_t> preprocess_batch(const Config& config, Stats& stats,
                                           std::string_view fastq, std::uint64_t first_ordinal) {
  if (stats.per_sample.size() != config.samples.size()) stats.per_sample.assign(config.samples.size(), 0);
  const auto adapter = conserved_adapter(config);
  const auto groups = primer_groups(config);
  std::vector<std::string> ids; std::vector<std::vector<TemplateSymbol>> patterns;
  ids.reserve(config.samples.size()); patterns.reserve(config.samples.size());
  for (const auto& sample : config.samples) { ids.push_back(sample_id(sample.cdna_primer)); patterns.push_back(tag_template(sample)); }
  auto records = parse_fastq(fastq, stats);
  std::vector<std::uint8_t> output;
  for (std::size_t record_index = 0; record_index < records.size(); ++record_index) {
    auto& record = records[record_index];
    ++stats.total_reads;
    const double error = mean_error(record.quality);
    if (error >= config.parameters.error_rate) { ++stats.bad_reads; continue; }
    if (record.sequence.size() <= config.parameters.min_length) { ++stats.short_reads; continue; }
    if (record.sequence.size() >= config.parameters.max_length) { ++stats.long_reads; continue; }
    ++stats.quality_reads;
    std::ostringstream name; name << "seq" << first_ordinal + record_index << "|ee=" << std::setprecision(16) << error;
    bool had_primers = false, assigned = false;
    for (const auto& group : groups) {
      const auto trimmed = primer_pair_trim(record, group.forward, adapter, config.parameters);
      if (!trimmed) continue;
      had_primers = true;
      std::string oriented_sequence = reverse_complement(trimmed->sequence);
      std::string oriented_quality = trimmed->quality; std::reverse(oriented_quality.begin(), oriented_quality.end());
      for (const auto sample_index : group.samples) {
        const auto start = exact_find(std::string_view(oriented_sequence).substr(0, std::min<std::size_t>(10, oriented_sequence.size())), ids[sample_index]);
        if (!start) continue;
        auto sample_sequence = oriented_sequence.substr(*start);
        auto sample_quality = oriented_quality.substr(*start);
        const auto suffix = patterns[sample_index].size() - std::min(patterns[sample_index].size(), ids[sample_index].size());
        const std::size_t tag_end = std::min(sample_sequence.size(), std::size_t{6} + suffix + 2);
        const auto tag = extract_tag(std::string_view(sample_sequence).substr(0, tag_end), patterns[sample_index]);
        ++stats.demultiplexed_reads; ++stats.per_sample[sample_index]; assigned = true;
        if (tag.errors > 4) ++stats.bpb_rejects;
        SpoolRecord spool;
        spool.sample = static_cast<std::uint16_t>(sample_index);
        spool.umi = tag.errors > 4 ? "REJECTS" : (tag.tag.empty() ? "NO_TAG" : tag.tag);
        spool.name = name.str(); spool.sequence = std::move(sample_sequence); spool.quality = std::move(sample_quality);
        spool.sampling_hash = stable_hash(config.parameters.deterministic_seed ^ sample_index, spool.name, spool.umi);
        std::vector<std::uint8_t> encoded; encode_spool_record(spool, encoded);
        const auto partition = stable_hash(sample_index, spool.umi) & (config.parameters.spool_partitions - 1);
        output.push_back(static_cast<std::uint8_t>(partition));
        binary::number(output, static_cast<std::uint32_t>(encoded.size()));
        output.insert(output.end(), encoded.begin(), encoded.end());
      }
    }
    if (!assigned) { if (had_primers) ++stats.id_rejects; else ++stats.primer_rejects; }
  }
  return output;
}

} // namespace webporpid
