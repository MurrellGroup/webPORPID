#include "webporpid/core.hpp"
#include "webporpid/binary.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <numeric>
#include <unordered_map>

#ifdef __wasm_simd128__
#include <wasm_simd128.h>
#endif

namespace webporpid {
namespace {

// Keep the dense scratch vectors off the small WASI stack.  A family centroid
// used to allocate and fill a 4,096-double vector twice per read.  The sparse
// representation below counts every read once, retains only observed bins,
// and evaluates exactly the same corrected squared distance from the dense
// mean.
struct SparseKmers {
  std::vector<std::pair<std::uint16_t, std::uint32_t>> bins;
  double total = 0.0;
};

int base_index(char value) {
  switch (std::toupper(static_cast<unsigned char>(value))) {
    case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1;
  }
}

SparseKmers kmer_counts(std::string_view sequence, std::vector<std::uint32_t>& scratch,
                        std::vector<std::uint16_t>& touched) {
  SparseKmers output; std::size_t code = 0, valid = 0;
  for (char value : sequence) {
    const int base = base_index(value);
    if (base < 0) { code = 0; valid = 0; continue; }
    code = ((code << 2) | static_cast<std::size_t>(base)) & 4095; ++valid;
    if (valid >= 6) {
      if (scratch[code]++ == 0) touched.push_back(static_cast<std::uint16_t>(code));
      output.total += 1.0;
    }
  }
  std::sort(touched.begin(), touched.end());
  output.bins.reserve(touched.size());
  for (const auto bin : touched) {
    output.bins.emplace_back(bin, scratch[bin]);
    scratch[bin] = 0;
  }
  touched.clear();
  return output;
}

std::string mode(const std::vector<std::string>& values) {
  if (values.empty()) return {};
  std::unordered_map<std::string, std::size_t> counts;
  std::string best = values.front(); std::size_t maximum = 1; counts[best] = 1;
  for (std::size_t index = 1; index < values.size(); ++index) {
    const auto count = ++counts[values[index]];
    if (count > maximum) { maximum = count; best = values[index]; }
  }
  return best;
}

std::string degap(std::string_view value) {
  std::string output; output.reserve(value.size());
  for (char base : value) if (base != '-') output.push_back(base);
  return output;
}

std::vector<std::size_t> reference_map(const Alignment& alignment) {
  std::vector<std::size_t> output; output.reserve(alignment.reference.size());
  for (std::size_t index = 0; index < alignment.reference.size(); ++index)
    if (alignment.reference[index] != '-') output.push_back(index);
  return output;
}

std::string centroid(const std::vector<std::string>& reads) {
  std::vector<double> mean(4096, 0.0);
  std::vector<std::uint32_t> scratch(4096, 0);
  std::vector<std::uint16_t> touched; touched.reserve(4096);
  std::vector<SparseKmers> profiles; profiles.reserve(reads.size());
  double mean_total = 0.0;
  for (const auto& read : reads) {
    auto counts = kmer_counts(read, scratch, touched);
    for (const auto& [bin, count] : counts.bins) mean[bin] += count;
    mean_total += counts.total;
    profiles.push_back(std::move(counts));
  }
  for (auto& value : mean) value /= reads.size();
  mean_total /= reads.size();
  double mean_squared = 0.0;
  for (const auto value : mean) mean_squared += value * value;
  std::size_t best = 0; double distance = std::numeric_limits<double>::infinity();
  for (std::size_t index = 0; index < profiles.size(); ++index) {
    double squared = mean_squared;
    for (const auto& [bin, count] : profiles[index].bins)
      squared += static_cast<double>(count) * count - 2.0 * count * mean[bin];
    const double total = profiles[index].total + mean_total;
    const double candidate = total == 0.0 ? 0.0 : squared / (6.0 * total);
    if (candidate < distance) { distance = candidate; best = index; }
  }
  return reads[best];
}

struct MappedAlignment { Alignment alignment; std::vector<std::size_t> map; };

std::vector<MappedAlignment> align_reads(std::string_view reference, const std::vector<std::string>& reads, bool seeded = true) {
  std::vector<MappedAlignment> output; output.reserve(reads.size());
  for (const auto& read : reads) {
    auto alignment = seeded ? seeded_global_align(reference, read) : needleman_wunsch(reference, read);
    auto map = reference_map(alignment); output.push_back({std::move(alignment), std::move(map)});
  }
  return output;
}

std::string extension_consensus(std::vector<std::string> values, bool front) {
  std::string consensus;
  const auto longest = [&] {
    std::size_t maximum = 0;
    for (const auto& value : values) maximum = std::max(maximum, value.size());
    return maximum;
  }();
  for (std::size_t iteration = 0; iteration <= longest
       && std::count_if(values.begin(), values.end(), [](const auto& value) { return !value.empty(); }) > values.size() / 2.0;
       ++iteration) {
    std::vector<std::string> active;
    for (const auto& value : values) if (!value.empty()) active.push_back(value);
    consensus = mode(active);
    auto alignments = align_reads(consensus, active, false);
    std::vector<std::string> next; next.reserve(alignments.size());
    for (const auto& mapped : alignments) {
      if (mapped.map.empty()) { next.emplace_back(); continue; }
      if (front) next.push_back(degap(std::string_view(mapped.alignment.query).substr(0, mapped.map.front())));
      else next.push_back(degap(std::string_view(mapped.alignment.query).substr(mapped.map.back() + 1)));
    }
    if (next == values) break;
    values = std::move(next);
  }
  return consensus;
}

std::string refine_reference(std::string_view candidate, const std::vector<std::string>& reads,
                             std::vector<MappedAlignment>* aligned_reads = nullptr) {
  if (candidate.size() < 2) return std::string(candidate);
  auto alignments = align_reads(candidate, reads);
  std::vector<bool> good(candidate.size(), true);
  for (std::size_t position = 0; position + 1 < candidate.size(); ++position) {
    std::size_t equal = 0;
    for (const auto& mapped : alignments) {
      if (mapped.map.size() <= position + 1) continue;
      char observed[2]{}; std::size_t observed_count = 0;
      for (std::size_t column = mapped.map[position]; column <= mapped.map[position + 1]; ++column) if (mapped.alignment.query[column] != '-') {
        if (observed_count < 2) observed[observed_count] = mapped.alignment.query[column];
        observed_count++;
      }
      if (observed_count == 2 && observed[0] == candidate[position] && observed[1] == candidate[position + 1]) ++equal;
    }
    if (static_cast<double>(equal) / reads.size() < 0.7) good[position] = good[position + 1] = false;
  }
  std::vector<std::string> fronts, ends;
  for (const auto& mapped : alignments) {
    fronts.push_back(mapped.map.empty() ? std::string{} : degap(std::string_view(mapped.alignment.query).substr(0, mapped.map.front())));
    ends.push_back(mapped.map.empty() ? std::string{} : degap(std::string_view(mapped.alignment.query).substr(mapped.map.back() + 1)));
  }
  const auto front = extension_consensus(std::move(fronts), true);
  const auto end = extension_consensus(std::move(ends), false);
  if (std::all_of(good.begin(), good.end(), [](bool value) { return value; })) {
    auto result = front + std::string(candidate) + end;
    if (aligned_reads) *aligned_reads = std::move(alignments);
    return result;
  }

  std::string rebuilt = front;
  std::size_t position = 0;
  while (position < candidate.size()) {
    if (good[position]) { rebuilt.push_back(candidate[position++]); continue; }
    const std::size_t start = position;
    while (position < candidate.size() && !good[position]) ++position;
    const std::size_t stop = position - 1;
    std::vector<std::string> alternatives; alternatives.reserve(alignments.size());
    for (const auto& mapped : alignments) {
      if (mapped.map.size() <= stop) alternatives.emplace_back();
      else alternatives.push_back(degap(std::string_view(mapped.alignment.query).substr(
          mapped.map[start], mapped.map[stop] - mapped.map[start] + 1)));
    }
    rebuilt += mode(alternatives);
  }
  rebuilt += end;
  if (aligned_reads) *aligned_reads = std::move(alignments);
  return rebuilt;
}

void agreement_summary(std::string_view candidate, const std::vector<std::string>& reads,
                       double& minimum_agreement, std::vector<LowAgreementSite>& low_sites,
                       const std::vector<MappedAlignment>* reusable = nullptr) {
  // RobustAmpliconDenoising.get_matches uses kmer_seeded_align here.  The old
  // port accidentally forced a full O(length^2) NW alignment for every read,
  // after already doing three seeded refinement passes.
  const auto owned_alignments = reusable ? std::vector<MappedAlignment>{} : align_reads(candidate, reads, true);
  const auto& alignments = reusable ? *reusable : owned_alignments;
  std::vector<double> agreements(candidate.size(), 0.0);
  std::vector<char> modal(candidate.size(), '-');
  for (std::size_t position = 0; position < candidate.size(); ++position) {
    // Symbol 0 represents the empty/gap observation; byte values occupy 1..256.
    std::array<std::uint32_t, 257> counts{}; std::size_t equal = 0;
    std::uint16_t best_symbol = 0; std::uint32_t best_count = 0;
    for (const auto& mapped : alignments) {
      std::uint16_t symbol = 0;
      if (mapped.map.size() > position && mapped.alignment.query[mapped.map[position]] != '-')
        symbol = static_cast<std::uint16_t>(std::toupper(static_cast<unsigned char>(mapped.alignment.query[mapped.map[position]]))) + 1;
      const auto count = ++counts[symbol]; if (count > best_count) { best_count = count; best_symbol = symbol; }
      if (symbol && symbol - 1 == std::toupper(static_cast<unsigned char>(candidate[position]))) ++equal;
    }
    agreements[position] = static_cast<double>(equal) / reads.size();
    if (best_symbol) modal[position] = static_cast<char>(best_symbol - 1);
  }
  if (agreements.size() > 2) minimum_agreement = *std::min_element(agreements.begin() + 1, agreements.end() - 1);
  minimum_agreement = std::round(minimum_agreement * 100.0) / 100.0;
  std::vector<std::size_t> runs(modal.size(), 0);
  for (std::size_t start = 0; start < modal.size();) {
    std::size_t end = start + 1; while (end < modal.size() && modal[end] == modal[start]) ++end;
    std::fill(runs.begin() + start, runs.begin() + end, end - start); start = end;
  }
  for (std::size_t position = 0; position < candidate.size(); ++position) {
    const auto rounded_agreement = std::round(agreements[position] * 100.0) / 100.0;
    if (rounded_agreement <= minimum_agreement && modal[position] != '-')
      low_sites.push_back({static_cast<std::uint32_t>(candidate.size() - position),
        static_cast<float>(rounded_agreement), complement(modal[position]),
        static_cast<std::uint32_t>(runs[position])});
  }
}

std::vector<std::uint64_t> cutoffs(std::span<const std::uint8_t> bytes, std::string& error) {
  binary::Reader reader(bytes); std::uint32_t count = 0;
  if (!reader.magic("WPT1") || !reader.number(count) || count > 65535) { error = "Invalid downsampling thresholds."; return {}; }
  std::vector<std::uint64_t> output(count);
  for (auto& value : output) if (!reader.number(value)) { error = "Truncated downsampling thresholds."; return {}; }
  if (!reader.done()) { error = "Downsampling thresholds have trailing bytes."; return {}; }
  return output;
}

bool selected(const SpoolRecord& record, const std::vector<std::uint64_t>& values) {
  return record.sample >= values.size() || record.sampling_hash <= values[record.sample];
}

std::string fixed_two(double value) {
  const auto hundredths = static_cast<unsigned>(std::lround(std::clamp(value, 0.0, 1.0) * 100.0));
  std::string output = std::to_string(hundredths / 100);
  output.push_back('.'); output.push_back(static_cast<char>('0' + (hundredths / 10) % 10));
  output.push_back(static_cast<char>('0' + hundredths % 10)); return output;
}

double mean(std::span<const double> values) {
  return values.empty() ? 0.0 : std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double ln_gamma(double value) {
  constexpr std::array<double, 9> coefficients{{0.9999999999998099, 676.5203681218851, -1259.1392167224028,
    771.32342877767613, -176.61502916214059, 12.507343278686905, -0.13857109526572012,
    9.9843695780195716e-6, 1.5056327351493116e-7}};
  if (value < 0.5) return std::log(std::numbers::pi) - std::log(std::sin(std::numbers::pi * value)) - ln_gamma(1.0 - value);
  const double z = value - 1.0; double x = coefficients[0];
  for (std::size_t index = 1; index < coefficients.size(); ++index) x += coefficients[index] / (z + index);
  const double t = z + 7.5;
  return 0.5 * std::log(2.0 * std::numbers::pi) + (z + 0.5) * std::log(t) - t + std::log(x);
}

double beta_fraction(double x, double a, double b) {
  double c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0); if (std::abs(d) < 1e-30) d = 1e-30;
  d = 1.0 / d; double h = d;
  for (int m = 1; m <= 200; ++m) {
    const double m2 = 2.0 * m;
    double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
    d = 1.0 + aa * d; if (std::abs(d) < 1e-30) d = 1e-30;
    c = 1.0 + aa / c; if (std::abs(c) < 1e-30) c = 1e-30;
    d = 1.0 / d; h *= d * c;
    aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
    d = 1.0 + aa * d; if (std::abs(d) < 1e-30) d = 1e-30;
    c = 1.0 + aa / c; if (std::abs(c) < 1e-30) c = 1e-30;
    d = 1.0 / d; const double delta = d * c; h *= delta;
    if (std::abs(delta - 1.0) < 3e-12) break;
  }
  return h;
}

double regularized_beta(double x, double a, double b) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  const double front = std::exp(ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * std::log(x) + b * std::log(1.0 - x));
  return x < (a + 1.0) / (a + b + 2.0) ? front * beta_fraction(x, a, b) / a
       : 1.0 - front * beta_fraction(1.0 - x, b, a) / b;
}

double t_test(std::span<const double> left, std::span<const double> right) {
  if (left.size() < 2 || right.size() < 2) return 1.0;
  const double lm = mean(left), rm = mean(right);
  double lv = 0.0, rv = 0.0;
  for (double value : left) lv += (value - lm) * (value - lm);
  for (double value : right) rv += (value - rm) * (value - rm);
  lv /= left.size() - 1; rv /= right.size() - 1;
  const double df = left.size() + right.size() - 2;
  const double pooled = std::sqrt(((left.size() - 1) * lv + (right.size() - 1) * rv) / df);
  if (pooled == 0.0) return lm == rm ? 1.0 : 0.0;
  const double t = std::abs(lm - rm) / (pooled * std::sqrt(1.0 / left.size() + 1.0 / right.size()));
  return regularized_beta(df / (df + t * t), df / 2.0, 0.5);
}

std::size_t sample_id_length(std::string_view primer) {
  std::size_t length = 0; bool started = false;
  for (char value : primer) {
    if (std::islower(static_cast<unsigned char>(value))) { started = true; ++length; }
    else if (started) break;
  }
  return length;
}

std::size_t umi_length(std::string_view primer) {
  const auto start = primer.find_first_of("Nn"); if (start == std::string_view::npos) return 0;
  const auto end = primer.find_first_not_of("Nn", start); return (end == std::string_view::npos ? primer.size() : end) - start;
}

bool heteroduplex(const std::vector<SpoolRecord>& reads, const Sample& sample) {
  if (reads.empty() || std::any_of(reads.begin(), reads.end(), [](const auto& read) { return read.quality.size() < 50; })) return false;
  std::array<double, 50> averages{};
  for (const auto& read : reads) for (std::size_t position = 0; position < 50; ++position)
    averages[position] += static_cast<unsigned char>(read.quality[position]) >= 33
      ? (static_cast<unsigned char>(read.quality[position]) - 33.0) / reads.size() : 0.0;
  const std::size_t start = sample_id_length(sample.cdna_primer), stop = start + umi_length(sample.cdna_primer);
  if (start >= stop || stop + 25 > 50) return false;
  const auto barcode = std::span(averages).subspan(start, stop - start);
  const auto downstream = std::span(averages).subspan(stop + 4, 21);
  if (mean(barcode) >= mean(downstream)) return false;
  return t_test(barcode, downstream) < 1.0 / (5.0 * reads.size() * reads.size())
      && *std::min_element(barcode.begin(), barcode.end()) < mean(downstream) / 2.0;
}

std::uint64_t read_ordinal(std::string_view name) {
  std::size_t position = name.starts_with("seq") ? 3 : 0;
  std::uint64_t value = 0; bool found = false;
  while (position < name.size() && std::isdigit(static_cast<unsigned char>(name[position]))) {
    found = true; value = value * 10 + static_cast<unsigned char>(name[position]) - '0'; ++position;
  }
  return found ? value : std::numeric_limits<std::uint64_t>::max();
}

Alignment iupac_primer_alignment(std::string_view reference, std::string_view query) {
  const auto mask = [](char value) {
    switch (std::toupper(static_cast<unsigned char>(value))) {
      case 'A': return 1; case 'C': return 2; case 'G': return 4; case 'T': return 8; case 'N': return 15;
      default: return 15;
    }
  };
  constexpr int gap = 100000, edge_gap = 50000, compatible_penalty = 1;
  const std::size_t rows = reference.size() + 1, columns = query.size() + 1;
  std::vector<int> previous(columns), current(columns); std::vector<std::uint8_t> trace(rows * columns);
  for (std::size_t column = 1; column < columns; ++column) { previous[column] = previous[column - 1] - edge_gap; trace[column] = 2; }
  for (std::size_t row = 1; row < rows; ++row) {
    current[0] = previous[0] - edge_gap; trace[row * columns] = 3;
    for (std::size_t column = 1; column < columns; ++column) {
      const auto left_base = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[row - 1])));
      const auto right_base = static_cast<char>(std::toupper(static_cast<unsigned char>(query[column - 1])));
      const int emission = left_base == right_base ? 0
        : ((mask(left_base) & mask(right_base)) ? -compatible_penalty : -gap);
      const int diagonal = previous[column - 1] + emission;
      const int left = current[column - 1] - (row == rows - 1 ? edge_gap : gap);
      const int up = previous[column] - (column == columns - 1 ? edge_gap : gap);
      if (diagonal >= left && diagonal >= up) { current[column] = diagonal; trace[row * columns + column] = 1; }
      else if (left >= up) { current[column] = left; trace[row * columns + column] = 2; }
      else { current[column] = up; trace[row * columns + column] = 3; }
    }
    previous.swap(current);
  }
  Alignment output; std::size_t row = reference.size(), column = query.size();
  while (row || column) {
    if (row == 0) { output.reference.push_back('-'); output.query.push_back(query[--column]); continue; }
    if (column == 0) { output.reference.push_back(reference[--row]); output.query.push_back('-'); continue; }
    const auto op = trace[row * columns + column];
    if (op == 1) { output.reference.push_back(reference[--row]); output.query.push_back(query[--column]); }
    else if (op == 2) { output.reference.push_back('-'); output.query.push_back(query[--column]); }
    else { output.reference.push_back(reference[--row]); output.query.push_back('-'); }
  }
  std::reverse(output.reference.begin(), output.reference.end()); std::reverse(output.query.begin(), output.query.end()); return output;
}

std::string primer_trim(std::string_view sequence, std::string_view primer) {
  if (sequence.empty() || primer.empty()) return std::string(sequence);
  if (sequence.size() >= primer.size()) {
    bool compatible = true;
    for (std::size_t index = 0; index < primer.size(); ++index) {
      const char expected = static_cast<char>(std::toupper(static_cast<unsigned char>(primer[index])));
      const char observed = static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[index])));
      if (expected != 'N' && expected != observed) { compatible = false; break; }
    }
    if (compatible) return std::string(sequence.substr(primer.size()));
  }
  const auto prefix = sequence.substr(0, std::min(sequence.size(), primer.size() + 3));
  const auto alignment = iupac_primer_alignment(primer, prefix);
  const auto last = alignment.reference.find_last_not_of('-');
  if (last == std::string::npos) return std::string(sequence);
  std::size_t consumed = 0; for (std::size_t index = 0; index <= last; ++index) if (alignment.query[index] != '-') ++consumed;
  return std::string(sequence.substr(std::min(consumed, sequence.size())));
}

} // namespace

std::string family_consensus(const std::vector<std::string>& reads, double& minimum_agreement,
                             std::vector<LowAgreementSite>& low_sites) {
  minimum_agreement = 0.0; low_sites.clear();
  if (reads.empty()) return {};
  if (reads.size() == 1) { minimum_agreement = 1.0; return reads.front(); }
  if (std::all_of(reads.begin() + 1, reads.end(), [&](const auto& read) { return read == reads.front(); })) {
    minimum_agreement = 1.0; return reads.front();
  }
  std::string candidate = centroid(reads);
  std::vector<MappedAlignment> reusable_alignments;
  for (int pass = 0; pass < 3; ++pass) {
    std::vector<MappedAlignment> pass_alignments;
    auto refined = refine_reference(candidate, reads, &pass_alignments);
    // Refinement is deterministic.  Once fixed, the remaining Julia-equivalent
    // passes can no longer alter either sequence or agreement behavior.
    if (refined == candidate) { reusable_alignments = std::move(pass_alignments); break; }
    candidate = std::move(refined);
  }
  agreement_summary(candidate, reads, minimum_agreement, low_sites,
    reusable_alignments.empty() ? nullptr : &reusable_alignments);
  return candidate;
}

std::vector<std::uint8_t> process_consensus_partition(std::span<const std::uint8_t> bytes,
                                                      std::span<const std::uint8_t> cutoff_bytes,
                                                      const std::vector<FamilyDecision>& model,
                                                      const Config& config, std::string& error) {
  const auto thresholds = cutoffs(cutoff_bytes, error); if (!error.empty()) return {};
  const auto records = decode_spool(bytes, error); if (!error.empty()) return {};
  std::vector<std::size_t> decisions_per_sample(config.samples.size());
  for (const auto& decision : model) if (decision.sample < decisions_per_sample.size()) decisions_per_sample[decision.sample]++;
  std::vector<std::unordered_map<std::string, const FamilyDecision*>> decisions(config.samples.size());
  for (std::size_t sample = 0; sample < decisions.size(); ++sample) decisions[sample].reserve(decisions_per_sample[sample]);
  for (const auto& decision : model) if (decision.sample < decisions.size()) decisions[decision.sample].emplace(decision.umi, &decision);
  // Callers sort decoded consensus records, so an ordered family map only adds
  // O(reads log families) work. Per-sample hash tables also avoid allocating a
  // concatenated sample/tag lookup key for every read. Within-family read order
  // is still restored explicitly below.
  std::vector<std::unordered_map<std::string, std::vector<SpoolRecord>>> grouped(config.samples.size());
  for (std::size_t sample = 0; sample < grouped.size(); ++sample) grouped[sample].reserve(decisions_per_sample[sample]);
  for (const auto& record : records) if (record.sample < grouped.size() && selected(record, thresholds)
      && decisions[record.sample].contains(record.umi)) grouped[record.sample][record.umi].push_back(record);
  std::vector<ConsensusRecord> consensuses;
  std::vector<std::pair<std::uint16_t, std::string>> heteroduplexes;
  for (std::size_t sample_index = 0; sample_index < grouped.size(); ++sample_index) for (auto& [umi, family_reads] : grouped[sample_index]) {
    const auto* decision = decisions[sample_index].at(umi);
    std::stable_sort(family_reads.begin(), family_reads.end(), [](const auto& left, const auto& right) {
      const auto left_ordinal = read_ordinal(left.name), right_ordinal = read_ordinal(right.name);
      return left_ordinal != right_ordinal ? left_ordinal < right_ordinal : left.name < right.name;
    });
    if (decision->disposition == FamilyDisposition::bpb_reject) continue;
    if (heteroduplex(family_reads, config.samples[sample_index])) {
      heteroduplexes.emplace_back(static_cast<std::uint16_t>(sample_index), umi); continue;
    }
    if (decision->disposition != FamilyDisposition::likely_real) continue;
    std::vector<std::string> sequences; sequences.reserve(family_reads.size());
    for (const auto& read : family_reads) sequences.push_back(read.sequence);
    ConsensusRecord record; record.sample = static_cast<std::uint16_t>(sample_index); record.umi = umi;
    record.family_size = family_reads.size();
    auto raw = family_consensus(sequences, record.minimum_agreement, record.low_sites);
    const auto& full_primer = config.samples[sample_index].cdna_primer;
    const auto start = std::find_if(full_primer.begin(), full_primer.end(), [](char value) { return std::islower(static_cast<unsigned char>(value)); });
    std::string trim_primer(start, full_primer.end());
    std::transform(trim_primer.begin(), trim_primer.end(), trim_primer.begin(), [](unsigned char value) { return std::toupper(value); });
    record.sequence = reverse_complement(primer_trim(raw, trim_primer));
    record.id = config.samples[sample_index].name + umi + " fs=" + std::to_string(record.family_size)
      + " minag=" + fixed_two(record.minimum_agreement);
    consensuses.push_back(std::move(record));
  }
  std::vector<std::uint8_t> output; binary::magic(output, "WPO1");
  binary::number(output, static_cast<std::uint32_t>(consensuses.size()));
  for (const auto& record : consensuses) {
    binary::number(output, record.sample); binary::string(output, record.id); binary::string(output, record.umi);
    binary::number(output, record.family_size); binary::number(output, record.minimum_agreement); binary::string(output, record.sequence);
    binary::number(output, static_cast<std::uint32_t>(record.low_sites.size()));
    for (const auto& site : record.low_sites) {
      binary::number(output, site.position); binary::number(output, site.agreement);
      binary::number(output, static_cast<std::uint8_t>(site.modal_base)); binary::number(output, site.run_length);
    }
  }
  binary::number(output, static_cast<std::uint32_t>(heteroduplexes.size()));
  for (const auto& [sample, umi] : heteroduplexes) { binary::number(output, sample); binary::string(output, umi); }
  return output;
}

} // namespace webporpid
