#include "webporpid/core.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace webporpid {
namespace {

void append(Alignment& target, Alignment source) {
  target.reference += source.reference; target.query += source.query;
}

std::vector<std::pair<std::size_t, std::size_t>> longest_increasing_subsequence(
    const std::vector<std::pair<std::size_t, std::size_t>>& matches) {
  if (matches.empty()) return {};
  std::vector<std::size_t> previous(matches.size(), 0), current(matches.size(), 0);
  std::size_t maximum = 0;
  for (std::size_t index = 0; index < matches.size(); ++index) {
    std::size_t low = 0, high = maximum;
    while (low < high) {
      const std::size_t middle = (low + high) / 2;
      if (matches[current[middle]].second < matches[index].second) low = middle + 1;
      else high = middle;
    }
    previous[index] = low ? current[low - 1] : 0;
    current[low] = index; maximum = std::max(maximum, low + 1);
  }
  std::vector<std::pair<std::size_t, std::size_t>> output(maximum);
  std::size_t index = current[maximum - 1];
  for (std::size_t position = maximum; position-- > 0;) { output[position] = matches[index]; index = previous[index]; }
  return output;
}

bool inconsistent(const std::vector<std::pair<std::size_t, std::size_t>>& matches) {
  for (std::size_t index = 1; index < matches.size(); ++index)
    if (matches[index].second <= matches[index - 1].second) return true;
  return false;
}

Alignment banded_needleman_wunsch(std::string_view reference, std::string_view query) {
  if (reference.empty()) return {std::string(query.size(), '-'), std::string(query)};
  if (query.empty()) return {std::string(reference), std::string(reference.size(), '-')};
  const std::size_t rows = reference.size(), columns = query.size();
  const auto difference = static_cast<std::ptrdiff_t>(rows) - static_cast<std::ptrdiff_t>(columns);
  const std::ptrdiff_t bandwidth = static_cast<std::ptrdiff_t>(std::ceil(std::sqrt((rows + columns) / 2.0)));
  const std::ptrdiff_t lower = -bandwidth - std::max<std::ptrdiff_t>(0, difference);
  const std::ptrdiff_t upper = bandwidth + std::max<std::ptrdiff_t>(0, -difference);
  const std::size_t band_columns = static_cast<std::size_t>(upper - lower + 1);
  constexpr std::int32_t impossible = std::numeric_limits<std::int32_t>::min() / 4;
  std::vector<std::uint8_t> trace((rows + 1) * band_columns, 0);
  std::vector<std::int32_t> previous(columns + 1, impossible), current(columns + 1, impossible);
  auto range = [&](std::size_t row) {
    const auto signed_row = static_cast<std::ptrdiff_t>(row);
    return std::pair<std::size_t, std::size_t>{
      static_cast<std::size_t>(std::max<std::ptrdiff_t>(0, signed_row + lower)),
      static_cast<std::size_t>(std::min<std::ptrdiff_t>(columns, signed_row + upper))};
  };
  auto trace_at = [&](std::size_t row, std::size_t column) -> std::uint8_t& {
    const auto band = static_cast<std::ptrdiff_t>(column) - static_cast<std::ptrdiff_t>(row) - lower;
    return trace[row * band_columns + static_cast<std::size_t>(band)];
  };
  auto [previous_min, previous_max] = range(0);
  previous[0] = 0;
  for (std::size_t column = 1; column <= previous_max; ++column) { previous[column] = previous[column - 1] - 99; trace_at(0, column) = 2; }
  for (std::size_t row = 1; row <= rows; ++row) {
    const auto [current_min, current_max] = range(row);
    std::fill(current.begin() + current_min, current.begin() + current_max + 1, impossible);
    if (current_min == 0) { current[0] = previous_min == 0 && previous[0] != impossible ? previous[0] - 99 : -99 * static_cast<std::int32_t>(row); trace_at(row, 0) = 3; }
    for (std::size_t column = std::max<std::size_t>(1, current_min); column <= current_max; ++column) {
      const bool equal = std::toupper(static_cast<unsigned char>(reference[row - 1])) == std::toupper(static_cast<unsigned char>(query[column - 1]));
      const auto diagonal = column - 1 >= previous_min && column - 1 <= previous_max ? previous[column - 1] + (equal ? 100 : -100) : impossible;
      const auto left = column > current_min && current[column - 1] != impossible ? current[column - 1] - (row == rows ? 99 : 100) : impossible;
      const auto up = column >= previous_min && column <= previous_max && previous[column] != impossible ? previous[column] - (column == columns ? 99 : 100) : impossible;
      if (diagonal >= left && diagonal >= up) { current[column] = diagonal; trace_at(row, column) = 1; }
      else if (left >= up) { current[column] = left; trace_at(row, column) = 2; }
      else { current[column] = up; trace_at(row, column) = 3; }
    }
    previous.swap(current); previous_min = current_min; previous_max = current_max;
  }
  if (previous[columns] == impossible) return needleman_wunsch(reference, query);
  Alignment output; output.reference.reserve(rows + columns); output.query.reserve(rows + columns);
  std::size_t row = rows, column = columns;
  while (row || column) {
    if (row == 0) { output.reference.push_back('-'); output.query.push_back(query[--column]); continue; }
    if (column == 0) { output.reference.push_back(reference[--row]); output.query.push_back('-'); continue; }
    const auto operation = trace_at(row, column);
    if (operation == 1) { output.reference.push_back(reference[--row]); output.query.push_back(query[--column]); }
    else if (operation == 2) { output.reference.push_back('-'); output.query.push_back(query[--column]); }
    else if (operation == 3) { output.reference.push_back(reference[--row]); output.query.push_back('-'); }
    else return needleman_wunsch(reference, query);
  }
  std::reverse(output.reference.begin(), output.reference.end()); std::reverse(output.query.begin(), output.query.end()); return output;
}

} // namespace

Alignment needleman_wunsch(std::string_view reference, std::string_view query) {
  const std::size_t rows = reference.size() + 1, columns = query.size() + 1;
  if (reference.empty()) return {std::string(query.size(), '-'), std::string(query)};
  if (query.empty()) return {std::string(reference), std::string(reference.size(), '-')};
  std::vector<std::uint8_t> trace(rows * columns, 0);
  std::vector<std::int32_t> previous(columns, 0), current(columns, 0);
  for (std::size_t column = 1; column < columns; ++column) { previous[column] = previous[column - 1] - 99; trace[column] = 2; }
  for (std::size_t row = 1; row < rows; ++row) {
    current[0] = previous[0] - 99; trace[row * columns] = 3;
    for (std::size_t column = 1; column < columns; ++column) {
      const bool equal = std::toupper(static_cast<unsigned char>(reference[row - 1]))
                       == std::toupper(static_cast<unsigned char>(query[column - 1]));
      const auto diagonal = previous[column - 1] + (equal ? 100 : -100);
      // Julia calls horizontal movement `del`, and applies the lower-edge
      // multiplier when the reference is at its final row.
      const auto left = current[column - 1] - (row == rows - 1 ? 99 : 100);
      const auto up = previous[column] - (column == columns - 1 ? 99 : 100);
      if (diagonal >= left && diagonal >= up) { current[column] = diagonal; trace[row * columns + column] = 1; }
      else if (left >= up) { current[column] = left; trace[row * columns + column] = 2; }
      else { current[column] = up; trace[row * columns + column] = 3; }
    }
    previous.swap(current);
  }
  Alignment output; output.reference.reserve(reference.size() + query.size()); output.query.reserve(reference.size() + query.size());
  std::size_t row = reference.size(), column = query.size();
  while (row || column) {
    if (row == 0) { output.reference.push_back('-'); output.query.push_back(query[--column]); continue; }
    if (column == 0) { output.reference.push_back(reference[--row]); output.query.push_back('-'); continue; }
    const auto operation = trace[row * columns + column];
    if (operation == 1) { output.reference.push_back(reference[--row]); output.query.push_back(query[--column]); }
    else if (operation == 2) { output.reference.push_back('-'); output.query.push_back(query[--column]); }
    else { output.reference.push_back(reference[--row]); output.query.push_back('-'); }
  }
  std::reverse(output.reference.begin(), output.reference.end()); std::reverse(output.query.begin(), output.query.end());
  return output;
}

Alignment seeded_global_align(std::string_view reference, std::string_view query) {
  constexpr std::size_t word = 30, skip = 10;
  if (reference.size() < word || query.size() < word) return needleman_wunsch(reference, query);
  struct Location { std::size_t position = 0; bool unique = true; };
  std::vector<std::pair<std::size_t, std::size_t>> matches;
  const auto canonical = [](std::string_view sequence) {
    return std::all_of(sequence.begin(), sequence.end(), [](char value) { return value == 'A' || value == 'C' || value == 'G' || value == 'T'; });
  };
  if (canonical(reference) && canonical(query)) {
    // Thirty canonical bases fit exactly in 60 bits.  This avoids hashing 30
    // bytes for every overlapping query word while remaining collision-free.
    constexpr std::uint64_t mask = (std::uint64_t{1} << (2 * word)) - 1;
    const auto base = [](char value) -> std::uint64_t { return value == 'A' ? 0 : value == 'C' ? 1 : value == 'G' ? 2 : 3; };
    const auto encoded = [&](std::string_view sequence, std::size_t position) {
      std::uint64_t code = 0; for (std::size_t index = 0; index < word; ++index) code = (code << 2) | base(sequence[position + index]); return code;
    };
    std::unordered_map<std::uint64_t, Location> reference_words, query_words;
    reference_words.reserve(reference.size() / skip + 1); query_words.reserve(query.size() - word + 1);
    for (std::size_t position = 0; position + word <= reference.size(); position += skip) {
      const auto value = encoded(reference, position);
      auto [iterator, inserted] = reference_words.emplace(value, Location{position, true}); if (!inserted) iterator->second.unique = false;
    }
    std::uint64_t code = encoded(query, 0);
    for (std::size_t position = 0; position + word <= query.size(); ++position) {
      if (position) code = ((code << 2) | base(query[position + word - 1])) & mask;
      auto [iterator, inserted] = query_words.emplace(code, Location{position, true}); if (!inserted) iterator->second.unique = false;
    }
    for (const auto& [value, left] : reference_words) {
      if (!left.unique) continue;
      const auto found = query_words.find(value);
      if (found != query_words.end() && found->second.unique) matches.emplace_back(left.position, found->second.position);
    }
  } else {
    // Preserve the original byte-exact seed behavior for ambiguous/lowercase
    // input; the canonical hot path above covers normal FASTQ consensus work.
    std::unordered_map<std::string_view, Location> reference_words, query_words;
    reference_words.reserve(reference.size() / skip + 1); query_words.reserve(query.size() - word + 1);
    for (std::size_t position = 0; position + word <= reference.size(); position += skip) {
      const auto value = reference.substr(position, word);
      auto [iterator, inserted] = reference_words.emplace(value, Location{position, true}); if (!inserted) iterator->second.unique = false;
    }
    for (std::size_t position = 0; position + word <= query.size(); ++position) {
      const auto value = query.substr(position, word);
      auto [iterator, inserted] = query_words.emplace(value, Location{position, true}); if (!inserted) iterator->second.unique = false;
    }
    for (const auto& [value, left] : reference_words) {
      if (!left.unique) continue;
      const auto found = query_words.find(value);
      if (found != query_words.end() && found->second.unique) matches.emplace_back(left.position, found->second.position);
    }
  }
  if (matches.empty()) return banded_needleman_wunsch(reference, query);
  std::sort(matches.begin(), matches.end());
  if (inconsistent(matches)) matches = longest_increasing_subsequence(matches);
  if (matches.empty() || inconsistent(matches)) return banded_needleman_wunsch(reference, query);

  std::vector<std::pair<std::size_t, std::size_t>> cleaned{matches.front()};
  for (std::size_t index = 1; index < matches.size(); ++index) {
    const auto dr = matches[index].first - cleaned.back().first;
    const auto dq = matches[index].second - cleaned.back().second;
    if ((skip < dr && dr < word) || (dq != skip && dq < word)) continue;
    cleaned.push_back(matches[index]);
  }
  struct Range { std::size_t first, last; };
  std::vector<Range> ranges{{0, 0}};
  for (std::size_t index = 1; index < cleaned.size(); ++index) {
    if (cleaned[index].first - cleaned[index - 1].first == skip
        && cleaned[index].second - cleaned[index - 1].second == skip) ranges.back().last = index;
    else ranges.push_back({index, index});
  }

  Alignment output;
  std::size_t previous_reference = 0, previous_query = 0;
  for (const auto& range : ranges) {
    const auto start_reference = cleaned[range.first].first;
    const auto start_query = cleaned[range.first].second;
    append(output, needleman_wunsch(reference.substr(previous_reference, start_reference - previous_reference),
                                     query.substr(previous_query, start_query - previous_query)));
    const auto end_reference = cleaned[range.last].first + word;
    const auto end_query = cleaned[range.last].second + word;
    output.reference.append(reference.substr(start_reference, end_reference - start_reference));
    output.query.append(query.substr(start_query, end_query - start_query));
    previous_reference = end_reference; previous_query = end_query;
  }
  append(output, needleman_wunsch(reference.substr(previous_reference), query.substr(previous_query)));
  return output;
}

} // namespace webporpid
