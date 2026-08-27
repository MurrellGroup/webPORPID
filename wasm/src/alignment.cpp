#include "webporpid/core.hpp"

#include <algorithm>
#include <cctype>
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
  std::unordered_map<std::string_view, Location> reference_words, query_words;
  for (std::size_t position = 0; position + word <= reference.size(); position += skip) {
    const auto value = reference.substr(position, word);
    auto [iterator, inserted] = reference_words.emplace(value, Location{position, true});
    if (!inserted) iterator->second.unique = false;
  }
  for (std::size_t position = 0; position + word <= query.size(); ++position) {
    const auto value = query.substr(position, word);
    auto [iterator, inserted] = query_words.emplace(value, Location{position, true});
    if (!inserted) iterator->second.unique = false;
  }
  std::vector<std::pair<std::size_t, std::size_t>> matches;
  for (const auto& [value, left] : reference_words) {
    if (!left.unique) continue;
    const auto found = query_words.find(value);
    if (found != query_words.end() && found->second.unique) matches.emplace_back(left.position, found->second.position);
  }
  if (matches.empty()) return needleman_wunsch(reference, query);
  std::sort(matches.begin(), matches.end());
  if (inconsistent(matches)) matches = longest_increasing_subsequence(matches);
  if (matches.empty() || inconsistent(matches)) return needleman_wunsch(reference, query);

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
