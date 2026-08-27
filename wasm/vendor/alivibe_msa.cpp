#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// A port of the refinedMSA path in MurrellGroup/WebWidgets revision
// cbcd02719dd0a5f1f05d3127666f00e8579f2423. The unusual scoring, order
// dependence, strict-greater-than tie breaking, POA column packing, and three
// refinement passes are deliberate compatibility constraints. Swig's
// nucleotide and amino-acid entry points each differ in one intentional way:
// N or X, respectively, is a wildcard substitution match. The legacy entry
// point retains literal scoring for differential compatibility tests.

namespace {

constexpr float kGraphNegative = -1.0e9F;
char active_ambiguity_wildcard = '\0';

bool substitution_match(char left, char right, bool nucleotide_n_wildcard) {
    return left == right || (
        nucleotide_n_wildcard
        && active_ambiguity_wildcard != '\0'
        && (left == active_ambiguity_wildcard || right == active_ambiguity_wildcard)
    );
}

struct PairAlignment {
    std::string first;
    std::string second;
};

struct Match {
    int i = 0;
    int j = 0;
    int length = 0;
};

std::string javascript_substring(const std::string& value, int start, int end) {
    const int length = static_cast<int>(value.size());
    start = std::max(0, std::min(start, length));
    end = std::max(0, std::min(end, length));
    if (start > end) std::swap(start, end);
    return value.substr(static_cast<std::size_t>(start), static_cast<std::size_t>(end - start));
}

std::pair<double, std::uint8_t> find_max_3(double a, double b, double c) {
    double value = a;
    std::uint8_t index = 1;
    if (b > value) {
        value = b;
        index = 2;
    }
    if (c > value) {
        value = c;
        index = 3;
    }
    return {value, index};
}

std::pair<double, std::uint8_t> find_max_2(double a, double b) {
    if (b > a) return {b, 2};
    return {a, 1};
}

PairAlignment affine_nw_align(
    const std::string& s1,
    const std::string& s2,
    double gap_open = -10.0,
    double gap_extend = -0.2,
    double match_cost = 1.0,
    double mismatch_cost = -0.7,
    double boundary_gap_factor = 10.0,
    bool nucleotide_n_wildcard = false
) {
    const int n = static_cast<int>(s1.size());
    const int m = static_cast<int>(s2.size());
    const std::size_t width = static_cast<std::size_t>(m + 1);
    const std::size_t cells = static_cast<std::size_t>(n + 1) * width;
    const auto at = [width](int i, int j) { return static_cast<std::size_t>(i) * width + static_cast<std::size_t>(j); };

    std::vector<double> matrix_m(cells, 0.0);
    std::vector<double> matrix_ix(cells, 0.0);
    std::vector<double> matrix_iy(cells, 0.0);
    for (int i = 0; i <= n; ++i) {
        matrix_ix[at(i, 0)] = gap_open + gap_extend * i;
        matrix_iy[at(i, 0)] = -std::numeric_limits<double>::infinity();
        matrix_m[at(i, 0)] = gap_open + gap_extend * i;
    }
    for (int j = 0; j <= m; ++j) {
        matrix_ix[at(0, j)] = -std::numeric_limits<double>::infinity();
        matrix_iy[at(0, j)] = gap_open + gap_extend * j;
        matrix_m[at(0, j)] = gap_open + gap_extend * j;
    }

    std::vector<std::uint8_t> trace_m(cells, 0);
    std::vector<std::uint8_t> trace_ix(cells, 0);
    std::vector<std::uint8_t> trace_iy(cells, 0);
    for (int j = 0; j <= m; ++j) {
        trace_m[at(0, j)] = 3;
        trace_ix[at(0, j)] = 3;
        trace_iy[at(0, j)] = 3;
    }
    for (int i = 1; i <= n; ++i) {
        trace_m[at(i, 0)] = 2;
        trace_ix[at(i, 0)] = 2;
        trace_iy[at(i, 0)] = 2;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const double substitution = substitution_match(
                s1[static_cast<std::size_t>(i - 1)],
                s2[static_cast<std::size_t>(j - 1)],
                nucleotide_n_wildcard
            )
                ? match_cost : mismatch_cost;
            auto [m_value, m_trace] = find_max_3(
                matrix_m[at(i - 1, j - 1)] + substitution,
                matrix_ix[at(i - 1, j - 1)] + substitution,
                matrix_iy[at(i - 1, j - 1)] + substitution
            );
            matrix_m[at(i, j)] = m_value;
            trace_m[at(i, j)] = m_trace;

            const double extend_x = j == m ? gap_extend / boundary_gap_factor : gap_extend;
            auto [ix_value, ix_trace] = find_max_2(
                matrix_m[at(i - 1, j)] + gap_open,
                matrix_ix[at(i - 1, j)] + extend_x
            );
            matrix_ix[at(i, j)] = ix_value;
            trace_ix[at(i, j)] = ix_trace;

            const double extend_y = i == n ? gap_extend / boundary_gap_factor : gap_extend;
            auto [iy_value, iy_compact_trace] = find_max_2(
                matrix_m[at(i, j - 1)] + gap_open,
                matrix_iy[at(i, j - 1)] + extend_y
            );
            matrix_iy[at(i, j)] = iy_value;
            trace_iy[at(i, j)] = iy_compact_trace == 1 ? 1 : 3;
        }
    }

    std::string reverse_1;
    std::string reverse_2;
    reverse_1.reserve(static_cast<std::size_t>(n + m));
    reverse_2.reserve(static_cast<std::size_t>(n + m));
    int x = n;
    int y = m;
    std::uint8_t state = find_max_3(matrix_m[at(x, y)], matrix_ix[at(x, y)], matrix_iy[at(x, y)]).second;
    while (x > 0 && y > 0) {
        std::uint8_t next = 0;
        if (state == 1) next = trace_m[at(x, y)];
        else if (state == 2) next = trace_ix[at(x, y)];
        else next = trace_iy[at(x, y)];
        if (state == 1) {
            reverse_1.push_back(s1[static_cast<std::size_t>(x - 1)]);
            reverse_2.push_back(s2[static_cast<std::size_t>(y - 1)]);
            --x;
            --y;
        } else if (state == 2) {
            reverse_1.push_back(s1[static_cast<std::size_t>(x - 1)]);
            reverse_2.push_back('-');
            --x;
        } else {
            reverse_1.push_back('-');
            reverse_2.push_back(s2[static_cast<std::size_t>(y - 1)]);
            --y;
        }
        state = next;
    }
    while (x > 0) {
        reverse_1.push_back(s1[static_cast<std::size_t>(--x)]);
        reverse_2.push_back('-');
    }
    while (y > 0) {
        reverse_1.push_back('-');
        reverse_2.push_back(s2[static_cast<std::size_t>(--y)]);
    }
    std::reverse(reverse_1.begin(), reverse_1.end());
    std::reverse(reverse_2.begin(), reverse_2.end());
    return {std::move(reverse_1), std::move(reverse_2)};
}

PairAlignment constrained_nw_align(
    const std::string& s1,
    const std::string& s2,
    double gap_open,
    double gap_extend,
    double match_cost,
    double mismatch_cost,
    double boundary_gap_factor,
    bool /* is_start */,
    bool is_end,
    bool nucleotide_n_wildcard
) {
    const int n = static_cast<int>(s1.size());
    const int m = static_cast<int>(s2.size());
    const double use_boundary_factor = is_end ? boundary_gap_factor : 1.0;
    const std::size_t width = static_cast<std::size_t>(m + 1);
    const std::size_t cells = static_cast<std::size_t>(n + 1) * width;
    const auto at = [width](int i, int j) { return static_cast<std::size_t>(i) * width + static_cast<std::size_t>(j); };

    std::vector<double> matrix_m(cells, 0.0);
    std::vector<double> matrix_ix(cells, 0.0);
    std::vector<double> matrix_iy(cells, 0.0);
    matrix_m[0] = 0.0;
    matrix_ix[0] = -std::numeric_limits<double>::infinity();
    matrix_iy[0] = -std::numeric_limits<double>::infinity();
    for (int i = 1; i <= n; ++i) {
        matrix_ix[at(i, 0)] = gap_open + gap_extend * i;
        matrix_iy[at(i, 0)] = -std::numeric_limits<double>::infinity();
        matrix_m[at(i, 0)] = gap_open + gap_extend * i;
    }
    for (int j = 1; j <= m; ++j) {
        matrix_ix[at(0, j)] = -std::numeric_limits<double>::infinity();
        matrix_iy[at(0, j)] = gap_open + gap_extend * j;
        matrix_m[at(0, j)] = gap_open + gap_extend * j;
    }

    std::vector<std::uint8_t> trace_m(cells, 0);
    std::vector<std::uint8_t> trace_ix(cells, 0);
    std::vector<std::uint8_t> trace_iy(cells, 0);
    for (int j = 0; j <= m; ++j) {
        trace_m[at(0, j)] = 3;
        trace_ix[at(0, j)] = 3;
        trace_iy[at(0, j)] = 3;
    }
    for (int i = 1; i <= n; ++i) {
        trace_m[at(i, 0)] = 2;
        trace_ix[at(i, 0)] = 2;
        trace_iy[at(i, 0)] = 2;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const double substitution = substitution_match(
                s1[static_cast<std::size_t>(i - 1)],
                s2[static_cast<std::size_t>(j - 1)],
                nucleotide_n_wildcard
            )
                ? match_cost : mismatch_cost;
            auto [m_value, m_trace] = find_max_3(
                matrix_m[at(i - 1, j - 1)] + substitution,
                matrix_ix[at(i - 1, j - 1)] + substitution,
                matrix_iy[at(i - 1, j - 1)] + substitution
            );
            matrix_m[at(i, j)] = m_value;
            trace_m[at(i, j)] = m_trace;

            const double extend_x = j == m ? gap_extend / use_boundary_factor : gap_extend;
            auto [ix_value, ix_trace] = find_max_2(
                matrix_m[at(i - 1, j)] + gap_open,
                matrix_ix[at(i - 1, j)] + extend_x
            );
            matrix_ix[at(i, j)] = ix_value;
            trace_ix[at(i, j)] = ix_trace;

            const double extend_y = i == n ? gap_extend / use_boundary_factor : gap_extend;
            auto [iy_value, iy_compact_trace] = find_max_2(
                matrix_m[at(i, j - 1)] + gap_open,
                matrix_iy[at(i, j - 1)] + extend_y
            );
            matrix_iy[at(i, j)] = iy_value;
            trace_iy[at(i, j)] = iy_compact_trace == 1 ? 1 : 3;
        }
    }

    std::string reverse_1;
    std::string reverse_2;
    reverse_1.reserve(static_cast<std::size_t>(n + m));
    reverse_2.reserve(static_cast<std::size_t>(n + m));
    int x = n;
    int y = m;
    std::uint8_t state = find_max_3(matrix_m[at(x, y)], matrix_ix[at(x, y)], matrix_iy[at(x, y)]).second;
    while (x > 0 && y > 0) {
        std::uint8_t next = state == 1 ? trace_m[at(x, y)] : state == 2 ? trace_ix[at(x, y)] : trace_iy[at(x, y)];
        if (state == 1) {
            reverse_1.push_back(s1[static_cast<std::size_t>(x - 1)]);
            reverse_2.push_back(s2[static_cast<std::size_t>(y - 1)]);
            --x;
            --y;
        } else if (state == 2) {
            reverse_1.push_back(s1[static_cast<std::size_t>(x - 1)]);
            reverse_2.push_back('-');
            --x;
        } else {
            reverse_1.push_back('-');
            reverse_2.push_back(s2[static_cast<std::size_t>(y - 1)]);
            --y;
        }
        state = next;
    }
    while (x > 0) {
        reverse_1.push_back(s1[static_cast<std::size_t>(--x)]);
        reverse_2.push_back('-');
    }
    while (y > 0) {
        reverse_1.push_back('-');
        reverse_2.push_back(s2[static_cast<std::size_t>(--y)]);
    }
    std::reverse(reverse_1.begin(), reverse_1.end());
    std::reverse(reverse_2.begin(), reverse_2.end());
    return {std::move(reverse_1), std::move(reverse_2)};
}

std::vector<Match> merge_and_erode_anchors(const std::vector<Match>& matches, int margin) {
    if (matches.empty()) return {};
    std::vector<Match> merged;
    Match current = matches.front();
    for (std::size_t k = 1; k < matches.size(); ++k) {
        const Match& next = matches[k];
        const int current_diagonal = current.j - current.i;
        const int next_diagonal = next.j - next.i;
        const bool connected = next.i <= current.i + current.length;
        if (current_diagonal == next_diagonal && connected) {
            const int new_end = std::max(current.i + current.length, next.i + next.length);
            current.length = new_end - current.i;
        } else {
            merged.push_back(current);
            current = next;
        }
    }
    merged.push_back(current);

    std::vector<Match> eroded;
    eroded.reserve(merged.size());
    for (const Match& match : merged) {
        const int new_length = match.length - 2 * margin;
        if (new_length > 0) eroded.push_back({match.i + margin, match.j + margin, new_length});
    }
    return eroded;
}

PairAlignment stitch_alignments(
    const std::string& s1,
    const std::string& s2,
    const std::vector<Match>& anchors,
    double gap_open,
    double gap_extend,
    double match_cost,
    double mismatch_cost,
    double boundary_gap_factor,
    bool nucleotide_n_wildcard
) {
    std::string final_1;
    std::string final_2;
    int index_1 = 0;
    int index_2 = 0;
    for (std::size_t k = 0; k < anchors.size(); ++k) {
        const Match& anchor = anchors[k];
        const std::string sub_1 = javascript_substring(s1, index_1, anchor.i);
        const std::string sub_2 = javascript_substring(s2, index_2, anchor.j);
        if (!sub_1.empty() || !sub_2.empty()) {
            PairAlignment gap = constrained_nw_align(
                sub_1, sub_2, gap_open, gap_extend, match_cost, mismatch_cost,
                boundary_gap_factor, k == 0 && index_1 == 0, false, nucleotide_n_wildcard
            );
            final_1 += gap.first;
            final_2 += gap.second;
        }
        const std::string anchor_sequence = s1.substr(static_cast<std::size_t>(anchor.i), static_cast<std::size_t>(anchor.length));
        final_1 += anchor_sequence;
        final_2 += anchor_sequence;
        index_1 = anchor.i + anchor.length;
        index_2 = anchor.j + anchor.length;
    }
    if (index_1 < static_cast<int>(s1.size()) || index_2 < static_cast<int>(s2.size())) {
        PairAlignment tail = constrained_nw_align(
            s1.substr(static_cast<std::size_t>(index_1)),
            s2.substr(static_cast<std::size_t>(index_2)),
            gap_open, gap_extend, match_cost, mismatch_cost, boundary_gap_factor, false, true,
            nucleotide_n_wildcard
        );
        final_1 += tail.first;
        final_2 += tail.second;
    }
    return {std::move(final_1), std::move(final_2)};
}

struct KmerEntry {
    std::string value;
    int count = 0;
    int last_index = 0;
};

std::vector<KmerEntry> ordered_kmers(const std::string& sequence, int kmer_size) {
    std::vector<KmerEntry> entries;
    if (static_cast<int>(sequence.size()) < kmer_size) return entries;
    std::unordered_map<std::string, std::size_t> lookup;
    lookup.reserve(sequence.size());
    for (int i = 0; i <= static_cast<int>(sequence.size()) - kmer_size; ++i) {
        std::string kmer = sequence.substr(static_cast<std::size_t>(i), static_cast<std::size_t>(kmer_size));
        auto found = lookup.find(kmer);
        if (found == lookup.end()) {
            const std::size_t index = entries.size();
            entries.push_back({std::move(kmer), 1, i});
            lookup.emplace(entries.back().value, index);
        } else {
            KmerEntry& entry = entries[found->second];
            ++entry.count;
            entry.last_index = i;
        }
    }
    return entries;
}

std::vector<Match> simple_lis(const std::vector<Match>& matches) {
    if (matches.empty()) return {};
    std::vector<int> tails;
    std::vector<int> parent(matches.size(), -1);
    for (int i = 0; i < static_cast<int>(matches.size()); ++i) {
        const int value = matches[static_cast<std::size_t>(i)].j;
        int left = 0;
        int right = static_cast<int>(tails.size());
        while (left < right) {
            const int middle = (left + right) >> 1;
            if (matches[static_cast<std::size_t>(tails[static_cast<std::size_t>(middle)])].j < value) left = middle + 1;
            else right = middle;
        }
        if (left < static_cast<int>(tails.size())) {
            tails[static_cast<std::size_t>(left)] = i;
            parent[static_cast<std::size_t>(i)] = left > 0 ? tails[static_cast<std::size_t>(left - 1)] : -1;
        } else {
            tails.push_back(i);
            parent[static_cast<std::size_t>(i)] = tails.size() > 1 ? tails[tails.size() - 2] : -1;
        }
    }
    std::vector<Match> result;
    for (int current = tails.back(); current != -1; current = parent[static_cast<std::size_t>(current)]) {
        result.push_back(matches[static_cast<std::size_t>(current)]);
    }
    std::reverse(result.begin(), result.end());
    return result;
}

PairAlignment double_dp_nw_align(
    const std::string& s1,
    const std::string& s2,
    double gap_open = -10.0,
    double gap_extend = -0.2,
    double match_cost = 1.0,
    double mismatch_cost = -0.7,
    double boundary_gap_factor = 10.0,
    bool nucleotide_n_wildcard = false
) {
    constexpr int kmer_size = 11;
    if (std::min(s1.size(), s2.size()) < static_cast<std::size_t>(3 * kmer_size)) {
        return affine_nw_align(
            s1, s2, gap_open, gap_extend, match_cost, mismatch_cost,
            boundary_gap_factor, nucleotide_n_wildcard
        );
    }

    std::unordered_map<std::string, std::vector<int>> positions;
    positions.reserve(s1.size());
    for (int i = 0; i <= static_cast<int>(s1.size()) - kmer_size; ++i) {
        positions[s1.substr(static_cast<std::size_t>(i), kmer_size)].push_back(i);
    }
    std::vector<Match> matches;
    for (int j = 0; j <= static_cast<int>(s2.size()) - kmer_size; ++j) {
        auto found = positions.find(s2.substr(static_cast<std::size_t>(j), kmer_size));
        if (found == positions.end() || found->second.size() > 25) continue;
        for (int i : found->second) matches.push_back({i, j, kmer_size});
    }
    if (matches.empty()) {
        return affine_nw_align(
            s1, s2, gap_open, gap_extend, match_cost, mismatch_cost,
            boundary_gap_factor, nucleotide_n_wildcard
        );
    }
    std::stable_sort(matches.begin(), matches.end(), [](const Match& left, const Match& right) {
        if (left.i != right.i) return left.i < right.i;
        return left.j < right.j;
    });

    std::vector<double> scores(matches.size(), 0.0);
    std::vector<int> parents(matches.size(), -1);
    for (int current = 0; current < static_cast<int>(matches.size()); ++current) {
        const Match& current_match = matches[static_cast<std::size_t>(current)];
        double maximum = current_match.length;
        const int start = std::max(0, current - 80);
        for (int previous = current - 1; previous >= start; --previous) {
            const Match& previous_match = matches[static_cast<std::size_t>(previous)];
            if (previous_match.i >= current_match.i || previous_match.j >= current_match.j) continue;
            const int diagonal_difference = std::abs(
                (current_match.j - current_match.i) - (previous_match.j - previous_match.i)
            );
            const int gap_i = current_match.i - (previous_match.i + previous_match.length);
            const int gap_j = current_match.j - (previous_match.j + previous_match.length);
            const int distance = std::max(0, gap_i) + std::max(0, gap_j);
            const double penalty = diagonal_difference * 3.0 + distance * 0.1;
            const double candidate = scores[static_cast<std::size_t>(previous)] + current_match.length - penalty;
            if (candidate > maximum) {
                maximum = candidate;
                parents[static_cast<std::size_t>(current)] = previous;
            }
        }
        scores[static_cast<std::size_t>(current)] = maximum;
    }
    int best = 0;
    double best_score = scores.front();
    for (int i = 1; i < static_cast<int>(scores.size()); ++i) {
        if (scores[static_cast<std::size_t>(i)] > best_score) {
            best_score = scores[static_cast<std::size_t>(i)];
            best = i;
        }
    }
    std::vector<Match> chain;
    for (int current = best; current != -1; current = parents[static_cast<std::size_t>(current)]) {
        chain.push_back(matches[static_cast<std::size_t>(current)]);
    }
    std::reverse(chain.begin(), chain.end());
    return stitch_alignments(
        s1, s2, merge_and_erode_anchors(chain, kmer_size),
        gap_open, gap_extend, match_cost, mismatch_cost, boundary_gap_factor,
        nucleotide_n_wildcard
    );
}

class SequenceSet {
public:
    void add(int sequence) {
        const std::size_t word = static_cast<std::size_t>(sequence) >> 6U;
        if (words_.size() <= word) words_.resize(word + 1, 0);
        const std::uint64_t bit = std::uint64_t{1} << (static_cast<unsigned>(sequence) & 63U);
        if ((words_[word] & bit) == 0) {
            words_[word] |= bit;
            ++size_;
        }
    }

    [[nodiscard]] bool contains(int sequence) const {
        const std::size_t word = static_cast<std::size_t>(sequence) >> 6U;
        if (word >= words_.size()) return false;
        return (words_[word] & (std::uint64_t{1} << (static_cast<unsigned>(sequence) & 63U))) != 0;
    }

    [[nodiscard]] bool intersects(const SequenceSet& other) const {
        const std::size_t limit = std::min(words_.size(), other.words_.size());
        for (std::size_t i = 0; i < limit; ++i) {
            if ((words_[i] & other.words_[i]) != 0) return true;
        }
        return false;
    }

    [[nodiscard]] int size() const { return size_; }

private:
    std::vector<std::uint64_t> words_;
    int size_ = 0;
};

struct PoaNode {
    char character = 0;
    std::vector<int> next;
    std::vector<int> previous;
    SequenceSet sequences;
};

class PoaGraph {
public:
    int create_node(char character) {
        nodes.push_back(PoaNode{character, {}, {}, {}});
        return static_cast<int>(nodes.size()) - 1;
    }

    void add_edge(int from, int to) {
        auto& next = nodes[static_cast<std::size_t>(from)].next;
        if (std::find(next.begin(), next.end(), to) == next.end()) next.push_back(to);
        auto& previous = nodes[static_cast<std::size_t>(to)].previous;
        if (std::find(previous.begin(), previous.end(), from) == previous.end()) previous.push_back(from);
    }

    void initialize_first(const std::string& sequence) {
        sequences.push_back(sequence);
        int previous = -1;
        for (char character : sequence) {
            const int node = create_node(character);
            nodes[static_cast<std::size_t>(node)].sequences.add(0);
            if (previous != -1) add_edge(previous, node);
            previous = node;
        }
    }

    [[nodiscard]] std::vector<int> topological_sort() const {
        std::vector<std::uint8_t> visited(nodes.size(), 0);
        std::vector<int> postorder;
        postorder.reserve(nodes.size());
        struct Frame { int node; std::size_t next_index; };
        std::vector<Frame> stack;
        for (int start = 0; start < static_cast<int>(nodes.size()); ++start) {
            if (visited[static_cast<std::size_t>(start)]) continue;
            visited[static_cast<std::size_t>(start)] = 1;
            stack.push_back({start, 0});
            while (!stack.empty()) {
                Frame& frame = stack.back();
                const auto& next = nodes[static_cast<std::size_t>(frame.node)].next;
                if (frame.next_index < next.size()) {
                    const int child = next[frame.next_index++];
                    if (!visited[static_cast<std::size_t>(child)]) {
                        visited[static_cast<std::size_t>(child)] = 1;
                        stack.push_back({child, 0});
                    }
                    continue;
                }
                postorder.push_back(frame.node);
                stack.pop_back();
            }
        }
        std::reverse(postorder.begin(), postorder.end());
        return postorder;
    }

    struct Consensus {
        std::string sequence;
        std::vector<int> mapping;
    };

    [[nodiscard]] Consensus consensus_path() const {
        if (nodes.empty()) return {};
        const std::vector<int> sorted = topological_sort();
        std::vector<int> path_score(nodes.size(), 0);
        std::vector<int> best_previous(nodes.size(), -1);
        for (int node_id : sorted) {
            const PoaNode& node = nodes[static_cast<std::size_t>(node_id)];
            int maximum_previous_score = 0;
            int best = -1;
            for (int previous : node.previous) {
                if (path_score[static_cast<std::size_t>(previous)] > maximum_previous_score) {
                    maximum_previous_score = path_score[static_cast<std::size_t>(previous)];
                    best = previous;
                }
            }
            path_score[static_cast<std::size_t>(node_id)] = node.sequences.size() + maximum_previous_score;
            best_previous[static_cast<std::size_t>(node_id)] = best;
        }
        int best_end = 0;
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            if (path_score[static_cast<std::size_t>(i)] > path_score[static_cast<std::size_t>(best_end)]) best_end = i;
        }
        std::vector<int> mapping;
        for (int current = best_end; current != -1; current = best_previous[static_cast<std::size_t>(current)]) {
            mapping.push_back(current);
        }
        std::reverse(mapping.begin(), mapping.end());
        std::string sequence;
        sequence.reserve(mapping.size());
        for (int node : mapping) sequence.push_back(nodes[static_cast<std::size_t>(node)].character);
        return {std::move(sequence), std::move(mapping)};
    }

    std::vector<PoaNode> nodes;
    std::vector<std::string> sequences;
};

std::vector<Match> sequence_to_graph_anchors(const std::string& sequence, const std::string& consensus) {
    constexpr int kmer_size = 15;
    if (sequence.size() < kmer_size || consensus.size() < kmer_size) return {};
    const std::vector<KmerEntry> first = ordered_kmers(sequence, kmer_size);
    const std::vector<KmerEntry> second = ordered_kmers(consensus, kmer_size);
    std::unordered_map<std::string, const KmerEntry*> second_lookup;
    second_lookup.reserve(second.size());
    for (const KmerEntry& entry : second) second_lookup.emplace(entry.value, &entry);
    std::vector<Match> matches;
    for (const KmerEntry& entry : first) {
        auto found = second_lookup.find(entry.value);
        if (entry.count == 1 && found != second_lookup.end() && found->second->count == 1) {
            matches.push_back({entry.last_index, found->second->last_index, kmer_size});
        }
    }
    std::stable_sort(matches.begin(), matches.end(), [](const Match& left, const Match& right) {
        return left.i < right.i;
    });
    return merge_and_erode_anchors(simple_lis(matches), 2);
}

std::vector<int> subgraph_between(
    const PoaGraph& graph,
    int start,
    int end,
    const std::vector<int>& topological,
    const std::vector<int>& topological_index
) {
    if (start >= static_cast<int>(topological_index.size())) return {};
    const int start_index = start == -1 ? -1 : topological_index[static_cast<std::size_t>(start)];
    const int end_index = end == -1 ? std::numeric_limits<int>::max() : topological_index[static_cast<std::size_t>(end)];
    std::vector<int> result;
    for (int node : topological) {
        const int index = topological_index[static_cast<std::size_t>(node)];
        if (index > start_index && index < end_index) result.push_back(node);
    }
    return result;
}

enum class GraphOperationType : std::uint8_t { Match, Insert, Delete };

struct GraphOperation {
    GraphOperationType type = GraphOperationType::Insert;
    char character = 0;
    int node = -1;
};

std::vector<GraphOperation> sequence_to_graph_dp(
    const std::string& sequence,
    const std::vector<int>& graph_nodes,
    const PoaGraph& graph,
    double gap_open,
    double gap_extend,
    double match_cost,
    double mismatch_cost,
    bool nucleotide_n_wildcard
) {
    const int sequence_length = static_cast<int>(sequence.size());
    const int node_count = static_cast<int>(graph_nodes.size());
    std::vector<int> real_to_sub(graph.nodes.size(), -1);
    for (int i = 0; i < node_count; ++i) real_to_sub[static_cast<std::size_t>(graph_nodes[static_cast<std::size_t>(i)])] = i;
    const std::size_t cells = static_cast<std::size_t>(sequence_length + 1) * static_cast<std::size_t>(node_count);
    std::vector<float> score(cells, kGraphNegative);
    std::vector<std::uint8_t> trace(cells, 0);
    std::vector<int> trace_previous(cells, -1);
    const auto at = [node_count](int i, int node) {
        return static_cast<std::size_t>(i) * static_cast<std::size_t>(node_count) + static_cast<std::size_t>(node);
    };

    for (int u = 0; u < node_count; ++u) {
        const int real = graph_nodes[static_cast<std::size_t>(u)];
        const PoaNode& node = graph.nodes[static_cast<std::size_t>(real)];
        double maximum = kGraphNegative;
        int best = -1;
        bool connected_to_start = false;
        for (int previous : node.previous) {
            if (previous >= static_cast<int>(real_to_sub.size()) || real_to_sub[static_cast<std::size_t>(previous)] == -1) {
                connected_to_start = true;
                break;
            }
        }
        if (connected_to_start) maximum = gap_open;
        for (int previous : node.previous) {
            if (previous >= static_cast<int>(real_to_sub.size())) continue;
            const int previous_sub = real_to_sub[static_cast<std::size_t>(previous)];
            if (previous_sub == -1) continue;
            const float previous_score = score[at(0, previous_sub)];
            if (previous_score > kGraphNegative) {
                const double candidate = static_cast<double>(previous_score) + gap_extend;
                if (candidate > maximum) {
                    maximum = candidate;
                    best = previous_sub;
                }
            }
        }
        score[at(0, u)] = static_cast<float>(maximum);
        trace[at(0, u)] = 3;
        trace_previous[at(0, u)] = best;
    }

    for (int i = 1; i <= sequence_length; ++i) {
        const char character = sequence[static_cast<std::size_t>(i - 1)];
        for (int u = 0; u < node_count; ++u) {
            const int real = graph_nodes[static_cast<std::size_t>(u)];
            const PoaNode& node = graph.nodes[static_cast<std::size_t>(real)];
            const double substitution = substitution_match(character, node.character, nucleotide_n_wildcard)
                ? match_cost : mismatch_cost;
            bool connected_to_start = false;
            for (int previous : node.previous) {
                if (previous >= static_cast<int>(real_to_sub.size()) || real_to_sub[static_cast<std::size_t>(previous)] == -1) {
                    connected_to_start = true;
                    break;
                }
            }

            double diagonal = kGraphNegative;
            int diagonal_previous = -1;
            if (connected_to_start && i == 1) diagonal = substitution;
            for (int previous : node.previous) {
                if (previous >= static_cast<int>(real_to_sub.size())) continue;
                const int previous_sub = real_to_sub[static_cast<std::size_t>(previous)];
                if (previous_sub == -1) continue;
                const float previous_score = score[at(i - 1, previous_sub)];
                if (previous_score > kGraphNegative) {
                    const double candidate = static_cast<double>(previous_score) + substitution;
                    if (candidate > diagonal) {
                        diagonal = candidate;
                        diagonal_previous = previous_sub;
                    }
                }
            }
            const double insertion = static_cast<double>(score[at(i - 1, u)]) + gap_extend;
            double deletion = kGraphNegative;
            int deletion_previous = -1;
            for (int previous : node.previous) {
                if (previous >= static_cast<int>(real_to_sub.size())) continue;
                const int previous_sub = real_to_sub[static_cast<std::size_t>(previous)];
                if (previous_sub == -1) continue;
                const float previous_score = score[at(i, previous_sub)];
                if (previous_score > kGraphNegative) {
                    const double candidate = static_cast<double>(previous_score) + gap_extend;
                    if (candidate > deletion) {
                        deletion = candidate;
                        deletion_previous = previous_sub;
                    }
                }
            }

            double best_score = diagonal;
            std::uint8_t type = 1;
            int previous = diagonal_previous;
            if (deletion > best_score) {
                best_score = deletion;
                type = 3;
                previous = deletion_previous;
            }
            if (insertion > best_score) {
                best_score = insertion;
                type = 2;
                previous = u;
            }
            score[at(i, u)] = static_cast<float>(best_score);
            trace[at(i, u)] = type;
            trace_previous[at(i, u)] = previous;
        }
    }

    double maximum_end = kGraphNegative;
    int current_node = -1;
    for (int u = 0; u < node_count; ++u) {
        if (score[at(sequence_length, u)] > maximum_end) {
            maximum_end = score[at(sequence_length, u)];
            current_node = u;
        }
    }
    std::vector<GraphOperation> path;
    int current_sequence = sequence_length;
    while (current_sequence > 0 || current_node != -1) {
        if (current_node == -1) {
            path.push_back({GraphOperationType::Insert, sequence[static_cast<std::size_t>(current_sequence - 1)], -1});
            --current_sequence;
            continue;
        }
        const std::uint8_t type = trace[at(current_sequence, current_node)];
        const int previous = trace_previous[at(current_sequence, current_node)];
        if (type == 1) {
            path.push_back({
                GraphOperationType::Match,
                sequence[static_cast<std::size_t>(current_sequence - 1)],
                graph_nodes[static_cast<std::size_t>(current_node)]
            });
            --current_sequence;
            current_node = previous;
        } else if (type == 2) {
            path.push_back({
                GraphOperationType::Insert,
                sequence[static_cast<std::size_t>(current_sequence - 1)],
                graph_nodes[static_cast<std::size_t>(current_node)]
            });
            --current_sequence;
        } else if (type == 3) {
            path.push_back({GraphOperationType::Delete, 0, graph_nodes[static_cast<std::size_t>(current_node)]});
            current_node = previous;
        } else if (current_sequence > 0) {
            path.push_back({GraphOperationType::Insert, sequence[static_cast<std::size_t>(current_sequence - 1)], -1});
            --current_sequence;
        } else {
            break;
        }
    }
    std::reverse(path.begin(), path.end());
    return path;
}

int apply_alignment_to_graph(PoaGraph& graph, const std::vector<GraphOperation>& path, int start_node, int sequence_index) {
    int current = start_node;
    for (const GraphOperation& operation : path) {
        if (operation.type == GraphOperationType::Match) {
            const bool same = graph.nodes[static_cast<std::size_t>(operation.node)].character == operation.character;
            if (same) {
                graph.nodes[static_cast<std::size_t>(operation.node)].sequences.add(sequence_index);
                if (current != -1 && current != operation.node) graph.add_edge(current, operation.node);
                current = operation.node;
            } else {
                const int created = graph.create_node(operation.character);
                graph.nodes[static_cast<std::size_t>(created)].sequences.add(sequence_index);
                if (current != -1) graph.add_edge(current, created);
                current = created;
            }
        } else if (operation.type == GraphOperationType::Insert) {
            const int created = graph.create_node(operation.character);
            graph.nodes[static_cast<std::size_t>(created)].sequences.add(sequence_index);
            if (current != -1) graph.add_edge(current, created);
            current = created;
        }
    }
    return current;
}

void add_sequence_to_graph(PoaGraph& graph, const std::string& sequence, bool nucleotide_n_wildcard) {
    const int sequence_index = static_cast<int>(graph.sequences.size());
    graph.sequences.push_back(sequence);
    const PoaGraph::Consensus consensus = graph.consensus_path();
    const std::vector<Match> anchors = sequence_to_graph_anchors(sequence, consensus.sequence);
    struct Constraint { int sequence_start; int sequence_end; int node_start; };
    std::vector<Constraint> constraints;
    constraints.reserve(anchors.size());
    for (const Match& anchor : anchors) {
        constraints.push_back({
            anchor.i,
            anchor.i + anchor.length,
            consensus.mapping[static_cast<std::size_t>(anchor.j)]
        });
    }

    int sequence_position = 0;
    int previous_node = -1;
    const std::vector<int> topological = graph.topological_sort();
    std::vector<int> topological_index(graph.nodes.size(), -1);
    for (int i = 0; i < static_cast<int>(topological.size()); ++i) {
        topological_index[static_cast<std::size_t>(topological[static_cast<std::size_t>(i)])] = i;
    }

    for (std::size_t k = 0; k <= constraints.size(); ++k) {
        const int block_sequence_end = k < constraints.size()
            ? constraints[k].sequence_start : static_cast<int>(sequence.size());
        const int block_node_target = k < constraints.size() ? constraints[k].node_start : -1;
        if (sequence_position < block_sequence_end || (previous_node != -1 && block_node_target != -1)) {
            const std::string sub_sequence = javascript_substring(sequence, sequence_position, block_sequence_end);
            const std::vector<int> subgraph = subgraph_between(
                graph, previous_node, block_node_target, topological, topological_index
            );
            if (!sub_sequence.empty() || !subgraph.empty()) {
                const std::vector<GraphOperation> path = sequence_to_graph_dp(
                    sub_sequence, subgraph, graph, -10.0, -1.0, 1.0, -1.0,
                    nucleotide_n_wildcard
                );
                previous_node = apply_alignment_to_graph(graph, path, previous_node, sequence_index);
            }
        }

        if (k < constraints.size()) {
            const Constraint& constraint = constraints[k];
            auto found = std::find(consensus.mapping.begin(), consensus.mapping.end(), constraint.node_start);
            const int consensus_start = static_cast<int>(std::distance(consensus.mapping.begin(), found));
            const int length = constraint.sequence_end - constraint.sequence_start;
            for (int i = 0; i < length; ++i) {
                const char character = sequence[static_cast<std::size_t>(constraint.sequence_start + i)];
                const int target = consensus.mapping[static_cast<std::size_t>(consensus_start + i)];
                graph.nodes[static_cast<std::size_t>(target)].sequences.add(sequence_index);
                if (previous_node != -1 && previous_node != target) graph.add_edge(previous_node, target);
                previous_node = target;
                (void)character; // The source assumes exact seed identity and fuses directly.
            }
            sequence_position = constraint.sequence_end;
        }
    }
}

std::vector<std::string> generate_msa_from_graph(const PoaGraph& graph) {
    const std::vector<int> topological = graph.topological_sort();
    std::vector<std::vector<int>> columns;
    std::vector<int> node_to_column(graph.nodes.size(), -1);
    for (int node_id : topological) {
        const PoaNode& node = graph.nodes[static_cast<std::size_t>(node_id)];
        int column = 0;
        for (int previous : node.previous) {
            if (node_to_column[static_cast<std::size_t>(previous)] >= column) {
                column = node_to_column[static_cast<std::size_t>(previous)] + 1;
            }
        }
        while (true) {
            if (column >= static_cast<int>(columns.size())) {
                columns.emplace_back();
                break;
            }
            bool conflict = false;
            for (int existing : columns[static_cast<std::size_t>(column)]) {
                if (node.sequences.intersects(graph.nodes[static_cast<std::size_t>(existing)].sequences)) {
                    conflict = true;
                    break;
                }
            }
            if (!conflict) break;
            ++column;
        }
        columns[static_cast<std::size_t>(column)].push_back(node_id);
        node_to_column[static_cast<std::size_t>(node_id)] = column;
    }

    std::vector<std::string> result(
        graph.sequences.size(),
        std::string(columns.size(), '-')
    );
    for (int node_id = 0; node_id < static_cast<int>(graph.nodes.size()); ++node_id) {
        const PoaNode& node = graph.nodes[static_cast<std::size_t>(node_id)];
        const int column = node_to_column[static_cast<std::size_t>(node_id)];
        for (int sequence = 0; sequence < static_cast<int>(graph.sequences.size()); ++sequence) {
            if (node.sequences.contains(sequence)) {
                result[static_cast<std::size_t>(sequence)][static_cast<std::size_t>(column)] = node.character;
            }
        }
    }
    return result;
}

std::vector<std::string> multiple_sequence_align(
    const std::vector<std::string>& sequences,
    bool nucleotide_n_wildcard
) {
    if (sequences.empty()) return {};
    if (sequences.size() == 1) return sequences;
    PoaGraph graph;
    graph.initialize_first(sequences.front());
    for (std::size_t i = 1; i < sequences.size(); ++i) {
        add_sequence_to_graph(graph, sequences[i], nucleotide_n_wildcard);
    }
    return generate_msa_from_graph(graph);
}

std::vector<std::string> remove_gap_only_columns(const std::vector<std::string>& msa) {
    if (msa.empty()) return {};
    const std::size_t length = msa.front().size();
    std::vector<std::uint8_t> keep(length, 0);
    for (std::size_t column = 0; column < length; ++column) {
        for (const std::string& sequence : msa) {
            if (sequence[column] != '-') {
                keep[column] = 1;
                break;
            }
        }
    }
    std::vector<std::string> result(msa.size());
    for (std::string& sequence : result) sequence.reserve(length);
    for (std::size_t column = 0; column < length; ++column) {
        if (!keep[column]) continue;
        for (std::size_t row = 0; row < msa.size(); ++row) {
            result[row].push_back(msa[row][column]);
        }
    }
    return result;
}

struct RefinementBlock {
    std::size_t start = 0;
    std::size_t end = 0;
    bool refine = false;
};

std::vector<RefinementBlock> identify_refinement_blocks(const std::vector<std::string>& msa, int flank_size) {
    const std::size_t length = msa.front().size();
    if (length == 0) return {{0, 0, false}};
    std::vector<std::uint8_t> unstable(length, 0);
    for (std::size_t column = 0; column < length; ++column) {
        for (const std::string& sequence : msa) {
            if (sequence[column] == '-') {
                unstable[column] = 1;
                break;
            }
        }
    }
    std::vector<std::uint8_t> mask(length, 0);
    for (std::size_t column = 0; column < length; ++column) {
        if (!unstable[column]) continue;
        const std::size_t start = column > static_cast<std::size_t>(flank_size)
            ? column - static_cast<std::size_t>(flank_size) : 0;
        const std::size_t end = std::min(length, column + static_cast<std::size_t>(flank_size) + 1);
        for (std::size_t k = start; k < end; ++k) mask[k] = 1;
    }
    std::vector<RefinementBlock> blocks;
    std::size_t current_start = 0;
    bool current_refine = mask[0] == 1;
    for (std::size_t column = 1; column < length; ++column) {
        const bool refine = mask[column] == 1;
        if (refine == current_refine) continue;
        blocks.push_back({current_start, column, current_refine});
        current_start = column;
        current_refine = refine;
    }
    blocks.push_back({current_start, length, current_refine});
    return blocks;
}

std::vector<std::string> refine_block_slice(
    const std::vector<std::string>& slice,
    bool nucleotide_n_wildcard
) {
    const int sequence_count = static_cast<int>(slice.size());
    std::vector<std::string> raw(slice.size());
    bool every_empty = true;
    for (std::size_t row = 0; row < slice.size(); ++row) {
        raw[row].reserve(slice[row].size());
        for (char character : slice[row]) {
            if (character != '-') raw[row].push_back(character);
        }
        if (!raw[row].empty()) every_empty = false;
    }
    if (every_empty) return slice;

    int best_index = 0;
    int maximum_length = -1;
    for (int i = 0; i < sequence_count; ++i) {
        const int length = static_cast<int>(raw[static_cast<std::size_t>(i)].size());
        if (length > maximum_length) {
            maximum_length = length;
            best_index = i;
        }
    }
    const std::string& reference = raw[static_cast<std::size_t>(best_index)];
    const int reference_length = static_cast<int>(reference.size());
    std::vector<std::vector<char>> matches(
        static_cast<std::size_t>(reference_length),
        std::vector<char>(slice.size(), '-')
    );
    std::vector<std::vector<std::string>> insertions(
        static_cast<std::size_t>(reference_length + 1),
        std::vector<std::string>(slice.size())
    );

    for (int row = 0; row < sequence_count; ++row) {
        if (row == best_index) {
            for (int k = 0; k < reference_length; ++k) {
                matches[static_cast<std::size_t>(k)][static_cast<std::size_t>(row)] = reference[static_cast<std::size_t>(k)];
            }
            continue;
        }
        const std::string& query = raw[static_cast<std::size_t>(row)];
        if (query.empty()) continue;
        const PairAlignment aligned = double_dp_nw_align(
            reference, query, -10.0, -0.2, 1.0, -0.7, 10.0,
            nucleotide_n_wildcard
        );
        int reference_position = 0;
        std::string active_insertion;
        for (std::size_t k = 0; k < aligned.first.size(); ++k) {
            const char reference_character = aligned.first[k];
            const char query_character = aligned.second[k];
            if (reference_character != '-') {
                if (!active_insertion.empty()) {
                    insertions[static_cast<std::size_t>(reference_position)][static_cast<std::size_t>(row)] = std::move(active_insertion);
                    active_insertion.clear();
                }
                matches[static_cast<std::size_t>(reference_position)][static_cast<std::size_t>(row)] = query_character;
                ++reference_position;
            } else if (query_character != '-') {
                active_insertion.push_back(query_character);
            }
        }
        if (!active_insertion.empty()) {
            insertions[static_cast<std::size_t>(reference_position)][static_cast<std::size_t>(row)] = std::move(active_insertion);
        }
    }

    std::vector<std::string> result(slice.size());
    for (int position = 0; position <= reference_length; ++position) {
        std::size_t maximum_insertion = 0;
        for (int row = 0; row < sequence_count; ++row) {
            maximum_insertion = std::max(
                maximum_insertion,
                insertions[static_cast<std::size_t>(position)][static_cast<std::size_t>(row)].size()
            );
        }
        for (std::size_t insertion_position = 0; insertion_position < maximum_insertion; ++insertion_position) {
            for (int row = 0; row < sequence_count; ++row) {
                const std::string& insertion = insertions[static_cast<std::size_t>(position)][static_cast<std::size_t>(row)];
                result[static_cast<std::size_t>(row)].push_back(
                    insertion_position < insertion.size() ? insertion[insertion_position] : '-'
                );
            }
        }
        if (position < reference_length) {
            for (int row = 0; row < sequence_count; ++row) {
                result[static_cast<std::size_t>(row)].push_back(
                    matches[static_cast<std::size_t>(position)][static_cast<std::size_t>(row)]
                );
            }
        }
    }
    return result;
}

std::vector<std::string> refined_msa(
    const std::vector<std::string>& sequences,
    int iterations,
    bool nucleotide_n_wildcard
) {
    std::vector<std::string> current = multiple_sequence_align(sequences, nucleotide_n_wildcard);
    if (current.size() <= 1) return current;
    for (int iteration = 0; iteration < iterations; ++iteration) {
        current = remove_gap_only_columns(current);
        const std::vector<RefinementBlock> blocks = identify_refinement_blocks(current, 20);
        std::vector<std::string> next(current.size());
        for (const RefinementBlock& block : blocks) {
            std::vector<std::string> block_slice(current.size());
            for (std::size_t row = 0; row < current.size(); ++row) {
                block_slice[row] = current[row].substr(block.start, block.end - block.start);
            }
            const std::vector<std::string> resolved = block.refine
                ? refine_block_slice(block_slice, nucleotide_n_wildcard)
                : block_slice;
            for (std::size_t row = 0; row < next.size(); ++row) next[row] += resolved[row];
        }
        current = std::move(next);
    }
    return remove_gap_only_columns(current);
}

std::vector<std::uint8_t> result_bytes;
std::string error_text;

bool read_u32(const std::uint8_t*& cursor, const std::uint8_t* end, std::uint32_t& value) {
    if (static_cast<std::size_t>(end - cursor) < 4) return false;
    value = static_cast<std::uint32_t>(cursor[0])
        | (static_cast<std::uint32_t>(cursor[1]) << 8U)
        | (static_cast<std::uint32_t>(cursor[2]) << 16U)
        | (static_cast<std::uint32_t>(cursor[3]) << 24U);
    cursor += 4;
    return true;
}

void write_u32(std::vector<std::uint8_t>& target, std::uint32_t value) {
    target.push_back(static_cast<std::uint8_t>(value));
    target.push_back(static_cast<std::uint8_t>(value >> 8U));
    target.push_back(static_cast<std::uint8_t>(value >> 16U));
    target.push_back(static_cast<std::uint8_t>(value >> 24U));
}

bool decode_sequences(const std::uint8_t* input, std::size_t length, std::vector<std::string>& sequences) {
    if (!input || length < 8 || std::memcmp(input, "AMSA", 4) != 0) {
        error_text = "Invalid Alivibe MSA request.";
        return false;
    }
    const std::uint8_t* cursor = input + 4;
    const std::uint8_t* end = input + length;
    std::uint32_t count = 0;
    if (!read_u32(cursor, end, count) || count > 10000) {
        error_text = "Invalid Alivibe MSA sequence count.";
        return false;
    }
    sequences.clear();
    sequences.reserve(count);
    std::uint64_t total = 0;
    for (std::uint32_t i = 0; i < count; ++i) {
        std::uint32_t sequence_length = 0;
        if (!read_u32(cursor, end, sequence_length) || sequence_length > static_cast<std::uint32_t>(end - cursor)) {
            error_text = "Truncated Alivibe MSA request.";
            return false;
        }
        total += sequence_length;
        if (total > 256ULL * 1024ULL * 1024ULL) {
            error_text = "Alivibe MSA input is too large.";
            return false;
        }
        sequences.emplace_back(reinterpret_cast<const char*>(cursor), sequence_length);
        cursor += sequence_length;
    }
    if (cursor != end) {
        error_text = "Unexpected bytes after the Alivibe MSA request.";
        return false;
    }
    if (sequences.size() > 1 && sequences.front().empty()) {
        error_text = "The first Alivibe MSA sequence is empty.";
        return false;
    }
    return true;
}

bool encode_sequences(const std::vector<std::string>& sequences) {
    std::uint64_t size = 8 + sequences.size() * 4ULL;
    for (const std::string& sequence : sequences) size += sequence.size();
    if (size > std::numeric_limits<std::uint32_t>::max()) {
        error_text = "Alivibe MSA output is too large.";
        return false;
    }
    result_bytes.clear();
    result_bytes.reserve(static_cast<std::size_t>(size));
    result_bytes.insert(result_bytes.end(), {'A', 'M', 'S', 'A'});
    write_u32(result_bytes, static_cast<std::uint32_t>(sequences.size()));
    for (const std::string& sequence : sequences) {
        write_u32(result_bytes, static_cast<std::uint32_t>(sequence.size()));
        result_bytes.insert(result_bytes.end(), sequence.begin(), sequence.end());
    }
    return true;
}

}  // namespace

#define ALIVIBE_EXPORT(name) __attribute__((export_name(name)))

extern "C" {

ALIVIBE_EXPORT("alivibe_msa_alloc") std::uint8_t* alivibe_msa_alloc(std::uint32_t size) {
    return static_cast<std::uint8_t*>(std::malloc(size));
}

ALIVIBE_EXPORT("alivibe_msa_free") void alivibe_msa_free(std::uint8_t* pointer) {
    std::free(pointer);
}

std::int32_t run_msa(
    const std::uint8_t* input,
    std::uint32_t length,
    std::int32_t iterations,
    char ambiguity_wildcard
) {
    error_text.clear();
    result_bytes.clear();
    if (iterations < 0 || iterations > 20) {
        error_text = "Invalid Alivibe MSA refinement count.";
        return -1;
    }
    std::vector<std::string> sequences;
    if (!decode_sequences(input, length, sequences)) return -1;
    active_ambiguity_wildcard = ambiguity_wildcard;
    const std::vector<std::string> result = refined_msa(sequences, iterations, ambiguity_wildcard != '\0');
    active_ambiguity_wildcard = '\0';
    if (result.size() != sequences.size()) {
        error_text = "Alivibe MSA returned the wrong number of rows.";
        return -1;
    }
    if (!encode_sequences(result)) return -1;
    return static_cast<std::int32_t>(result.size());
}

// Literal scoring retained for ABI compatibility and exact oracle tests.
ALIVIBE_EXPORT("alivibe_msa_run") std::int32_t alivibe_msa_run(
    const std::uint8_t* input,
    std::uint32_t length,
    std::int32_t iterations
) {
    return run_msa(input, length, iterations, '\0');
}

// Nucleotide-mode Swig route: N receives the ordinary match score against any
// non-gap residue. Gap costs are unchanged because gaps are DP transitions,
// not substitution characters.
ALIVIBE_EXPORT("alivibe_msa_run_nucleotide") std::int32_t alivibe_msa_run_nucleotide(
    const std::uint8_t* input,
    std::uint32_t length,
    std::int32_t iterations
) {
    return run_msa(input, length, iterations, 'N');
}

// Amino-acid Swig route: X is an unknown-residue wildcard. Asparagine N
// remains an ordinary literal amino acid, and gaps retain the same affine and
// POA transition costs as the pinned Alivibe implementation.
ALIVIBE_EXPORT("alivibe_msa_run_amino_acid") std::int32_t alivibe_msa_run_amino_acid(
    const std::uint8_t* input,
    std::uint32_t length,
    std::int32_t iterations
) {
    return run_msa(input, length, iterations, 'X');
}

ALIVIBE_EXPORT("alivibe_msa_result_ptr") const std::uint8_t* alivibe_msa_result_ptr() {
    return result_bytes.empty() ? nullptr : result_bytes.data();
}

ALIVIBE_EXPORT("alivibe_msa_result_len") std::uint32_t alivibe_msa_result_len() {
    return static_cast<std::uint32_t>(result_bytes.size());
}

ALIVIBE_EXPORT("alivibe_msa_error_ptr") const char* alivibe_msa_error_ptr() {
    return error_text.empty() ? nullptr : error_text.data();
}

ALIVIBE_EXPORT("alivibe_msa_error_len") std::uint32_t alivibe_msa_error_len() {
    return static_cast<std::uint32_t>(error_text.size());
}

}  // extern "C"
