#include "StarAligner.hpp"
#include "../Utils/Pseudo.hpp"
#include "../Utils/Utils.hpp"
#include "../Utils/Graph.hpp"
#include "../PairwiseAlignment/NeedlemanWunshReusable.hpp"

// Function to align sequences using star alignment
std::vector<std::vector<unsigned char>> star_alignment::StarAligner::align(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center) {
    return StarAligner(insertions, sequences, thresh, center)._align();
}

// Function to get gaps in sequences using star alignment
void star_alignment::StarAligner::get_gaps(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center) {
    StarAligner(insertions, sequences, thresh, center)._get_gaps();
}

// Constructor for StarAligner class
star_alignment::StarAligner::StarAligner(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center)
    : thresh1(thresh)
    , Insertions(insertions)
    , _sequences(sequences)
    , _row(_sequences.size())
    , _lengths(_set_lengths())
    , _centre(_set_centre())
    , _centre_len(_sequences[_centre].size()) {
    if (center != -1) {
        _centre = center;
        _centre_len = _sequences[_centre].size();
    }
}

// Function to set the lengths of sequences
std::vector<size_t> star_alignment::StarAligner::_set_lengths() const {
    std::vector<size_t> lengths(_row);
    for (size_t i = 0; i != _row; ++i) lengths[i] = _sequences[i].size();
    return lengths;
}

// Function to set the central sequence
size_t star_alignment::StarAligner::_set_centre() const {
    size_t centre_index = 0;
    for (size_t i = 1; i != _row; ++i)
        if (_lengths[i] > _lengths[centre_index])
            centre_index = i;
    return centre_index;
}

// Function to perform alignment
std::vector<std::vector<unsigned char>> star_alignment::StarAligner::_align() const {
    return _insert_gaps(_merge_results(_pairwise_align()));
}

// Function to get gaps for alignment
void star_alignment::StarAligner::_get_gaps() const {
    mul_pairwise_align();
}

// Function to perform pairwise alignment
auto star_alignment::StarAligner::_pairwise_align() const -> std::vector<std::array<std::vector<utils::Insertion>, 2>> {
    suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER> st(_sequences[_centre].cbegin(), _sequences[_centre].cend(), nucleic_acid_pseudo::end_mark); // Create suffix array
    std::vector<std::array<std::vector<utils::Insertion>, 2>> all_pairwise_gaps;
    wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh); // Create WFA aligner
    
    for (size_t i = 0; i != _row; ++i) {
        auto common_substrings = _optimal_path(st.get_common_substrings(_sequences[i].cbegin(), _sequences[i].cend(), thresh1));
        
        // Define alignment intervals
        std::vector<quadra> intervals;
        intervals.reserve(common_substrings.size() + 1);
        if (common_substrings.empty()) {
            intervals.emplace_back(quadra({0, _centre_len, 0, _sequences[i].size()}));
        } else {
            if (common_substrings[0][0] || common_substrings[0][1])
                intervals.emplace_back(quadra({0, common_substrings[0][0], 0, common_substrings[0][1]}));
            for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
                if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                    common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                    intervals.emplace_back(quadra({
                        common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                        common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                    }));
            if (common_substrings.back()[0] + common_substrings.back()[2] != _centre_len ||
                common_substrings.back()[1] + common_substrings.back()[2] != _lengths[i])
                intervals.emplace_back(quadra({
                    common_substrings.back()[0] + common_substrings.back()[2], _centre_len,
                    common_substrings.back()[1] + common_substrings.back()[2], _lengths[i]
                }));
        }

        // Perform pairwise alignment for each interval
        std::array<std::vector<utils::Insertion>, 2> pairwise_gaps;
        for (size_t j = 0; j != intervals.size(); ++j) {
            const size_t centre_begin = intervals[j][0];
            const size_t centre_end = intervals[j][1];
            const size_t sequence_begin = intervals[j][2];
            const size_t sequence_end = intervals[j][3];

            auto [lhs_gaps, rhs_gaps] = mywfa(aligner, _sequences[_centre], centre_begin, centre_end,
                _sequences[i], sequence_begin, sequence_end); // Perform alignment using WFA
            
            // Collect gaps for alignment
            for (int ii = 0; ii < lhs_gaps.size(); ii++) {
                if ((!pairwise_gaps[0].empty()) && pairwise_gaps[0].back().index == std::get<0>(lhs_gaps[ii]))
                    pairwise_gaps[0].back().number += std::get<1>(lhs_gaps[ii]);
                else
                    pairwise_gaps[0].emplace_back(utils::Insertion({(size_t)std::get<0>(lhs_gaps[ii]), (size_t)std::get<1>(lhs_gaps[ii])}));
            }
            for (int ii = 0; ii < rhs_gaps.size(); ii++) {
                if ((!pairwise_gaps[1].empty()) && pairwise_gaps[1].back().index == std::get<0>(rhs_gaps[ii]))
                    pairwise_gaps[1].back().number += std::get<1>(rhs_gaps[ii]);
                else
                    pairwise_gaps[1].emplace_back(utils::Insertion({(size_t)std::get<0>(rhs_gaps[ii]), (size_t)std::get<1>(rhs_gaps[ii])}));
            }
        }
        all_pairwise_gaps.emplace_back(pairwise_gaps);
        _sequences[i].clear();
    }
    return all_pairwise_gaps;
}

// Helper function to find optimal path in common substrings
auto star_alignment::StarAligner::_optimal_path(const std::vector<triple>& common_substrings) -> std::vector<triple> {
    std::vector<triple> optimal_common_substrings;
    if (common_substrings.empty()) return optimal_common_substrings;

    const size_t pair_num = common_substrings.size();
    utils::AdjacencyList graph(pair_num + 1);

    // Build graph of common substrings
    for (size_t i = 0; i != pair_num; ++i)
        for (size_t j = 0; j != pair_num; ++j)
            if (i != j && common_substrings[i][0] + common_substrings[i][2] < common_substrings[j][0] + common_substrings[j][2]
                && common_substrings[i][1] + common_substrings[i][2] < common_substrings[j][1] + common_substrings[j][2]) {
                const int possible_overlap = std::max(
                    static_cast<int>(common_substrings[i][0] + common_substrings[i][2]) - static_cast<int>(common_substrings[j][0]),
                    static_cast<int>(common_substrings[i][1] + common_substrings[i][2]) - static_cast<int>(common_substrings[j][1])
                );
                unsigned weight = common_substrings[j][2];
                if (possible_overlap > 0) weight -= possible_overlap;
                graph.add_edge(i + 1, j + 1, weight);
            }
    for (size_t i = 0; i != pair_num; ++i)
        graph.add_edge(0, i + 1, common_substrings[i][2]);

    // Find the longest path
    const auto optimal_path = graph.get_longest_path();
    optimal_common_substrings.reserve(optimal_path.size());
    optimal_common_substrings.emplace_back(triple({common_substrings[optimal_path[0] - 1][0],
                                                   common_substrings[optimal_path[0] - 1][1],
                                                   common_substrings[optimal_path[0] - 1][2]}));

    for (size_t i = 0; i < optimal_path.size() - 1; ++i) {
        size_t new_len = graph.get_weight(optimal_path[i], optimal_path[i + 1]);
        size_t old_len = common_substrings[optimal_path[i + 1] - 1][2];
        int difference = static_cast<int>(old_len) - static_cast<int>(new_len);
        size_t lhs_first = common_substrings[optimal_path[i + 1] - 1][0];
        size_t rhs_first = common_substrings[optimal_path[i + 1] - 1][1];
        if (difference > 0) {
            lhs_first += difference; rhs_first += difference;
        }
        optimal_common_substrings.emplace_back(triple({lhs_first, rhs_first, new_len}));
    }

    return optimal_common_substrings;
}

// Function to append gaps
void star_alignment::StarAligner::_append(const std::vector<size_t>& src_gaps, std::vector<utils::Insertion>& des_gaps, size_t start) {
    for (size_t i = 0; i != src_gaps.size(); ++i)
        if (src_gaps[i]) {
            if (!des_gaps.empty() && des_gaps.back().index == start + i)
                des_gaps.back().number += src_gaps[i];
            else
                des_gaps.emplace_back(utils::Insertion({start + i, src_gaps[i]}));
        }
}

// Function to merge pairwise alignment results
auto star_alignment::StarAligner::_merge_results(const std::vector<std::array<std::vector<utils::Insertion>, 2>>& pairwise_gaps) const -> std::vector<std::vector<utils::Insertion>> {
    // Implementation of merging results from pairwise alignments
}

// Function to insert gaps into sequences
auto star_alignment::StarAligner::_insert_gaps(const std::vector<std::vector<utils::Insertion>>& gaps) const -> std::vector<sequence_type> {
    // Insert gaps into sequences based on pairwise gaps
}

// Function to perform multi-threaded pairwise alignment
void star_alignment::StarAligner::mul_pairwise_align() const {
    // Perform pairwise alignment using multiple threads
}

// Helper function for multi-threaded alignment of sequences
void star_alignment::StarAligner::mul_fasta_func(int i, const suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER>& st,
    std::vector<std::array<std::vector<utils::Insertion>, 2>>& all_pairwise_gaps, int threshold1) const {
    // Function to align a sequence with the center sequence
}

// Helper function for backtracking to find optimal path
std::vector<int> star_alignment::StarAligner::_trace_back_bp(const std::vector<triple>& common_substrings, int* p) {
    // Backtrack to find the optimal alignment path
}

// Function to find optimal path with backtracking
auto star_alignment::StarAligner::_optimal_path_bp(const std::vector<triple>& optimal_common_substrings) -> std::vector<triple> {
    // Use backtracking to find the optimal alignment path
}
