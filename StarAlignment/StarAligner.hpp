#pragma once
#include "../SuffixArray/SuffixArray.hpp"  // Include Suffix Array utility
#include "../Utils/Utils.hpp"  // Include general utilities
#include "../multi-thread/multi.hpp"  // Include multi-threading utilities

#include <vector>
#include <array>
#include <string>
#include <algorithm>

namespace star_alignment // Namespace for star alignment
{

    // Class for performing star alignment
    class StarAligner
    {
    private:
        using triple = std::array<size_t, 3>; // Define a triple array with 3 elements
        using quadra = std::array<size_t, 4>; // Define a quadra array with 4 elements
        using sequence_type = std::vector<unsigned char>; // Sequence type definition

    public:
        // Static function to align sequences based on insertions and threshold
        static std::vector<sequence_type> align(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center);

        // Static function to obtain gaps in sequences
        static void get_gaps(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center);

        // Functions to get optimal alignment paths and trace back the alignment
        static std::vector<triple> _optimal_path(const std::vector<triple>& common_substrings);
        static std::vector<int> _trace_back_bp(const std::vector<triple>& common_substrings, int* p);
        static std::vector<triple> _optimal_path_bp(const std::vector<triple>& optimal_common_substrings);

    private:
        // Constructor for StarAligner class
        StarAligner(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center);

        // Main alignment function
        std::vector<sequence_type> _align() const;

        // Function to get gaps for alignment
        void _get_gaps() const;

        // Set lengths of sequences
        std::vector<size_t> _set_lengths() const;

        // Set the center sequence
        size_t _set_centre() const;

        // Main steps of the star alignment
        std::vector<std::array<std::vector<utils::Insertion>, 2>> _pairwise_align() const; // Perform pairwise alignment
        std::vector<std::vector<utils::Insertion>> _merge_results(const std::vector<std::array<std::vector<utils::Insertion>, 2>>& pairwise_gaps) const; // Merge pairwise alignment results
        std::vector<sequence_type> _insert_gaps(const std::vector<std::vector<utils::Insertion>>& gaps) const; // Insert gaps into sequences

        // Multi-threaded pairwise alignment
        void mul_pairwise_align() const;

        // Helper function for multi-threaded alignment
        void mul_fasta_func(int i, const suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER>& st,
            std::vector<std::array<std::vector<utils::Insertion>, 2>>& all_pairwise_gaps, int threshold1) const;

        // Support function for appending gaps to the sequences
        static void _append(const std::vector<size_t>& src_gaps, std::vector<utils::Insertion>& des_gaps, size_t start);

        // Data members
        std::vector<std::vector<utils::Insertion>>& Insertions; // Reference to vector of insertions
        std::vector<sequence_type>& _sequences; // Reference to vector of sequences
        const size_t _row; // Number of sequences
        std::vector<size_t> _lengths; // Lengths of sequences
        size_t thresh1; // Threshold value for alignment
        size_t _centre; // Index of the center sequence
        size_t _centre_len; // Length of the center sequence
    };

}
