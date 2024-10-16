#pragma once
// Utilities for handling insertions (gaps) in sequences
#include <cstddef>
#include <iterator>

namespace utils
{
    // Struct for representing an insertion with nucleotide and gap counts
    struct Insertion2
    {
        size_t index; // The position where insertion starts
        size_t n_num; // Number of nucleotides inserted at the given position
        size_t gap_num; // Number of gaps inserted at the given position
    };

    // Struct for representing an insertion (gap)
    struct Insertion
    {
        size_t index; // The position where insertion starts
        size_t number; // Number of gaps inserted at the given position

        // Comparison operator to check equality between insertions
        bool operator==(const Insertion &rhs) const noexcept;

        // Static function to merge two sorted sequences of insertions (with addition of gap counts)
        template<typename InIt1, typename InIt2, typename OutIt>
        static void plus(InIt1 lhs_first, InIt1 lhs_last, InIt2 rhs_first, InIt2 rhs_last, OutIt dest)
        {
            // Iterate through both sequences and merge them
            for (; lhs_first != lhs_last && rhs_first != rhs_last;)
                if (lhs_first->index < rhs_first->index)
                {
                    *dest++ = *lhs_first++;
                }
                else if (lhs_first->index > rhs_first->index)
                {
                    *dest++ = *rhs_first++;
                }
                else // If indices are equal, add the gap numbers
                {
                    *dest++ = utils::Insertion({ lhs_first->index, lhs_first->number + rhs_first->number });
                    ++lhs_first;
                    ++rhs_first;
                }

            // Copy remaining elements from either sequence
            std::copy(lhs_first, lhs_last, dest);
            std::copy(rhs_first, rhs_last, dest);
        }

        // Static function to subtract gaps from one sequence from another (lhs >= rhs)
        template<typename InIt1, typename InIt2, typename OutIt>
        static void minus(InIt1 lhs_first, InIt1 lhs_last, InIt2 rhs_first, InIt2 rhs_last, OutIt dest)
        {
            // Iterate through both sequences and subtract gaps in rhs from lhs
            for (; rhs_first != rhs_last; ++lhs_first, ++rhs_first)
            {
                while (lhs_first->index != rhs_first->index)
                    *dest++ = *lhs_first++;

                const size_t difference = lhs_first->number - rhs_first->number;
                if (difference) *dest++ = utils::Insertion({ lhs_first->index, difference });
            }

            // Copy remaining elements from lhs
            std::copy(lhs_first, lhs_last, dest);
        }

        // Static function to insert gaps into a sequence
        template<typename InIt1, typename InIt2, typename OutIt>
        static void insert_gaps(InIt1 sequence_first, InIt1 sequence_last,
                                InIt2 insertion_first, InIt2 insertion_last, OutIt dest,
                                typename std::iterator_traits<InIt1>::value_type gap_symbol)
        {
            // Insert gaps at specified positions into the sequence
            for (unsigned last_index = 0; insertion_first != insertion_last; ++insertion_first)
            {
                auto sequence_stop = sequence_first;
                std::advance(sequence_stop, insertion_first->index - last_index);

                // Copy the sequence up to the point where gaps need to be inserted
                dest = std::copy(sequence_first, sequence_stop, dest);
                // Fill with the specified number of gap symbols
                dest = std::fill_n(dest, insertion_first->number, gap_symbol);

                sequence_first = sequence_stop;
                last_index = insertion_first->index;
            }

            // Copy the remainder of the sequence
            std::copy(sequence_first, sequence_last, dest);
        }

    };

}
