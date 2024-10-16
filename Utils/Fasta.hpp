#pragma once
// Header file for reading and writing .fasta files
#include <string>
#include <vector>
#include <iostream>

namespace utils
{
    // Class for handling fasta file operations
    class Fasta
    {
    private:
        // Private function to read sequences from input stream
        void _read(std::istream &is);

    public:
        static constexpr unsigned max_line_length = 80; // Maximum line length for fasta file output

        // Data members for storing sequences and their identifiers
        std::vector<std::string> sequences;
        std::vector<std::string> identifications;

        // Constructor: reads fasta sequences from input stream
        explicit Fasta(std::istream &is);

        // Function to write sequences to an output stream
        void write_to(std::ostream &os, bool with_idification = true) const;

        // Static function to cut a sequence into lines of max length and write to output stream
        static void cut_and_write(std::ostream &os, const std::string &sequence);

        // Static function to write sequences to output stream without identifiers
        template<typename InputIterator>
        static void write_to(std::ostream &os, InputIterator sequence_first, InputIterator sequence_last)
        {
            if (sequence_first == sequence_last) return; // Return if no sequences to write

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len; ++sequence_first, ++i)
            {
                os << *sequence_first; // Write each sequence
                if (i != len - 1) os << '\n'; // Add newline after each sequence except the last one
            }
        }

        // Static function to write sequences along with identifiers to output stream
        template<typename InputIterator1, typename InputIterator2>
        static void write_to(std::ostream &os, InputIterator1 sequence_first, InputIterator1 sequence_last,
                             InputIterator2 identification_first)
        {
            if (sequence_first == sequence_last) return; // Return if no sequences to write

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len; ++sequence_first, ++identification_first, ++i)
            {
                os << '>' << *identification_first << '\n'; // Write identifier line
                cut_and_write(os, *sequence_first); // Write sequence with line breaks
                if (i != len - 1) os << '\n'; // Add newline after each sequence except the last one
            }
        }
    };
}
