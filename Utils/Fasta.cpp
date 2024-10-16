#include "Fasta.hpp"
// Functions for reading and writing .fasta files
#include <cstring>

// Constructor: reads fasta sequences from input stream
utils::Fasta::Fasta(std::istream &is)
{
    _read(is);
}

// Function to write fasta sequences to an output stream
void utils::Fasta::write_to(std::ostream &os, bool with_identification) const
{
    if (with_identification)
        write_to(os, sequences.cbegin(), sequences.cend(), identifications.cbegin());
    else
        write_to(os, sequences.cbegin(), sequences.cend());
}

// Private function to read fasta file from input stream
void utils::Fasta::_read(std::istream &is)
{
    std::string each_line; // Stores each line of input
    std::string each_sequence; // Stores the current sequence
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0)
            continue; // Skip empty lines

        if (each_line[0] == '>') // Header line
        {
            identifications.emplace_back(each_line.substr(1)); // Extract sequence identifier (without '>')
            if (flag)
                sequences.emplace_back(std::move(each_sequence)); // Store the previous sequence if any
            flag = true;
            each_sequence.clear(); // Prepare for a new sequence
        }
        else if (flag) // Sequence lines
        {
            each_sequence += each_line; // Concatenate the sequence lines
        }
    }
    sequences.emplace_back(each_sequence); // Store the last sequence
}

// Function to cut sequence into multiple lines and write to output stream
void utils::Fasta::cut_and_write(std::ostream &os, const std::string &sequence)
{
    const size_t sequence_length = sequence.size();

    // Allocate buffer for the sequence with line breaks added
    char *cut_sequence = new char[sequence_length + sequence_length / max_line_length + 1];
    size_t des_index = 0;
    for (size_t src_index = 0; src_index < sequence_length; src_index += max_line_length)
    {
        if (src_index) cut_sequence[des_index++] = '\n'; // Add a newline after each segment

        // Determine length to write (either max line length or remaining length)
        size_t write_length = sequence_length - src_index;
        if (write_length > max_line_length) write_length = max_line_length;

        memcpy(cut_sequence + des_index, sequence.data() + src_index, write_length);
        des_index += write_length;
    }
    cut_sequence[des_index] = 0; // Null terminate the string

    os << cut_sequence; // Write the sequence to output stream
    delete[] cut_sequence; // Free the allocated memory
}
