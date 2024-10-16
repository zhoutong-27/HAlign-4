#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <tuple>
#include <cmath>
#include "../SuffixArray/SuffixArray.hpp"
#include "../Utils/Utils.hpp"
#include <iostream>
#include <string>
#include "WFA2-lib/bindings/cpp/WFAligner.hpp"

// Constant definition
const int thresh0 = 10001;

// Type aliases for convenience
using RandomAccessIterator = std::vector<unsigned char>::const_iterator;
using gap_vector_type = std::vector<int>;
using char_vector_type = std::vector<unsigned char>;
using triple = std::array<int, 4>;
using quadra = std::array<int, 5>;
using insert = std::vector<std::tuple<int, int>>;
using in = std::tuple<int, int>;

// Class for K-band alignment implementation
class Kband {
public:
    // Alignment scoring parameters
    int match;
    int mismatch;
    int d;    // Gap open penalty
    int e;    // Gap extension penalty
    int my_INT_MIN;

    // Pointers for sequences
    unsigned char* A;
    unsigned char* B;
    unsigned char* Aq;
    unsigned char* Bq;

    // Score matrices
    int** pm;
    int** pm2;
    int** pmt1;
    int** pmt2;
    int** pmt;

    // Backtrace pointers
    unsigned char** bt;

    // Vectors for storing sequences
    std::vector<unsigned char> seq_A;
    std::vector<unsigned char> seq_B;

public:
    // Constructor
    Kband();

    // Destructor
    ~Kband();

    // Scoring function for comparing two characters
    inline int score(unsigned char xi, unsigned char yi);

    // Function to check if a given point is inside the strip defined by parameters
    inline bool InsiderStrip(int i, int j, int k, int diff = 0);

    // Helper function for indexing
    inline int index(int i, int l);

    // Function to get the maximum of two integers
    inline int maxi(int a, int b);

    // Function to get the maximum of three integers
    inline int maxi(int a, int b, int c);

    // Function to initialize alignment parameters
    void Init(int m, int k, int diff);

    // Function to initialize two sets of alignment parameters
    void InitTwo(int ii, int k, int diff);

    // Function to choose the alignment path based on scores
    int ChooseWay(int p0, int p1, int p2, bool state = true);

    // Function to parse a given value
    inline int parse(int b, int s);

    // Function to perform sequence alignment using the K-band algorithm
    std::tuple<insert, insert> PSA_AGP_Kband3(const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end,
                                              const std::vector<unsigned char>& sequence2, int b_begin, int b_end,
                                              int cmatch = 1, int cmismatch = -2, int cd = 3, int ce = 1);
};

// Function to perform sequence alignment using WFA (Wavefront Alignment)
std::tuple<insert, insert> mywfa(wfa::WFAlignerGapAffine& aligner, const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end,
                                 const std::vector<unsigned char>& sequence2, int b_begin, int b_end);

// Function to parse CIGAR string representation of alignment
std::tuple<insert, insert> parseCigar(const std::string& cigar);

// Overloaded function to parse CIGAR with specified beginning positions
std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin);

// Overloaded function to parse CIGAR with a tag option
std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin, bool tag);

// Function to insert gaps into a string based on given insertions
void insertGaps(std::string& str, const insert& insertions);
