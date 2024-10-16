#include "NeedlemanWunshReusable.hpp"

// Constructor for Kband class
Kband::Kband() {
    match = 1; // Match score
    mismatch = -2; // Mismatch penalty
    d = 3; // Gap open penalty
    e = 1; // Gap extension penalty
    my_INT_MIN = 0; // Initialize with a minimum value

    // Allocate memory for sequences and score matrices
    A = new unsigned char[thresh0 + 2];
    B = new unsigned char[thresh0 + 2];
    pm = new int* [3];
    pm2 = new int* [3];
    for (int i = 0; i < 3; ++i) {
        pm[i] = new int[thresh0];
        pm2[i] = new int[thresh0];
    }
    pmt1 = pm;
    pmt2 = pm2;
    pmt = pm;

    // Allocate memory for backtrace matrix
    bt = new unsigned char* [thresh0];
    for (int i = 0; i < thresh0; ++i) {
        bt[i] = new unsigned char[thresh0];
    }

    // Resize vectors to store sequences
    seq_A.resize(thresh0);
    seq_B.resize(thresh0);
}

// Destructor for Kband class
Kband::~Kband() {
    // Free allocated memory
    delete[] A;
    delete[] B;
    for (int i = 0; i < 3; ++i) {
        delete[] pmt1[i];
        delete[] pmt2[i];
    }
    delete[] pmt1;
    delete[] pmt2;
    for (int i = 0; i < thresh0; ++i) 
        delete[] bt[i];
    delete[] bt;

    // Clear sequences
    char_vector_type().swap(seq_A);
    char_vector_type().swap(seq_B);
}

// Function to calculate score for matching two characters
inline int Kband::score(unsigned char xi, unsigned char yi) {
    return (xi == yi) ? match : mismatch;
}

// Check if the given point is inside the k-band strip
inline bool Kband::InsiderStrip(int i, int j, int k, int diff) {
    return ((-k <= (j - i)) && ((j - i) <= (k + diff)));
}

// Helper function for indexing
inline int Kband::index(int i, int l) {
    return (i + l) % l;
}

// Function to get the maximum of two integers
inline int Kband::maxi(int a, int b) {
    return (a > b) ? a : b;
}

// Function to get the maximum of three integers
inline int Kband::maxi(int a, int b, int c) {
    int max = (a > b) ? a : b;
    return (max > c) ? max : c;
}

// Function to initialize matrices for alignment
void Kband::Init(int m, int k, int diff) {
    for (int i = 0; i < (m + 1); i++) {
        for (int j = 0; j < (diff + 2 * k + 1); j++)
            bt[i][j] = '\0';
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < (diff + 2 * k + 1); j++)
            pm[i][j] = my_INT_MIN;
    }
    pm[0][k] = 0;
    bt[0][k] = '\16';
    for (int j = 0; j < (diff + k + 1); j++) {
        pm[1][j + k] = -d - e * (j - 1);
        bt[0][j + k] = (char)8;
    }
    for (int i = 0; i < (k + 1); i++)
        bt[i][k - i] = '\3';
}

// Function to initialize secondary matrices
void Kband::InitTwo(int ii, int k, int diff) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < (diff + 2 * k + 1); j++)
            pm2[i][j] = my_INT_MIN;
    }
    if (ii < k + 1)
        pm2[2][index(k - ii, diff + 2 * k + 1)] = -d - e * (ii - 1);
}

// Function to choose the alignment path based on scores
int Kband::ChooseWay(int p0, int p1, int p2, bool state) {
    if (p0 >= p1) {
        return (p0 >= p2) ? (state ? 16 : 0) : (state ? 48 : 2);
    } else {
        return (p1 >= p2) ? (state ? 32 : 1) : (state ? 48 : 2);
    }
}

// Function to parse a value for backtrace
inline int Kband::parse(int b, int s) {
    b = (b >> (4 - s * 2));
    return (b & 3) - 1;
}

// Function to perform sequence alignment using K-band method
std::tuple<insert, insert> Kband::PSA_AGP_Kband3(const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end, const std::vector<unsigned char>& sequence2, int b_begin, int b_end,
                                                 int cmatch, int cmismatch, int cd, int ce) {
    // Initialize scoring parameters
    match = cmatch;
    mismatch = cmismatch;
    d = cd;
    e = ce;

    // Other alignment logic follows...
}

// Function to parse the CIGAR string representation of an alignment
std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin, bool tag) {
    // Parses the CIGAR string and returns the gaps in two sequences
}

// Function to insert gaps into a string based on given insertions
void insertGaps(std::string& str, const insert& insertions) {
    for (auto it = insertions.rbegin(); it != insertions.rend(); ++it) {
        int index = std::get<0>(*it);
        int num = std::get<1>(*it);
        str.insert(index, num, '-');
    }
}

// Function to perform sequence alignment using WFA (Wavefront Alignment)
std::tuple<insert, insert> mywfa(wfa::WFAlignerGapAffine& aligner, const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end, const std::vector<unsigned char>& sequence2, int b_begin, int b_end) {
    // Performs alignment using WFA and returns the parsed CIGAR
}
