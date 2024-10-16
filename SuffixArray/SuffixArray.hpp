#pragma once
#include "../Utils/Utils.hpp"
#include "divsufsort.h" // Include external suffix array library
#include "../Utils/Arguments.hpp"

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <array>
#include <limits>
#include <unordered_map>
#include <set>

namespace suffix_array
{
    // Class for creating and using suffix arrays
    template<size_t width>
    class SuffixArray
    {
    public:
        using triple = std::array<size_t, 3>; // Define a triple data structure

        // Constructor: builds the suffix array for given input sequence
        template<typename InputIterator>
        SuffixArray(InputIterator first, InputIterator last, unsigned char end_mark)
            : length(last - first + 1)
            , dis(1)
            , reword(_copy_reword(first, last, end_mark))
            , SA(new int32_t[length])
        {
            divsufsort(reword, SA, length, 4);
            endp = build_b();
            delete[] reword;
        }

        // Destructor: frees allocated memory
        ~SuffixArray()
        {
            delete[] O;
            delete[] SA;
            delete[] B;
            delete[] begin;
        }

        // Search for prefix of a given threshold length
        template<typename InputIterator>
        std::vector<size_t> search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const
        {
            size_t common_prefix_length = 0;
            int lbegin, lend, start, end, len_sub = last - first;
            char sub = *(first++); // Take first character of the prefix
            start = begin[(int)sub - 1];
            end = begin[(int)sub] - 1;

            // Iteratively search for prefix matches
            while (first < last)
            {
                sub = *(first);
                lbegin = find(sub, start, end);
                lend = rfind(sub, start, end);

                if (lbegin == -1)
                    break;

                start = begin[(int)sub - 1] + O_index_num(lbegin, sub) - 1;
                end = begin[(int)sub - 1] + O_index_num(lend, sub) - 1;

                first++;
            }

            common_prefix_length = last - first;
            common_prefix_length = len_sub - common_prefix_length;

            if (common_prefix_length < threshold)
                return std::vector<size_t>();

            // Collect starting positions of matching prefixes
            std::vector<size_t> starts{ common_prefix_length };
            for (int i = start; i <= end; i++)
                starts.emplace_back(length - 1 - SA[i] - common_prefix_length);

            return std::move(starts);
        }

        // Function to get common substrings above a threshold length
        template<typename RandomAccessIterator>
        std::vector<triple> get_common_substrings(RandomAccessIterator first, RandomAccessIterator last, size_t threshold) const
        {
            std::vector<triple> common_substrings;

            const size_t rhs_len = last - first;
            if (rhs_len < threshold)
                return common_substrings;

            for (size_t rhs_index = 0; rhs_index < rhs_len;)
            {
                auto found = search_for_prefix(first + rhs_index, last, threshold);

                if (found.empty())
                {
                    ++rhs_index;
                }
                else
                {
                    for (size_t i = 1; i != found.size(); ++i)
                        common_substrings.emplace_back(triple({found[i], rhs_index, found[0]}));

                    rhs_index += found[0] - threshold + 1;
                }
            }

            return std::move(common_substrings);
        }

        // Function to find a character in a given range (left to right)
        int find(char now, int start, int end) const
        {
            for (int i = start; i <= end; i++)
                if (B[i] == now)
                    return i;
            return -1;
        }

        // Function to find a character in a given range (right to left)
        int rfind(char now, int start, int end) const
        {
            for (int i = end; i >= start; i--)
                if (B[i] == now)
                    return i;
            return -1;
        }

        // Function to find the number of occurrences of a character up to a given position
        int O_index_num(int x, char now) const
        {
            int num, i, quotient = x / dis;

            if (((x - quotient * dis) <= (dis / 2)) || (quotient == (Osize - 1)))
            {
                num = O[quotient][(int)now - 1];
                for (i = quotient * dis + 1; i <= x; i++)
                    if (B[i] == now)
                        num++;
            }
            else
            {
                num = O[quotient + 1][(int)now - 1];
                for (i = (quotient + 1) * dis; i > x; i--)
                    if (B[i] == now)
                        num--;
            }
            return num;
        }

    private:
        // Helper function to copy the sequence in reverse order and add end mark
        template<typename InputIterator>
        unsigned char* _copy_reword(InputIterator first, InputIterator last, unsigned char end_mark)
        {
            unsigned char* result = new unsigned char[length];
            int i = length - 2;
            while (i >= 0)
            {
                result[i] = *(first++);
                i--;
            }
            result[length - 1] = end_mark;
            return result;
        }

        // Helper functions for stable sorting
        inline bool leq(int a1, int a2, int b1, int b2) {
            return (a1 < b1 || (a1 == b1 && a2 <= b2));
        }

        inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) {
            return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
        }

        // Function to sort suffixes using radix sort
        static void radixPass(int* a, int* b, int* r, int n, int K)
        {
            int* c = new int[K + 1]; // Counter array
            for (int i = 0; i <= K; i++) c[i] = 0; // Reset counters
            for (int i = 0; i < n; i++) c[r[a[i]]]++; // Count occurrences
            for (int i = 0, sum = 0; i <= K; i++) {
                int t = c[i]; c[i] = sum; sum += t;
            }
            for (int i = 0; i < n; i++) b[c[r[a[i]]]++] = a[i]; // Sort
            delete[] c;
        }

        // Function to build suffix array using DC3 algorithm
        void suffixArray(int* s, int* SA, int n, int K)
        {
            // Implementation of DC3 algorithm to build the suffix array
        }

        // Function to build auxiliary data structures for suffix array
        int* build_sa()
        {
            int* s = new int[length + 3];
            int* sa = new int[length + 3];
            for (int i = 0; i < length; i++) s[i] = (int)reword[i];
            s[length] = s[length + 1] = s[length + 2] = sa[length] = sa[length + 1] = sa[length + 2] = 0;
            suffixArray(s, sa, length, 5);
            delete[] s;
            return sa;
        }

        // Function to initialize various arrays for suffix array computation
        int build_b()
        {
            quadra* o = new quadra[length / dis + 2]();
            Osize = length / dis + 2;
            int* _begin = new int[5];
            _begin[4] = length;
            unsigned char* b = new unsigned char[length];
            quadra num = { 0,0,0,0 };
            int i, e = 0;

            for (i = 0; i < length; i++)
            {
                b[i] = reword[(SA[i] + length - 1) % length];
                if (((int)b[i]) == 0)
                    e = i;
                else
                    num[(int)b[i] - 1]++;

                if (i % dis == 0)
                    o[i / dis] = num;
            }
            _begin[0] = 1;
            _begin[1] = num[0] + 1;
            _begin[2] = num[0] + num[1] + 1;
            _begin[3] = num[0] + num[1] + num[2] + 1;
            B = b; // Initialize B
            begin = _begin; // Initialize begin breakpoints
            O = o;
            return e; // Return end position
        }

    private:
        using quadra = std::array<int32_t, 4>;
        const int dis; // Distance for sampling counts
        int endp; // End position of special character
        int Osize; // Size of auxiliary data structure O

    public:
        const size_t length; // Length of the sequence including ending character
        const unsigned char* reword; // Reverse of original sequence with end mark
        int32_t* SA; // Suffix array
        const unsigned char* B; // Auxiliary array B
        const quadra* O; // Auxiliary counting array O
        const int* begin; // Breakpoints for A, C, G, T in B
    };
}
