#if defined(_WIN32)
    #pragma once
#endif

#include "PairwiseAlignment/NeedlemanWunshReusable.hpp"
#include "StarAlignment/StarAligner.hpp"
#include "Utils/Fasta.hpp"
#include "Utils/Arguments.hpp"
#include "Utils/CommandLine.hpp"

#include <tuple>
#include <fstream>
#include <filesystem>
#include <string>

int main(int argc, char* argv[])
{
    SmpCommandLine userCommands(argc, argv);
    threadPool0 = NULL; // Initialize thread pool to NULL
    int center = -1; // Default reference index -1
    std::string center_name = ""; // Name of the reference sequence
    const auto start_point = std::chrono::high_resolution_clock::now(); // Record start time
    int thresh1 = 15;  
    int numThreads = 1;
    
    // Extract flagged arguments from the command line
    center_name = userCommands.getString("r", "reference", "[Longest]", "The reference sequence name [Please delete all whitespace]");
    numThreads = userCommands.getInteger("t", "threads", 1, "The number of threads");
    thresh1 = userCommands.getInteger("sa", "sa", 15, "The global sa threshold");
    arguments::in_file_name = userCommands.getString(1, "", " Input file/folder path[Please use .fasta as the file suffix or a forder]");
    arguments::out_file_name = userCommands.getString(2, "", " Output file path[Please use .fasta as the file suffix]");

    // Show help message if required
    if (userCommands.helpMessageWanted() || argc == 1)
    {
        userCommands.showHelpMessage();
        std::cout << "http://lab.malab.cn/soft/halign/\n";
        exit(0);
    }

    // Check if the minimum number of arguments is provided
    if (argc < 3)
    {
        userCommands.showHelpMessage();
        std::cout << "\n";
        std::cout << "http://lab.malab.cn/soft/halign/\n";
        std::cout << "Insufficient input parameter!\n";
        exit(1);
    }

    // Resolve absolute path of the input file/folder
    std::filesystem::path absolutePath = std::filesystem::absolute(arguments::in_file_name);
    if (std::filesystem::exists(absolutePath)) {
        arguments::in_file_name = absolutePath.generic_string();
        std::replace(arguments::in_file_name.begin(), arguments::in_file_name.end(), '\\', '/');
        if (std::filesystem::is_directory(absolutePath)) {
            if (arguments::in_file_name.back() != '/')
                arguments::in_file_name += '/';
        }
        else if (std::filesystem::is_regular_file(absolutePath) && absolutePath.extension().string().substr(1, 2) == "fa") {
            // Input is a fasta file
        }
        else {
            std::cout << "The input file/folder path does not represent a .fasta file or a directory." << std::endl;
            exit(1);
        }
    }
    else {
        std::cout << "The input file/folder path does not exist." << std::endl;
        exit(1);
    }

    // Resolve absolute path of the output file
    absolutePath = std::filesystem::absolute(arguments::out_file_name);
    arguments::out_file_name = absolutePath.generic_string();
    std::replace(arguments::out_file_name.begin(), arguments::out_file_name.end(), '\\', '/');

    // Print input/output information
    std::cout << "[  Input_name  ] : " << arguments::in_file_name << std::endl;
    std::cout << "[  Output_name ] : " << arguments::out_file_name << std::endl;
    std::cout << "[   Reference  ] = " << center_name << std::endl;
    std::cout << "[    Threads   ] = " << numThreads << std::endl;
    std::cout << "[      SA      ] = " << thresh1 << std::endl;
    
    threadPool0 = new ThreadPool(numThreads); // Create a thread pool

    std::vector<std::vector<unsigned char>> pseudo_sequences;
    cout_cur_time();
    std::cout << "Start: Read and data preprocessing: ";
    int II = 0;

    const auto read_T = std::chrono::high_resolution_clock::now(); // Record read start time
    if (arguments::in_file_name[arguments::in_file_name.size() - 1] == '/')
    {
        std::vector<std::string> files;
        // Get files in input directory
#if defined(_WIN32)
        getFiles_win(arguments::in_file_name, files);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
        getFiles_linux(arguments::in_file_name, files);
#endif  
        if (files.size() == 0)
        {
            std::cout << "\nThe input folder is empty!" << std::endl;
            exit(-1);
        }
        std::sort(files.begin(), files.end());
        std::cout << files.size() << " files\n";
        for (int i = 0; i < files.size(); i++)
        {
            std::ifstream ifs(files[i].c_str(), std::ios::binary | std::ios::in);
            for (auto x : utils::read_to_pseudo(ifs, center_name, II, center))
                pseudo_sequences.emplace_back(x);
            ifs.clear();
        }
    }
    else if (arguments::in_file_name.substr(arguments::in_file_name.find_last_of('.') + 1, 2) == "fa")
    {
        // Read from single fasta file
        std::ifstream ifs(arguments::in_file_name, std::ios::binary | std::ios::in);
        if (!ifs)
        {
            std::cout << "cannot access file " << arguments::in_file_name << '\n';
            exit(0);
        }
        std::cout << "1 files\n";
        pseudo_sequences = utils::read_to_pseudo(ifs, center_name, II, center);
        ifs.clear();
    }
    std::cout << "                    | Info : read consumes : " << (std::chrono::high_resolution_clock::now() - read_T) << "\n";
    cout_cur_time();
    std::cout << "End  : " << pseudo_sequences.size() << " sequences were discovered\n";
    if (pseudo_sequences.size() < 2)
    {
        std::cout << "The number of input sequences is less than two!\n";
        exit(1); // Exit if fewer than two sequences
    }

    // Start alignment process
    const auto align_start = std::chrono::high_resolution_clock::now(); // Record alignment start time
    std::vector<std::vector<utils::Insertion>> insertions(pseudo_sequences.size());
    star_alignment::StarAligner::get_gaps(insertions, pseudo_sequences, thresh1, center); // Perform MSA
    std::cout << "                    | Info : align time consumes : " << (std::chrono::high_resolution_clock::now() - align_start) << "\n";
    std::cout << "                    | Info : align memory peak   : " << getPeakRSS() << " B\n"; // Output memory usage
    
    const auto INSERT_T = std::chrono::high_resolution_clock::now(); // Record insertion start time
    size_t JJ = 0;
    if (arguments::out_file_name.substr(arguments::out_file_name.find_last_of('.') + 1, 2) == "fa")
    {
        std::ofstream ofs(arguments::out_file_name, std::ios::binary | std::ios::out); // Open output file
        if (!ofs)
        {
            std::cout << "cannot write file " << arguments::out_file_name << '\n';
            exit(0);
        }

        // Write to output fasta file
        if (arguments::in_file_name[arguments::in_file_name.size() - 1] == '/')
        {
            std::vector<std::string> files;
#if defined(_WIN32)
            getFiles_win(arguments::in_file_name, files);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
            getFiles_linux(arguments::in_file_name, files);
#endif  
            std::sort(files.begin(), files.end());
            for (int i = 0; i < files.size(); i++)
            {
                std::ifstream ifs(files[i].c_str(), std::ios::binary | std::ios::in);
                utils::write_to_fasta(ofs, ifs, insertions, JJ);
                ifs.clear();
            }
        }
        else if (arguments::in_file_name.substr(arguments::in_file_name.find_last_of('.') + 1, 2) == "fa")
        {
            std::ifstream ifs(arguments::in_file_name, std::ios::binary | std::ios::in);
            utils::write_to_fasta(ofs, ifs, insertions, JJ);
            ifs.clear();
        }
        ofs.close();
    }
    std::cout << "                    | Info : write consumes: " << (std::chrono::high_resolution_clock::now() - INSERT_T) << "\n";

    // Print process statistics
    std::cout << "                    | Info : Current pid   : " << getpid() << std::endl;
    std::cout << "                    | Info : Time consumes : " << (std::chrono::high_resolution_clock::now() - start_point) << "\n";
    std::cout << "                    | Info : Memory usage  : " << getPeakRSS() << " B" << std::endl;
    std::cout << "                    | Info : Finished\n";
    std::cout << "http://lab.malab.cn/soft/halign/\n";
    return 0;
}
