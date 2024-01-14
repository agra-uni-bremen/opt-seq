// *********************************************************************************************************************
// Includes
// *********************************************************************************************************************

#include "instruction-sequence.hpp"  // path_node

#include <nlohmann/json.hpp>  // nlohmann::json

#include <algorithm>  // std::find
#include <cassert>    // assert
#include <cctype>     // std::isdigit
#include <cstdint>    // std::int32_t
#include <cstdlib>    // EXIT_SUCCESS
#include <fstream>    // std::ifstream
#include <iostream>   // std::cerr
#include <string>     // std::string::npos
#include <utility>    // std::move
#include <vector>     // std::vector

// *********************************************************************************************************************
// Types
// *********************************************************************************************************************

enum class mode
{
    STD,
    SUBSEQ,
    VARIANT,
    FULL
};

// *********************************************************************************************************************
// Functions
// *********************************************************************************************************************

// =====================================================================================================================
// Parsing
// =====================================================================================================================

auto static parse(std::ifstream& json, mode const m)
{
    assert(json.is_open());

    std::vector<std::vector<path_node>> seqs;

    nlohmann::json data;
    json >> data;

    seqs.reserve(data.size());

    for (auto const& data_seq : data.items())
    {
        assert(data_seq.value().is_array());

        switch (m)
        {
            case mode::STD:
                if (data_seq.key().find('v') != std::string::npos || data_seq.key().find('-') != std::string::npos)
                {
                    continue;
                }
                break;
            case mode::SUBSEQ:
                if (data_seq.key().find('v') != std::string::npos)
                {  // variant was found
                    continue;
                }
                break;
            case mode::VARIANT:
                if (data_seq.key().find('-') != std::string::npos)
                {  // subsequence was found
                    continue;
                }
                break;
            case mode::FULL:  // fall-through
            default: break;
        }

        std::vector<path_node> seq;
        seq.reserve(data_seq.value().size());

        for (auto const& data_node : data_seq.value())
        {
            path_node v;

            v.instr = instruction::UNDEF;
            for (auto i = 1; i < static_cast<std::int32_t>(instruction::COUNT); ++i)
            {
                if (data_node["instruction"] == instr2str(static_cast<instruction>(i)))
                {
                    v.instr = static_cast<instruction>(i);
                    break;
                }
            }

            v.w = data_node["weight"];

            auto fill = [&](auto& dependencies, auto const& data_dependencies) {
                for (auto const offset : data_dependencies)
                {
                    // dependencies are treated in the same manner
                    if (std::find(dependencies.begin(), dependencies.end(), offset) == dependencies.end())
                    {
                        dependencies.push_back(offset);
                    }
                }
            };
            v.dependencies.resize(1);
            fill(v.dependencies[0], data_node["dependencies_true"]);
            fill(v.dependencies[0], data_node["dependencies_anti"]);
            fill(v.dependencies[0], data_node["dependencies_output"]);

            seq.push_back(std::move(v));
        }

        seqs.push_back(std::move(seq));
    }

    return seqs;
}

// =====================================================================================================================
// Startup
// =====================================================================================================================

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "optseq> Pass a file path and mode" << std::endl;
        return EXIT_FAILURE;
    }

    if (!std::isdigit(*argv[2]) || std::atoi(argv[2]) > 3)
    {
        std::cerr << "optseq> The mode must be an integer between 0 and 3" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream file{argv[1]};
    if (!file.is_open())
    {
        std::cerr << "optseq> Failed to open " << argv[1] << std::endl;
        return EXIT_FAILURE;
    }

    // read sequences
    std::vector<std::vector<path_node>> parsed_seqs;
    try
    {
        parsed_seqs = parse(file, static_cast<mode>(std::atoi(argv[2])));
    }
    catch (nlohmann::json::exception const& e)
    {
        std::cerr << "optseq> JSON parsing error: " << e.what() << std::endl;
        file.close();
        return EXIT_FAILURE;
    }
    file.close();

    instruction_sequence instr_seq{parsed_seqs};

    instr_seq.optimize();  // merge as many sequences as possible

    std::cout << "optseq> " << instr_seq << std::endl;

    return EXIT_SUCCESS;
}
