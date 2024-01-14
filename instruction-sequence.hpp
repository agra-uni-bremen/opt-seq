// *********************************************************************************************************************
// Includes
// *********************************************************************************************************************

#include <algorithm>  // std::count
#include <cassert>    // assert
#include <cstdint>    // std::int32_t
#include <iostream>   // std::cout
#include <numeric>    // std::accumulate
#include <ostream>    // std::ostream
#include <utility>    // std::move
#include <vector>     // std::vector

// *********************************************************************************************************************
// Types
// *********************************************************************************************************************

enum class instruction
{
    UNDEF,
    // RV32I
    LUI,
    AUIPC,
    JAL,
    JALR,
    BEQ,
    BNE,
    BLT,
    BGE,
    BLTU,
    BGEU,
    LB,
    LH,
    LW,
    LBU,
    LHU,
    SB,
    SH,
    SW,
    ADDI,
    SLTI,
    SLTIU,
    XORI,
    ORI,
    ANDI,
    SLLI,
    SRLI,
    SRAI,
    ADD,
    SUB,
    SLL,
    SLT,
    SLTU,
    XOR,
    SRL,
    SRA,
    OR,
    AND,
    FENCE,
    ECALL,
    EBREAK,
    // Zifencei
    FENCE_I,
    // Zicsr
    CSRRW,
    CSRRS,
    CSRRC,
    CSRRWI,
    CSRRSI,
    CSRRCI,
    // RV32M
    MUL,
    MULH,
    MULHSU,
    MULHU,
    DIV,
    DIVU,
    REM,
    REMU,
    // RV32A
    LR_W,
    SC_W,
    AMOSWAP_W,
    AMOADD_W,
    AMOXOR_W,
    AMOAND_W,
    AMOOR_W,
    AMOMIN_W,
    AMOMAX_W,
    AMOMINU_W,
    AMOMAXU_W,
    // RV64I
    LWU,
    LD,
    SD,
    ADDIW,
    SLLIW,
    SRLIW,
    SRAIW,
    ADDW,
    SUBW,
    SLLW,
    SRLW,
    SRAW,
    // RV64M
    MULW,
    DIVW,
    DIVUW,
    REMW,
    REMUW,
    // RV64A
    LR_D,
    SC_D,
    AMOSWAP_D,
    AMOADD_D,
    AMOXOR_D,
    AMOAND_D,
    AMOOR_D,
    AMOMIN_D,
    AMOMAX_D,
    AMOMINU_D,
    AMOMAXU_D,
    // RV32F
    FLW,
    FSW,
    FMADD_S,
    FMSUB_S,
    FNMADD_S,
    FNMSUB_S,
    FADD_S,
    FSUB_S,
    FMUL_S,
    FDIV_S,
    FSQRT_S,
    FSGNJ_S,
    FSGNJN_S,
    FSGNJX_S,
    FMIN_S,
    FMAX_S,
    FCVT_W_S,
    FCVT_WU_S,
    FMV_X_W,
    FEQ_S,
    FLT_S,
    FLE_S,
    FCLASS_S,
    FCVT_S_W,
    FCVT_S_WU,
    FMV_W_X,
    // RV64F
    FCVT_L_S,
    FCVT_LU_S,
    FCVT_S_L,
    FCVT_S_LU,
    // RV32D
    FLD,
    FSD,
    FMADD_D,
    FMSUB_D,
    FNMSUB_D,
    FNMADD_D,
    FADD_D,
    FSUB_D,
    FMUL_D,
    FDIV_D,
    FSQRT_D,
    FSGNJ_D,
    FSGNJN_D,
    FSGNJX_D,
    FMIN_D,
    FMAX_D,
    FCVT_S_D,
    FCVT_D_S,
    FEQ_D,
    FLT_D,
    FLE_D,
    FCLASS_D,
    FCVT_W_D,
    FCVT_WU_D,
    FCVT_D_W,
    FCVT_D_WU,
    // RV64D
    FCVT_L_D,
    FCVT_LU_D,
    FMV_X_D,
    FCVT_D_L,
    FCVT_D_LU,
    FMV_D_X,
    // privileged
    URET,
    SRET,
    MRET,
    WFI,
    SFENCE_VMA,

    COUNT  // invalid argument
};

auto constexpr instr2str(instruction const instr) noexcept -> char const*
{
    switch (instr)
    {
        case instruction::UNDEF: return "ZERO-INVALID";
        case instruction::LUI: return "LUI";
        case instruction::AUIPC: return "AUIPC";
        case instruction::JAL: return "JAL";
        case instruction::JALR: return "JALR";
        case instruction::BEQ: return "BEQ";
        case instruction::BNE: return "BNE";
        case instruction::BLT: return "BLT";
        case instruction::BGE: return "BGE";
        case instruction::BLTU: return "BLTU";
        case instruction::BGEU: return "BGEU";
        case instruction::LB: return "LB";
        case instruction::LH: return "LH";
        case instruction::LW: return "LW";
        case instruction::LBU: return "LBU";
        case instruction::LHU: return "LHU";
        case instruction::SB: return "SB";
        case instruction::SH: return "SH";
        case instruction::SW: return "SW";
        case instruction::ADDI: return "ADDI";
        case instruction::SLTI: return "SLTI";
        case instruction::SLTIU: return "SLTIU";
        case instruction::XORI: return "XORI";
        case instruction::ORI: return "ORI";
        case instruction::ANDI: return "ANDI";
        case instruction::SLLI: return "SLLI";
        case instruction::SRLI: return "SRLI";
        case instruction::SRAI: return "SRAI";
        case instruction::ADD: return "ADD";
        case instruction::SUB: return "SUB";
        case instruction::SLL: return "SLL";
        case instruction::SLT: return "SLT";
        case instruction::SLTU: return "SLTU";
        case instruction::XOR: return "XOR";
        case instruction::SRL: return "SRL";
        case instruction::SRA: return "SRA";
        case instruction::OR: return "OR";
        case instruction::AND: return "AND";
        case instruction::FENCE: return "FENCE";
        case instruction::ECALL: return "ECALL";
        case instruction::EBREAK: return "EBREAK";
        case instruction::FENCE_I: return "FENCE_I";
        case instruction::CSRRW: return "CSRRW";
        case instruction::CSRRS: return "CSRRS";
        case instruction::CSRRC: return "CSRRC";
        case instruction::CSRRWI: return "CSRRWI";
        case instruction::CSRRSI: return "CSRRSI";
        case instruction::CSRRCI: return "CSRRCI";
        case instruction::MUL: return "MUL";
        case instruction::MULH: return "MULH";
        case instruction::MULHSU: return "MULHSU";
        case instruction::MULHU: return "MULHU";
        case instruction::DIV: return "DIV";
        case instruction::DIVU: return "DIVU";
        case instruction::REM: return "REM";
        case instruction::REMU: return "REMU";
        case instruction::LR_W: return "LR_W";
        case instruction::SC_W: return "SC_W";
        case instruction::AMOSWAP_W: return "AMOSWAP_W";
        case instruction::AMOADD_W: return "AMOADD_W";
        case instruction::AMOXOR_W: return "AMOXOR_W";
        case instruction::AMOAND_W: return "AMOAND_W";
        case instruction::AMOOR_W: return "AMOOR_W";
        case instruction::AMOMIN_W: return "AMOMIN_W";
        case instruction::AMOMAX_W: return "AMOMAX_W";
        case instruction::AMOMINU_W: return "AMOMINU_W";
        case instruction::AMOMAXU_W: return "AMOMAXU_W";
        case instruction::LWU: return "LWU";
        case instruction::LD: return "LD";
        case instruction::SD: return "SD";
        case instruction::ADDIW: return "ADDIW";
        case instruction::SLLIW: return "SLLIW";
        case instruction::SRLIW: return "SRLIW";
        case instruction::SRAIW: return "SRAIW";
        case instruction::ADDW: return "ADDW";
        case instruction::SUBW: return "SUBW";
        case instruction::SLLW: return "SLLW";
        case instruction::SRLW: return "SRLW";
        case instruction::SRAW: return "SRAW";
        case instruction::MULW: return "MULW";
        case instruction::DIVW: return "DIVW";
        case instruction::DIVUW: return "DIVUW";
        case instruction::REMW: return "REMW";
        case instruction::REMUW: return "REMUW";
        case instruction::LR_D: return "LR_D";
        case instruction::SC_D: return "SC_D";
        case instruction::AMOSWAP_D: return "AMOSWAP_D";
        case instruction::AMOADD_D: return "AMOADD_D";
        case instruction::AMOXOR_D: return "AMOXOR_D";
        case instruction::AMOAND_D: return "AMOAND_D";
        case instruction::AMOOR_D: return "AMOOR_D";
        case instruction::AMOMIN_D: return "AMOMIN_D";
        case instruction::AMOMAX_D: return "AMOMAX_D";
        case instruction::AMOMINU_D: return "AMOMINU_D";
        case instruction::AMOMAXU_D: return "AMOMAXU_D";
        case instruction::FLW: return "FLW";
        case instruction::FSW: return "FSW";
        case instruction::FMADD_S: return "FMADD_S";
        case instruction::FMSUB_S: return "FMSUB_S";
        case instruction::FNMADD_S: return "FNMADD_S";
        case instruction::FNMSUB_S: return "FNMSUB_S";
        case instruction::FADD_S: return "FADD_S";
        case instruction::FSUB_S: return "FSUB_S";
        case instruction::FMUL_S: return "FMUL_S";
        case instruction::FDIV_S: return "FDIV_S";
        case instruction::FSQRT_S: return "FSQRT_S";
        case instruction::FSGNJ_S: return "FSGNJ_S";
        case instruction::FSGNJN_S: return "FSGNJN_S";
        case instruction::FSGNJX_S: return "FSGNJX_S";
        case instruction::FMIN_S: return "FMIN_S";
        case instruction::FMAX_S: return "FMAX_S";
        case instruction::FCVT_W_S: return "FCVT_W_S";
        case instruction::FCVT_WU_S: return "FCVT_WU_S";
        case instruction::FMV_X_W: return "FMV_X_W";
        case instruction::FEQ_S: return "FEQ_S";
        case instruction::FLT_S: return "FLT_S";
        case instruction::FLE_S: return "FLE_S";
        case instruction::FCLASS_S: return "FCLASS_S";
        case instruction::FCVT_S_W: return "FCVT_S_W";
        case instruction::FCVT_S_WU: return "FCVT_S_WU";
        case instruction::FMV_W_X: return "FMV_W_X";
        case instruction::FCVT_L_S: return "FCVT_L_S";
        case instruction::FCVT_LU_S: return "FCVT_LU_S";
        case instruction::FCVT_S_L: return "FCVT_S_L";
        case instruction::FCVT_S_LU: return "FCVT_S_LU";
        case instruction::FLD: return "FLD";
        case instruction::FSD: return "FSD";
        case instruction::FMADD_D: return "FMADD_D";
        case instruction::FMSUB_D: return "FMSUB_D";
        case instruction::FNMSUB_D: return "FNMSUB_D";
        case instruction::FNMADD_D: return "FNMADD_D";
        case instruction::FADD_D: return "FADD_D";
        case instruction::FSUB_D: return "FSUB_D";
        case instruction::FMUL_D: return "FMUL_D";
        case instruction::FDIV_D: return "FDIV_D";
        case instruction::FSQRT_D: return "FSQRT_D";
        case instruction::FSGNJ_D: return "FSGNJ_D";
        case instruction::FSGNJN_D: return "FSGNJN_D";
        case instruction::FSGNJX_D: return "FSGNJX_D";
        case instruction::FMIN_D: return "FMIN_D";
        case instruction::FMAX_D: return "FMAX_D";
        case instruction::FCVT_S_D: return "FCVT_S_D";
        case instruction::FCVT_D_S: return "FCVT_D_S";
        case instruction::FEQ_D: return "FEQ_D";
        case instruction::FLT_D: return "FLT_D";
        case instruction::FLE_D: return "FLE_D";
        case instruction::FCLASS_D: return "FCLASS_D";
        case instruction::FCVT_W_D: return "FCVT_W_D";
        case instruction::FCVT_WU_D: return "FCVT_WU_D";
        case instruction::FCVT_D_W: return "FCVT_D_W";
        case instruction::FCVT_D_WU: return "FCVT_D_WU";
        case instruction::FCVT_L_D: return "FCVT_L_D";
        case instruction::FCVT_LU_D: return "FCVT_LU_D";
        case instruction::FMV_X_D: return "FMV_X_D";
        case instruction::FCVT_D_L: return "FCVT_D_L";
        case instruction::FCVT_D_LU: return "FCVT_D_LU";
        case instruction::FMV_D_X: return "FMV_D_X";
        case instruction::URET: return "URET";
        case instruction::SRET: return "SRET";
        case instruction::MRET: return "MRET";
        case instruction::WFI: return "WFI";
        case instruction::SFENCE_VMA: return "SFENCE_VMA";
        case instruction::COUNT:  // fall-through
        default: assert(false); return "";
    }
}

struct path_node
{
    auto friend operator==(path_node const& lhs, path_node const& rhs) noexcept
    {
        return (lhs.instr == rhs.instr);
    }

    instruction instr;

    std::int32_t w;  // number of calls

    std::vector<std::vector<std::int32_t>> dependencies;  // sequences of dependencies (tokens)

    // extensions
    std::vector<std::vector<std::int32_t>> dependencies_from;

    std::vector<std::int32_t> shared_seqs;

    std::vector<std::int32_t> tokens;  // indices
};

class instruction_sequence
{
  public:
    explicit instruction_sequence(std::vector<std::vector<path_node>> const& seqs)
    {
        sort(seqs);  // based on mappings without dependencies

        // prepare sequences for token management
        auto idx = 0;
        auto ofs = 0;
        for (auto i = 0; i < static_cast<std::int32_t>(base_seqs.size()); ++i)
        {
            for (auto& v : base_seqs[i])
            {
                assert(v.dependencies.size() == 1);

                auto const offsets = v.dependencies[0];
                v.dependencies.clear();
                v.dependencies.resize(1);
                v.dependencies_from.resize(1);

                for (auto const offset : offsets)
                {
                    v.dependencies[0].push_back(idx - offset);                             // absolute value
                    base_seqs[i][idx - offset - ofs].dependencies_from[0].push_back(idx);  // reverse dependencies
                }

                v.shared_seqs.push_back(i);
                v.tokens.push_back(idx);

                data.push_back(v);

                tokens.push_back(idx++);
            }
            ofs += static_cast<std::int32_t>(base_seqs[i].size());
        }

        assert(std::accumulate(data.begin(), data.end(), 0, [](auto const sum, auto const& v) {
                   return (sum + v.tokens.size());
               }) == static_cast<std::int32_t>(tokens.size()));
    }

    auto friend operator<<(std::ostream& os, instruction_sequence const& instr_seq) -> std::ostream&
    {
        auto output = false;  // helper
        for (auto i = 0; i < static_cast<std::int32_t>(instr_seq.data.size()); ++i)
        {
            if (!instr_seq.data[i].tokens.empty())
            {  // not rearranged (possibly merged)
                os << instr2str(instr_seq.data[i].instr);
                output = true;

                if (instr_seq.data[i].tokens.size() == 1)
                {  // NOP bypass (not merged)
                    os << "/NOP";
                }
            }

            // lookahead
            if (i + 1 != static_cast<std::int32_t>(instr_seq.data.size()) && !instr_seq.data[i + 1].tokens.empty() &&
                output)
            {  // not the last iteration
                os << "->";
            }
        }
        return os;
    }

    auto optimize()
    {  // bruteforce
        auto cur_score = 0.0f;
        auto ofs = 0;  // first node of the second sequence to be merged
        for (auto k = 0; k < static_cast<std::int32_t>(base_seqs.size()) - 1; ++k)
        {
            auto const old_score = cur_score;
            ofs += static_cast<std::int32_t>(base_seqs[k].size());

            std::vector<path_node> seq1(data.begin(), data.begin() + ofs);
            std::vector<path_node> seq2(data.begin() + ofs,
                                        data.begin() + ofs + static_cast<std::int32_t>(base_seqs[k + 1].size()));
#ifndef NDEBUG
            std::cout << "seq1 = ";
            print(seq1);
            std::cout << std::endl;
            std::cout << "seq2 = ";
            print(seq2);
            std::cout << std::endl;
#endif
            // bottom-up to avoid blocking optimizations
            for (auto i = static_cast<std::int32_t>(seq1.size()) - 1; i >= 0; --i)
            {
                // consider a node with dependencies
                if (seq1[i].dependencies[0].empty() || !seq1[i].dependencies_from[0].empty())
                {
                    continue;
                }

                for (auto j = static_cast<std::int32_t>(seq2.size()) - 1; j >= 0; --j)
                {
                    if (seq1[i] == seq2[j] && std::find(seq2[j].shared_seqs.begin(), seq2[j].shared_seqs.end(), k) ==
                                                  seq2[j].shared_seqs.end())
                    {  // merge into seq2
                        merge(j + ofs, i);
                        seq1 = std::vector(data.begin(), data.begin() + ofs);
                        seq2 = std::vector(data.begin() + ofs,
                                           data.begin() + ofs + static_cast<std::int32_t>(base_seqs[k + 1].size()));
                        break;
                    }
                }
            }

            // top-down
            for (auto i = 0; i < static_cast<std::int32_t>(seq2.size()); ++i)
            {
                // consider a node with reverse dependencies
                if (seq2[i].dependencies_from[0].empty() || !seq2[i].dependencies[0].empty())
                {
                    continue;
                }

                for (auto j = 0; j < static_cast<std::int32_t>(seq1.size()); ++j)
                {
                    if (seq2[i] == seq1[j] && std::find(seq1[j].shared_seqs.begin(), seq1[j].shared_seqs.end(),
                                                        k + 1) == seq1[j].shared_seqs.end())
                    {  // merge into seq1
                        merge(j, i + ofs);
                        seq1 = std::vector(data.begin(), data.begin() + ofs);
                        seq2 = std::vector(data.begin() + ofs,
                                           data.begin() + ofs + static_cast<std::int32_t>(base_seqs[k + 1].size()));
                        break;
                    }
                }
            }

            // map the rest (bottom-up)
            for (auto i = static_cast<std::int32_t>(seq1.size()) - 1; i >= 0; --i)
            {
                for (auto j = static_cast<std::int32_t>(seq2.size()) - 1; j >= 0; --j)
                {
                    if (seq1[i] == seq2[j] && std::find(seq2[j].shared_seqs.begin(), seq2[j].shared_seqs.end(), k) ==
                                                  seq2[j].shared_seqs.end())
                    {
                        merge(j + ofs, i);
                        seq1 = std::vector(data.begin(), data.begin() + ofs);
                        seq2 = std::vector(data.begin() + ofs,
                                           data.begin() + ofs + static_cast<std::int32_t>(base_seqs[k + 1].size()));
                        break;
                    }
                }
            }

            // map the rest (top-down)
            for (auto i = 0; i < static_cast<std::int32_t>(seq2.size()); ++i)
            {
                for (auto j = 0; j < static_cast<std::int32_t>(seq1.size()); ++j)
                {
                    if (seq2[i] == seq1[j] && std::find(seq1[j].shared_seqs.begin(), seq1[j].shared_seqs.end(),
                                                        k + 1) == seq1[j].shared_seqs.end())
                    {
                        merge(j, i + ofs);
                        seq1 = std::vector(data.begin(), data.begin() + ofs);
                        seq2 = std::vector(data.begin() + ofs,
                                           data.begin() + ofs + static_cast<std::int32_t>(base_seqs[k + 1].size()));
                        break;
                    }
                }
            }

            cur_score = score(k + 2);

            assert(cur_score >= 0.0f);

            // terminating case
            if (cur_score <= old_score)
            {
                break;
            }
        }
    }
  private:
    auto size(std::int32_t const k) const noexcept
    {  // depending on the number of sequences
        assert(k <= static_cast<std::int32_t>(base_seqs.size()));

        return std::accumulate(base_seqs.begin(), base_seqs.begin() + k, 0,
                               [](auto const sum, auto const& seq) { return (sum + seq.size()); });
    }

    auto length(std::int32_t const k) const noexcept
    {
        return std::count_if(data.begin(), data.begin() + size(k), [](auto const& v) { return v.tokens.empty(); });
    }

    auto merge(std::int32_t const dst, std::int32_t const src) -> void
    {
        // check for dependency problems
        for (auto const& seq : data[src].dependencies)
        {
            for (auto const idx : seq)
            {
                if (dst < tokens[idx])
                {  // dependency cannot be resolved
                    return;
                }
            }
        }
        for (auto const& seq : data[src].dependencies_from)
        {
            for (auto const idx : seq)
            {
                if (dst > tokens[idx])
                {
                    return;
                }
            }
        }

        // rearrange the source node
        data[dst].dependencies.insert(data[dst].dependencies.end(), data[src].dependencies.begin(),
                                      data[src].dependencies.end());
        data[dst].dependencies_from.insert(data[dst].dependencies_from.end(), data[src].dependencies_from.begin(),
                                           data[src].dependencies_from.end());
        data[dst].shared_seqs.insert(data[dst].shared_seqs.end(), data[src].shared_seqs.begin(),
                                     data[src].shared_seqs.end());
        data[dst].tokens.insert(data[dst].tokens.end(), data[src].tokens.begin(), data[src].tokens.end());
        while (!data[src].tokens.empty())
        {
            tokens[data[src].tokens.back()] = dst;
            data[src].tokens.pop_back();
        }
    }

    auto print(std::vector<path_node> const& seq) const -> void
    {
        auto output = false;
        for (auto i = 0; i < static_cast<std::int32_t>(seq.size()); ++i)
        {
            if (!seq[i].tokens.empty())
            {
                std::cout << instr2str(seq[i].instr);
                output = true;

                if (seq[i].tokens.size() == 1)
                {
                    std::cout << "/NOP";
                }
            }

            if (i + 1 != static_cast<std::int32_t>(seq.size()) && !seq[i + 1].tokens.empty() && output)
            {
                std::cout << "->";
            }
        }
    }

    auto score(std::int32_t const k) const noexcept -> float
    {
        assert(k <= static_cast<std::int32_t>(base_seqs.size()));

        auto score = 0.0f;
        for (auto i = 0; i < k; ++i)
        {
            score +=
                base_seqs[i].back().w * static_cast<std::int32_t>(base_seqs[i].size());  // weights are already summed
            // for execution time (NOP factor)
            score -= (length(i) - static_cast<std::int32_t>(base_seqs[i].size())) * base_seqs[i].back().w * 0.1f;
        }
        return score;
    }

    auto sort(std::vector<std::vector<path_node>> const& seqs) -> void
    {
        // find start sequence
        std::vector<path_node> seq0;
        auto max_score = 0;
        for (auto const& seq : seqs)
        {
            auto const score = seq.back().w * static_cast<std::int32_t>(seq.size());
            if (score >= max_score)
            {
                seq0 = seq;
                max_score = score;
            }
        }

        base_seqs.reserve(seqs.size());
        base_seqs.push_back(std::move(seq0));

        // maximize mapping
        for (auto k = 1; k < static_cast<std::int32_t>(seqs.size()); ++k)
        {
            std::vector<path_node> tmp_seq;
            auto max_count = 0;
            for (auto const& seq : seqs)
            {
                if (std::find(base_seqs.begin(), base_seqs.end(), seq) != base_seqs.end())
                {  // this sequence is already a base sequence
                    continue;
                }

                auto count = 0;
                for (auto const& base_seq : base_seqs)
                {
                    if (std::search(seq.begin(), seq.end(), base_seq.begin(), base_seq.end()) == seq.end() &&
                        std::search(base_seq.begin(), base_seq.end(), seq.begin(), seq.end()) == base_seq.end())
                    {  // there is no subsequence
                        for (auto const& v : base_seq)
                        {
                            count += std::count(seq.begin(), seq.end(), v);
                        }
                    }
                }
                if (count >= max_count)
                {
                    tmp_seq = seq;
                    max_count = count;
                }
            }
            base_seqs.push_back(std::move(tmp_seq));
        }

        assert(base_seqs.size() == seqs.size());
    }

    std::vector<std::vector<path_node>> base_seqs;  // to determine the boundaries for the instruction sequence (data)

    std::vector<path_node> data;  // merged

    std::vector<std::int32_t> tokens;
};
