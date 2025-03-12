#ifndef COMMON_HPP
#define COMMON_HPP

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <string>

#include "TCanvas.h"
#include "TVectorD.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "easylogging++.h"
#include "nlohmann/json.hpp"
#include "argparse/argparse.hpp"

#ifndef FPGA_CHANNEL_NUMBER
#define FPGA_CHANNEL_NUMBER 152
#endif

#ifndef FPGA_CHANNEL_NUMBER_VALID
#define FPGA_CHANNEL_NUMBER_VALID 144
#endif

using json = nlohmann::json;

void set_easylogger();

struct ScriptOptions {
    std::string input_file;
    std::string output_file;
    std::string output_folder;
    int n_events;
    bool verbose;
    std::string script_name;
    std::string script_version;
};

ScriptOptions parse_arguments_single_root(int argc, char **argv, const std::string& version = "0.1");

ScriptOptions parse_arguments_single_json(int argc, char **argv, const std::string& version = "0.1");

int get_valid_fpga_channel(int fpga_channel);
int get_total_fpga_channel(int fpga_channel);
inline int get_unified_fpga_channel(int fpga_id, int fpga_channel){
    return fpga_id * FPGA_CHANNEL_NUMBER + fpga_channel;
}
inline int get_unified_valid_fpga_channel(int fpga_id, int fpga_channel){
    return fpga_id * FPGA_CHANNEL_NUMBER_VALID + fpga_channel;
}

inline UInt_t decode_tot_value(UInt_t val1) {
    UInt_t mask = -(val1 >= 512);  // mask is 0xFFFFFFFF if val1 >= 512, else 0
    return (val1 & ~mask) | ((val1 - 512) * 8 & mask);
}

inline double decode_toa_value_ns(UInt_t val2) {
    constexpr double scale0 = 0.025;
    constexpr double scale1 = 0.2;
    constexpr double scale2 = 6.25;

    UInt_t part0 = val2 & 0x07;
    UInt_t part1 = (val2 >> 3) & 0x1F;
    UInt_t part2 = (val2 >> 8) & 0x1F;

    return part0 * scale0 + part1 * scale1 + part2 * scale2;
}

#endif // COMMON_HPP