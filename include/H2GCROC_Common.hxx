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
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "easylogging++.h"
#include "argparse/argparse.hpp"

#ifndef FPGA_CHANNEL_NUMBER
#define FPGA_CHANNEL_NUMBER 152
#endif

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

ScriptOptions parse_arguments(int argc, char **argv, const std::string& version = "0.1");

#endif // COMMON_HPP