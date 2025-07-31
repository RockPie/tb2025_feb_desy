#pragma once

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <regex>

#include "TCanvas.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TH2.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TVirtualFitter.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadExecutor.hxx>


#include "easylogging++.h"
#include "nlohmann/json.hpp"
#include "argparse/argparse.hpp"

#include "H2GCROC_Common.hxx"

bool readRootMetaData(const char* _root_file_path, TFile*& _input_root, TTree*& _input_tree, int& _fpga_count, int& _machine_gun_samples, int& _max_entry, std::vector<UShort_t>& _legal_fpga_id_list);

bool drawTwoASIC(TCanvas *_canvas, TLegend *_legend, int _bin, std::vector<std::vector<double>> &_fpgas_data, std::vector<unsigned short> _fpga_id, std::vector<Color_t> &_fpga_colors, std::string _y_title);