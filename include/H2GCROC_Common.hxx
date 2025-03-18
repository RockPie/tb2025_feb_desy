#ifndef COMMON_HPP
#define COMMON_HPP

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <string>

#include "TCanvas.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "TPaveText.h"


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
    bool focal;
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

class GlobalChannelPainter {
public:
    GlobalChannelPainter(const std::string& mapping_file);
    GlobalChannelPainter(const std::string& mapping_file, const std::string& channel_mapping_file);
    ~GlobalChannelPainter();

    TCanvas* get_canvas() { return painter_canvas; }
    void draw_global_channel_hists2D(std::vector <TH2D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_global_channel_hists1D(std::vector <TH1D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_global_channel_hists1D_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels);
    void draw_global_channel_hists1D_run_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels);
    void draw_global_channel_canvas(std::vector <TCanvas*> canvas_list,  std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_canvas_components(TCanvas* source_canvas, TPad* target_pad);

    TCanvas* draw_module_channel_canvas(std::vector <TCanvas*> canvas_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, int module_index);

private:
    void clear_canvas();

private:
    TCanvas *painter_canvas;
    json mapping_json;

    json focal_mapping_json;
    json focal_channel_mapping_json;

    bool is_EEEMCal_mapping;
    bool is_FoCal_mapping;

    std::vector <int> module_fpga_list;
    std::vector <int> module_asic_list;
    std::vector <int> module_connector_list;
    std::vector <std::vector <int>> connector_list_list;

    std::vector <std::vector <int>> focal_module_board_list;
    std::vector <std::vector <int>> focal_module_channel_list;
    std::unordered_map <int, int> focal_channel_map;
    std::unordered_map <int, int> focal_fpga_map;

    std::vector <TCanvas*> sub_canvas_list;
};

#endif // COMMON_HPP