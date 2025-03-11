#include <iostream>
#include <unistd.h>
#include "TCanvas.h" 
#include "TVectorD.h"
#include "TVector.h"
#include "TF1.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "easylogging++.h"
#include "argparse/argparse.hpp"
#include <vector>
#include <unordered_map>
#include <sstream>
#include <string>
#include <cmath>

#ifndef FPGA_CHANNEL_NUMBER
#define FPGA_CHANNEL_NUMBER 152
#endif

using namespace std;

INITIALIZE_EASYLOGGINGPP

void set_easylogger();

struct LCSResult {
    vector<Long64_t> sequence;
    // indices[i] holds the positions from arrays[i] that form the subsequence.
    vector<vector<int>> indices;
};

// Helper function to create a unique key from the current state.
// The state consists of the index in the reference array and the current starting indices in each of the other arrays.
string makeKey(int pos0, const vector<int>& positions) {
    stringstream ss;
    ss << pos0;
    for (int p : positions) {
        ss << ',' << p;
    }
    return ss.str();
}

LCSResult longestMatchingSequenceRec(std::vector<Long64_t> arrays[], int n, int pos0, const vector<int>& pos, int variance, unordered_map<string, LCSResult>& memo);

LCSResult longestMatchingSequence(std::vector<Long64_t> arrays[], int n, int variance);

int main(int argc, char **argv){

    // * --- Initial settings -----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string script_input_file, script_output_file;
    int script_n_events;
    bool script_verbose = false;
    std::string script_name = __FILE__;
    std::string script_version = "0.1";
    std::string script_output_folder;

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    argparse::ArgumentParser program(script_name, script_version);

    program.add_argument("-f", "--file").help("Input .root file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    program.add_argument("-v", "--verbose").help("Verbose mode").default_value(false).implicit_value(true);
    try {
        program.parse_args(argc, argv);
        script_input_file  = program.get<std::string>("--file");
        script_output_file = program.get<std::string>("--output");
        auto script_n_events_str = program.get<std::string>("--events");
        script_n_events    = std::stoi(script_n_events_str);
        script_verbose     = program.get<bool>("--verbose");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        return 1;
    }

    LOG(INFO) << "Input file: " << script_input_file;

    if (access(script_input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << script_input_file << " does not exist!";
        return 1;
    }
    if (script_input_file.substr(script_input_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Input file " << script_input_file << " should end with .root!";
        return 1;
    }
    if (script_output_file.substr(script_output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << script_output_file << " should end with .root!";
        return 1;
    }
    script_output_folder = script_output_file.substr(0, script_output_file.find_last_of("/\\"));
    if (script_output_folder.empty()) {
        script_output_folder = "./dump/" + script_name;
    }
    if (access(script_output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << script_output_folder;
        if (mkdir(script_output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << script_output_folder;
            return 1;
        }
    }
    if (access(script_output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << script_output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << script_name;
    LOG(INFO) << "Input file: " << script_input_file;
    LOG(INFO) << "Output file: " << script_output_file << " in " << script_output_folder;
    LOG(INFO) << "Number of events: " << script_n_events;

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int machine_gun_samples = 0;
    int fpga_number = 0;

    TFile *input_root = new TFile(script_input_file.c_str(), "READ");
    if (input_root->IsZombie()) {
        LOG(ERROR) << "Failed to open input file " << script_input_file;
        return 1;
    }
    TTree *input_tree = (TTree*) input_root->Get("data_tree");
    if (input_tree == nullptr) {
        LOG(ERROR) << "Failed to get data tree from input file " << script_input_file;
        return 1;
    }
    TNamed *legal_fpga_id_list_tnamed = (TNamed*) input_root->Get("Rootifier_legal_fpga_id_list");
    if (legal_fpga_id_list_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get legal fpga id list from input file " << script_input_file;
        return 1;
    }
    std::string legal_fpga_id_list_str = legal_fpga_id_list_tnamed->GetTitle();
    std::vector <UShort_t> legal_fpga_id_list;
    std::istringstream legal_fpga_id_list_stream(legal_fpga_id_list_str);
    UShort_t legal_fpga_id;
    while (legal_fpga_id_list_stream >> legal_fpga_id) {
        legal_fpga_id_list.push_back(legal_fpga_id);
    }
    fpga_number = legal_fpga_id_list.size();
    TNamed *input_machine_gun_samples_tnamed = (TNamed*) input_root->Get("EventRecon_machine_gun_samples");
    if (input_machine_gun_samples_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get machine gun samples from input file " << script_input_file;
        return 1;
    }
    machine_gun_samples = std::stoi(input_machine_gun_samples_tnamed->GetTitle());
    LOG(INFO) << "Number of FPGAs: " << fpga_number;
    LOG(INFO) << "Machine gun samples: " << machine_gun_samples;

    UShort_t input_fpga_id;
    ULong64_t input_timestamps[machine_gun_samples];
    UInt_t input_daqh_list[4 * machine_gun_samples];
    Bool_t input_tc_list[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    Bool_t input_tp_list[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    UInt_t input_val0_list[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    UInt_t input_val1_list[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    UInt_t input_val2_list[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    UInt_t input_crc32_list[4 * machine_gun_samples];
    UInt_t input_last_heartbeat[machine_gun_samples];

    input_tree->SetBranchAddress("fpga_id", &input_fpga_id);
    input_tree->SetBranchAddress("timestamps", input_timestamps);
    input_tree->SetBranchAddress("daqh_list", input_daqh_list);
    input_tree->SetBranchAddress("tc_list", input_tc_list);
    input_tree->SetBranchAddress("tp_list", input_tp_list);
    input_tree->SetBranchAddress("val0_list", input_val0_list);
    input_tree->SetBranchAddress("val1_list", input_val1_list);
    input_tree->SetBranchAddress("val2_list", input_val2_list);
    input_tree->SetBranchAddress("crc32_list", input_crc32_list);
    input_tree->SetBranchAddress("last_heartbeat", input_last_heartbeat);

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << script_output_file;
        return 1;
    }
    TTree* output_tree = new TTree("data_tree", "Data tree");
    output_tree->SetDirectory(output_root);

    std::vector <ULong64_t*> branch_timestamps_list;
    std::vector <UInt_t*> branch_daqh_list_list;
    std::vector <Bool_t*> branch_tc_list_list;
    std::vector <Bool_t*> branch_tp_list_list;
    std::vector <UInt_t*> branch_val0_list_list;
    std::vector <UInt_t*> branch_val1_list_list;
    std::vector <UInt_t*> branch_val2_list_list;
    std::vector <UInt_t*> branch_crc32_list_list;
    std::vector <UInt_t*> branch_last_heartbeat_list;

    for (int i = 0; i < fpga_number; i++) {
        auto _fpga_id = legal_fpga_id_list[i];
        auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
        auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
        auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
        auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];

        output_tree->Branch(("timestamps_" + std::to_string(_fpga_id)).c_str(), branch_timestamps, ("timestamps_" + std::to_string(_fpga_id) + "[" + std::to_string(machine_gun_samples) + "]/l").c_str());
        output_tree->Branch(("daqh_list_" + std::to_string(_fpga_id)).c_str(), branch_daqh_list, ("daqh_list_" + std::to_string(_fpga_id) + "[" + std::to_string(4 * machine_gun_samples) + "]/i").c_str());
        output_tree->Branch(("tc_list_" + std::to_string(_fpga_id)).c_str(), branch_tc_list, ("tc_list_" + std::to_string(_fpga_id) + "[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/O").c_str());
        output_tree->Branch(("tp_list_" + std::to_string(_fpga_id)).c_str(), branch_tp_list, ("tp_list_" + std::to_string(_fpga_id) + "[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/O").c_str());
        output_tree->Branch(("val0_list_" + std::to_string(_fpga_id)).c_str(), branch_val0_list, ("val0_list_" + std::to_string(_fpga_id) + "[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/i").c_str());
        output_tree->Branch(("val1_list_" + std::to_string(_fpga_id)).c_str(), branch_val1_list, ("val1_list_" + std::to_string(_fpga_id) + "[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/i").c_str());
        output_tree->Branch(("val2_list_" + std::to_string(_fpga_id)).c_str(), branch_val2_list, ("val2_list_" + std::to_string(_fpga_id) + "[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/i").c_str());
        output_tree->Branch(("crc32_list_" + std::to_string(_fpga_id)).c_str(), branch_crc32_list, ("crc32_list_" + std::to_string(_fpga_id) + "[" + std::to_string(4 * machine_gun_samples) + "]/i").c_str());
        output_tree->Branch(("last_heartbeat_" + std::to_string(_fpga_id)).c_str(), branch_last_heartbeat, ("last_heartbeat_" + std::to_string(_fpga_id) + "[" + std::to_string(machine_gun_samples) + "]/i").c_str());

        branch_timestamps_list.push_back(branch_timestamps);
        branch_daqh_list_list.push_back(branch_daqh_list);
        branch_tc_list_list.push_back(branch_tc_list);
        branch_tp_list_list.push_back(branch_tp_list);
        branch_val0_list_list.push_back(branch_val0_list);
        branch_val1_list_list.push_back(branch_val1_list);
        branch_val2_list_list.push_back(branch_val2_list);
        branch_crc32_list_list.push_back(branch_crc32_list);
        branch_last_heartbeat_list.push_back(branch_last_heartbeat);
    }

    // ! --- Do the event matching between FPGAs here -----------------------------------
    // ! --------------------------------------------------------------------------------
    const int SWMA_window_size = 50;
    const int SWMA_window_max = 200; // if one pool reaches this size, all pools are cleared
    const int SWMA_variance = 200;
    int entry_max = input_tree->GetEntries();
    if (script_n_events > 0 && script_n_events < entry_max) {
        entry_max = script_n_events;
    }
    bool flag_last_entry = false;
    std::vector <ULong64_t> timestamp_pools[fpga_number];
    std::vector <ULong64_t*> timestamp_pools_original[fpga_number];
    std::vector <UInt_t*> daqh_list_pools[fpga_number];
    std::vector <Bool_t*> tc_list_pools[fpga_number];
    std::vector <Bool_t*> tp_list_pools[fpga_number];
    std::vector <UInt_t*> val0_list_pools[fpga_number];
    std::vector <UInt_t*> val1_list_pools[fpga_number];
    std::vector <UInt_t*> val2_list_pools[fpga_number];
    std::vector <UInt_t*> crc32_list_pools[fpga_number];
    std::vector <UInt_t*> last_heartbeat_pools[fpga_number];

    std::vector <int> matched_length;
    UShort_t legal_fpga_ids[fpga_number];

    for (int i = 0; i < fpga_number; i++) {
        legal_fpga_ids[i] = legal_fpga_id_list[i];
        timestamp_pools[i].reserve(SWMA_window_size);
        timestamp_pools_original[i].reserve(SWMA_window_size);
        daqh_list_pools[i].reserve(SWMA_window_size);
        tc_list_pools[i].reserve(SWMA_window_size);
        tp_list_pools[i].reserve(SWMA_window_size);
        val0_list_pools[i].reserve(SWMA_window_size);
        val1_list_pools[i].reserve(SWMA_window_size);
        val2_list_pools[i].reserve(SWMA_window_size);
        crc32_list_pools[i].reserve(SWMA_window_size);
        last_heartbeat_pools[i].reserve(SWMA_window_size);

        matched_length.reserve(int(entry_max / fpga_number / SWMA_window_size));
    }
    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        if (_entry == entry_max - 1) {
            flag_last_entry = true;
        }
        auto _fpga_id = input_fpga_id;
        auto _fpga_index = -1;
        for (int i = 0; i < fpga_number; i++) {
            if (_fpga_id == legal_fpga_ids[i]) {
                _fpga_index = i;
                break;
            }
        }
        auto _timestamp = input_timestamps[0];
        // check if the fpga_id is legal
        bool flag_legal_fpga_id = false;
        for (int i = 0; i < fpga_number; i++) {
            if (_fpga_id == legal_fpga_ids[i]) {
                flag_legal_fpga_id = true;
                break;
            }
        }
        if (!flag_legal_fpga_id) {
            LOG(ERROR) << "Illegal FPGA ID: " << _fpga_id;
            continue;
        }

        auto _timestamp_original = new ULong64_t[machine_gun_samples];
        auto _daqh_list = new UInt_t[4 * machine_gun_samples];
        auto _tc_list = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto _tp_list = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto _val0_list = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto _val1_list = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto _val2_list = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto _crc32_list = new UInt_t[4 * machine_gun_samples];
        auto _last_heartbeat = new UInt_t[machine_gun_samples];

        std::copy(input_timestamps, input_timestamps + machine_gun_samples, _timestamp_original);
        std::copy(input_daqh_list, input_daqh_list + 4 * machine_gun_samples, _daqh_list);
        std::copy(input_tc_list, input_tc_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, _tc_list);
        std::copy(input_tp_list, input_tp_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, _tp_list);
        std::copy(input_val0_list, input_val0_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, _val0_list);
        std::copy(input_val1_list, input_val1_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, _val1_list);
        std::copy(input_val2_list, input_val2_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, _val2_list);
        std::copy(input_crc32_list, input_crc32_list + 4 * machine_gun_samples, _crc32_list);
        std::copy(input_last_heartbeat, input_last_heartbeat + machine_gun_samples, _last_heartbeat);

        timestamp_pools[_fpga_index].push_back(_timestamp);
        timestamp_pools_original[_fpga_index].push_back(_timestamp_original);
        daqh_list_pools[_fpga_index].push_back(_daqh_list);
        tc_list_pools[_fpga_index].push_back(_tc_list);
        tp_list_pools[_fpga_index].push_back(_tp_list);
        val0_list_pools[_fpga_index].push_back(_val0_list);
        val1_list_pools[_fpga_index].push_back(_val1_list);
        val2_list_pools[_fpga_index].push_back(_val2_list);
        crc32_list_pools[_fpga_index].push_back(_crc32_list);
        last_heartbeat_pools[_fpga_index].push_back(_last_heartbeat);

        bool flag_pools_full = true;
        for (int i = 0; i < fpga_number; i++) {
            if (timestamp_pools[i].size() < SWMA_window_size) {
                flag_pools_full = false;
                break;
            }
        }
        bool flag_last_entry_valid = flag_last_entry;
        if (flag_last_entry){
            // check if there are at least 2 entries in each pool
            for (int i = 0; i < fpga_number; i++) {
                if (timestamp_pools[i].size() < 2) {
                    flag_last_entry_valid = false;
                    break;
                }
            }
        }
        if (!flag_pools_full && !flag_last_entry_valid) {
            bool flag_clear_pools = false;
            for (int i = 0; i < fpga_number; i++) {
                if (timestamp_pools[i].size() > SWMA_window_max) {
                    flag_clear_pools = true;
                    break;
                }
            }
            if (flag_clear_pools) {
                for (int i = 0; i < fpga_number; i++) {
                    if (script_verbose) {
                        LOG(DEBUG) << "Clearing pool " << i;
                    }
                    for (int j = 0; j < timestamp_pools[i].size(); j++) {
                        if (timestamp_pools_original[i][j] != nullptr) {
                            delete[] timestamp_pools_original[i][j];
                        }
                    }
                    for (int j = 0; j < daqh_list_pools[i].size(); j++) {
                        if (daqh_list_pools[i][j] != nullptr) {
                            delete[] daqh_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < tc_list_pools[i].size(); j++) {
                        if (tc_list_pools[i][j] != nullptr) {
                            delete[] tc_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < tp_list_pools[i].size(); j++) {
                        if (tp_list_pools[i][j] != nullptr) {
                            delete[] tp_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < val0_list_pools[i].size(); j++) {
                        if (val0_list_pools[i][j] != nullptr) {
                            delete[] val0_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < val1_list_pools[i].size(); j++) {
                        if (val1_list_pools[i][j] != nullptr) {
                            delete[] val1_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < val2_list_pools[i].size(); j++) {
                        if (val2_list_pools[i][j] != nullptr) {
                            delete[] val2_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < crc32_list_pools[i].size(); j++) {
                        if (crc32_list_pools[i][j] != nullptr) {
                            delete[] crc32_list_pools[i][j];
                        }
                    }
                    for (int j = 0; j < last_heartbeat_pools[i].size(); j++) {
                        if (last_heartbeat_pools[i][j] != nullptr) {
                            delete[] last_heartbeat_pools[i][j];
                        }
                    }
                    timestamp_pools[i].clear();
                    timestamp_pools_original[i].clear();
                    daqh_list_pools[i].clear();
                    tc_list_pools[i].clear();
                    tp_list_pools[i].clear();
                    val0_list_pools[i].clear();
                    val1_list_pools[i].clear();
                    val2_list_pools[i].clear();
                    crc32_list_pools[i].clear();
                    last_heartbeat_pools[i].clear();
                }
            }
            continue;
        }
        if (flag_pools_full || flag_last_entry_valid) {
            // Calculate the difference between the timestamps
            std::vector <Long64_t> timestamp_diffs[fpga_number];
            for (int i = 0; i < fpga_number; i++) {
                int _diff_size = SWMA_window_size - 1;
                if (flag_last_entry) {
                    _diff_size = timestamp_pools[i].size() - 1;
                }
                timestamp_diffs[i].reserve(_diff_size);
                for (int j = 0; j < _diff_size; j++) {
                    timestamp_diffs[i].push_back(Long64_t(timestamp_pools[i][j + 1] - timestamp_pools[i][j]));
                }
                // print the timestamp_diffs
                if (script_verbose) {
                    std::ostringstream oss;
                    oss << "FPGA " << i << " ";
                    for (int j = 0; j < _diff_size; j++) {
                        oss << timestamp_diffs[i][j] << " ";
                    }
                    LOG(INFO) << oss.str();
                }
            }

            LCSResult _result = longestMatchingSequence(timestamp_diffs, fpga_number, SWMA_variance);
            if (script_verbose) {
                cout << "Longest matching sequence: ";
                for (Long64_t val : _result.sequence) {
                    cout << val << " ";
                }
                cout << endl;
                cout << "Indices: " << endl;
                for (int j = 0; j < fpga_number; j++) {
                    cout << "FPGA " << j << ": ";
                    for (int index : _result.indices[j]) {
                        cout << index << " ";
                    }
                    cout << endl;
                }
            }
            auto _result_sequence_size = _result.sequence.size();
    
            if (_result_sequence_size > 0) {
                matched_length.push_back(_result_sequence_size);
                for (int _matched_index = 0; _matched_index < _result_sequence_size; _matched_index++) {
                    for (int i = 0; i < fpga_number; i++) {
                        auto _index          = _result.indices[i][_matched_index];
                        auto _timestamp      = timestamp_pools_original[i][_index];
                        auto _daqh_list      = daqh_list_pools[i][_index];
                        auto _tc_list        = tc_list_pools[i][_index];
                        auto _tp_list        = tp_list_pools[i][_index];
                        auto _val0_list      = val0_list_pools[i][_index];
                        auto _val1_list      = val1_list_pools[i][_index];
                        auto _val2_list      = val2_list_pools[i][_index];
                        auto _crc32_list     = crc32_list_pools[i][_index];
                        auto _last_heartbeat = last_heartbeat_pools[i][_index];

                        auto *branch_timestamps     = branch_timestamps_list[i];
                        auto *branch_daqh_list      = branch_daqh_list_list[i];
                        auto *branch_tc_list        = branch_tc_list_list[i];
                        auto *branch_tp_list        = branch_tp_list_list[i];
                        auto *branch_val0_list      = branch_val0_list_list[i];
                        auto *branch_val1_list      = branch_val1_list_list[i];
                        auto *branch_val2_list      = branch_val2_list_list[i];
                        auto *branch_crc32_list     = branch_crc32_list_list[i];
                        auto *branch_last_heartbeat = branch_last_heartbeat_list[i];

                        std::copy(_timestamp, _timestamp + machine_gun_samples, branch_timestamps);
                        std::copy(_daqh_list, _daqh_list + 4 * machine_gun_samples, branch_daqh_list);
                        std::copy(_tc_list, _tc_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, branch_tc_list);
                        std::copy(_tp_list, _tp_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, branch_tp_list);
                        std::copy(_val0_list, _val0_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, branch_val0_list);
                        std::copy(_val1_list, _val1_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, branch_val1_list);
                        std::copy(_val2_list, _val2_list + FPGA_CHANNEL_NUMBER * machine_gun_samples, branch_val2_list);
                        std::copy(_crc32_list, _crc32_list + 4 * machine_gun_samples, branch_crc32_list);
                        std::copy(_last_heartbeat, _last_heartbeat + machine_gun_samples, branch_last_heartbeat);
                    }
                    output_tree->Fill();
                }
            }
            // Erase the part til the last matched indices (fix to delete memory before erase)
            for (int i = 0; i < fpga_number; i++) {
                if (_result.indices[i].size() == 0) {
                    continue;
                }
                int eraseCount = _result.indices[i].back() + 1; // number of elements to remove
                // Free memory for each pointer being removed.
                for (int j = 0; j < eraseCount; j++) {
                    if (j < timestamp_pools_original[i].size()) {
                        delete[] timestamp_pools_original[i][j];
                    }
                    if (j < daqh_list_pools[i].size()) {
                        delete[] daqh_list_pools[i][j];
                    }
                    if (j < tc_list_pools[i].size()) {
                        delete[] tc_list_pools[i][j];
                    }
                    if (j < tp_list_pools[i].size()) {
                        delete[] tp_list_pools[i][j];
                    }
                    if (j < val0_list_pools[i].size()) {
                        delete[] val0_list_pools[i][j];
                    }
                    if (j < val1_list_pools[i].size()) {
                        delete[] val1_list_pools[i][j];
                    }
                    if (j < val2_list_pools[i].size()) {
                        delete[] val2_list_pools[i][j];
                    }
                    if (j < crc32_list_pools[i].size()) {
                        delete[] crc32_list_pools[i][j];
                    }
                    if (j < last_heartbeat_pools[i].size()) {
                        delete[] last_heartbeat_pools[i][j];
                    }
                }
                // Erase the pointers from the vectors now that the memory is freed.
                timestamp_pools[i].erase(timestamp_pools[i].begin(), timestamp_pools[i].begin() + eraseCount);
                timestamp_pools_original[i].erase(timestamp_pools_original[i].begin(), timestamp_pools_original[i].begin() + eraseCount);
                daqh_list_pools[i].erase(daqh_list_pools[i].begin(), daqh_list_pools[i].begin() + eraseCount);
                tc_list_pools[i].erase(tc_list_pools[i].begin(), tc_list_pools[i].begin() + eraseCount);
                tp_list_pools[i].erase(tp_list_pools[i].begin(), tp_list_pools[i].begin() + eraseCount);
                val0_list_pools[i].erase(val0_list_pools[i].begin(), val0_list_pools[i].begin() + eraseCount);
                val1_list_pools[i].erase(val1_list_pools[i].begin(), val1_list_pools[i].begin() + eraseCount);
                val2_list_pools[i].erase(val2_list_pools[i].begin(), val2_list_pools[i].begin() + eraseCount);
                crc32_list_pools[i].erase(crc32_list_pools[i].begin(), crc32_list_pools[i].begin() + eraseCount);
                last_heartbeat_pools[i].erase(last_heartbeat_pools[i].begin(), last_heartbeat_pools[i].begin() + eraseCount);
            }
        }
    }

    // print the matched_length
    if (script_verbose) {
        std::ostringstream oss;
        oss << "Matched length: ";
        for (int i = 0; i < matched_length.size(); i++) {
            oss << matched_length[i] << " ";
        }
        LOG(INFO) << oss.str();
    }

    Long64_t total_matched_events = 0;
    for (int i = 0; i < matched_length.size(); i++) {
        total_matched_events += matched_length[i] + 1;
    }

    if (entry_max == 0) {
        LOG(ERROR) << "No events in the input file!";
    } else {
        if (fpga_number == 0) {
            LOG(ERROR) << "No legal FPGA ID in the input file!";
        }
        else {
            LOG(INFO) << "Total matched events: " << total_matched_events << " (" << (float) total_matched_events / entry_max * 100 * fpga_number << "%)";
        }
    }

    // Read the meta data from the input file
    TNamed *input_script_rootifier_name_tnamed = (TNamed*) input_root->Get("Rootifier_script_name");
    TNamed *input_script_rootifier_version_tnamed = (TNamed*) input_root->Get("Rootifier_script_version");
    TNamed *input_script_rootifier_input_file_tnamed = (TNamed*) input_root->Get("Rootifier_script_input_file");
    TNamed *input_script_rootifier_output_file_tnamed = (TNamed*) input_root->Get("Rootifier_script_output_file");
    TNamed *input_script_rootifier_n_events_tnamed = (TNamed*) input_root->Get("Rootifier_script_n_events");
    TNamed *input_script_rootifier_output_folder_tnamed = (TNamed*) input_root->Get("Rootifier_script_output_folder");

    TNamed *input_script_eventrecon_name_tnamed = (TNamed*) input_root->Get("EventRecon_script_name");
    TNamed *input_script_eventrecon_version_tnamed = (TNamed*) input_root->Get("EventRecon_script_version");
    TNamed *input_script_eventrecon_input_file_tnamed = (TNamed*) input_root->Get("EventRecon_script_input_file");
    TNamed *input_script_eventrecon_output_file_tnamed = (TNamed*) input_root->Get("EventRecon_script_output_file");
    TNamed *input_script_eventrecon_n_events_tnamed = (TNamed*) input_root->Get("EventRecon_script_n_events");
    TNamed *input_script_eventrecon_output_folder_tnamed = (TNamed*) input_root->Get("EventRecon_script_output_folder");

    // Read the running info from the input file
    TNamed *input_counter_header_line_tnamed = (TNamed*) input_root->Get("Rootifier_counter_header_line");
    TNamed *input_counter_heartbeat_packet_tnamed = (TNamed*) input_root->Get("Rootifier_counter_heartbeat_packet");
    TNamed *input_counter_complete_20_lines_tnamed = (TNamed*) input_root->Get("Rootifier_counter_complete_20_lines");
    TNamed *input_counter_not_complete_20_lines_tnamed = (TNamed*) input_root->Get("Rootifier_counter_not_complete_20_lines");
    TNamed *input_counter_not_complete_20_lines_total_tnamed = (TNamed*) input_root->Get("Rootifier_counter_not_complete_20_lines_total");
    TNamed *input_counter_invalid_line_tnamed = (TNamed*) input_root->Get("Rootifier_counter_invalid_line");
    TNamed *input_counter_valid_line_tnamed = (TNamed*) input_root->Get("Rootifier_counter_valid_line");
    TNamed *input_legal_fpga_id_list_tnamed = (TNamed*) input_root->Get("Rootifier_legal_fpga_id_list");

    TNamed *input_eventrecon_good_machine_gun_events_tnamed = (TNamed*) input_root->Get("EventRecon_good_machine_gun_events");
    TNamed *input_eventrecon_bad_machine_gun_events_tnamed = (TNamed*) input_root->Get("EventRecon_bad_machine_gun_events");
    TNamed *input_eventrecon_machine_gun_samples_tnamed = (TNamed*) input_root->Get("EventRecon_machine_gun_samples");
    input_root->Close();

    output_root->cd();

    output_tree->Write();

    // -- Write the meta data from the input file
    if (input_script_rootifier_name_tnamed != nullptr) {
        input_script_rootifier_name_tnamed->Write();
    }
    if (input_script_rootifier_version_tnamed != nullptr) {
        input_script_rootifier_version_tnamed->Write();
    }
    if (input_script_rootifier_input_file_tnamed != nullptr) {
        input_script_rootifier_input_file_tnamed->Write();
    }
    if (input_script_rootifier_output_file_tnamed != nullptr) {
        input_script_rootifier_output_file_tnamed->Write();
    }
    if (input_script_rootifier_n_events_tnamed != nullptr) {
        input_script_rootifier_n_events_tnamed->Write();
    }
    if (input_script_rootifier_output_folder_tnamed != nullptr) {
        input_script_rootifier_output_folder_tnamed->Write();
    }

    if (input_script_eventrecon_name_tnamed != nullptr) {
        input_script_eventrecon_name_tnamed->Write();
    }
    if (input_script_eventrecon_version_tnamed != nullptr) {
        input_script_eventrecon_version_tnamed->Write();
    }
    if (input_script_eventrecon_input_file_tnamed != nullptr) {
        input_script_eventrecon_input_file_tnamed->Write();
    }
    if (input_script_eventrecon_output_file_tnamed != nullptr) {
        input_script_eventrecon_output_file_tnamed->Write();
    }
    if (input_script_eventrecon_n_events_tnamed != nullptr) {
        input_script_eventrecon_n_events_tnamed->Write();
    }
    if (input_script_eventrecon_output_folder_tnamed != nullptr) {
        input_script_eventrecon_output_folder_tnamed->Write();
    }

    // -- Write the running info from the input file
    if (input_counter_header_line_tnamed != nullptr) {
        input_counter_header_line_tnamed->Write();
    }
    if (input_counter_heartbeat_packet_tnamed != nullptr) {
        input_counter_heartbeat_packet_tnamed->Write();
    }
    if (input_counter_complete_20_lines_tnamed != nullptr) {
        input_counter_complete_20_lines_tnamed->Write();
    }
    if (input_counter_not_complete_20_lines_tnamed != nullptr) {
        input_counter_not_complete_20_lines_tnamed->Write();
    }
    if (input_counter_not_complete_20_lines_total_tnamed != nullptr) {
        input_counter_not_complete_20_lines_total_tnamed->Write();
    }
    if (input_counter_invalid_line_tnamed != nullptr) {
        input_counter_invalid_line_tnamed->Write();
    }
    if (input_counter_valid_line_tnamed != nullptr) {
        input_counter_valid_line_tnamed->Write();
    }
    if (input_legal_fpga_id_list_tnamed != nullptr) {
        input_legal_fpga_id_list_tnamed->Write();
    }

    if (input_eventrecon_good_machine_gun_events_tnamed != nullptr) {
        input_eventrecon_good_machine_gun_events_tnamed->Write();
    }
    if (input_eventrecon_bad_machine_gun_events_tnamed != nullptr) {
        input_eventrecon_bad_machine_gun_events_tnamed->Write();
    }
    if (input_eventrecon_machine_gun_samples_tnamed != nullptr) {
        input_eventrecon_machine_gun_samples_tnamed->Write();
    }


    // -- Write the meta data
    TNamed("EventMatch_script_name", script_name.c_str()).Write();
    TNamed("EventMatch_script_version", script_version.c_str()).Write();
    TNamed("EventMatch_script_input_file", script_input_file.c_str()).Write();
    TNamed("EventMatch_script_output_file", script_output_file.c_str()).Write();
    TNamed("EventMatch_script_n_events", std::to_string(script_n_events).c_str()).Write();
    TNamed("EventMatch_script_output_folder", script_output_folder.c_str()).Write();
    // -- Write the running info
    TNamed("EventMatch_total_matched_events", std::to_string(total_matched_events).c_str()).Write();
    
    for (int i = 0; i < fpga_number; i++) {
        delete[] branch_timestamps_list[i];
        delete[] branch_daqh_list_list[i];
        delete[] branch_tc_list_list[i];
        delete[] branch_tp_list_list[i];
        delete[] branch_val0_list_list[i];
        delete[] branch_val1_list_list[i];
        delete[] branch_val2_list_list[i];
        delete[] branch_crc32_list_list[i];
        delete[] branch_last_heartbeat_list[i];
    }

    for (int i = 0; i < fpga_number; i++) {
        for (auto ptr : timestamp_pools_original[i]) {
            delete[] ptr;
        }
        for (auto ptr : daqh_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : tc_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : tp_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : val0_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : val1_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : val2_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : crc32_list_pools[i]) {
            delete[] ptr;
        }
        for (auto ptr : last_heartbeat_pools[i]) {
            delete[] ptr;
        }
    }

    output_root->Close();
    return 0;
}

void set_easylogger(){
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s}[%levshort] (%fbase) %msg");
    defaultConf.set(el::Level::Info,    el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;34m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;33m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Error,   el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;31m%levshort\033[0m] (%fbase) %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

LCSResult longestMatchingSequenceRec(std::vector<Long64_t> arrays[], int n, int pos0, const vector<int>& pos, int variance, unordered_map<string, LCSResult>& memo) {
    // Build a key for memoization from the current state.
    string key = makeKey(pos0, pos);
    if (memo.find(key) != memo.end())
        return memo[key];
    
    // Initialize best result as an empty sequence.
    LCSResult best;
    best.sequence = {};
    best.indices = vector<vector<int>>(n, vector<int>());
    
    // Iterate over the reference array (arrays[0]) starting from pos0.
    int m0 = arrays[0].size();
    for (int i = pos0; i < m0; i++) {
        Long64_t candidate = arrays[0][i];
        vector<int> newPos(pos.size());
        bool valid = true;
        // For each other array, try to find the next element (starting from pos[j]) that matches candidate within tolerance.
        for (int j = 0; j < static_cast<int>(pos.size()); j++) {
            int start = pos[j];
            const auto& arr = arrays[j + 1]; // pos[0] corresponds to arrays[1], etc.
            int found = -1;
            for (int k = start; k < static_cast<int>(arr.size()); k++) {
                if (std::abs(arr[k] - candidate) <= variance) {
                    found = k;
                    break;
                }
            }
            if (found == -1) {
                valid = false;
                break; // Candidate not found in one of the arrays.
            }
            // Save the next starting position (one past the matching index) for the j+1-th array.
            newPos[j] = found + 1;
        }
        if (!valid)
            continue;
        
        // Log that the candidate was found in all arrays.
        // LOG(INFO) << "Candidate " << candidate << " found in all arrays at reference index " << i;
        
        // Recursively search for the rest of the subsequence.
        LCSResult subRes = longestMatchingSequenceRec(arrays, n, i + 1, newPos, variance, memo);
        
        // Build the current candidate result: candidate followed by the subsequence found.
        LCSResult current;
        current.sequence.push_back(candidate);
        current.sequence.insert(current.sequence.end(), subRes.sequence.begin(), subRes.sequence.end());
        
        // For indices, record the index from the reference array...
        current.indices = vector<vector<int>>(n, vector<int>());
        current.indices[0].push_back(i);
        // ...and for each of the other arrays, record the index where candidate was found.
        for (int j = 0; j < static_cast<int>(pos.size()); j++) {
            current.indices[j + 1].push_back(newPos[j] - 1);
        }
        // Append the indices from the recursive call.
        for (int r = 0; r < n; r++) {
            current.indices[r].insert(current.indices[r].end(), subRes.indices[r].begin(), subRes.indices[r].end());
        }
        
        // Update best result if the current sequence is longer.
        if (current.sequence.size() > best.sequence.size()) {
            best = current;
        }
    }
    
    // Cache and return the best result from the current state.
    memo[key] = best;
    return best;
}

// Main function to compute the longest common (approximate) subsequence.
// - arrays: a C-style array of vectors (each vector<Long64_t>), where arrays[0] is the reference.
// - n: number of arrays in the input (_n)
// - variance: allowed tolerance (0-100) for matching elements.
LCSResult longestMatchingSequence(std::vector<Long64_t> arrays[], int n, int variance) {
    if (n < 2) {
        LCSResult result;
        if (n == 1) {
            result.sequence = arrays[0];
            result.indices = vector<vector<int>>(1, vector<int>());
            for (int i = 0; i < static_cast<int>(arrays[0].size()); i++)
                result.indices[0].push_back(i);
        }
        return result;
    }
    // For arrays[1]..arrays[n-1], the initial search starts at index 0.
    vector<int> initPos(n - 1, 0);
    unordered_map<string, LCSResult> memo;
    return longestMatchingSequenceRec(arrays, n, 0, initPos, variance, memo);
}