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

#define FPGA_CHANNEL_NUMBER 152
#define HEARTBEAT_HEADER_START 0x23
#define HEARTBEAT_HEADER_END   0x23
#define DATA_PACKET_HEADER_START 0x23
#define DATA_PACKET_HEADER_END   0x23
#define LINE_NUMBER_PER_FPGA  20

INITIALIZE_EASYLOGGINGPP

void set_easylogger();

int main(int argc, char **argv){

    // * --- Initial settings -----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string script_input_file, script_output_file;
    int script_n_events;
    std::string script_name = __FILE__;
    std::string script_version = "0.1";
    std::string script_output_folder;

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    argparse::ArgumentParser program(script_name, script_version);

    program.add_argument("-f", "--file").help("Input .h2g file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    try {
        program.parse_args(argc, argv);
        script_input_file  = program.get<std::string>("--file");
        script_output_file = program.get<std::string>("--output");
        auto script_n_events_str = program.get<std::string>("--events");
        script_n_events    = std::stoi(script_n_events_str);
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

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << script_output_file;
        return 1;
    }
    TTree* output_tree = new TTree("data_tree", "Data tree");
    output_tree->SetDirectory(output_root);

    // create branches
    // fpga_id is a 16-bit unsigned integer
    // timestamp is a 64-bit unsigned integer
    // daqh_list is 4 32-bit unsigned integers
    // tc and tp are boolen values
    auto *branch_fpga_id    = new UShort_t();   // 0 - 8
    auto *branch_timestamp  = new ULong64_t();  // 64 bits
    auto *branch_daqh_list  = new UInt_t[4];    // 32 bits
    auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER];
    auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_crc32_list = new UInt_t[4];
    auto *branch_last_heartbeat = new UInt_t();

    output_tree->Branch("fpga_id", branch_fpga_id, "fpga_id/s");
    output_tree->Branch("timestamp", branch_timestamp, "timestamp/l");
    output_tree->Branch("daqh_list", branch_daqh_list, "daqh_list[4]/i");
    output_tree->Branch("tc_list", branch_tc_list, ("tc_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    output_tree->Branch("tp_list", branch_tp_list, ("tp_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    output_tree->Branch("val0_list", branch_val0_list, ("val0_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    output_tree->Branch("val1_list", branch_val1_list, ("val1_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    output_tree->Branch("val2_list", branch_val2_list, ("val2_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    output_tree->Branch("crc32_list", branch_crc32_list, "crc32_list[4]/i");
    output_tree->Branch("last_heartbeat", branch_last_heartbeat, "last_heartbeat/i");

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int legal_line_id_list[] = {0x00, 0x01, 0x02, 0x03, 0x04};
    int legal_asic_id_list[] = {0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7};
    int legal_half_id_list[] = {0x24, 0x25};
    std::vector <int> legal_fpga_id_list;
    legal_fpga_id_list.reserve(32);

    const int input_file_read_chunk_size = 4096; // in bytes
    const int input_file_bytes_thershold = 1452; // in bytes
    const int SWTFC_pool_size_base = 200;
    const int SWTFC_window_size_base = 150;
    std::string input_file_buffer;
    char input_file_chunk[input_file_read_chunk_size];
    std::ifstream input_file(script_input_file, std::ios::binary);
    if (!input_file.is_open()) {
        LOG(ERROR) << "Failed to open input file " << script_input_file;
        return 1;
    }

    const int line_pool_size = 2000;
    std::vector <std::string> valid_line_pool;
    valid_line_pool.reserve(line_pool_size);
    std::vector <ULong64_t> timestamp_pool;
    timestamp_pool.reserve(line_pool_size);
    std::vector <UShort_t> fpga_id_pool;
    fpga_id_pool.reserve(line_pool_size);
    std::vector <bool> line_matched_flags;
    line_matched_flags.reserve(line_pool_size);

    bool flag_continue_reading = true;
    bool input_file_header_found = false;
    bool input_file_header_end_found = false;
    bool input_file_is_last_chunk = false;

    int counter_heartbeat_packet = 0;
    int counter_header_line = 0;
    int counter_complete_20_lines = 0;
    int counter_not_complete_20_lines = 0;
    int counter_not_complete_20_lines_total = 0;
    int counter_invalid_line = 0;
    int counter_valid_line = 0;

    UInt_t current_hear_beat_value = 0;

    while ((input_file.read(input_file_chunk, input_file_read_chunk_size) || input_file.gcount() > 0) && flag_continue_reading) {
        input_file_buffer.append(input_file_chunk, input_file.gcount());
        input_file_is_last_chunk = input_file.eof(); // Check if it is the last chunk
        if (input_file_is_last_chunk) {
            LOG(INFO) << "Last chunk found!";
        }
        while ((input_file_buffer.size() > input_file_bytes_thershold || (input_file_is_last_chunk && input_file_buffer.size() >= input_file_bytes_thershold)) && flag_continue_reading) {
            if (input_file_buffer.at(0) == DATA_PACKET_HEADER_START && !(input_file_header_found && input_file_header_end_found)) {
                if (!input_file_header_found) {
                    if (input_file_buffer.at(10) == DATA_PACKET_HEADER_END) {
                        input_file_header_found = true;
                        LOG(INFO) << "Header start found!";
                        size_t pos = input_file_buffer.find('\n');
                        if (pos != std::string::npos) {
                            input_file_buffer.erase(0, pos + 1);
                        }
                        counter_header_line++;
                        continue;
                    }
                }
                if (input_file_header_found && !input_file_header_end_found) {
                    if (input_file_buffer.at(10) == DATA_PACKET_HEADER_END) {
                        input_file_header_end_found = true;
                        LOG(INFO) << "Header end found!";
                        size_t pos = input_file_buffer.find('\n');
                        if (pos != std::string::npos) {
                            input_file_buffer.erase(0, pos + 1);
                        }
                        counter_header_line++;
                        continue;
                    } else {
                        size_t pos = input_file_buffer.find('\n');
                        if (pos == std::string::npos) {
                            break;
                        }
                        counter_header_line++;
                        std::string line = input_file_buffer.substr(0, pos);
                        input_file_buffer.erase(0, pos + 1); // Remove processed line
                        // print out the line in
                        LOG(INFO) << line;
                    }
                }
            } else {
                if (input_file_header_found && input_file_header_end_found) {
                    if (input_file_buffer.at(0) == HEARTBEAT_HEADER_START && input_file_buffer.at(10) == HEARTBEAT_HEADER_END) {
                        input_file_buffer.erase(0, 1452); // Adjust the size to remove the processed frame
                        counter_heartbeat_packet++;
                    } else {
                         // ! -- data packet processing --------------------------------
                        // remove the first 12
                        int _index_base = 0;
                        int _index_end = 1452;
                        if (input_file_is_last_chunk) {
                            if (input_file_buffer.size() <= 1452) {
                                _index_end = input_file_buffer.size() - 40;
                            }
                        }
                        while (_index_base < _index_end) {
                            try {
                                unsigned char _asic_id = input_file_buffer.at(_index_base);
                                if (std::find(std::begin(legal_asic_id_list), std::end(legal_asic_id_list), _asic_id) == std::end(legal_asic_id_list)) {
                                    throw std::runtime_error("Bad ASIC ID: " + std::to_string((int)_asic_id));
                                }
                                unsigned char _half_id = input_file_buffer.at(_index_base + 2);
                                if (std::find(std::begin(legal_half_id_list), std::end(legal_half_id_list), _half_id) == std::end(legal_half_id_list)) {
                                    throw std::runtime_error("Bad HALF ID: " + std::to_string((int)_half_id));
                                }
                                unsigned char _line_id = input_file_buffer.at(_index_base + 3);
                                if (std::find(std::begin(legal_line_id_list), std::end(legal_line_id_list), _line_id) == std::end(legal_line_id_list)) {
                                    throw std::runtime_error("Bad LINE ID: " + std::to_string((int)_line_id));
                                }
                                unsigned char _fpga_id = input_file_buffer.at(_index_base + 1);
                                if (_fpga_id > 31) {
                                    throw std::runtime_error("Bad FPGA ID: " + std::to_string((int)_fpga_id));
                                }
                            } catch (const std::runtime_error& e) {
                                auto _start = _index_base - 40;
                                if (_start < 0) {
                                    _start = 0;
                                }
                                auto _end = _index_base + 40;
                                if (_end > 1452) {
                                    _end = 1452;
                                }
                                _index_base += 12;
                                counter_invalid_line++;
                                continue;
                            }

                            unsigned char _fpga_id = (unsigned char) input_file_buffer.at(_index_base + 1);
                            if (std::find(legal_fpga_id_list.begin(), legal_fpga_id_list.end(), _fpga_id) == legal_fpga_id_list.end()) {
                                    legal_fpga_id_list.push_back(_fpga_id);
                            }
                            counter_valid_line++;
                            auto _new_line = input_file_buffer.substr(_index_base, 40);
                            if (_new_line.size() != 40) {
                                LOG(ERROR) << "Invalid line size: " << _new_line.size();
                                break;
                            }
                            valid_line_pool.push_back(_new_line);
                            ULong64_t _timestamp = (unsigned char) input_file_buffer.at(_index_base + 4) << 24 | (unsigned char) input_file_buffer.at(_index_base + 5) << 16 | (unsigned char) input_file_buffer.at(_index_base + 6) << 8 | (unsigned char) input_file_buffer.at(_index_base + 7);
                            timestamp_pool.push_back(_timestamp);
                            fpga_id_pool.push_back(_fpga_id);
                            line_matched_flags.push_back(false);
                            _index_base += 40;
                        }
                            
                        int target_SWTFC_pool_size = SWTFC_pool_size_base * legal_fpga_id_list.size();
                        int _SWTFC_window_size = SWTFC_window_size_base * legal_fpga_id_list.size();
                        if (target_SWTFC_pool_size > line_pool_size) {
                            LOG(WARNING) << "Target pool size is too large: " << target_SWTFC_pool_size;
                            target_SWTFC_pool_size = line_pool_size;
                        }

                        if (valid_line_pool.size() >= target_SWTFC_pool_size || input_file_is_last_chunk) {
                            for (int i = 0; i < valid_line_pool.size() - _SWTFC_window_size; i++) {
                                // LOG(DEBUG) << "Checking line: " << i;
                                if (line_matched_flags.at(i)) {
                                    continue;
                                }
                                auto _timestamp_ref = timestamp_pool.at(i);
                                auto _fpga_id_ref = fpga_id_pool.at(i);
                                int _match_line_counter = 1;
                                std::vector <int> _match_line_index;
                                _match_line_index.push_back(i);
                                line_matched_flags.at(i) = true;
                                _match_line_index.reserve(_SWTFC_window_size);
                                int _match_searching_end = valid_line_pool.size() - 1;
                                if (i + _SWTFC_window_size < _match_searching_end) {
                                    _match_searching_end = i + _SWTFC_window_size;
                                }
                                for (int j = i + 1; j < _match_searching_end; j++) {
                                    if (line_matched_flags.at(j)) {
                                        continue;
                                    }
                                    auto _timestamp = timestamp_pool.at(j);
                                    auto _fpga_id = fpga_id_pool.at(j);
                                    if (_timestamp == _timestamp_ref && _fpga_id == _fpga_id_ref) {
                                        _match_line_counter++;
                                        _match_line_index.push_back(j);
                                        line_matched_flags.at(j) = true;
                                    }
                                }
                                if (_match_line_counter == LINE_NUMBER_PER_FPGA){
                                    counter_complete_20_lines += LINE_NUMBER_PER_FPGA;
                                    if (counter_complete_20_lines / 20 >= script_n_events && script_n_events > 0) {
                                        flag_continue_reading = false;
                                        break;
                                    }
                                    for (auto &index : _match_line_index) {
                                        // LOG(INFO) << "Matched line: " << index;
                                        auto _line_id = (unsigned char) valid_line_pool.at(index).at(3);
                                        auto _asic_id = (unsigned char) valid_line_pool.at(index).at(0);
                                        auto _half_id = (unsigned char) valid_line_pool.at(index).at(2);
                                        // LOG(INFO) << "Matched line: " << index << " ASIC ID: " << std::hex << (int)_asic_id << " HALF ID: " << (int)_half_id << " LINE ID: " << (int)_line_id << std::dec;
                                        int _half_index = (_asic_id - 0xa0) * 2 + (_half_id - 0x24);
                                        int _half_channel_offset = _half_index * 38;
                                        switch (_line_id) {
                                            case 0x00:
                                            {
                                                if (_asic_id == 0xa0 && _half_id == 0x24) {
                                                    *branch_fpga_id = _fpga_id_ref;
                                                    *branch_timestamp = _timestamp_ref;
                                                    *branch_last_heartbeat = current_hear_beat_value;
                                                }
                                                UInt_t _daqh = (unsigned char) valid_line_pool.at(index).at(8) << 24 | (unsigned char) valid_line_pool.at(index).at(9) << 16 | (unsigned char) valid_line_pool.at(index).at(10) << 8 | (unsigned char) valid_line_pool.at(index).at(11);
                                                branch_daqh_list[_half_index] = _daqh;
                                                for (int _word_index = 1; _word_index < 8; _word_index++) {
                                                    ULong64_t _word = (unsigned char) valid_line_pool.at(index).at(8 + 4 * _word_index) << 24 | (unsigned char) valid_line_pool.at(index).at(9 + 4 * _word_index) << 16 | (unsigned char) valid_line_pool.at(index).at(10 + 4 * _word_index) << 8 | (unsigned char) valid_line_pool.at(index).at(11 + 4 * _word_index);
                                                    Bool_t _tc = (_word >> 31) & 0x1;
                                                    Bool_t _tp = (_word >> 30) & 0x1;
                                                    UInt_t _val0 = (_word >> 20) & 0x3ff;
                                                    UInt_t _val1 = (_word >> 10) & 0x3ff;
                                                    UInt_t _val2 = _word & 0x3ff;
                                                    branch_tc_list[_word_index - 1 + _half_channel_offset] = _tc;
                                                    branch_tp_list[_word_index - 1 + _half_channel_offset] = _tp;
                                                    branch_val0_list[_word_index - 1 + _half_channel_offset] = _val0;
                                                    branch_val1_list[_word_index - 1 + _half_channel_offset] = _val1;
                                                    branch_val2_list[_word_index - 1 + _half_channel_offset] = _val2;
                                                }
                                                break;
                                            }
                                            case 0x01:
                                            {
                                                for (int _word_index = 0; _word_index < 8; _word_index++) {
                                                    ULong64_t _word = (unsigned char) valid_line_pool.at(index).at(8 + 4 * _word_index) << 24 | (unsigned char) valid_line_pool.at(index).at(9 + 4 * _word_index) << 16 | (unsigned char) valid_line_pool.at(index).at(10 + 4 * _word_index) << 8 | (unsigned char) valid_line_pool.at(index).at(11 + 4 * _word_index);
                                                    Bool_t _tc = (_word >> 31) & 0x1;
                                                    Bool_t _tp = (_word >> 30) & 0x1;
                                                    UInt_t _val0 = (_word >> 20) & 0x3ff;
                                                    UInt_t _val1 = (_word >> 10) & 0x3ff;
                                                    UInt_t _val2 = _word & 0x3ff;
                                                    branch_tc_list[_word_index + 7 + _half_channel_offset] = _tc;
                                                    branch_tp_list[_word_index + 7 + _half_channel_offset] = _tp;
                                                    branch_val0_list[_word_index + 7 + _half_channel_offset] = _val0;
                                                    branch_val1_list[_word_index + 7 + _half_channel_offset] = _val1;
                                                    branch_val2_list[_word_index + 7 + _half_channel_offset] = _val2;
                                                }
                                                break;
                                            }
                                            case 0x02:
                                            {
                                                for (int _word_index = 0; _word_index < 8; _word_index++) {
                                                    ULong64_t _word = (unsigned char) valid_line_pool.at(index).at(8 + 4 * _word_index) << 24 | (unsigned char) valid_line_pool.at(index).at(9 + 4 * _word_index) << 16 | (unsigned char) valid_line_pool.at(index).at(10 + 4 * _word_index) << 8 | (unsigned char) valid_line_pool.at(index).at(11 + 4 * _word_index);
                                                    Bool_t _tc = (_word >> 31) & 0x1;
                                                    Bool_t _tp = (_word >> 30) & 0x1;
                                                    UInt_t _val0 = (_word >> 20) & 0x3ff;
                                                    UInt_t _val1 = (_word >> 10) & 0x3ff;
                                                    UInt_t _val2 = _word & 0x3ff;
                                                    branch_tc_list[_word_index + 15 + _half_channel_offset] = _tc;
                                                    branch_tp_list[_word_index + 15 + _half_channel_offset] = _tp;
                                                    branch_val0_list[_word_index + 15 + _half_channel_offset] = _val0;
                                                    branch_val1_list[_word_index + 15 + _half_channel_offset] = _val1;
                                                    branch_val2_list[_word_index + 15 + _half_channel_offset] = _val2;
                                                }
                                                break;
                                            }
                                            case 0x03:
                                            {
                                                for (int _word_index = 0; _word_index < 8; _word_index++) {
                                                    ULong64_t _word = (unsigned char) valid_line_pool.at(index).at(8 + 4 * _word_index) << 24 | (unsigned char) valid_line_pool.at(index).at(9 + 4 * _word_index) << 16 | (unsigned char) valid_line_pool.at(index).at(10 + 4 * _word_index) << 8 | (unsigned char) valid_line_pool.at(index).at(11 + 4 * _word_index);
                                                    Bool_t _tc = (_word >> 31) & 0x1;
                                                    Bool_t _tp = (_word >> 30) & 0x1;
                                                    UInt_t _val0 = (_word >> 20) & 0x3ff;
                                                    UInt_t _val1 = (_word >> 10) & 0x3ff;
                                                    UInt_t _val2 = _word & 0x3ff;
                                                    branch_tc_list[_word_index + 23 + _half_channel_offset] = _tc;
                                                    branch_tp_list[_word_index + 23 + _half_channel_offset] = _tp;
                                                    branch_val0_list[_word_index + 23 + _half_channel_offset] = _val0;
                                                    branch_val1_list[_word_index + 23 + _half_channel_offset] = _val1;
                                                    branch_val2_list[_word_index + 23 + _half_channel_offset] = _val2;
                                                }
                                                break;
                                            }
                                            case 0x04:
                                            {
                                                for (int _word_index = 0; _word_index < 7; _word_index++) {
                                                    // LOG(DEBUG) << "Word index: " << _word_index;
                                                    ULong64_t _word = (unsigned char) valid_line_pool.at(index).at(8 + 4 * _word_index) << 24 | (unsigned char) valid_line_pool.at(index).at(9 + 4 * _word_index) << 16 | (unsigned char) valid_line_pool.at(index).at(10 + 4 * _word_index) << 8 | (unsigned char) valid_line_pool.at(index).at(11 + 4 * _word_index);
                                                    Bool_t _tc = (_word >> 31) & 0x1;
                                                    Bool_t _tp = (_word >> 30) & 0x1;
                                                    UInt_t _val0 = (_word >> 20) & 0x3ff;
                                                    UInt_t _val1 = (_word >> 10) & 0x3ff;
                                                    UInt_t _val2 = _word & 0x3ff;
                                                    branch_tc_list[_word_index + 31 + _half_channel_offset] = _tc;
                                                    branch_tp_list[_word_index + 31 + _half_channel_offset] = _tp;
                                                    branch_val0_list[_word_index + 31 + _half_channel_offset] = _val0;
                                                    branch_val1_list[_word_index + 31 + _half_channel_offset] = _val1;
                                                    branch_val2_list[_word_index + 31 + _half_channel_offset] = _val2;
                                                }
                                                UInt_t _crc32 = (unsigned char) valid_line_pool.at(index).at(8 + 4 * 7) << 24 | (unsigned char) valid_line_pool.at(index).at(9 + 4 * 7) << 16 | (unsigned char) valid_line_pool.at(index).at(10 + 4 * 7) << 8 | (unsigned char) valid_line_pool.at(index).at(11 + 4 * 7);
                                                branch_crc32_list[_half_index] = _crc32;
                                                break;
                                            }
                                            default:
                                                break;
                                        }
                                    }
                                    output_tree->Fill();
                                } else {
                                    counter_not_complete_20_lines += _match_line_counter;
                                    counter_not_complete_20_lines_total += 1;
                                }
                            }
                            valid_line_pool.erase(valid_line_pool.begin(), valid_line_pool.end() - _SWTFC_window_size);
                            timestamp_pool.erase(timestamp_pool.begin(), timestamp_pool.end() - _SWTFC_window_size);
                            fpga_id_pool.erase(fpga_id_pool.begin(), fpga_id_pool.end() - _SWTFC_window_size);
                            line_matched_flags.erase(line_matched_flags.begin(), line_matched_flags.end() - _SWTFC_window_size);
                        }
                        input_file_buffer.erase(0, 1452);
                    }
                }
            }
        }
    }

    LOG(INFO) << "Total header lines: " << counter_header_line;
    LOG(INFO) << "Total heartbeat packets: " << counter_heartbeat_packet;
    if (legal_fpga_id_list.size() == 0) {
        LOG(ERROR) << "No legal FPGA ID found!";
    } else {
        LOG(INFO) << "Total complete 20 lines: " << counter_complete_20_lines << " (" << counter_complete_20_lines / legal_fpga_id_list.size() / 20 << " samples per FPGA)";
        LOG(INFO) << "Total not complete 20 lines: " << counter_not_complete_20_lines << " (" << counter_not_complete_20_lines_total / legal_fpga_id_list.size() << " samples per FPGA)";
    }
    if ((counter_complete_20_lines + counter_not_complete_20_lines) == 0) {
        LOG(ERROR) << "No valid line found!";
    } else {
        LOG(INFO) << "Complete 20 lines percentage: " << (float)counter_complete_20_lines / (counter_complete_20_lines + counter_not_complete_20_lines) * 100 << "%";
        LOG(INFO) << "Total valid line percentage: " << (float)counter_valid_line / (counter_valid_line + counter_invalid_line / 40) * 100 << "%";
    }
    std::string _legal_fpga_id_list_str = "";
    for (auto &fpga_id : legal_fpga_id_list) {
        _legal_fpga_id_list_str += std::to_string(fpga_id) + " ";
    }
    LOG(INFO) << "Legal FPGA ID list: " << _legal_fpga_id_list_str;

    input_file.close();

    // * --- Write output file ----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    // -- Write the tree
    output_tree->Write();
    // -- Write the meta data
    TNamed("Rootifier_script_name", script_name.c_str()).Write();
    TNamed("Rootifier_script_version", script_version.c_str()).Write();
    TNamed("Rootifier_script_input_file", script_input_file.c_str()).Write();
    TNamed("Rootifier_script_output_file", script_output_file.c_str()).Write();
    TNamed("Rootifier_script_n_events", std::to_string(script_n_events).c_str()).Write();
    TNamed("Rootifier_script_output_folder", script_output_folder.c_str()).Write();
    // -- Write the runnning info
    TNamed("Rootifier_counter_header_line", std::to_string(counter_header_line).c_str()).Write();
    TNamed("Rootifier_counter_heartbeat_packet", std::to_string(counter_heartbeat_packet).c_str()).Write();
    TNamed("Rootifier_counter_complete_20_lines", std::to_string(counter_complete_20_lines).c_str()).Write();
    TNamed("Rootifier_counter_not_complete_20_lines", std::to_string(counter_not_complete_20_lines).c_str()).Write();
    TNamed("Rootifier_counter_not_complete_20_lines_total", std::to_string(counter_not_complete_20_lines_total).c_str()).Write();
    TNamed("Rootifier_counter_invalid_line", std::to_string(counter_invalid_line).c_str()).Write();
    TNamed("Rootifier_counter_valid_line", std::to_string(counter_valid_line).c_str()).Write();
    // -- Write the legal fpga id list
    TNamed("Rootifier_legal_fpga_id_list", _legal_fpga_id_list_str.c_str()).Write();
    
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