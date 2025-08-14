#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_root(argc, argv, "1.0");

    std::vector<Color_t> fpga_colors = {
        kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kPink, kViolet
    };

    gROOT->SetBatch(kTRUE);
    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 1.5625; // unit: ns
    // get run number from the input file name
    size_t run_number_pos = opts.input_file.find("Run");
    std::string run_info_str = "Run ";
    if (run_number_pos != std::string::npos) {
        run_info_str += opts.input_file.substr(run_number_pos + 3, 3);
    } else {
        run_info_str += "Unknown";
    }

    // * --- Read the root file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int fpga_count = -1;
    int machine_gun_samples = -1;
    int entry_max = -1;

    TFile *input_root;
    TTree *input_tree;
    std::vector<UShort_t> legal_fpga_id_list;

    if (!readRootMetaData(opts.input_file.c_str(), input_root, input_tree, fpga_count, machine_gun_samples, entry_max, legal_fpga_id_list)) {
        LOG(ERROR) << "Failed to read metadata from input file " << opts.input_file;
        return 1;
    }

    if (entry_max == 0) {
        LOG(ERROR) << "No events in the input file!";
        return 1;
    }
    if (opts.n_events > 0 && opts.n_events < entry_max) {
        entry_max = opts.n_events;
    } else {
        if (opts.n_events > entry_max) {
            LOG(WARNING) << "Requested number of events " << opts.n_events << " is larger than the number of events in the input file " << entry_max;
        }
    }

    // * --- Create the output file -----------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    std::unordered_map <int, int> unifiedToHistIndex;
    std::unordered_map <int, int> histIndexToUnified;

    int channel_array_index = 0;
    int sample_time_bins = int((machine_gun_samples * sample_time) / phase_shift_time);
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
            auto _channel_valid = get_valid_fpga_channel(_channel_index);
            if (_channel_valid == -1){
                continue;
            }
            auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
            unifiedToHistIndex[_unified_valid_channel_number] = channel_array_index;
            histIndexToUnified[channel_array_index] = _unified_valid_channel_number;

            channel_array_index++;
        } // end of channel loop
    } // end of fpga loop
   
    LOG(INFO) << "FPGA count: " << fpga_count << " Machine gun samples: " << machine_gun_samples << " Entry max: " << entry_max;

    std::vector <ULong64_t*> branch_timestamps_list;
    std::vector <UInt_t*> branch_daqh_list_list;
    std::vector <Bool_t*> branch_tc_list_list;
    std::vector <Bool_t*> branch_tp_list_list;
    std::vector <UInt_t*> branch_val0_list_list;
    std::vector <UInt_t*> branch_val1_list_list;
    std::vector <UInt_t*> branch_val2_list_list;
    std::vector <UInt_t*> branch_crc32_list_list;
    std::vector <UInt_t*> branch_last_heartbeat_list;

    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
        auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
        auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
        auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];

        input_tree->SetBranchAddress(("timestamps_"    + std::to_string(_fpga_id)).c_str(), branch_timestamps);
        input_tree->SetBranchAddress(("daqh_list_"     + std::to_string(_fpga_id)).c_str(), branch_daqh_list);
        input_tree->SetBranchAddress(("tc_list_"       + std::to_string(_fpga_id)).c_str(), branch_tc_list);
        input_tree->SetBranchAddress(("tp_list_"       + std::to_string(_fpga_id)).c_str(), branch_tp_list);
        input_tree->SetBranchAddress(("val0_list_"     + std::to_string(_fpga_id)).c_str(), branch_val0_list);
        input_tree->SetBranchAddress(("val1_list_"     + std::to_string(_fpga_id)).c_str(), branch_val1_list);
        input_tree->SetBranchAddress(("val2_list_"     + std::to_string(_fpga_id)).c_str(), branch_val2_list);
        input_tree->SetBranchAddress(("crc32_list_"    + std::to_string(_fpga_id)).c_str(), branch_crc32_list);
        input_tree->SetBranchAddress(("last_heartbeat_"+ std::to_string(_fpga_id)).c_str(), branch_last_heartbeat);

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

    long counter_total_events = 0;
    long counter_hamming_code_error = 0;
    long counter_bad_daqh_start_end = 0;

    long counter_double_tot = 0;
    long counter_double_toa = 0;
    long counter_valid_tot  = 0;
    long counter_valid_toa  = 0;

    // * --- Create data structures for plotting ----------------------------------------
    // * --------------------------------------------------------------------------------
    // std::vector <std::vector<double>> fpgas_adc_sum_list;
    // std::vector <std::vector<double>> fpgas_tot_sum_list;
    // std::vector <std::vector<double>> fpgas_tot_count_list;
    // std::vector <std::vector<double>> fpgas_toa_count_list;
    // std::vector <std::vector<double>> fpgas_hamming_code_error_rate;

    // std::vector <std::vector<double>> fpgas_hamming_code_error_count;
    // std::vector <std::vector<double>> fpgas_daqh_error_count;

    // double* fpgas_hamming_code_counter = new double[fpga_count];
    // for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //     fpgas_hamming_code_counter[_fpga_index] = 0.0;
    // }
    // fpgas_adc_sum_list.resize(fpga_count*2);
    // fpgas_tot_sum_list.resize(fpga_count*2);
    // fpgas_tot_count_list.resize(fpga_count*2);
    // fpgas_toa_count_list.resize(fpga_count*2);
    // fpgas_hamming_code_error_count.resize(fpga_count*2);
    // fpgas_hamming_code_error_rate.resize(fpga_count*2);
    // fpgas_daqh_error_count.resize(fpga_count*2);
    // for (int _fpga_index = 0; _fpga_index < fpga_count*2; _fpga_index++) {
    //     fpgas_adc_sum_list[_fpga_index].resize(entry_max, -1.0);
    //     fpgas_tot_sum_list[_fpga_index].resize(entry_max, -1.0);
    //     fpgas_tot_count_list[_fpga_index].resize(entry_max, -1.0);
    //     fpgas_toa_count_list[_fpga_index].resize(entry_max, -1.0);
    //     fpgas_hamming_code_error_count[_fpga_index].resize(entry_max, -1.0);
    //     fpgas_hamming_code_error_rate[_fpga_index].resize(entry_max, -1.0);
    //     fpgas_daqh_error_count[_fpga_index].resize(entry_max, -1.0);
    // }

    // std::vector <std::vector<double>> fpgas_timestamp_list;
    // fpgas_timestamp_list.resize(fpga_count);
    // for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //     fpgas_timestamp_list[_fpga_index].reserve(entry_max);
    // }

    // std::vector <std::vector<double>> fpgas_event_index_list;
    // fpgas_event_index_list.resize(fpga_count);
    // for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //     fpgas_event_index_list[_fpga_index].reserve(entry_max);
    // }

    std::vector <std::vector<double>> fpgas_adc_sum_list(fpga_count * 2, std::vector<double>(entry_max, -1.0));
    std::vector <std::vector<double>> fpgas_hamming_code_error_count(fpga_count * 2, std::vector<double>(entry_max, -1.0));
    std::vector <std::vector<double>> fpgas_daqh_error_count(fpga_count * 2, std::vector<double>(entry_max, -1.0));
    std::vector <std::vector<int>> fpgas_signal_shape_pass(fpga_count * 2, std::vector<int>(entry_max, 0));

    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        counter_total_events++;

        if (_entry % 5000 == 0) {
            LOG(INFO) << "Processing entry " << _entry << " / " << entry_max;
        }

        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id    = legal_fpga_id_list[_fpga_index];
            auto _timestamp  = branch_timestamps_list[_fpga_index][0];
            auto _daqh_list  = branch_daqh_list_list[_fpga_index];
            auto _tc_list    = branch_tc_list_list[_fpga_index];
            auto _tp_list    = branch_tp_list_list[_fpga_index];
            auto _val0_list  = branch_val0_list_list[_fpga_index];
            auto _val1_list  = branch_val1_list_list[_fpga_index];
            auto _val2_list  = branch_val2_list_list[_fpga_index];

            double _fpga_event_adc_sum_a0 = -1.0;
            double _fpga_event_adc_sum_a1 = -1.0;
            double _fpga_event_tot_sum_a0 = -1.0;
            double _fpga_event_tot_sum_a1 = -1.0;
            double _fpga_event_tot_count_a0 = -1.0;
            double _fpga_event_tot_count_a1 = -1.0;
            double _fpga_event_toa_count_a0 = -1.0;
            double _fpga_event_toa_count_a1 = -1.0;

            double _hamming_code_error_count_a0 = 0.0;
            double _hamming_code_error_count_a1 = 0.0;

            double _daqh_error_count_a0 = 0.0;
            double _daqh_error_count_a1 = 0.0;

            // -- Check Hamming code and header bytes --
            // -----------------------------------------
            std::vector <bool> _hamming_code_pass_list;
            bool skip_event_hamming_code = false;
            bool good_header_bytes = true;
            for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                bool _hamming_code_pass = true;
                for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
                    auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
                    auto _h1h2h3 = (_daqh >> 4) & 0x7;
                    if (_h1h2h3 != 0x00){
                        _hamming_code_pass = false;
                        if (_daqh_index < 2) {
                            _hamming_code_error_count_a0++;
                        } else {
                            _hamming_code_error_count_a1++;
                        }
                    }
                    auto _daqh_first_half_byte = (_daqh >> 28) & 0xf;
                    auto _daqh_last_half_byte  = _daqh & 0xf;
                    if (_daqh_first_half_byte != 5 ||  _daqh_last_half_byte != 5) {
                        good_header_bytes = false;
                        if (_daqh_index < 2) {
                            _daqh_error_count_a0++;
                        } else {
                            _daqh_error_count_a1++;
                        }
                    }
                }
                _hamming_code_pass_list.push_back(_hamming_code_pass);
                if (!_hamming_code_pass){
                    skip_event_hamming_code = true;
                    break;
                }
            }

            fpgas_hamming_code_error_count[_fpga_index * 2][_entry] = _hamming_code_error_count_a0;
            fpgas_hamming_code_error_count[_fpga_index * 2 + 1][_entry] = _hamming_code_error_count_a1;

            fpgas_daqh_error_count[_fpga_index * 2][_entry] = _daqh_error_count_a0;
            fpgas_daqh_error_count[_fpga_index * 2 + 1][_entry] = _daqh_error_count_a1;

            // fpgas_timestamp_list[_fpga_index].push_back(_timestamp);
            // fpgas_event_index_list[_fpga_index].push_back(counter_total_events);

            // if (skip_event_hamming_code){
            //     counter_hamming_code_error++;
            //     fpgas_hamming_code_counter[_fpga_index]++;
            //     for (int _fpga_index2 = 0; _fpga_index2 < fpga_count; _fpga_index2++) {
            //         fpgas_hamming_code_error_rate[_fpga_index2][_entry] = fpgas_hamming_code_counter[_fpga_index2];
            //     }
            //     continue;
            // }
            // if (!good_header_bytes){
            //     counter_bad_daqh_start_end++;
            //     for (int _fpga_index2 = 0; _fpga_index2 < fpga_count; _fpga_index2++) {
            //         fpgas_hamming_code_error_rate[_fpga_index2][_entry] = fpgas_hamming_code_counter[_fpga_index2];
            //     }
            //     continue;
            // }

            // fpgas_hamming_code_error_rate[_fpga_index][_entry] = double(fpgas_hamming_code_counter[_fpga_index]);

            std::vector <int> _channel_val0_max_list;
            std::vector <int> _channel_val1_max_list;
            std::vector <int> _channel_val2_max_list;
            std::vector <int> _channel_val0_max_index_list;
            std::vector <int> _channel_val1_max_index_list;
            std::vector <int> _channel_val2_max_index_list;

            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                auto _channel_valid = get_valid_fpga_channel(_channel_index);
                if (_channel_valid == -1){
                    
                    _channel_val0_max_list.push_back(-1);
                    _channel_val1_max_list.push_back(-1);
                    _channel_val2_max_list.push_back(-1);
                    _channel_val0_max_index_list.push_back(-1);
                    _channel_val1_max_index_list.push_back(-1);
                    _channel_val2_max_index_list.push_back(-1);
                    continue;
                }
                auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                // check if the channel is in the target channels
                // if (std::find(target_channels.begin(), target_channels.end(), _unified_valid_channel_number) == target_channels.end()) {
                //     _channel_val0_max_list.push_back(-1);
                //     _channel_val1_max_list.push_back(-1);
                //     _channel_val2_max_list.push_back(-1);
                //     _channel_val0_max_index_list.push_back(-1);
                //     _channel_val1_max_index_list.push_back(-1);
                //     _channel_val2_max_index_list.push_back(-1);
                //     continue;
                // }

                int _val0_max = -1;
                int _val0_max_index = -1;
                int _val1_max = -1;
                int _val1_max_index = -1;
                int _val2_max = -1;
                int _val2_max_index = -1;

                int _adc_pedestal = 0; // TODO: get the pedestal value from the metadata

                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    auto _val0 = int(_val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER]);
                    auto _val1 = int(_val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER]);
                    auto _val2 = int(_val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER]);
                    if (_sample_index == 0) {
                        _adc_pedestal = _val0; // assuming the first sample is the pedestal
                    }
                    // if (_sample_index == 3) {
                    //     if (_channel_index < 76) { // for the first ASIC
                    //         if (_fpga_event_adc_sum_a0 < 0) {
                    //             _fpga_event_adc_sum_a0 = 0;
                    //         }
                    //         _fpga_event_adc_sum_a0 += _val0 - _adc_pedestal;
                    //     } else {
                    //         if (_fpga_event_adc_sum_a1 < 0) {
                    //             _fpga_event_adc_sum_a1 = 0;
                    //         }
                    //         _fpga_event_adc_sum_a1 += _val0 - _adc_pedestal;
                    //     }
                    // }

                    if (_val0 > 0 && _val0 > _val0_max){
                        _val0_max = _val0;
                        _val0_max_index = _sample_index;
                    }
                    if (_val1 > 0){
                        if (_val1_max == -1) {
                            _val1_max = _val1;
                            _val1_max_index = _sample_index;
                            counter_valid_tot++;
                            if (_channel_index < 76) { // for the first ASIC
                                if (_fpga_event_tot_sum_a0 < 0) {
                                    _fpga_event_tot_sum_a0 = 0;
                                }
                                if (_fpga_event_tot_count_a0 < 0) {
                                    _fpga_event_tot_count_a0 = 0;
                                }
                            } else {
                                if (_fpga_event_tot_sum_a1 < 0) {
                                    _fpga_event_tot_sum_a1 = 0;
                                }
                                if (_fpga_event_tot_count_a1 < 0) {
                                    _fpga_event_tot_count_a1 = 0;
                                }
                            }
                            double _val1_decoded = _val1;
                            if (_val1_decoded > 512) {
                                _val1_decoded -= 512;
                                _val1_decoded *= 8.0;
                            }
                            if (_channel_index < 76) { // for the first ASIC
                                _fpga_event_tot_sum_a0 += _val1_decoded;
                                _fpga_event_tot_count_a0 += 1.0;
                            } else {
                                _fpga_event_tot_sum_a1 += _val1_decoded;
                                _fpga_event_tot_count_a1 += 1.0;
                            }
                            // _fpga_event_tot_sum += _val1_decoded;
                            // if (_fpga_event_tot_count < 0) {
                            //     _fpga_event_tot_count = 0;
                            // }
                            // _fpga_event_tot_count += 1.0;
                        } else
                            counter_double_tot++;
                    }
                    if (_val2 > 0) {
                        if (_val2_max == -1) {
                            _val2_max = _val2;
                            _val2_max_index = _sample_index;
                            counter_valid_toa++;
                            if (_channel_index < 76) { // for the first ASIC
                                if (_fpga_event_toa_count_a0 < 0) {
                                    _fpga_event_toa_count_a0 = 0;
                                }
                                _fpga_event_toa_count_a0 += 1.0;
                            } else {
                                if (_fpga_event_toa_count_a1 < 0) {
                                    _fpga_event_toa_count_a1 = 0;
                                }
                                _fpga_event_toa_count_a1 += 1.0;
                            }
                        } else
                            counter_double_toa++;
                    }
                } // end of sample loop
                _channel_val0_max_list.push_back(_val0_max);
                _channel_val1_max_list.push_back(_val1_max);
                _channel_val2_max_list.push_back(_val2_max);
                // LOG(DEBUG) << "toa code raw: " << _val2_max;
                _channel_val0_max_index_list.push_back(_val0_max_index);
                _channel_val1_max_index_list.push_back(_val1_max_index);
                _channel_val2_max_index_list.push_back(_val2_max_index);

                if (_val0_max_index < 2 || _val0_max_index > 15) {
                    fpgas_signal_shape_pass[_fpga_index * 2][_entry] += 1; // signal shape pass
                }

                if (_channel_index < 76) { // for the first ASIC
                    if (_fpga_event_adc_sum_a0 < 0) {
                        _fpga_event_adc_sum_a0 = 0;
                    }
                    _fpga_event_adc_sum_a0 += _val0_max - _adc_pedestal;
                } else {
                    if (_fpga_event_adc_sum_a1 < 0) {
                        _fpga_event_adc_sum_a1 = 0;
                    }
                    _fpga_event_adc_sum_a1 += _val0_max - _adc_pedestal;
                }

            } // end of channel loop
            fpgas_adc_sum_list[_fpga_index * 2][_entry] = _fpga_event_adc_sum_a0;
            fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] = _fpga_event_adc_sum_a1;
            // fpgas_tot_sum_list[_fpga_index * 2][_entry] = _fpga_event_tot_sum_a0;
            // fpgas_tot_sum_list[_fpga_index * 2 + 1][_entry] = _fpga_event_tot_sum_a1;
            // fpgas_tot_count_list[_fpga_index * 2][_entry] = _fpga_event_tot_count_a0;
            // fpgas_tot_count_list[_fpga_index * 2 + 1][_entry] = _fpga_event_tot_count_a1;
            // fpgas_toa_count_list[_fpga_index * 2][_entry] = _fpga_event_toa_count_a0;
            // fpgas_toa_count_list[_fpga_index * 2 + 1][_entry] = _fpga_event_toa_count_a1;
            // fpgas_adc_sum_list[_fpga_index][_entry] = _fpga_event_adc_sum;
            // fpgas_tot_sum_list[_fpga_index][_entry] = _fpga_event_tot_sum;
            // fpgas_tot_count_list[_fpga_index][_entry] = _fpga_event_tot_count;
            // fpgas_toa_count_list[_fpga_index][_entry] = _fpga_event_toa_count;
        } // end of fpga loop
    } // end of entry loop

    LOG(INFO) << "=== Data Reading Summary ===============================";
    LOG(INFO) << "Total events: " << counter_total_events;
    LOG(INFO) << "Hamming code errors: " << counter_hamming_code_error << " (" << (double(counter_hamming_code_error) / counter_total_events * 100.0) << "%)";
    LOG(INFO) << "Bad DAQH start/end: " << counter_bad_daqh_start_end << " (" << (double(counter_bad_daqh_start_end) / counter_total_events * 100.0) << "%)";
    LOG(INFO) << "Valid TOT: " << counter_valid_tot << " (" << (double(counter_valid_tot) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Valid TOA: " << counter_valid_toa << " (" << (double(counter_valid_toa) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Double TOT: " << counter_double_tot << " (" << (double(counter_double_tot) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Double TOA: " << counter_double_toa << " (" << (double(counter_double_toa) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "=== End of Data Reading Summary ========================";

    input_root->Close();

    output_root->cd();

    // * -- Create and save histograms here ----------------------------------------
    // * ---------------------------------------------------------------------------

    // * -- Draw histograms of ADC sums for each FPGA and all FPGAs combined -------
    auto *canvas_adc_sum = new TCanvas("canvas_adc_sum", "ADC Sum", 1200, 800);
    int canvas_hist_bin_num = 256;
    double canvas_fpga_max = 0.0;
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        for (int _entry = 0; _entry < entry_max; _entry++) {
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                canvas_fpga_max = std::max(canvas_fpga_max, fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry]);
            }
        }
    }
    canvas_adc_sum->Divide(1,2);
    canvas_adc_sum->cd(1);

    auto hist_legend = new TLegend(0.7, 0.6, 0.89, 0.89);
    hist_legend->SetFillStyle(0);
    hist_legend->SetBorderSize(0);
    hist_legend->SetTextSize(0.03);
    
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        auto *hist_adc_sum = new TH1D(("hist_adc_sum_" + std::to_string(_fpga_id)).c_str(), "ADC Sum for individual FPGAs", canvas_hist_bin_num, 0, canvas_fpga_max * 0.7);
        for (int _entry = 0; _entry < entry_max; _entry++) {
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                hist_adc_sum->Fill(fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry]);
            }
        }
        hist_adc_sum->SetLineColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        hist_adc_sum->SetMarkerColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        hist_adc_sum->SetMarkerStyle(20);
        hist_adc_sum->SetMarkerSize(0.5);
        // set the maximum value of the histogram to 1.1 times the maximum value
        hist_adc_sum->SetMaximum(hist_adc_sum->GetMaximum() * 6.0);
        hist_adc_sum->SetStats(0);
        hist_adc_sum->GetXaxis()->SetTitle("ADC Sum");
        hist_adc_sum->GetYaxis()->SetTitle("Counts");
        hist_adc_sum->Draw("HIST SAME");
        hist_legend->AddEntry(hist_adc_sum, ("FPGA " + std::to_string(_fpga_id)).c_str(), "l");
    }

    // Draw each FPGA with hamming code and daqh error filtering
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        auto *hist_adc_sum_filtered = new TH1D(("hist_adc_sum_filtered_" + std::to_string(_fpga_id)).c_str(), "ADC Sum for individual FPGAs (filtered)", canvas_hist_bin_num, 0, canvas_fpga_max * 0.7);
        for (int _entry = 0; _entry < entry_max; _entry++) {
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                if (fpgas_hamming_code_error_count[_fpga_index * 2][_entry] == 0 && fpgas_hamming_code_error_count[_fpga_index * 2 + 1][_entry] == 0 &&
                    fpgas_daqh_error_count[_fpga_index * 2][_entry] == 0 && fpgas_daqh_error_count[_fpga_index * 2 + 1][_entry] == 0) {
                    hist_adc_sum_filtered->Fill(fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry]);
                }
            }
        }
        hist_adc_sum_filtered->SetLineColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        hist_adc_sum_filtered->SetLineStyle(2); // dashed line for filtered
        hist_adc_sum_filtered->SetMarkerColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        hist_adc_sum_filtered->SetMarkerStyle(20);
        hist_adc_sum_filtered->SetMarkerSize(0.5);
        hist_adc_sum_filtered->SetStats(0);
        hist_adc_sum_filtered->GetXaxis()->SetTitle("ADC Sum");
        hist_adc_sum_filtered->GetYaxis()->SetTitle("Counts");
        hist_adc_sum_filtered->Draw("HIST SAME");
        hist_legend->AddEntry(hist_adc_sum_filtered, ("FPGA " + std::to_string(_fpga_id) + " (filtered)").c_str(), "l");
    }
    hist_legend->Draw("SAME");

    canvas_adc_sum->cd(2);

    auto *legend_fpga = new TLegend(0.5, 0.7, 0.89, 0.89);
    legend_fpga->SetFillStyle(0);
    legend_fpga->SetBorderSize(0);
    legend_fpga->SetTextSize(0.03);

    canvas_fpga_max = 0.0;
    for (int _entry = 0; _entry < entry_max; _entry++) {
        double _adc_sum = 0.0;
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                _adc_sum += fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry];
            }
        }
        canvas_fpga_max = std::max(canvas_fpga_max, _adc_sum);
    }

    auto *hist_adc_sum_all = new TH1D("hist_adc_sum_all", "ADC Sum for All FPGAs Combined", canvas_hist_bin_num, 0, canvas_fpga_max * 0.7);
    for (int _entry = 0; _entry < entry_max; _entry++) {
        double _adc_sum = 0.0;
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                _adc_sum += fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry];
            }
        }
        hist_adc_sum_all->Fill(_adc_sum);
    }
    hist_adc_sum_all->SetMaximum(hist_adc_sum_all->GetMaximum() * 1.4);
    hist_adc_sum_all->SetLineColor(kBlack);
    hist_adc_sum_all->SetMarkerColor(kBlack);
    hist_adc_sum_all->SetMarkerStyle(20);
    hist_adc_sum_all->SetMarkerSize(0.5);
    hist_adc_sum_all->SetStats(0);
    hist_adc_sum_all->GetXaxis()->SetTitle("ADC Sum");
    hist_adc_sum_all->GetYaxis()->SetTitle("Counts");
    hist_adc_sum_all->Draw("HIST SAME");
    legend_fpga->AddEntry(hist_adc_sum_all, ("All Combined, " + std::to_string(int(hist_adc_sum_all->GetEntries())) + " entries").c_str(), "l");

    // draw the sum without FPGA 3
    auto *hist_adc_sum_no_fpga3 = new TH1D("hist_adc_sum_no_fpga3", "ADC Sum for All FPGAs Combined (without FPGA 3)", canvas_hist_bin_num, 0, canvas_fpga_max * 0.7);
    for (int _entry = 0; _entry < entry_max; _entry++) {
        double _adc_sum = 0.0;
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id = legal_fpga_id_list[_fpga_index];
            if (_fpga_id == 3) {
                continue; // skip FPGA 3
            }
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                _adc_sum += fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry];
            }
        }
        hist_adc_sum_no_fpga3->Fill(_adc_sum);
    }
    hist_adc_sum_no_fpga3->SetLineColor(kGray);
    hist_adc_sum_no_fpga3->SetMarkerColor(kGray);
    hist_adc_sum_no_fpga3->SetMarkerStyle(20);
    hist_adc_sum_no_fpga3->SetMarkerSize(0.5);
    hist_adc_sum_no_fpga3->SetStats(0);
    hist_adc_sum_no_fpga3->Draw("HIST SAME");
    legend_fpga->AddEntry(hist_adc_sum_no_fpga3, ("All Combined (without FPGA 3), " + std::to_string(int(hist_adc_sum_no_fpga3->GetEntries())) + " entries").c_str(), "l");

    // draw the sum with hamming code and daqh error filtering
    auto *hist_adc_sum_all_filtered = new TH1D("hist_adc_sum_all_filtered", "ADC Sum for All FPGAs Combined (filtered)", canvas_hist_bin_num, 0, canvas_fpga_max * 0.7);
    for (int _entry = 0; _entry < entry_max; _entry++) {
        double _adc_sum = 0.0;
        bool _all_good_flags = true;
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id = legal_fpga_id_list[_fpga_index];
            if (_fpga_id == 3) {
                continue; // skip FPGA 3
            }
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                if (fpgas_hamming_code_error_count[_fpga_index * 2][_entry] == 0 && fpgas_hamming_code_error_count[_fpga_index * 2 + 1][_entry] == 0 && fpgas_daqh_error_count[_fpga_index * 2][_entry] == 0 && fpgas_daqh_error_count[_fpga_index * 2 + 1][_entry] == 0) {
                // if (fpgas_hamming_code_error_count[_fpga_index * 2][_entry] == 0 && fpgas_hamming_code_error_count[_fpga_index * 2 + 1][_entry] == 0) {
                    _adc_sum += fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry];
                }
                else {
                    _all_good_flags = false; // if any FPGA has hamming code or daqh error, skip this entry
                    break;
                }
            }
        }
        if (_adc_sum > 0 && _all_good_flags) {
            hist_adc_sum_all_filtered->Fill(_adc_sum);
        }
    }
    hist_adc_sum_all_filtered->SetLineColor(kRed);
    hist_adc_sum_all_filtered->SetMarkerColor(kRed);
    hist_adc_sum_all_filtered->SetMarkerStyle(20);
    hist_adc_sum_all_filtered->SetMarkerSize(0.5);
    hist_adc_sum_all_filtered->SetStats(0);
    hist_adc_sum_all_filtered->GetXaxis()->SetTitle("ADC Sum");
    hist_adc_sum_all_filtered->GetYaxis()->SetTitle("Counts");
    hist_adc_sum_all_filtered->Draw("HIST SAME");

    // do gaussian fit over 7000 - 9000 range
    // TF1 *gaus_fit = new TF1("gaus", "gaus", 7000, 13000);
    // hist_adc_sum_all_filtered->Fit(gaus_fit, "RQ");
    // hist_adc_sum_all_filtered->GetFunction("gaus")->SetLineColor(kRed);
    // hist_adc_sum_all_filtered->GetFunction("gaus")->SetLineWidth(2);
    // hist_adc_sum_all_filtered->GetFunction("gaus")->SetLineStyle(2);

    double hist_max_x = canvas_fpga_max * 0.7;
    
    
    int max_bin_first_half_x = 0;
    int max_bin_second_half_x = 0;
    double max_count_first_half = 0.0;
    double max_count_second_half = 0.0;
    for (int bin = 1; bin <= hist_adc_sum_all_filtered->GetNbinsX(); bin++) {
        double bin_center = hist_adc_sum_all_filtered->GetBinCenter(bin);
        double bin_content = hist_adc_sum_all_filtered->GetBinContent(bin);
        if (bin_center < hist_max_x / 2.0) {
            if (bin_content > max_count_first_half) {
                max_count_first_half = bin_content;
                max_bin_first_half_x = bin_center;
            }
        } else {
            if (bin_content > max_count_second_half) {
                max_count_second_half = bin_content;
                max_bin_second_half_x = bin_center;
            }
        }
    }
    LOG(INFO) << "Max Bin First Half X: " << max_bin_first_half_x;
    LOG(INFO) << "Max Bin Second Half X: " << max_bin_second_half_x;

    // TF1 *double_gaus_fit = new TF1("double_gaus", "gaus(0) + gaus(3)", 0, hist_max_x);
    // double_gaus_fit->SetParameters(500.0, max_bin_first_half_x, 1000.0, 500.0, max_bin_second_half_x, 1000.0);
    // // double_gaus_fit->SetParameters(500.0, 4000, 2000, 500.0, 9000, 2000);
    // hist_adc_sum_all_filtered->Fit(double_gaus_fit, "RQ");
    // hist_adc_sum_all_filtered->GetFunction("double_gaus")->SetLineColor(kBlue);
    // hist_adc_sum_all_filtered->GetFunction("double_gaus")->SetLineWidth(2);
    // hist_adc_sum_all_filtered->GetFunction("double_gaus")->SetLineStyle(2);

    // auto hist_adc_sum_fit_result_gaussian_0_mean = double_gaus_fit->GetParameter(1);
    // auto hist_adc_sum_fit_result_gaussian_0_sigma = abs(double_gaus_fit->GetParameter(2));
    // auto hist_adc_sum_fit_result_gaussian_1_mean = double_gaus_fit->GetParameter(4);
    // auto hist_adc_sum_fit_result_gaussian_1_sigma = abs(double_gaus_fit->GetParameter(5));
    // LOG(INFO) << "ADC Sum Fit Result Gaussian 0 Mean: " << hist_adc_sum_fit_result_gaussian_0_mean;
    // LOG(INFO) << "ADC Sum Fit Result Gaussian 0 Sigma: " << hist_adc_sum_fit_result_gaussian_0_sigma;
    // LOG(INFO) << "ADC Sum Fit Result Gaussian 1 Mean: " << hist_adc_sum_fit_result_gaussian_1_mean;
    // LOG(INFO) << "ADC Sum Fit Result Gaussian 1 Sigma: " << hist_adc_sum_fit_result_gaussian_1_sigma;
    // auto hist_adc_sum_fit_result_resolution = (hist_adc_sum_fit_result_gaussian_1_sigma / hist_adc_sum_fit_result_gaussian_1_mean) * 100.0;

    // LOG(INFO) << "ADC Sum Fit Result Resolution: " << hist_adc_sum_fit_result_resolution << "%";
    // double_gaus_fit->SetLineColor(kPink);
    // double_gaus_fit->SetLineWidth(2);
    // double_gaus_fit->Draw("SAME");

    TF1 *gaus_fit = new TF1("gaus", "gaus", max_bin_second_half_x - 500, max_bin_second_half_x + 1500);
    hist_adc_sum_all_filtered->Fit(gaus_fit, "RQ", "", max_bin_second_half_x - 600, max_bin_second_half_x + 3000);
    hist_adc_sum_all_filtered->GetFunction("gaus")->SetLineColor(kRed);
    hist_adc_sum_all_filtered->GetFunction("gaus")->SetLineWidth(2);
    hist_adc_sum_all_filtered->GetFunction("gaus")->SetLineStyle(2);
    
    auto hist_adc_sum_fit_result_gaussian_mean = gaus_fit->GetParameter(1);
    auto hist_adc_sum_fit_result_gaussian_sigma = abs(gaus_fit->GetParameter(2));
    LOG(INFO) << "ADC Sum Fit Result Gaussian Mean: " << hist_adc_sum_fit_result_gaussian_mean;
    LOG(INFO) << "ADC Sum Fit Result Gaussian Sigma: " << hist_adc_sum_fit_result_gaussian_sigma;
    auto hist_adc_sum_fit_result_resolution = (hist_adc_sum_fit_result_gaussian_sigma / hist_adc_sum_fit_result_gaussian_mean) * 100.0;

    LOG(INFO) << "ADC Sum Fit Result Resolution: " << hist_adc_sum_fit_result_resolution << "%";

    // Add the fit result to the legend
    


    // Draw the histogram
    // print the shape pass
    LOG(DEBUG) << "Signal Shape Pass:";
    for (int _entry = 0; _entry < 20; _entry++) {
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id = legal_fpga_id_list[_fpga_index];
            LOG(DEBUG) << "FPGA " << _fpga_id << ": " << fpgas_signal_shape_pass[_fpga_index * 2][_entry] << " (A0), " << fpgas_signal_shape_pass[_fpga_index * 2 + 1][_entry] << " (A1)";
        }
    }
    legend_fpga->AddEntry(hist_adc_sum_all_filtered, ("All Combined (filtered), " + std::to_string(int(hist_adc_sum_all_filtered->GetEntries())) + " entries").c_str(), "l");
    legend_fpga->AddEntry(gaus_fit, ("Mean: " + std::to_string(hist_adc_sum_fit_result_gaussian_mean) + ", Sigma: " + std::to_string(hist_adc_sum_fit_result_gaussian_sigma)).c_str(), "l");
    legend_fpga->AddEntry(hist_adc_sum_all_filtered, ("Resolution: " + std::to_string(hist_adc_sum_fit_result_resolution) + "%").c_str(), "l");

    // histogram which also filters out the signal shape pass
    auto *hist_adc_sum_filtered_signal_shape = new TH1D("hist_adc_sum_filtered_signal_shape", "ADC Sum for All FPGAs Combined (filtered, signal shape pass)", canvas_hist_bin_num, 0, canvas_fpga_max * 0.7);
    for (int _entry = 0; _entry < entry_max; _entry++) {
        double _adc_sum = 0.0;
        bool _all_good_flags = true;
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id = legal_fpga_id_list[_fpga_index];
            if (_fpga_id == 3) {
                continue; // skip FPGA 3
            }
            if (fpgas_adc_sum_list[_fpga_index * 2][_entry] >= 0 && fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] >= 0) {
                if (fpgas_hamming_code_error_count[_fpga_index * 2][_entry] == 0 && fpgas_hamming_code_error_count[_fpga_index * 2 + 1][_entry] == 0 && fpgas_daqh_error_count[_fpga_index * 2][_entry] == 0 && fpgas_daqh_error_count[_fpga_index * 2 + 1][_entry] == 0 && fpgas_signal_shape_pass[_fpga_index * 2][_entry] < 70) {
                    _adc_sum += fpgas_adc_sum_list[_fpga_index * 2][_entry] + fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry];
                }
                else {
                    _all_good_flags = false; // if any FPGA has hamming code or daqh error, skip this entry
                    break;
                }
            }
        }
        if (_adc_sum > 0 && _all_good_flags) {
            hist_adc_sum_filtered_signal_shape->Fill(_adc_sum);
        }
    }
    hist_adc_sum_filtered_signal_shape->SetLineColor(kGreen);
    hist_adc_sum_filtered_signal_shape->SetMarkerColor(kGreen);
    hist_adc_sum_filtered_signal_shape->SetMarkerStyle(20);
    hist_adc_sum_filtered_signal_shape->SetMarkerSize(0.5);
    hist_adc_sum_filtered_signal_shape->SetStats(0);
    hist_adc_sum_filtered_signal_shape->GetXaxis()->SetTitle("ADC Sum");
    hist_adc_sum_filtered_signal_shape->GetYaxis()->SetTitle("Counts");
    hist_adc_sum_filtered_signal_shape->Draw("HIST SAME");

    TF1 *shape_fit = new TF1("shape_fit", "gaus", max_bin_second_half_x - 500, max_bin_second_half_x + 1500);
    hist_adc_sum_filtered_signal_shape->Fit(shape_fit, "RQ", "", max_bin_second_half_x - 600, max_bin_second_half_x + 3000);
    hist_adc_sum_filtered_signal_shape->GetFunction("shape_fit")->SetLineColor(kGreen);
    hist_adc_sum_filtered_signal_shape->GetFunction("shape_fit")->SetLineWidth(2);
    hist_adc_sum_filtered_signal_shape->GetFunction("shape_fit")->SetLineStyle(2);

    auto hist_adc_sum_filtered_shape_fit_result_gaussian_mean = shape_fit->GetParameter(1);
    auto hist_adc_sum_filtered_shape_fit_result_gaussian_sigma = abs(shape_fit->GetParameter(2));
    LOG(INFO) << "ADC Sum Filtered Shape Fit Result Gaussian Mean: " << hist_adc_sum_filtered_shape_fit_result_gaussian_mean;
    LOG(INFO) << "ADC Sum Filtered Shape Fit Result Gaussian Sigma: " << hist_adc_sum_filtered_shape_fit_result_gaussian_sigma;
    auto hist_adc_sum_filtered_shape_fit_result_resolution = (hist_adc_sum_filtered_shape_fit_result_gaussian_sigma / hist_adc_sum_filtered_shape_fit_result_gaussian_mean) * 100.0;

    LOG(INFO) << "ADC Sum Filtered Shape Fit Result Resolution: " << hist_adc_sum_filtered_shape_fit_result_resolution << "%";
    legend_fpga->AddEntry(hist_adc_sum_filtered_signal_shape, ("All Combined (shape filtered), " + std::to_string(int(hist_adc_sum_filtered_signal_shape->GetEntries())) + " entries").c_str(), "l");
    legend_fpga->AddEntry(shape_fit, ("Mean: " + std::to_string(hist_adc_sum_filtered_shape_fit_result_gaussian_mean) + ", Sigma: " + std::to_string(hist_adc_sum_filtered_shape_fit_result_gaussian_sigma)).c_str(), "l");
    legend_fpga->AddEntry(hist_adc_sum_filtered_signal_shape, ("Resolution: " + std::to_string(hist_adc_sum_filtered_shape_fit_result_resolution) + "%").c_str(), "l");
    
    legend_fpga->Draw("SAME");

    canvas_adc_sum->Update();
    canvas_adc_sum->Write();

    output_root->Close();

    return 0;
}