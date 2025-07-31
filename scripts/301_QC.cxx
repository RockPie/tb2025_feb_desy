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

    // for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //     auto _fpga_id = legal_fpga_id_list[_fpga_index];
    //     for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
    //         auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_index);
    //         auto _channel_valid = get_valid_fpga_channel(_channel_index);
            
    //         std::vector <std::vector <double> > _input_adc_time_time_list_temp;
    //         std::vector <std::vector <double> > _input_adc_time_adc_list_temp;
    //         std::vector <bool> _input_adc_time_list_temp_valid;
    //         std::vector <double> _input_adc_time_toa_time_list_temp;
    //         std::vector <double> _input_adc_time_tot_time_list_temp;
            
    //         _input_adc_time_time_list_temp.resize(entry_max);
    //         _input_adc_time_adc_list_temp.resize(entry_max);
    //         _input_adc_time_list_temp_valid.resize(entry_max);
    //         _input_adc_time_toa_time_list_temp.resize(entry_max);
    //         _input_adc_time_tot_time_list_temp.resize(entry_max);
            
    //     }
    // }
    long counter_total_events = 0;
    long counter_hamming_code_error = 0;
    long counter_bad_daqh_start_end = 0;

    long counter_double_tot = 0;
    long counter_double_toa = 0;
    long counter_valid_tot  = 0;
    long counter_valid_toa  = 0;

    // * --- Create data structures for plotting ----------------------------------------
    // * --------------------------------------------------------------------------------
    std::vector <std::vector<double>> fpgas_adc_sum_list;
    std::vector <std::vector<double>> fpgas_tot_sum_list;
    std::vector <std::vector<double>> fpgas_tot_count_list;
    std::vector <std::vector<double>> fpgas_toa_count_list;
    std::vector <std::vector<double>> fpgas_hamming_code_error_rate;

    std::vector <std::vector<double>> fpgas_hamming_code_error_count;
    std::vector <std::vector<double>> fpgas_daqh_error_count;

    double* fpgas_hamming_code_counter = new double[fpga_count];
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        fpgas_hamming_code_counter[_fpga_index] = 0.0;
    }
    fpgas_adc_sum_list.resize(fpga_count*2);
    fpgas_tot_sum_list.resize(fpga_count*2);
    fpgas_tot_count_list.resize(fpga_count*2);
    fpgas_toa_count_list.resize(fpga_count*2);
    fpgas_hamming_code_error_count.resize(fpga_count*2);
    fpgas_hamming_code_error_rate.resize(fpga_count*2);
    fpgas_daqh_error_count.resize(fpga_count*2);
    for (int _fpga_index = 0; _fpga_index < fpga_count*2; _fpga_index++) {
        fpgas_adc_sum_list[_fpga_index].resize(entry_max, -1.0);
        fpgas_tot_sum_list[_fpga_index].resize(entry_max, -1.0);
        fpgas_tot_count_list[_fpga_index].resize(entry_max, -1.0);
        fpgas_toa_count_list[_fpga_index].resize(entry_max, -1.0);
        fpgas_hamming_code_error_count[_fpga_index].resize(entry_max, -1.0);
        fpgas_hamming_code_error_rate[_fpga_index].resize(entry_max, -1.0);
        fpgas_daqh_error_count[_fpga_index].resize(entry_max, -1.0);
    }

    std::vector <std::vector<double>> fpgas_timestamp_list;
    fpgas_timestamp_list.resize(fpga_count);
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        fpgas_timestamp_list[_fpga_index].reserve(entry_max);
    }

    std::vector <std::vector<double>> fpgas_event_index_list;
    fpgas_event_index_list.resize(fpga_count);
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        fpgas_event_index_list[_fpga_index].reserve(entry_max);
    }

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

            fpgas_timestamp_list[_fpga_index].push_back(_timestamp);
            fpgas_event_index_list[_fpga_index].push_back(counter_total_events);

            if (skip_event_hamming_code){
                counter_hamming_code_error++;
                fpgas_hamming_code_counter[_fpga_index]++;
                for (int _fpga_index2 = 0; _fpga_index2 < fpga_count; _fpga_index2++) {
                    fpgas_hamming_code_error_rate[_fpga_index2][_entry] = fpgas_hamming_code_counter[_fpga_index2];
                }
                continue;
            }
            if (!good_header_bytes){
                counter_bad_daqh_start_end++;
                for (int _fpga_index2 = 0; _fpga_index2 < fpga_count; _fpga_index2++) {
                    fpgas_hamming_code_error_rate[_fpga_index2][_entry] = fpgas_hamming_code_counter[_fpga_index2];
                }
                continue;
            }

            fpgas_hamming_code_error_rate[_fpga_index][_entry] = double(fpgas_hamming_code_counter[_fpga_index]);

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

                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    auto _val0 = int(_val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER]);
                    auto _val1 = int(_val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER]);
                    auto _val2 = int(_val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER]);

                    if (_val0 > 0 && _val0 > _val0_max){
                        _val0_max = _val0;
                        _val0_max_index = _sample_index;
                        if (_channel_index < 76) { // for the first ASIC
                             if (_fpga_event_adc_sum_a0 < 0) {
                                _fpga_event_adc_sum_a0 = 0;
                            }
                            _fpga_event_adc_sum_a0 += _val0;
                        } else {
                            if (_fpga_event_adc_sum_a1 < 0) {
                                _fpga_event_adc_sum_a1 = 0;
                            }
                            _fpga_event_adc_sum_a1 += _val0;
                        }
                       
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
            } // end of channel loop
            fpgas_adc_sum_list[_fpga_index * 2][_entry] = _fpga_event_adc_sum_a0;
            fpgas_adc_sum_list[_fpga_index * 2 + 1][_entry] = _fpga_event_adc_sum_a1;
            fpgas_tot_sum_list[_fpga_index * 2][_entry] = _fpga_event_tot_sum_a0;
            fpgas_tot_sum_list[_fpga_index * 2 + 1][_entry] = _fpga_event_tot_sum_a1;
            fpgas_tot_count_list[_fpga_index * 2][_entry] = _fpga_event_tot_count_a0;
            fpgas_tot_count_list[_fpga_index * 2 + 1][_entry] = _fpga_event_tot_count_a1;
            fpgas_toa_count_list[_fpga_index * 2][_entry] = _fpga_event_toa_count_a0;
            fpgas_toa_count_list[_fpga_index * 2 + 1][_entry] = _fpga_event_toa_count_a1;
            // fpgas_adc_sum_list[_fpga_index][_entry] = _fpga_event_adc_sum;
            // fpgas_tot_sum_list[_fpga_index][_entry] = _fpga_event_tot_sum;
            // fpgas_tot_count_list[_fpga_index][_entry] = _fpga_event_tot_count;
            // fpgas_toa_count_list[_fpga_index][_entry] = _fpga_event_toa_count;
        } // end of fpga loop
    } // end of entry loop

    LOG(INFO) << "=== Data Reading Summary ===================================================";
    LOG(INFO) << "Total events: " << counter_total_events;
    LOG(INFO) << "Hamming code errors: " << counter_hamming_code_error << " (" << (double(counter_hamming_code_error) / counter_total_events * 100.0) << "%)";
    LOG(INFO) << "Bad DAQH start/end: " << counter_bad_daqh_start_end << " (" << (double(counter_bad_daqh_start_end) / counter_total_events * 100.0) << "%)";
    LOG(INFO) << "Valid TOT: " << counter_valid_tot << " (" << (double(counter_valid_tot) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Valid TOA: " << counter_valid_toa << " (" << (double(counter_valid_toa) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Double TOT: " << counter_double_tot << " (" << (double(counter_double_tot) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Double TOA: " << counter_double_toa << " (" << (double(counter_double_toa) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "=== End of Data Reading Summary ============================================";

    input_root->Close();

    output_root->cd();
    int draw_bin_number = 40;
    // * --- Draw the distributions of ADC sums for each FPGA channel with time ---------
    // * --------------------------------------------------------------------------------
    // auto time_cluster_size = int(entry_max / 40); // bin size for the histogram
    // if (time_cluster_size < 1) {
    //     time_cluster_size = 1;
    // }
    TCanvas *adc_sum_canvas = new TCanvas("adc_sum_canvas", "ADC Sum vs Time", 800, 600);
    TLegend *adc_sum_legend = new TLegend(0.15, 0.15, 0.5, 0.3);

    drawTwoASIC(adc_sum_canvas, adc_sum_legend, draw_bin_number, fpgas_adc_sum_list, legal_fpga_id_list, fpga_colors, "ADC Sum");

    TLatex *adc_sum_annotation = new TLatex();
    adc_sum_annotation->SetTextAlign(12); // center
    adc_sum_annotation->SetTextSize(0.04);
    adc_sum_annotation->SetTextFont(62);
    adc_sum_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    adc_sum_annotation->SetTextSize(0.03);
    adc_sum_annotation->SetTextFont(42);
    adc_sum_annotation->DrawLatexNDC(0.12, 0.82, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    adc_sum_annotation->DrawLatexNDC(0.12, 0.78, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    adc_sum_annotation->DrawLatexNDC(0.12, 0.74, ("Machine Gun Samples: " + std::to_string(machine_gun_samples)).c_str());
    adc_sum_canvas->Update();
    adc_sum_canvas->Write();

    // * --- Draw the distributions of TOT sums for each FPGA channel with time ---------
    // * --------------------------------------------------------------------------------
    TCanvas *tot_sum_canvas = new TCanvas("tot_sum_canvas", "TOT Sum vs Time", 800, 600);
    TLegend *tot_sum_legend = new TLegend(0.15, 0.15, 0.5, 0.3);

    drawTwoASIC(tot_sum_canvas, tot_sum_legend, draw_bin_number, fpgas_tot_sum_list, legal_fpga_id_list, fpga_colors, "ToT Sum");

    TLatex *tot_sum_annotation = new TLatex();
    tot_sum_annotation->SetTextAlign(12); // center
    tot_sum_annotation->SetTextSize(0.04);
    tot_sum_annotation->SetTextFont(62);
    tot_sum_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    tot_sum_annotation->SetTextSize(0.03);
    tot_sum_annotation->SetTextFont(42);
    tot_sum_annotation->DrawLatexNDC(0.12, 0.82, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    tot_sum_annotation->DrawLatexNDC(0.12, 0.78, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    tot_sum_annotation->DrawLatexNDC(0.12, 0.74, ("Machine Gun Samples: " + std::to_string(machine_gun_samples)).c_str());
    tot_sum_canvas->Update();
    tot_sum_canvas->Write();

    // * --- Draw the distributions of TOT counts for each FPGA channel with time --------
    // * ---------------------------------------------------------------------------------
    TCanvas *tot_count_canvas = new TCanvas("tot_count_canvas", "TOT Count vs Time", 800, 600);
    TLegend *tot_count_legend = new TLegend(0.15, 0.15, 0.5, 0.3);

    drawTwoASIC(tot_count_canvas, tot_count_legend, draw_bin_number, fpgas_tot_count_list, legal_fpga_id_list, fpga_colors, "ToT Count");

    TLatex *tot_count_annotation = new TLatex();
    tot_count_annotation->SetTextAlign(12); // center
    tot_count_annotation->SetTextSize(0.04);
    tot_count_annotation->SetTextFont(62);
    tot_count_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    tot_count_annotation->SetTextSize(0.03);
    tot_count_annotation->SetTextFont(42);
    tot_count_annotation->DrawLatexNDC(0.12, 0.82, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    tot_count_annotation->DrawLatexNDC(0.12, 0.78, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    tot_count_annotation->DrawLatexNDC(0.12, 0.74, ("Machine Gun Samples: " + std::to_string(machine_gun_samples)).c_str());
    tot_count_canvas->Update();
    tot_count_canvas->Write();

    // * --- Draw the distributions of TOA counts for each FPGA channel with time --------
    // * ---------------------------------------------------------------------------------
    TCanvas *toa_count_canvas = new TCanvas("toa_count_canvas", "TOA Count vs Time", 800, 600);
    TLegend *toa_count_legend = new TLegend(0.15, 0.15, 0.5, 0.3);

    drawTwoASIC(toa_count_canvas, toa_count_legend, draw_bin_number, fpgas_toa_count_list, legal_fpga_id_list, fpga_colors, "ToA Count");

    TLatex *toa_count_annotation = new TLatex();
    toa_count_annotation->SetTextAlign(12); // center
    toa_count_annotation->SetTextSize(0.04);
    toa_count_annotation->SetTextFont(62);
    toa_count_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    toa_count_annotation->SetTextSize(0.03);
    toa_count_annotation->SetTextFont(42);
    toa_count_annotation->DrawLatexNDC(0.12, 0.82, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    toa_count_annotation->DrawLatexNDC(0.12, 0.78, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    toa_count_annotation->DrawLatexNDC(0.12, 0.74, ("Machine Gun Samples: " + std::to_string(machine_gun_samples)).c_str());
    toa_count_canvas->Update();
    toa_count_canvas->Write();

    // * --- Draw the hamming code error count for each FPGA channel with time ---------
    // * --------------------------------------------------------------------------------
    TCanvas *hamming_code_count_canvas = new TCanvas("hamming_code_count_canvas", "Hamming Code Error Count vs Time", 800, 600);
    TLegend *hamming_code_count_legend = new TLegend(0.65, 0.6, 0.9, 0.9);
    
    drawTwoASIC(hamming_code_count_canvas, hamming_code_count_legend, draw_bin_number, fpgas_hamming_code_error_count, legal_fpga_id_list, fpga_colors, "Hamming Code Error Count");
    
    TLatex *hamming_code_count_annotation = new TLatex();
    hamming_code_count_annotation->SetTextAlign(12);
    hamming_code_count_annotation->SetTextSize(0.04);
    hamming_code_count_annotation->SetTextFont(62);
    hamming_code_count_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    hamming_code_count_annotation->SetTextSize(0.03);
    hamming_code_count_annotation->SetTextFont(42);
    hamming_code_count_annotation->DrawLatexNDC(0.12, 0.82, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    hamming_code_count_annotation->DrawLatexNDC(0.12, 0.78, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    hamming_code_count_annotation->DrawLatexNDC(0.12, 0.74, ("Machine Gun Samples: " + std::to_string(machine_gun_samples)).c_str());
    
    hamming_code_count_canvas->Update();
    hamming_code_count_canvas->Write();

    // * --- Draw the DAQH error count for each FPGA channel with time -----------------
    // * --------------------------------------------------------------------------------
    TCanvas *daqh_error_count_canvas = new TCanvas("daqh_error_count_canvas", "DAQH Error Count vs Time", 800, 600);
    TLegend *daqh_error_count_legend = new TLegend(0.65, 0.6, 0.9, 0.9);
    
    drawTwoASIC(daqh_error_count_canvas, daqh_error_count_legend, draw_bin_number, fpgas_daqh_error_count, legal_fpga_id_list, fpga_colors, "DAQH Error Count");
    
    TLatex *daqh_error_count_annotation = new TLatex();
    daqh_error_count_annotation->SetTextAlign(12);
    daqh_error_count_annotation->SetTextSize(0.04);
    daqh_error_count_annotation->SetTextFont(62);
    daqh_error_count_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    daqh_error_count_annotation->SetTextSize(0.03);
    daqh_error_count_annotation->SetTextFont(42);
    daqh_error_count_annotation->DrawLatexNDC(0.12, 0.82, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    daqh_error_count_annotation->DrawLatexNDC(0.12, 0.78, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    daqh_error_count_annotation->DrawLatexNDC(0.12, 0.74, ("Machine Gun Samples: " + std::to_string(machine_gun_samples)).c_str());
    
    daqh_error_count_canvas->Update();
    daqh_error_count_canvas->Write();

    // * --- Draw the incremental Hamming code error rate for each FPGA channel ---------
    // * --------------------------------------------------------------------------------
    TCanvas *hamming_code_error_canvas = new TCanvas("hamming_code_error_canvas", "Hamming Code Error Rate vs Time", 800, 600);
    TLegend *hamming_code_error_legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    // no boarder, no fill color, no fill style
    hamming_code_error_legend->SetBorderSize(0);
    hamming_code_error_legend->SetFillColor(0);
    hamming_code_error_legend->SetFillStyle(0);
    hamming_code_error_legend->SetTextSize(0.03);
    hamming_code_error_canvas->cd(); 

    double hamming_code_error_maximum = 0.0;
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_hamming_code_error_max = *std::max_element(fpgas_hamming_code_error_rate[_fpga_index].begin(), fpgas_hamming_code_error_rate[_fpga_index].end());
        _fpga_hamming_code_error_max /= double(counter_total_events) * double(fpga_count); // normalize by the number of events
        if (_fpga_hamming_code_error_max > hamming_code_error_maximum) {
            hamming_code_error_maximum = _fpga_hamming_code_error_max;
        }
    }

    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        auto _hamming_code_error_rate_list = fpgas_hamming_code_error_rate[_fpga_index];
        // divided by the number of events to get the error rate
        for (int _time_index = 1; _time_index < entry_max; _time_index++) {
            _hamming_code_error_rate_list[_time_index] /= double(_time_index) * double(fpga_count);
        }
        TGraph *hamming_code_error_graph = new TGraph(entry_max);
        for (int _time_index = 0; _time_index < entry_max; _time_index++) {
            hamming_code_error_graph->SetPoint(_time_index, _time_index, _hamming_code_error_rate_list[_time_index]);
        }
        std::string legend_name = "FPGA " + std::to_string(_fpga_id);
        hamming_code_error_legend->AddEntry(hamming_code_error_graph, legend_name.c_str(), "l");
        hamming_code_error_graph->SetTitle("");
        hamming_code_error_graph->GetXaxis()->SetTitle("Event Index");
        // hamming_code_error_graph->GetXaxis()->SetRangeUser(-0.05 * entry_max, 1.05 * entry_max);
        hamming_code_error_graph->GetYaxis()->SetTitle("Hamming Code Error Rate");
        hamming_code_error_graph->GetYaxis()->SetTitleOffset(1.2);
        hamming_code_error_graph->GetYaxis()->SetRangeUser(0, 1.2 * hamming_code_error_maximum);
        hamming_code_error_graph->SetMarkerStyle(20);
        hamming_code_error_graph->SetMarkerColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        hamming_code_error_graph->SetLineColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        hamming_code_error_graph->SetLineWidth(2);
        if (_fpga_index == 0) {
            hamming_code_error_graph->Draw("AL");
        } else {
            hamming_code_error_graph->Draw("L");
        }
    }
    hamming_code_error_legend->Draw();
    // add annotations
    TLatex *hamming_code_error_annotation = new TLatex();
    
    hamming_code_error_annotation->SetTextAlign(12); // center
    
    hamming_code_error_annotation->SetTextSize(0.04);
    hamming_code_error_annotation->SetTextFont(62);
    hamming_code_error_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());

    hamming_code_error_annotation->SetTextSize(0.03);
    hamming_code_error_annotation->SetTextFont(42);
    hamming_code_error_annotation->DrawLatexNDC(0.12, 0.82, "Hamming Code Error Rate");
    hamming_code_error_annotation->DrawLatexNDC(0.12, 0.78, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    hamming_code_error_annotation->DrawLatexNDC(0.12, 0.74, ("FPGA Count: " + std::to_string(fpga_count)).c_str());

    hamming_code_error_canvas->Update();
    hamming_code_error_canvas->Write();

    // * --- Draw the timestamps for each FPGA ------------------------------------------
    // * --------------------------------------------------------------------------------
    TCanvas *timestamp_canvas = new TCanvas("timestamp_canvas", "Timestamps vs event index", 800, 600);
    TLegend *timestamp_legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    // no boarder, no fill color, no fill style
    timestamp_legend->SetBorderSize(0);
    timestamp_legend->SetFillColor(0);  
    timestamp_legend->SetFillStyle(0);
    timestamp_legend->SetTextSize(0.03);
    timestamp_canvas->cd();
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        auto _timestamp_list = fpgas_timestamp_list[_fpga_index];
        auto _event_index_list = fpgas_event_index_list[_fpga_index];

        TGraph *timestamp_graph = new TGraph(_timestamp_list.size(), 
                                              _event_index_list.data(), 
                                              _timestamp_list.data());
        std::string legend_name = "FPGA " + std::to_string(_fpga_id);
        timestamp_legend->AddEntry(timestamp_graph, legend_name.c_str(), "l");
        timestamp_graph->SetTitle("");
        timestamp_graph->GetXaxis()->SetTitle("Event Index");
        // timestamp_graph->GetXaxis()->SetRangeUser(-0.05 * entry_max, 1.05 * entry_max);
        timestamp_graph->GetYaxis()->SetTitle("Timestamp [25 ns]");
        timestamp_graph->GetYaxis()->SetTitleOffset(1.2);
        timestamp_graph->SetMarkerStyle(20);
        timestamp_graph->SetMarkerColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        timestamp_graph->SetLineColor(fpga_colors[_fpga_id % fpga_colors.size()]);
        timestamp_graph->SetLineWidth(2);
        if (_fpga_index == 0) {
            timestamp_graph->Draw("AL");
        } else {
            timestamp_graph->Draw("L");
        }
    }

    timestamp_legend->Draw();
    // add annotations
    TLatex *timestamp_annotation = new TLatex();
    timestamp_annotation->SetTextAlign(12); // center
    timestamp_annotation->SetTextSize(0.04);
    timestamp_annotation->SetTextFont(62);
    timestamp_annotation->DrawLatexNDC(0.12, 0.87, run_info_str.c_str());
    timestamp_annotation->SetTextSize(0.03);
    timestamp_annotation->SetTextFont(42);
    timestamp_annotation->DrawLatexNDC(0.12, 0.82, "Timestamps for each FPGA");
    timestamp_annotation->DrawLatexNDC(0.12, 0.78, ("Total Events: " + std::to_string(counter_total_events)).c_str());
    timestamp_annotation->DrawLatexNDC(0.12, 0.74, ("FPGA Count: " + std::to_string(fpga_count)).c_str());
    timestamp_canvas->Update();
    timestamp_canvas->Write();

    
    output_root->Close();

    return 0;
}