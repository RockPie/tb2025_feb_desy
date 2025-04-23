#include "H2GCROC_Common.hxx"
#include "csv.hpp"

INITIALIZE_EASYLOGGINGPP

double ExpShifted(double *x, double *par) {
    double A = par[0];
    double x0 = par[1];
    double B = par[2];
    double C = par[3];
    if (x[0] < x0) {
        return 0;
    } else {
        return A * exp(-B * (x[0] - x0)) + C;
    }
}

double ExpCorrection(double *x, double *par) {
    double A = par[0];
    double x0 = par[1];
    double B = par[2];
    double C = par[3];
    if (x[0] < x0) {
        return 0;
    } else {
        return A * exp(-B * (x[0] - x0)) * -1.0;
    }
}

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_root_single_csv(argc, argv, "1.1");

    gROOT->SetBatch(kTRUE);

    bool enable_focal_mapping = opts.focal;
    bool enable_timewalk = opts.timewalk;
    std::string timewalk_file;
    std::vector <int> timewalk_correction_uni_channel;
    std::vector <double> timewalk_correction_parA;
    std::vector <double> timewalk_correction_parB;
    std::vector <double> timewalk_correction_parC;
    std::vector <double> timewalk_correction_parx0;
    if (!enable_timewalk) {
        timewalk_file = opts.timewalk_file;
        json timewalk_correction_json;
        std::ifstream timewalk_correction_json_file(timewalk_file);
        if (!timewalk_correction_json_file.is_open()) {
            LOG(ERROR) << "Failed to open timewalk correction json file " << timewalk_file;
            return 1;
        }
        timewalk_correction_json_file >> timewalk_correction_json;
        timewalk_correction_json_file.close();

        timewalk_correction_uni_channel = timewalk_correction_json["unified_channel_numbers"].get<std::vector<int>>();
        timewalk_correction_parA = timewalk_correction_json["A"].get<std::vector<double>>();
        timewalk_correction_parB = timewalk_correction_json["B"].get<std::vector<double>>();
        timewalk_correction_parC = timewalk_correction_json["C"].get<std::vector<double>>();
        timewalk_correction_parx0 = timewalk_correction_json["x0"].get<std::vector<double>>();
    }

    // * --- Constants ------------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 25.0 / 16.0; // unit: ns

    const int channel_adc_hist_bins   = 256;
    const double channel_adc_hist_min = 0;
    const double channel_adc_hist_max = 1024;

    const double toa_time_hist_min  = 0.0;
    const double toa_time_hist_max  = 200.0; // unit: ns
    const int toa_time_hist_bins    = int((toa_time_hist_max - toa_time_hist_min)/phase_shift_time);

    const double toa_code_hist_min  = 0.0;
    const double toa_code_hist_max  = 1024.0;
    const int toa_code_hist_bins    = 256;

    const int example_channel_number = 147;
    const int module_to_check = 5;

    // * --- Global channel painter -----------------------------------------------------
    // * --------------------------------------------------------------------------------
    GlobalChannelPainter *global_painter = nullptr;
    if (enable_focal_mapping){
        global_painter = new GlobalChannelPainter("data/SPS_2024/config/focalh_mapping.json", "data/SPS_2024/config/h2gcroc_mapping.json");
    } else {
        global_painter = new GlobalChannelPainter("data/DESY_2025/config/EEEMCal_Mapping_DESY_2025.json");
    }

    // * --- Read the csv file ----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string template_time_column_header = "Time[ns]";
    std::vector<double> template_time_column_values;
    std::vector<int> template_DAC_values;
    std::vector<std::vector<double>> template_DAC_phase_results;
    std::vector<std::vector<double>> template_DAC_phase_results_errors;

    try {
        csv::CSVReader reader(opts.csv_file);
        auto csv_headers = reader.get_col_names();
        LOG(INFO) << "CSV file " << opts.csv_file << " has " << csv_headers.size() << " columns";

        if (std::find(csv_headers.begin(), csv_headers.end(), template_time_column_header) == csv_headers.end()) {
            LOG(ERROR) << "CSV file " << opts.csv_file << " does not have column " << template_time_column_header;
            return -1;
        }

        std::vector<std::string> dac_headers, dac_error_headers;
        std::regex dac_regex("^DAC(\\d+)$");
        std::regex dac_err_regex("^DAC(\\d+)Err$");
        std::smatch match;

        for (const auto& _header: csv_headers) {
            if (std::regex_match(_header, match, dac_regex)) {
                dac_headers.push_back(_header);
            } else if (std::regex_match(_header, match, dac_err_regex)) {
                dac_error_headers.push_back(_header);
            }
        }

        if (dac_headers.size() != dac_error_headers.size()) {
            LOG(ERROR) << "DAC and DAC error headers do not match!";
            return -1;
        }

        for (const auto& _header: dac_headers) {
            std::regex_match(_header, match, dac_regex);
            template_DAC_values.push_back(std::stoi(match[1]));
        }

        template_DAC_phase_results.resize(dac_headers.size());
        template_DAC_phase_results_errors.resize(dac_headers.size());

        for (csv::CSVRow& row: reader) {
            double time_val = row[template_time_column_header].get<double>();
            template_time_column_values.push_back(time_val);

            // Read DAC phase results
            for (size_t i = 0; i < dac_headers.size(); ++i) {
                double dac_phase_val = row[dac_headers[i]].get<double>();
                template_DAC_phase_results[i].push_back(dac_phase_val);
            }

            // Read DAC errors
            for (size_t i = 0; i < dac_error_headers.size(); ++i) {
                double dac_err_val = row[dac_error_headers[i]].get<double>();
                template_DAC_phase_results_errors[i].push_back(dac_err_val);
            }
        }

    }
    catch (std::exception& e) {
        LOG(ERROR) << "CSV file " << opts.csv_file << " has error: " << e.what();
        return -1;
    }

    // * --- Read the root file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int fpga_count = -1;
    int machine_gun_samples = -1;

    TFile *input_root = new TFile(opts.input_file.c_str(), "READ");
    if (input_root->IsZombie()) {
        LOG(ERROR) << "Failed to open input file " << opts.input_file;
        return 1;
    }

    TTree *input_tree = (TTree*) input_root->Get("data_tree");
    if (input_tree == nullptr) {
        LOG(ERROR) << "Failed to get data tree from input file " << opts.input_file;
        return 1;
    }

    int entry_max = input_tree->GetEntries();
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

    TNamed *legal_fpga_id_list_tnamed = (TNamed*) input_root->Get("Rootifier_legal_fpga_id_list");
    if (legal_fpga_id_list_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get legal fpga id list from input file " << opts.input_file;
        return 1;
    }
    std::string legal_fpga_id_list_str = legal_fpga_id_list_tnamed->GetTitle();
    std::vector <UShort_t> legal_fpga_id_list;
    std::istringstream legal_fpga_id_list_stream(legal_fpga_id_list_str);
    UShort_t legal_fpga_id;
    while (legal_fpga_id_list_stream >> legal_fpga_id) {
        legal_fpga_id_list.push_back(legal_fpga_id);
    }
    if (fpga_count == -1) {
        fpga_count = legal_fpga_id_list.size();
    } else {
        if (fpga_count != legal_fpga_id_list.size()) {
            LOG(ERROR) << "Number of FPGAs in the input files do not match!";
            return 1;
        }
    }

    TNamed *_input_machine_gun_samples_tnamed = (TNamed*) input_root->Get("EventRecon_machine_gun_samples");
    if (_input_machine_gun_samples_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get machine gun samples from input file " << opts.input_file;
        return 1;
    }
    int _machine_gun_samples = std::stoi(_input_machine_gun_samples_tnamed->GetTitle());
    if (machine_gun_samples == -1) {
        machine_gun_samples = _machine_gun_samples;
    } else {
        if (machine_gun_samples != _machine_gun_samples) {
            LOG(ERROR) << "Machine gun samples in the input files do not match!";
            return 1;
        }
    }

    LOG(INFO) << "FPGA count: " << fpga_count << " Machine gun samples: " << machine_gun_samples << " Entry max: " << entry_max;

    // * --- Create the output file -----------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    std::unordered_map <int, int> channel_unified_channel_number_to_index_map;
    std::unordered_map <int, int> channel_index_to_unified_channel_number_map;

    std::vector <TH2D*> channel_adc_samples_hist_list;
    TDirectory *channel_adc_smaples_folder = output_root->mkdir("ChannelADCSamples");

    std::vector <TH2D*> channel_toa_time_toa_code_hist_list;
    TDirectory *channel_toa_time_toa_code_folder = output_root->mkdir("ChannelToaTimeToaCode");

    std::vector <TH2D*> channel_toa_time_adc_max_hist_list;
    TDirectory *channel_toa_time_adc_max_folder = output_root->mkdir("ChannelToaTimeADCMax");

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
            channel_unified_channel_number_to_index_map[channel_array_index] = _unified_valid_channel_number;
            channel_index_to_unified_channel_number_map[_unified_valid_channel_number] = channel_array_index;

            // -- Create the channel ADC samples histograms
            // ----------------------------------------------------------------
            auto *_adc_samples_hist = new TH2D(("channel_adc_samples_" + std::to_string(_unified_valid_channel_number)).c_str(), ("Channel ADC Samples " + std::to_string(_unified_valid_channel_number)).c_str(), sample_time_bins, 0, machine_gun_samples * sample_time, channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max);
            channel_adc_samples_hist_list.push_back(_adc_samples_hist);
            _adc_samples_hist->SetDirectory(channel_adc_smaples_folder);

            // -- Create the channel TOA time and TOA code histograms
            // ----------------------------------------------------------------
            auto *_toa_time_toa_code_hist = new TH2D(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa_time_toa_code").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA Time and TOA Code").c_str(), toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max, toa_code_hist_bins, toa_code_hist_min, toa_code_hist_max);
            channel_toa_time_toa_code_hist_list.push_back(_toa_time_toa_code_hist);
            _toa_time_toa_code_hist->SetDirectory(channel_toa_time_toa_code_folder);

            // -- Create the channel TOA time and ADC max histograms
            // ----------------------------------------------------------------
            auto *_toa_time_adc_max_hist = new TH2D(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa_time_adc_max").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA Time and ADC Max").c_str(), channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max, toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max);
            channel_toa_time_adc_max_hist_list.push_back(_toa_time_adc_max_hist);
            _toa_time_adc_max_hist->SetDirectory(channel_toa_time_adc_max_folder);

            channel_array_index++;
        }
    }

    // * --- Read the events ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
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

    int _input_double_tot_counter = 0;
    int _input_double_toa_counter = 0;

    int _input_valid_tot_counter = 0;
    int _input_valid_toa_counter = 0;

    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id    = legal_fpga_id_list[_fpga_index];
            auto _timestamp  = branch_timestamps_list[_fpga_index][0];
            auto _daqh_list  = branch_daqh_list_list[_fpga_index];
            auto _tc_list    = branch_tc_list_list[_fpga_index];
            auto _tp_list    = branch_tp_list_list[_fpga_index];
            auto _val0_list  = branch_val0_list_list[_fpga_index];
            auto _val1_list  = branch_val1_list_list[_fpga_index];
            auto _val2_list  = branch_val2_list_list[_fpga_index];

            std::vector <bool> _hamming_code_pass_list;
            for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                bool _hamming_code_pass = true;
                for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
                    auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
                    auto _h1h2h3 = (_daqh >> 4) & 0x7;
                    if (_h1h2h3 != 0x00){
                        _hamming_code_pass = false;
                    }
                }
                _hamming_code_pass_list.push_back(_hamming_code_pass);
            }

            std::vector <int> _channel_val0_max_list;
            std::vector <int> _channel_val1_max_list;
            std::vector <int> _channel_val2_max_list;
            std::vector <int> _channel_val0_max_index_list;
            std::vector <int> _channel_val1_max_index_list;
            std::vector <int> _channel_val2_max_index_list;

            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                auto _channel_valid = get_valid_fpga_channel(_channel_index);
                if (_channel_valid == -1){
                    continue;
                }
                auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);

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
                    if (_val0 > 0){
                        if (_val0 > _val0_max){
                            _val0_max = _val0;
                            _val0_max_index = _sample_index;
                        }
                    }
                    if (_val1 > 0){
                        if (_val1_max == -1) {
                            _val1_max = _val1;
                            _val1_max_index = _sample_index;
                            _input_valid_tot_counter++;
                        } else {
                            _input_double_tot_counter++;
                        }
                    }
                    if (_val2 > 0) {
                        if (_val2_max == -1) {
                            _val2_max = _val2;
                            _val2_max_index = _sample_index;
                            _input_valid_toa_counter++;
                        } else {
                            _input_double_toa_counter++;
                        }
                    }
                } // end of sample loop
                _channel_val0_max_list.push_back(_val0_max);
                _channel_val1_max_list.push_back(_val1_max);
                _channel_val2_max_list.push_back(_val2_max);
                _channel_val0_max_index_list.push_back(_val0_max_index);
                _channel_val1_max_index_list.push_back(_val1_max_index);
                _channel_val2_max_index_list.push_back(_val2_max_index);
            } // end of channel loop
            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                auto _channel_valid = get_valid_fpga_channel(_channel_index);
                auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                auto _asic_index = int(_unified_valid_channel_number / (FPGA_CHANNEL_NUMBER_VALID / 2));
                auto _asic_half_index = int(_unified_valid_channel_number / (FPGA_CHANNEL_NUMBER_VALID / 4)) % 2;
                auto _channel_hist_index = channel_unified_channel_number_to_index_map[_unified_valid_channel_number];
                if (_channel_valid == -1){
                    continue;
                }
                auto _val0_max = _channel_val0_max_list[_channel_valid];
                auto _val1_max = _channel_val1_max_list[_channel_valid];
                auto _val2_max = _channel_val2_max_list[_channel_valid];
                auto _val0_max_index = _channel_val0_max_index_list[_channel_valid];
                auto _val1_max_index = _channel_val1_max_index_list[_channel_valid];
                auto _val2_max_index = _channel_val2_max_index_list[_channel_valid];
                
                if (_val2_max > 0) {
                    auto _val2_time = decode_toa_value_ns(_val2_max);
                    double _toa_value_ns = _val2_max_index * sample_time + _val2_time;
                    switch (_asic_index) {
                        case 0:
                            if (_asic_half_index == 0) {
                                if (_val2_max > 734.0) {
                                    _toa_value_ns -= 25.0;
                                }
                                break;
                            } else {
                                if (_val2_max > 984.0) {
                                    _toa_value_ns -= 25.0;
                                }
                                break;
                            }
                        break;
                        case 1:
                            if (_asic_half_index == 0) {
                                if (_val2_max < 234.0) {
                                    _toa_value_ns += 25.0;
                                }
                                break;
                            } else {
                                break;
                            }
                        case 2:
                            if (_asic_half_index == 0) {
                                if (_val2_max < 234.0) {
                                    _toa_value_ns += 25.0;
                                }
                                break;
                            } else {
                                if (_val2_max > 734.0) {
                                    _toa_value_ns -= 25.0;
                                }
                                break;
                            }
                        case 3:
                            if (_asic_half_index == 0) {
                                if (_val2_max > 734.0) {
                                    _toa_value_ns -= 25.0;
                                }
                                break;
                            } else {
                                if (_val2_max > 734.0) {
                                    _toa_value_ns -= 25.0;
                                }
                                break;
                            }
                        default:
                            break;
                    }
                    // ! -- Manually fix the TOA value
                    // switch (_asic_index) {
                    //     case 0:
                    //         if (_val2_max > 734.0) {
                    //             _toa_value_ns -= 25.0;
                    //         }
                    //         break;
                    //     case 1:
                    //         if (_val2_max < 224.0) {
                    //             _toa_value_ns += 25.0;
                    //         }
                    //         break;
                    //     case 2:
                    //         if (_val2_max < 224.0) {
                    //             _toa_value_ns += 25.0;
                    //         }
                    //         break;
                    //     case 3:
                    //         if (_val2_max > 734.0) {
                    //             _toa_value_ns -= 25.0;
                    //         }
                    //         break;
                    //     default:
                    //         break;
                    // }
                    // ! -- Apply timewalk correction
                    if (!enable_timewalk) {
                        double _val0_max_double = double(_val0_max);
                        auto _timewalk_correction_index = std::find(timewalk_correction_uni_channel.begin(), timewalk_correction_uni_channel.end(), _unified_valid_channel_number);
                        if (_timewalk_correction_index != timewalk_correction_uni_channel.end()) {
                            auto _timewalk_correction_index_int = std::distance(timewalk_correction_uni_channel.begin(), _timewalk_correction_index);
                            double _timewalk_correction_parA = timewalk_correction_parA[_timewalk_correction_index_int];
                            double _timewalk_correction_parB = timewalk_correction_parB[_timewalk_correction_index_int];
                            double _timewalk_correction_parC = timewalk_correction_parC[_timewalk_correction_index_int];
                            double _timewalk_correction_parx0 = timewalk_correction_parx0[_timewalk_correction_index_int];
                            double _param[4] = {_timewalk_correction_parA, _timewalk_correction_parx0, _timewalk_correction_parB, _timewalk_correction_parC};
                            double _timewalk_correction = ExpCorrection(&_val0_max_double, _param);
                            _toa_value_ns += _timewalk_correction;
                        }
                    }
                    channel_toa_time_toa_code_hist_list[_channel_hist_index]->Fill(_toa_value_ns, _val2_max);

                    channel_toa_time_adc_max_hist_list[_channel_hist_index]->Fill(_val0_max, _toa_value_ns);

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        double _sample_time = _sample_index * sample_time;
                        _sample_time -= (_toa_value_ns - 88.0);
                        channel_adc_samples_hist_list[_channel_hist_index]->Fill(_sample_time, _val0);
                    } // end of the 2nd sample loop
                }
            } // end of the 2nd channel loop
        } // end of fpga loop
    } // end of event loop

    input_root->Close();


    TColor *softBlue   = new TColor(4001, 0.35, 0.55, 0.75);  // Steel Blue
    TColor *softRed    = new TColor(4002, 0.75, 0.35, 0.35);  // Soft Coral
    TColor *softGreen  = new TColor(4003, 0.45, 0.65, 0.45);  // Sage Green
    TColor *softPurple = new TColor(4004, 0.55, 0.45, 0.65);  // Lavender Gray
    TColor *softTeal   = new TColor(4005, 0.35, 0.65, 0.65);  // Light Teal
    TColor *softOrange = new TColor(4006, 0.85, 0.55, 0.35);  // Warm Orange
    TColor *softBrown  = new TColor(4007, 0.65, 0.55, 0.45);  // Sand Brown
    TColor *softGray   = new TColor(4008, 0.55, 0.55, 0.55);  // Mid Gray
    TColor *softOlive  = new TColor(4009, 0.65, 0.65, 0.45);  // Olive Green
    std::vector <Color_t> color_list = {kBlack, 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009};

    // * --- Save the histograms --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".pdf";

    // * --- Save the toa_time and toa_code histograms -----------------------------------
    // * --------------------------------------------------------------------------------
    std::vector <TCanvas*> canvas_channel_toa_time_toa_code_list;
    std::vector <TCanvas*> canvas_channel_toa_time_projection_list;
    for (int i = 0; i < channel_toa_time_toa_code_hist_list.size(); i++) {
        auto _canvas = new TCanvas(("ChannelTOATimeTOACode_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time and TOA Code " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        channel_toa_time_toa_code_hist_list[i]->SetStats(0);
        channel_toa_time_toa_code_hist_list[i]->Draw("COLZ");
        _canvas->Update();
        channel_toa_time_toa_code_folder->cd();
        channel_toa_time_toa_code_hist_list[i]->Write();
        
        canvas_channel_toa_time_toa_code_list.push_back(_canvas);

        auto _canvas_projection = new TCanvas(("ChannelTOATimeProjection_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time Projection " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);
        channel_toa_time_toa_code_hist_list[i]->ProjectionX()->Draw();
        _canvas_projection->Update();
        
        canvas_channel_toa_time_projection_list.push_back(_canvas_projection);
    }
    // Save only the central module histograms
    auto painter_canvas_toa_time_toa_code = global_painter->draw_module_channel_canvas(canvas_channel_toa_time_toa_code_list, channel_unified_channel_number_to_index_map, "ChannelTOATimeTOACode", "Channel TOA Time and TOA Code", module_to_check);
    if (painter_canvas_toa_time_toa_code == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_toa_time_toa_code->Write();
    painter_canvas_toa_time_toa_code->Print((pdf_file_name + "(").c_str());

    auto painter_canvas_toa_time_projection = global_painter->draw_module_channel_canvas(canvas_channel_toa_time_projection_list, channel_unified_channel_number_to_index_map, "ChannelTOATimeProjection", "Channel TOA Time Projection", module_to_check);
    if (painter_canvas_toa_time_projection == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_toa_time_projection->Write();
    painter_canvas_toa_time_projection->Print(pdf_file_name.c_str());

    // * --- Save the toa_time and adc_max histograms ------------------------------------
    // * --------------------------------------------------------------------------------
    std::vector <TCanvas*> canvas_channel_toa_time_adc_max_list;
    std::vector <TCanvas*> canvas_channel_toa_time_adc_max_fitted_list;
    std::vector <double> channel_toa_time_adc_max_threshold_channel_number_list;
    std::vector <double> channel_toa_time_adc_max_threshold_list;
    std::vector <double> channel_toa_time_adc_max_threshold_error_list;

    std::vector <double> fitting_results_A;
    std::vector <double> fitting_results_x0;
    std::vector <double> fitting_results_B;
    std::vector <double> fitting_results_C;
    std::vector <int> fitting_unified_channel_numbers;

    TGraphErrors *example_channel_toa_time_adc_max_graph = nullptr;
    TF1 *example_channel_toa_time_adc_max_fitted = nullptr;
    double example_channel_toa_time_adc_max_threshold = -1.0;

    LOG(INFO) << "Hist Array Size: " << channel_toa_time_adc_max_hist_list.size();
    for (int i = 0; i < channel_toa_time_adc_max_hist_list.size(); i++) {
        auto _canvas = new TCanvas(("ChannelTOATimeADCMax_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time and ADC Max " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        channel_toa_time_adc_max_hist_list[i]->SetStats(0);
        channel_toa_time_adc_max_hist_list[i]->Draw("COLZ");
        _canvas->Update();
        channel_toa_time_adc_max_folder->cd();
        channel_toa_time_adc_max_hist_list[i]->Write();
        
        canvas_channel_toa_time_adc_max_list.push_back(_canvas);

        auto _canvas_fitted = new TCanvas(("ChannelTOATimeADCMaxFitted_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time and ADC Max Fitted " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        std::vector <double> _bin_mean_list;
        std::vector <double> _bin_error_y_list;
        std::vector <double> _bin_x_center_list;
        std::vector <double> _bin_error_x_list;

        std::vector <double> _bin_mean_list_unfiltered;
        std::vector <double>  _bin_x_center_list_unfiltered;

        double _threshold_adc = -1.0;
        double _threshold_toa = -1.0;
        double _constant_toa  = 0.0;
        int _target_entries = int(double(channel_toa_time_adc_max_hist_list[i]->GetEntries()) / double(channel_toa_time_adc_max_hist_list[i]->GetNbinsX()) / 5.0);

        std::vector <bool> _bin_crossed_threshold_list;

        auto _hist_nbinsX = channel_toa_time_adc_max_hist_list[i]->GetNbinsX();
        for (int _bin_index = 1; _bin_index <= _hist_nbinsX; _bin_index++) {
            TH1D *_hist_projection = channel_toa_time_adc_max_hist_list[i]->ProjectionY("_py", _bin_index, _bin_index);
            double _bin_mean = _hist_projection->GetMean();
            double _bin_rms = _hist_projection->GetRMS();
            double _bin_entries = _hist_projection->GetEntries();
            _bin_crossed_threshold_list.push_back(_bin_entries > _target_entries);
            _bin_mean_list_unfiltered.push_back(_bin_mean);
            _bin_x_center_list_unfiltered.push_back(channel_toa_time_adc_max_hist_list[i]->GetXaxis()->GetBinCenter(_bin_index));
            if (_bin_entries > 0) {
                double _xCentre = channel_toa_time_adc_max_hist_list[i]->GetXaxis()->GetBinCenter(_bin_index);
                _bin_x_center_list.push_back(_xCentre);
                _bin_mean_list.push_back(_bin_mean);
                _bin_error_y_list.push_back(_bin_rms / std::sqrt(_bin_entries));
                _bin_error_x_list.push_back(channel_toa_time_adc_max_hist_list[i]->GetXaxis()->GetBinWidth(_bin_index) / 2.0);
            }
            delete _hist_projection;
        }

        // search for the first consecutive 10 bins that cross the threshold
        int _threshold_counter = 0;
        int _starting_bin_index = -1;
        for (int _bin_index = 0; _bin_index < _bin_crossed_threshold_list.size(); _bin_index++) {
            if (_bin_crossed_threshold_list[_bin_index]) {
                _threshold_counter++;
                if (_threshold_counter == 20) {
                    _starting_bin_index = _bin_index - 19;
                    break;
                }
            } else {
                _threshold_counter = 0;
            }
        }

        // search for the largest toa value in the first 10 bins that cross the threshold
        for (int _bin_index = _starting_bin_index; _bin_index < _starting_bin_index + 10; _bin_index++) {
            if (_threshold_toa < 0) {
                _threshold_toa = _bin_mean_list_unfiltered[_bin_index];
                _threshold_adc = _bin_x_center_list_unfiltered[_bin_index];
            } else {
                if (_bin_mean_list_unfiltered[_bin_index] > _threshold_toa) {
                    _threshold_toa = _bin_mean_list_unfiltered[_bin_index];
                    _threshold_adc = _bin_x_center_list_unfiltered[_bin_index];
                }
            }
        }

        // get the last 5 toa values as constant toa value
        for (int _bin_index = (_bin_mean_list.size() / 2) - 5 ; _bin_index < (_bin_mean_list.size() / 2); _bin_index++) {
            _constant_toa += _bin_mean_list[_bin_index];
        }
        _constant_toa /= 5.0;

        auto _canvas_legend = new TLegend(0.2, 0.6, 0.89, 0.89);
        _canvas_legend->SetFillColor(kWhite);
        _canvas_legend->SetLineColor(kWhite);
        _canvas_legend->SetShadowColor(kWhite);
        _canvas_legend->SetBorderSize(0);
        _canvas_legend->SetTextSize(0.05);

        TGraphErrors *_graph = new TGraphErrors(_bin_x_center_list.size(), &_bin_x_center_list[0], &_bin_mean_list[0], &_bin_error_x_list[0], &_bin_error_y_list[0]);
        _graph->SetTitle(("Channel " + std::to_string(channel_index_to_unified_channel_number_map[i]) + " TOA Time and ADC Max").c_str());
        _graph->SetMarkerStyle(20);
        _graph->SetMarkerSize(0.2);
        _graph->SetMarkerColor(kBlack);
        _graph->SetLineWidth(0);
        _graph->SetLineColorAlpha(kBlack, 0.5);
        _graph->GetXaxis()->SetRangeUser(channel_adc_hist_min, channel_adc_hist_max);
        _graph->GetYaxis()->SetRangeUser(50, 170);

        _canvas_fitted->cd();
        _graph->Draw("AP");

        if (channel_index_to_unified_channel_number_map[i] == example_channel_number) {
            example_channel_toa_time_adc_max_graph = (TGraphErrors*) _graph->Clone();
            example_channel_toa_time_adc_max_threshold = _threshold_adc;
        }

        // Draw the threshold line
        TLine *_threshold_line = new TLine(_threshold_adc, toa_time_hist_min, _threshold_adc, toa_time_hist_max);
        _threshold_line->SetLineColor(kRed);
        _threshold_line->SetLineWidth(1);
        _threshold_line->SetLineStyle(2);
        _threshold_line->Draw();
        _canvas_fitted->Update();

        // Draw the constant line
        TLine *_constant_line = new TLine(channel_adc_hist_min, _constant_toa, channel_adc_hist_max, _constant_toa);
        _constant_line->SetLineColor(kBlue);
        _constant_line->SetLineWidth(1);
        _constant_line->SetLineStyle(2);
        _constant_line->Draw();
        _canvas_fitted->Update();
        
        channel_toa_time_adc_max_folder->cd();
        _graph->Write();

        auto _exp_fit = new TF1(("exp_fit_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ExpShifted, _threshold_adc, channel_adc_hist_max - 400, 4);
        double _exp_fit_A = _threshold_toa - _constant_toa;
        if (_exp_fit_A <= 0) {
            _exp_fit_A = 50.0;
        }
        double _exp_fit_B = 0.015;
        double _exp_fit_C = _constant_toa;
        double _exp_fit_X0 = _threshold_adc;
        _exp_fit->SetParameters(_exp_fit_A, _exp_fit_X0, _exp_fit_B, _exp_fit_C);
        _exp_fit->SetParLimits(0, 0.5*_exp_fit_A, 2.0*_exp_fit_A);
        // _exp_fit->FixParameter(0, _exp_fit_A);
        // _exp_fit->SetParLimits(1, _threshold_adc-8, _threshold_adc+8);
        _exp_fit->FixParameter(1, _threshold_adc);
        _exp_fit->SetParLimits(2, 0.010, 0.022);
        _exp_fit->SetParLimits(3, _constant_toa-50.0, _constant_toa+50.0);

        _exp_fit->SetLineColor(kCyan + 3);
        _exp_fit->SetLineWidth(1);
        _exp_fit->SetLineStyle(1);

        if (enable_timewalk) {
            _graph->Fit(_exp_fit, "RNQ");
            _exp_fit->Draw("same");
            if (channel_index_to_unified_channel_number_map[i] == example_channel_number) {
                example_channel_toa_time_adc_max_fitted = (TF1*) _exp_fit->Clone();
            }
        }

        auto _exp_fit_A_value = _exp_fit->GetParameter(0);
        auto _exp_fit_B_value = _exp_fit->GetParameter(2);
        auto _exp_fit_C_value = _exp_fit->GetParameter(3);
        auto _exp_fit_X0_value = _exp_fit->GetParameter(1);

        fitting_results_A.push_back(_exp_fit_A_value);
        fitting_results_x0.push_back(_exp_fit_X0_value);
        fitting_results_B.push_back(_exp_fit_B_value);
        fitting_results_C.push_back(_exp_fit_C_value);
        fitting_unified_channel_numbers.push_back(channel_index_to_unified_channel_number_map[i]);

        LOG(DEBUG) << "For index: " << i << " Channel: " << channel_index_to_unified_channel_number_map[i] << " A: " << _exp_fit_A_value << " B: " << _exp_fit_B_value << " C: " << _exp_fit_C_value << " X0: " << _exp_fit_X0_value;
        channel_toa_time_adc_max_threshold_channel_number_list.push_back(channel_index_to_unified_channel_number_map[i]);
        channel_toa_time_adc_max_threshold_list.push_back(_threshold_adc);
        channel_toa_time_adc_max_threshold_error_list.push_back(channel_toa_time_adc_max_hist_list[i]->GetXaxis()->GetBinWidth(1) / 2.0);

        _canvas_legend->AddEntry(_graph, "Data", "p");
        _canvas_legend->AddEntry(_exp_fit, ("Fit: A = " + std::to_string(_exp_fit_A_value).substr(0, 5) + ", B = " + std::to_string(_exp_fit_B_value).substr(0, 5)).c_str(), "l");
        _canvas_legend->AddEntry(_exp_fit, ("C = " + std::to_string(_exp_fit_C_value).substr(0, 5) + ", X0 = " + std::to_string(_exp_fit_X0_value).substr(0, 5)).c_str(), "l");
        _canvas_legend->Draw();
        

        canvas_channel_toa_time_adc_max_fitted_list.push_back(_canvas_fitted);
    }
    // Save only the central module histograms
    auto painter_canvas_toa_time_adc_max = global_painter->draw_module_channel_canvas(canvas_channel_toa_time_adc_max_list, channel_unified_channel_number_to_index_map, "ChannelTOATimeADCMax", "Channel TOA Time and ADC Max", module_to_check);
    if (painter_canvas_toa_time_adc_max == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_toa_time_adc_max->Write();
    painter_canvas_toa_time_adc_max->Print(pdf_file_name.c_str());

    auto painter_canvas_toa_time_adc_max_fitted = global_painter->draw_module_channel_canvas(canvas_channel_toa_time_adc_max_fitted_list, channel_unified_channel_number_to_index_map, "ChannelTOATimeADCMaxFitted", "Channel TOA Time and ADC Max Fitted", module_to_check);
    if (painter_canvas_toa_time_adc_max_fitted == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_toa_time_adc_max_fitted->Write();
    painter_canvas_toa_time_adc_max_fitted->Print(pdf_file_name.c_str());

    // * --- Draw channel wise TOA Time and ADC Max threshold histograms -----------------
    // * --------------------------------------------------------------------------------
    auto _canvas_threshold = new TCanvas("ChannelTOATimeADCMaxThreshold", "Channel TOA Time and ADC Max Threshold", 800, 600);
    // LOG(INFO) << "Channel TOA Time and ADC Max Threshold Size: " << channel_toa_time_adc_max_threshold_channel_number_list.size();
    // LOG(INFO) << "The x-axis range is: ";
    // for (auto _channel_number : channel_toa_time_adc_max_threshold_channel_number_list) {
    //     LOG(INFO) << _channel_number;
    // }
    auto _graph_threshold = new TGraphErrors(channel_toa_time_adc_max_threshold_channel_number_list.size(), &channel_toa_time_adc_max_threshold_channel_number_list[0], &channel_toa_time_adc_max_threshold_list[0], 0, &channel_toa_time_adc_max_threshold_error_list[0]);
    _graph_threshold->SetTitle("Channel TOA Time and ADC Max Threshold");
    _graph_threshold->SetMarkerStyle(20);
    _graph_threshold->SetMarkerSize(0.2);
    _graph_threshold->SetMarkerColor(kBlack);
    _graph_threshold->SetLineWidth(0);
    _graph_threshold->SetLineColorAlpha(kBlack, 0.5);
    _graph_threshold->GetXaxis()->SetRangeUser(0, channel_toa_time_adc_max_threshold_channel_number_list.size());
    _graph_threshold->GetYaxis()->SetRangeUser(channel_adc_hist_min, channel_adc_hist_max);
    _graph_threshold->GetXaxis()->SetTitle("Channel Number");
    _graph_threshold->GetYaxis()->SetTitle("ADC Value");
    _canvas_threshold->cd();
    _graph_threshold->Draw("AP");
    _canvas_threshold->Update();

    auto _toa_theshold_latex = new TLatex();
    _toa_theshold_latex->SetNDC();
    _toa_theshold_latex->SetTextSize(0.04);
    _toa_theshold_latex->SetTextFont(62);
    const double _text_line_height = 0.04;
    const double _text_line_start = 0.85;
    const double _text_line_left = 0.13;
    _toa_theshold_latex->DrawLatex(_text_line_left, _text_line_start, "FoCal-H Prototype II");
    _toa_theshold_latex->SetTextSize(0.03);
    _toa_theshold_latex->SetTextFont(42);
    _toa_theshold_latex->DrawLatex(_text_line_left, _text_line_start - _text_line_height, "Channel TOA Time and ADC Max Threshold");
    _toa_theshold_latex->DrawLatex(_text_line_left, _text_line_start - 2 * _text_line_height, "SPS H2 Beam Line");
    _toa_theshold_latex->DrawLatex(_text_line_left, _text_line_start - 3 * _text_line_height, "September 2024");
    _toa_theshold_latex->SetTextFont(52);
    _toa_theshold_latex->SetTextColor(kGray + 3);
    _toa_theshold_latex->DrawLatex(_text_line_left, _text_line_start - 54* _text_line_height, "Work in Progress");

    _canvas_threshold->Update();
    output_root->cd();
    _canvas_threshold->Write();

    _canvas_threshold->Print(pdf_file_name.c_str());

    // * --- Save the adc_samples histograms --------------------------------------------
    // * --------------------------------------------------------------------------------
    std::vector <TCanvas*> canvas_channel_adc_samples_list;
    std::vector <TCanvas*> canvas_channel_adc_samples_fitted_list;
    for (int i = 0; i < channel_adc_samples_hist_list.size(); i++) {
        auto _canvas = new TCanvas(("ChannelADCSamples_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel ADC Samples " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        channel_adc_samples_hist_list[i]->SetStats(0);
        channel_adc_samples_hist_list[i]->Draw("COLZ");
        _canvas->Update();
        channel_adc_smaples_folder->cd();
        channel_adc_samples_hist_list[i]->Write();
        
        canvas_channel_adc_samples_list.push_back(_canvas);

        auto _canvas_fit = new TCanvas(("ChannelADCSamplesFitted_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel ADC Samples Fitted " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        std::vector <double> _bin_mean_list;
        std::vector <double> _bin_error_y_list;
        std::vector <double> _bin_x_center_list;
        std::vector <double> _bin_error_x_list;

        auto _hist_nbinsX = channel_adc_samples_hist_list[i]->GetNbinsX();
        for (int _bin_index = 1; _bin_index <= _hist_nbinsX; _bin_index++) {
            TH1D *_hist_projection = channel_adc_samples_hist_list[i]->ProjectionY("_py", _bin_index, _bin_index);
            double _bin_mean = _hist_projection->GetMean();
            double _bin_rms = _hist_projection->GetRMS();
            double _bin_entries = _hist_projection->GetEntries();
            if (_bin_entries > 0) {
                double _xCentre = channel_adc_samples_hist_list[i]->GetXaxis()->GetBinCenter(_bin_index);
                _bin_x_center_list.push_back(_xCentre);
                _bin_mean_list.push_back(_bin_mean);
                _bin_error_y_list.push_back(_bin_rms / std::sqrt(_bin_entries));
                _bin_error_x_list.push_back(channel_adc_samples_hist_list[i]->GetXaxis()->GetBinWidth(_bin_index) / 2.0);
            }
            delete _hist_projection;
        }

        TGraphErrors *_graph = new TGraphErrors(_bin_x_center_list.size(), &_bin_x_center_list[0], &_bin_mean_list[0], &_bin_error_x_list[0], &_bin_error_y_list[0]);
        _graph->SetTitle(("Channel " + std::to_string(channel_index_to_unified_channel_number_map[i]) + " ADC Samples").c_str());
        _graph->SetMarkerStyle(20);
        _graph->SetMarkerSize(0.2);
        _graph->SetMarkerColor(kBlack);
        _graph->SetLineWidth(0);
        _graph->SetLineColorAlpha(kBlack, 0.5);
        _graph->GetXaxis()->SetRangeUser(0, machine_gun_samples * sample_time);
        _graph->GetYaxis()->SetRangeUser(channel_adc_hist_min, channel_adc_hist_max);
        _graph->GetXaxis()->SetTitle("Time [ns]");
        _graph->GetYaxis()->SetTitle("ADC Value");
        _canvas_fit->cd();
        _graph->Draw("AP");
        _canvas_fit->Update();

        channel_adc_smaples_folder->cd();
        _graph->Write();

        canvas_channel_adc_samples_fitted_list.push_back(_canvas_fit);
    }
    // Save only the central module histograms
    auto painter_canvas = global_painter->draw_module_channel_canvas(canvas_channel_adc_samples_list, channel_unified_channel_number_to_index_map, "ChannelADCSamples", "Channel ADC Samples", module_to_check);
    if (painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas->Write();
    painter_canvas->Print(pdf_file_name.c_str());

    auto painter_canvas_fitted = global_painter->draw_module_channel_canvas(canvas_channel_adc_samples_fitted_list, channel_unified_channel_number_to_index_map, "ChannelADCSamplesFitted", "Channel ADC Samples Fitted", module_to_check);
    if (painter_canvas_fitted == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_fitted->Write();
    painter_canvas_fitted->Print(pdf_file_name.c_str());

    // * --- Save the example channel histograms -----------------------------------------
    // * --------------------------------------------------------------------------------
    auto canvas_example_channel_toa_time_adc_max = new TCanvas("ExampleChannelTOATimeADCMax", "Example Channel TOA Time and ADC Max", 800, 600);
    auto canvas_example_channel_toa_time_adc_max_legend = new TLegend(0.5, 0.7, 0.89, 0.89);
    canvas_example_channel_toa_time_adc_max_legend->SetFillColor(0);
    canvas_example_channel_toa_time_adc_max_legend->SetLineColor(0);
    canvas_example_channel_toa_time_adc_max_legend->SetShadowColor(0);
    canvas_example_channel_toa_time_adc_max_legend->SetBorderSize(0);
    canvas_example_channel_toa_time_adc_max_legend->SetTextSize(0.02);
    
    example_channel_toa_time_adc_max_graph->GetXaxis()->SetRangeUser(channel_adc_hist_min, 620);
    example_channel_toa_time_adc_max_graph->GetYaxis()->SetRangeUser(50, 170);
    example_channel_toa_time_adc_max_graph->GetXaxis()->SetTitle("ADC Value");
    example_channel_toa_time_adc_max_graph->GetYaxis()->SetTitle("TOA Time [ns]");
    example_channel_toa_time_adc_max_graph->SetTitle("Example Channel TOA Time and ADC Max");
    example_channel_toa_time_adc_max_graph->SetMarkerStyle(20);
    example_channel_toa_time_adc_max_graph->SetMarkerSize(0.2);
    example_channel_toa_time_adc_max_graph->SetMarkerColor(kBlack);
    example_channel_toa_time_adc_max_graph->SetLineWidth(1);
    example_channel_toa_time_adc_max_graph->SetLineColor(kBlack);

    example_channel_toa_time_adc_max_graph->Draw("APE");
    canvas_example_channel_toa_time_adc_max_legend->AddEntry(example_channel_toa_time_adc_max_graph, "Data", "pe");

    // Draw the threshold line
    TLine *threshold_line = new TLine(example_channel_toa_time_adc_max_threshold, 50, example_channel_toa_time_adc_max_threshold, 170);
    threshold_line->SetLineColor(kCyan + 3);
    threshold_line->SetLineWidth(1);
    threshold_line->SetLineStyle(2);
    threshold_line->Draw("same");
    canvas_example_channel_toa_time_adc_max_legend->AddEntry(threshold_line, ("ToA Threshold: " + std::to_string(int(example_channel_toa_time_adc_max_threshold)) + " ADC").c_str(), "l");

    // Draw the fitted line
    example_channel_toa_time_adc_max_fitted->SetLineColor(kPink + 7);
    example_channel_toa_time_adc_max_fitted->SetLineWidth(2);
    example_channel_toa_time_adc_max_fitted->SetLineStyle(1);
    example_channel_toa_time_adc_max_fitted->Draw("same");

    auto example_fitting_A = example_channel_toa_time_adc_max_fitted->GetParameter(0);
    auto example_fitting_B = example_channel_toa_time_adc_max_fitted->GetParameter(2);
    auto example_fitting_C = example_channel_toa_time_adc_max_fitted->GetParameter(3);
    auto example_fitting_x0 = example_channel_toa_time_adc_max_fitted->GetParameter(1);
    
    canvas_example_channel_toa_time_adc_max_legend->AddEntry(example_channel_toa_time_adc_max_fitted, ("Fit: " + std::to_string(example_fitting_A).substr(0, 5) + " * exp(" + std::to_string(example_fitting_B).substr(0, 5) + " * (x - " + std::to_string(example_fitting_x0).substr(0, 5) + ")) + " + std::to_string(example_fitting_C).substr(0, 5)).c_str(), "l");

    // // Draw the confident band
    // auto _fit_confidence_band = new TH1F("resolution_caen_fit_confidence_band", "resolution_caen_fit_confidence_band", 200, channel_adc_hist_min, channel_adc_hist_max);
    // TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    // if (fitter) {
    //     fitter->GetConfidenceIntervals(_fit_confidence_band, 0.68);
    // } else {
    //     LOG(WARNING) << "Fitter is not initialized. Confidence intervals will not be drawn.";
    // }
    // _fit_confidence_band->SetFillColorAlpha(kPink + 7, 0.5);
    // _fit_confidence_band->SetFillStyle(1001);
    // _fit_confidence_band->SetLineColorAlpha(kPink + 7, 0.0);
    // _fit_confidence_band->SetMarkerColorAlpha(kPink + 7, 0.0);
    // _fit_confidence_band->Draw("e3 same");

    canvas_example_channel_toa_time_adc_max->cd();
    canvas_example_channel_toa_time_adc_max_legend->Draw();

    auto example_channel_toa_time_adc_max_latex = new TLatex();
    example_channel_toa_time_adc_max_latex->SetNDC();
    example_channel_toa_time_adc_max_latex->SetTextSize(0.04);
    example_channel_toa_time_adc_max_latex->SetTextFont(62);
    example_channel_toa_time_adc_max_latex->DrawLatex(0.13, 0.85, "FoCal-H Prototype II");
    example_channel_toa_time_adc_max_latex->SetTextSize(0.03);
    example_channel_toa_time_adc_max_latex->SetTextFont(42);
    example_channel_toa_time_adc_max_latex->DrawLatex(0.13, 0.81, ("Channel " + std::to_string(example_channel_number) + ", ToA Time vs. ADC Max").c_str());
    example_channel_toa_time_adc_max_latex->DrawLatex(0.13, 0.77, "SPS H2 Beam Line");
    example_channel_toa_time_adc_max_latex->DrawLatex(0.13, 0.73, "September 2024");
    example_channel_toa_time_adc_max_latex->SetTextFont(52);
    example_channel_toa_time_adc_max_latex->SetTextColor(kGray + 3);
    example_channel_toa_time_adc_max_latex->DrawLatex(0.13, 0.69, "Work in Progress");

    canvas_example_channel_toa_time_adc_max->Update();
    output_root->cd();
    canvas_example_channel_toa_time_adc_max->Write();
    canvas_example_channel_toa_time_adc_max->Print(pdf_file_name.c_str());


    // Close the pdf file by saving a dummy canvas
    auto dummy_canvas = new TCanvas("dummy_canvas", "dummy_canvas", 800, 600);
    dummy_canvas->Print((pdf_file_name + ")").c_str());
    dummy_canvas->Close();

    output_root->Close();

    if (enable_timewalk){
        json timewalk_json;
        timewalk_json["unified_channel_numbers"] = fitting_unified_channel_numbers;
        timewalk_json["A"] = fitting_results_A;
        timewalk_json["x0"] = fitting_results_x0;
        timewalk_json["B"] = fitting_results_B;
        timewalk_json["C"] = fitting_results_C;

        std::string timewalk_json_file = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".json";
        std::ofstream timewalk_json_stream(timewalk_json_file);
        timewalk_json_stream << timewalk_json.dump(4);
        timewalk_json_stream.close();
    }
    return 0;
}