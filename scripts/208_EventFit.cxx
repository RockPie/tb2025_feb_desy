#include "H2GCROC_Common.hxx"
#include "csv.hpp"

INITIALIZE_EASYLOGGINGPP

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

double fitParameterizedFunction(double *x, double *p, std::vector<double> para_time, std::vector<double> para_values) {
    double _xx = x[0] - p[1];
    if (_xx < para_time.front() || _xx > para_time.back()) {
        return 0;
    }
    double _template_pedestal = (para_values[0] + para_values[1]) * 0.5;
    for (size_t i = 0; i < para_time.size() - 1; ++i) {
        if (_xx >= para_time[i] && _xx < para_time[i+1]) {
            double _slope = (para_values[i+1] - para_values[i]) / (para_time[i+1] - para_time[i]);
            double _prediction = para_values[i] + _slope * (_xx - para_time[i]);
            double _res = (p[0] * (_prediction - _template_pedestal)) + _template_pedestal;
            if (_res > 1023) {
                return 1023;
            } else {
                return _res;
            }
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_root_single_csv(argc, argv, "1.0");

    ROOT::EnableImplicitMT(30);
    gROOT->SetBatch(kTRUE);

    // * --- Read ToA Fine TDC DNL correction file -----------------------------------
    // * -----------------------------------------------------------------------------
    std::string DNL_fine_TDC_json_file = "dump/207_ToACalib/FoCal/Run770_toa_calib_fineDNL.json";
    json DNL_fine_TDC_json;
    std::ifstream DNL_fine_TDC_json_file_stream(DNL_fine_TDC_json_file);
    if (!DNL_fine_TDC_json_file_stream.is_open()) {
        LOG(ERROR) << "Failed to open DNL fine TDC json file " << DNL_fine_TDC_json_file;
        return 1;
    }
    DNL_fine_TDC_json_file_stream >> DNL_fine_TDC_json;

    std::vector<std::vector<double>> fine_DNL_factors;
    std::unordered_map<int, int> channel_to_fine_DNL_index;

    int fine_DNL_index = 0;
    for (auto& [key_str, value] : DNL_fine_TDC_json.items()) {
        //  DNL index " << fine_DNL_index;
        int channel_number = std::stoi(key_str);
        channel_to_fine_DNL_index[channel_number] = fine_DNL_index;
        // LOG(INFO) << "Channel " << channel_number << " has fine DNL index " << fine_DNL_index;
    
        std::vector<double> factors;
    
        // Check if value is an array and all elements are null
        bool all_null = value.is_array() && std::all_of(value.begin(), value.end(), [](const json& x) {
            return x.is_null();
        });
    
        if (all_null) {
            // No valid DNL — use default neutral correction
            factors = std::vector<double>(8, 1.0);  // assuming 3-bit TDC = 8 bins
            // LOG(WARNING) << "Channel " << channel_number << " has no valid DNL. Using neutral factors (1.0).";
        } else {
            factors = value.get<std::vector<double>>();
        }
    
        fine_DNL_factors.push_back(factors);
        ++fine_DNL_index;
    }

    DNL_fine_TDC_json_file_stream.close();

    LOG(INFO) << "Fine DNL factors size: " << fine_DNL_factors.size();

    // * --- Read ToA coarse TDC DNL correction file ---------------------------------
    // * -----------------------------------------------------------------------------
    std::string DNL_corase_TDC_json_file = "dump/207_ToACalib/FoCal/Run770_toa_calib_coraseDNL.json";
    json DNL_corase_TDC_json;
    std::ifstream DNL_corase_TDC_json_file_stream(DNL_corase_TDC_json_file);
    if (!DNL_corase_TDC_json_file_stream.is_open()) {
        LOG(ERROR) << "Failed to open DNL coarse TDC json file " << DNL_corase_TDC_json_file;
        return 1;
    }
    DNL_corase_TDC_json_file_stream >> DNL_corase_TDC_json;

    std::vector<std::vector<double>> corase_DNL_factors;
    std::unordered_map<int, int> channel_to_corase_DNL_index;

    int corase_DNL_index = 0;
    for (auto& [key_str, value] : DNL_corase_TDC_json.items()) {
        // LOG(INFO) << "Channel " << key_str << " has coarse DNL index " << corase_DNL_index;
        int channel_number = std::stoi(key_str);
        channel_to_corase_DNL_index[channel_number] = corase_DNL_index;
    
        std::vector<double> factors;
    
        // Check if value is an array and all elements are null
        bool all_null = value.is_array() && std::all_of(value.begin(), value.end(), [](const json& x) {
            return x.is_null();
        });
    
        if (all_null) {
            // No valid DNL — use default neutral correction
            factors = std::vector<double>(32, 1.0);  // assuming 5-bit TDC = 32 bins
           //  LOG(WARNING) << "Channel " << channel_number << " has no valid DNL in corase. Using neutral factors (1.0).";
        } else {
            factors = value.get<std::vector<double>>();
        }
    
        corase_DNL_factors.push_back(factors);
        ++corase_DNL_index;
    }

    DNL_corase_TDC_json_file_stream.close();

    // * --- Read timewalk correction file -------------------------------------------
    // * -----------------------------------------------------------------------------
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

    // * --- Constants and variables ------------------------------------------------------
    // * ----------------------------------------------------------------------------------
    const double sample_time      = 25.0;           // unit: ns
    const double phase_shift_time = 25.0 / 16.0;    // unit: ns

    const double time_peak_value    = 72.0;
    const double time_rising_ratio  = 1.60;
    const double time_falling_ratio = 1.75;

    const double toa_shift_offset = 88.0; // unit: ns

    const int fitted_events_to_print = 10;

    const double adc_sum_slice_min = 0;
    const double adc_sum_slice_max = 20000;

    std::vector<int> target_channels = {145};

    const double event_waveform_x_min = -25.0;
    const double event_waveform_x_max = 250.0;

    const double _text_line_height = 0.04;
    const double _text_line_start  = 0.85;
    const double _text_line_left   = 0.13;

    const double _dac_step_size = 16.0;

    const double injection_template_pedestal = 60.0;
    
    bool enable_fine_DNL_correction = true;

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
            if (time_val < time_peak_value) {
                time_val *= time_rising_ratio;
            } else {
                time_val = time_peak_value * time_rising_ratio + (time_val - time_peak_value) * time_falling_ratio;
            }
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

    // * --- Read the events ---------------------------------------------------------------
    // * ----------------------------------------------------------------------------------
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

    // ! -- Loop over the events --
    long counter_total_events           = 0;
    long counter_hamming_code_error     = 0;
    long counter_bad_daqh_start_end     = 0;
    long counter_hamming_code_passed    = 0;
    long counter_waveform_shape_passed  = 0;
    long counter_adc_sum_passed         = 0;

    long counter_double_tot = 0;
    long counter_double_toa = 0;
    long counter_valid_tot  = 0;
    long counter_valid_toa  = 0;

    std::vector <std::vector< std::vector <double> >> _input_adc_time_time_list_list;   // channel < event < sample
    std::vector <std::vector< std::vector <double> >> _input_adc_time_adc_list_list;    // channel < event < sample
    std::vector <std::vector< bool >> _input_adc_time_list_list_valid;                  // channel < event
    std::vector <std::vector< double >> _input_adc_time_toa_time_list;                  // channel < event
    std::vector <std::vector< double >> _input_adc_time_tot_time_list;                  // channel < event

    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
            auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_index);
            auto _channel_valid = get_valid_fpga_channel(_channel_index);
            
            std::vector <std::vector <double> > _input_adc_time_time_list_temp;
            std::vector <std::vector <double> > _input_adc_time_adc_list_temp;
            std::vector <bool> _input_adc_time_list_temp_valid;
            std::vector <double> _input_adc_time_toa_time_list_temp;
            std::vector <double> _input_adc_time_tot_time_list_temp;
            
            _input_adc_time_time_list_temp.resize(entry_max);
            _input_adc_time_adc_list_temp.resize(entry_max);
            _input_adc_time_list_temp_valid.resize(entry_max);
            _input_adc_time_toa_time_list_temp.resize(entry_max);
            _input_adc_time_tot_time_list_temp.resize(entry_max);
            
            _input_adc_time_time_list_list.push_back(_input_adc_time_time_list_temp);
            _input_adc_time_adc_list_list.push_back(_input_adc_time_adc_list_temp);
            _input_adc_time_list_list_valid.push_back(_input_adc_time_list_temp_valid);
            _input_adc_time_toa_time_list.push_back(_input_adc_time_toa_time_list_temp);
            _input_adc_time_tot_time_list.push_back(_input_adc_time_tot_time_list_temp);
        }
    }
    
    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        counter_total_events++;

        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id    = legal_fpga_id_list[_fpga_index];
            auto _timestamp  = branch_timestamps_list[_fpga_index][0];
            auto _daqh_list  = branch_daqh_list_list[_fpga_index];
            auto _tc_list    = branch_tc_list_list[_fpga_index];
            auto _tp_list    = branch_tp_list_list[_fpga_index];
            auto _val0_list  = branch_val0_list_list[_fpga_index];
            auto _val1_list  = branch_val1_list_list[_fpga_index];
            auto _val2_list  = branch_val2_list_list[_fpga_index];

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
                    }
                    auto _daqh_first_half_byte = (_daqh >> 28) & 0xf;
                    auto _daqh_last_half_byte  = _daqh & 0xf;
                    if (_daqh_first_half_byte != 5 ||  _daqh_last_half_byte != 5) {
                        good_header_bytes = false;
                    }
                }
                _hamming_code_pass_list.push_back(_hamming_code_pass);
                if (!_hamming_code_pass){
                    skip_event_hamming_code = true;
                    break;
                }
            }
            if (skip_event_hamming_code){
                counter_hamming_code_error++;
                continue;
            }
            if (!good_header_bytes){
                counter_bad_daqh_start_end++;
                continue;
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
                if (std::find(target_channels.begin(), target_channels.end(), _unified_valid_channel_number) == target_channels.end()) {
                    _channel_val0_max_list.push_back(-1);
                    _channel_val1_max_list.push_back(-1);
                    _channel_val2_max_list.push_back(-1);
                    _channel_val0_max_index_list.push_back(-1);
                    _channel_val1_max_index_list.push_back(-1);
                    _channel_val2_max_index_list.push_back(-1);
                    continue;
                }

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
                    }
                    if (_val1 > 0){
                        if (_val1_max == -1) {
                            _val1_max = _val1;
                            _val1_max_index = _sample_index;
                            counter_valid_tot++;
                        } else
                            counter_double_tot++;
                    }
                    if (_val2 > 0) {
                        if (_val2_max == -1) {
                            _val2_max = _val2;
                            _val2_max_index = _sample_index;
                            counter_valid_toa++;
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
            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                auto _channel_valid = get_valid_fpga_channel(_channel_index);
                auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                auto _asic_index = int(_unified_valid_channel_number / (FPGA_CHANNEL_NUMBER_VALID / 2));
                auto _asic_half_index = int(_unified_valid_channel_number / (FPGA_CHANNEL_NUMBER_VALID / 4)) % 2;
                auto _channel_hist_index = unifiedToHistIndex[_unified_valid_channel_number];
                if (_channel_valid == -1){
                    continue;
                }
                if (std::find(target_channels.begin(), target_channels.end(), _unified_valid_channel_number) == target_channels.end()) {
                    continue;
                }
                auto _val0_max = _channel_val0_max_list[_channel_index];
                auto _val1_max = _channel_val1_max_list[_channel_index];
                auto _val2_max = _channel_val2_max_list[_channel_index];
                // LOG(DEBUG) << "Channel " << _unified_valid_channel_number << " val0 max: " << _val0_max << " val1 max: " << _val1_max << " val2 max: " << _val2_max;
                auto _val0_max_index = _channel_val0_max_index_list[_channel_index];
                auto _val1_max_index = _channel_val1_max_index_list[_channel_index];
                auto _val2_max_index = _channel_val2_max_index_list[_channel_index];

                // ! --- user logic        ---
                if (_val2_max > 0 && _val0_max > 0) {
                    auto _fine_DNL_index        = channel_to_fine_DNL_index[_unified_valid_channel_number];
                    auto _fine_DNL_factor       = fine_DNL_factors[_fine_DNL_index];
                    auto _val2_max_fineTDC      = _val2_max % 8;

                    auto coarseDnlIdx    = channel_to_corase_DNL_index[_unified_valid_channel_number];
                    auto coarseDnlFactor = corase_DNL_factors[coarseDnlIdx];
                    auto _val2_max_coraseTDC    = (_val2_max / 8) % 32;

                    auto _fine_DNL_factor_value = _fine_DNL_factor[_val2_max_fineTDC] * coarseDnlFactor[_val2_max_coraseTDC];

                    auto _val2_time = decode_toa_value_ns(_val2_max);
                    double _toa_value_ns = _val2_max_index * sample_time + _val2_time;
                    double _tot_value = double(_val1_max);
                    if (_tot_value >= 512) {
                        _tot_value -= 512;
                        _tot_value = _tot_value * 8.0;
                    }
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
                    } // end of toa value fix
                    // if not used for timewalk correction, use the existing correction file
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
                    } // end of timewalk correction
                    double _val0_sum = 0;
                    double _val0_pedestal;
                    std::vector <double> _val0_samples;
                    std::vector <double> _val0_time;
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        double _sample_time = _sample_index * sample_time;
                        _sample_time -= _toa_value_ns;
                        _sample_time += toa_shift_offset;
                        _val0_samples.push_back(_val0);
                        _val0_time.push_back(_sample_time);
                        _val0_sum += _val0 - _val0_list[_channel_index];
                        // LOG(DEBUG) << "Sample index: " << _sample_index << " Sample time: " << _sample_time << " toa value: " << _toa_value_ns << " toa code: " << _val2_max;
                    } // end of sample loop
                    _val0_pedestal = std::min(_val0_samples[0], _val0_samples[1]);

                    counter_hamming_code_passed++;

                    bool _event_valid_for_adc_time_analysis = true;
                    // * filter on the ADC max value
                    if (_val0_sum < adc_sum_slice_min || _val0_sum > adc_sum_slice_max){
                        _event_valid_for_adc_time_analysis = false;
                    } else {
                        counter_adc_sum_passed++;
                    }
                    // * filter on the waveform shape
                    if (_val0_samples[2] < _val0_pedestal || _val0_samples[3] < _val0_pedestal || _val0_samples[3] < _val0_samples[2] ||  _val0_samples[4] < _val0_samples[2] ||_val0_samples[4] < _val0_samples[5] || _val0_samples[5] < _val0_samples[6] || _val0_samples[6] < _val0_samples[7] || _val0_samples[7] < _val0_samples[8] || _val0_samples[8] < _val0_samples[9]) {
                        _event_valid_for_adc_time_analysis = false;
                    } else {
                        counter_waveform_shape_passed++;
                    }

                    // * Save sample data to array
                    _input_adc_time_list_list_valid[_channel_hist_index].push_back(_event_valid_for_adc_time_analysis);
                    _input_adc_time_toa_time_list[_channel_hist_index].push_back(_toa_value_ns);
                    _input_adc_time_tot_time_list[_channel_hist_index].push_back(_tot_value);

                    auto _val0_time_temp = std::vector<double>(machine_gun_samples, 0);
                    auto _val0_samples_temp = std::vector<double>(machine_gun_samples, 0);

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        _val0_time_temp[_sample_index] = _val0_time[_sample_index];
                        double _pedestal_aligned_val0 = (_val0_samples[_sample_index] - _val0_pedestal) + injection_template_pedestal;
                        if (_pedestal_aligned_val0 < 0) {
                            _pedestal_aligned_val0 = 0;
                        }
                        if (_pedestal_aligned_val0 > 1023) {
                            _pedestal_aligned_val0 = 1023;
                        }
                        _val0_samples_temp[_sample_index] = _pedestal_aligned_val0;
                    } // end of the 2nd sample loop

                    _input_adc_time_time_list_list[_channel_hist_index].push_back(_val0_time_temp);
                    _input_adc_time_adc_list_list[_channel_hist_index].push_back(_val0_samples_temp);

                } // end of valid toa and adc value
                // ! --- end of user logic ---
            } // end of 2nd channel loop
        } // end of fpga loop
    } // end of entry loop

    LOG(INFO) << "--- Event Loop Finished ----------------------------------------";
    LOG(INFO) << "Total events: " << counter_total_events;
    LOG(INFO) << "Hamming code error: " << double(counter_hamming_code_error) / double(counter_total_events) * 100.0 << "%";
    LOG(INFO) << "Bad daqh start end: " << double(counter_bad_daqh_start_end) / double(counter_total_events) * 100.0 << "%";
    LOG(INFO) << "Hamming code passed: " << double(counter_hamming_code_passed) / double(counter_total_events) * 100.0 << "%";
    LOG(INFO) << "-- Waveform shape passed: " << double(counter_waveform_shape_passed) / double(counter_hamming_code_passed) * 100.0 << "%";
    LOG(INFO) << "-- ADC sum passed: " << double(counter_adc_sum_passed) / double(counter_hamming_code_passed) * 100.0 << "%";
    LOG(INFO) << "----------------------------------------------------------------";

    input_root->Close();

    // create pdf file
    std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".pdf";
    bool pdf_first_page = true;

    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
            auto _channel_valid = get_valid_fpga_channel(_channel_index);
            auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
            auto _channel_hist_index = unifiedToHistIndex[_unified_valid_channel_number];
            if (_channel_valid == -1){
                continue;
            }
            if (std::find(target_channels.begin(), target_channels.end(), _unified_valid_channel_number) == target_channels.end()) {
                continue;
            }

            auto _events_time_array   = _input_adc_time_time_list_list[_channel_hist_index];
            auto _events_adc_array    = _input_adc_time_adc_list_list[_channel_hist_index];
            auto _events_valid        = _input_adc_time_list_list_valid[_channel_hist_index];
            auto _events_toa_time     = _input_adc_time_toa_time_list[_channel_hist_index];
            auto _events_tot_time     = _input_adc_time_tot_time_list[_channel_hist_index];
            int _counter_event_printed = 0;

            std::vector <double> _events_fit_amp_array;
            std::vector <double> _events_fit_t0_array;
            std::vector <double> _events_fit_dac_array;
            std::vector <double> _events_fit_tot_array;
            std::vector <double> _events_fit_adc_sum_array;

            _events_fit_amp_array.reserve(_events_valid.size());
            _events_fit_t0_array.reserve(_events_valid.size());
            _events_fit_dac_array.reserve(_events_valid.size());
            _events_fit_tot_array.reserve(_events_valid.size());
            _events_fit_adc_sum_array.reserve(_events_valid.size());

            for (int _event_index = 0; _event_index < _events_valid.size(); _event_index++) {
                auto _event_valid = _events_valid[_event_index];
                if (_event_valid) {
                    auto _event_time_array = _events_time_array[_event_index];
                    auto _event_adc_array  = _events_adc_array[_event_index];
                    auto _event_toa_time   = _events_toa_time[_event_index];
                    auto _event_tot_time   = _events_tot_time[_event_index];

                    auto _event_graph_error = new TGraphErrors(machine_gun_samples, _event_time_array.data(), _event_adc_array.data(), nullptr, nullptr);
                    _event_graph_error->SetTitle("");

                    // ! DO THE FITTING
                    double _event_fit_range_min = event_waveform_x_min;
                    double _event_fit_range_max = event_waveform_x_max;

                    std::vector<double> _inj_dac_scan_chi2_array;
                    for (int _inj_dac_index=0; _inj_dac_index<template_DAC_phase_results.size(); _inj_dac_index++) {
                        auto _inj_dac_samples  = template_DAC_phase_results[_inj_dac_index];
                        auto _wrapped_fit_func = [=](double *x, double *par) {
                            return fitParameterizedFunction(x, par, template_time_column_values, _inj_dac_samples);
                        };
                        auto _inj_fit_func = new TF1(("fit_" + std::to_string(_unified_valid_channel_number)).c_str(),
                            _wrapped_fit_func, _event_time_array[0], _event_time_array[machine_gun_samples - 1], 2);

                        _inj_fit_func->SetParNames("A", "t0");
                        _inj_fit_func->SetParameters(1.0, -46.0);
                        _inj_fit_func->SetParLimits(0, (double(_inj_dac_index) - 1.0) / double(_inj_dac_index), (double(_inj_dac_index) + 1.0) / double(_inj_dac_index));
                        _inj_fit_func->SetParLimits(1, -100.0, 0.0);
                        _inj_fit_func->FixParameter(1, -46.0);

                        _event_graph_error->Fit(_inj_fit_func, "RNQ");
                        double _inj_fit_chi2 = _inj_fit_func->GetChisquare();

                        _inj_dac_scan_chi2_array.push_back(_inj_fit_chi2);
                    } // end of DAC scan

                    // find the minimum chi2
                    auto _min_chi2_index = std::min_element(_inj_dac_scan_chi2_array.begin(), _inj_dac_scan_chi2_array.end()) - _inj_dac_scan_chi2_array.begin();
                    auto _min_chi2 = _inj_dac_scan_chi2_array[_min_chi2_index];
                    // LOG(INFO) << "FPGA " << _fpga_id << " Channel " << _unified_valid_channel_number << " Event " << _event_index << " Min chi2: " << _min_chi2;

                    // do the fitting with the minimum chi2
                    auto _inj_dac_samples  = template_DAC_phase_results[_min_chi2_index];
                    auto _wrapped_fit_func = [=](double *x, double *par) {
                        return fitParameterizedFunction(x, par, template_time_column_values, _inj_dac_samples);
                    };
                    auto _inj_fit_func = new TF1(("fit_" + std::to_string(_unified_valid_channel_number)).c_str(),
                        _wrapped_fit_func, _event_time_array[0], _event_time_array[machine_gun_samples - 1], 2);

                    _inj_fit_func->SetParNames("A", "t0");
                    _inj_fit_func->SetParameters(1.0, -46.0);
                    _inj_fit_func->SetParLimits(0, (double(_min_chi2_index) - 1.0) / double(_min_chi2_index), (double(_min_chi2_index) + 1.0) / double(_min_chi2_index));

                    _inj_fit_func->SetParLimits(1, -100, 0.0);
                    _inj_fit_func->FixParameter(1, -46.0);
                    _event_graph_error->Fit(_inj_fit_func, "RNQ");

                    auto _event_fit_result_A = _inj_fit_func->GetParameter(0);
                    auto _event_fit_result_t0 = _inj_fit_func->GetParameter(1);
                    auto _event_fit_result_chi2 = _inj_fit_func->GetChisquare();
                    auto _event_fit_result_ndf = _inj_fit_func->GetNDF();

                    _events_fit_amp_array.push_back(_event_fit_result_A);
                    _events_fit_t0_array.push_back(_event_fit_result_t0);
                    _events_fit_dac_array.push_back(_min_chi2_index * _dac_step_size);
                    _events_fit_tot_array.push_back(_event_tot_time);
                    double _adc_sum = 0;
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        _adc_sum += _event_adc_array[_sample_index] - injection_template_pedestal;
                    }
                    _events_fit_adc_sum_array.push_back(_adc_sum);
                    
                    if (_counter_event_printed < fitted_events_to_print) {
                        auto _event_canvas = new TCanvas(Form("FPGA %d Channel %d Event %d", _fpga_id, _unified_valid_channel_number, _event_index), Form("FPGA %d Channel %d Event %d", _fpga_id, _unified_valid_channel_number, _event_index), 800, 600);

                        auto _event_canvas_axis_dummy_hist = new TH1D(Form("FPGA %d Channel %d Event %d", _fpga_id, _unified_valid_channel_number, _event_index), Form("FPGA %d Channel %d Event %d", _fpga_id, _unified_valid_channel_number, _event_index), 100, event_waveform_x_min, event_waveform_x_max);
                        _event_canvas_axis_dummy_hist->SetTitle("");
                        _event_canvas_axis_dummy_hist->GetXaxis()->SetTitle("Time [ns]");
                        _event_canvas_axis_dummy_hist->GetYaxis()->SetTitle("ADC");
                        _event_canvas_axis_dummy_hist->GetXaxis()->SetRangeUser(-25, 250);
                        _event_canvas_axis_dummy_hist->GetYaxis()->SetRangeUser(0, 1023 + 100);
                        _event_canvas_axis_dummy_hist->SetStats(0);
                        _event_canvas_axis_dummy_hist->Draw("AXIS");

                        _event_graph_error->SetMarkerStyle(20);
                        _event_graph_error->SetMarkerSize(0.5);
                        _event_graph_error->SetMarkerColor(kRed);
                        _event_graph_error->SetLineColor(kRed);

                        _event_graph_error->Draw("P");

                        _inj_fit_func->SetLineColor(kBlue);
                        _inj_fit_func->SetLineWidth(2);
                        _inj_fit_func->Draw("same");

                        auto _event_1023_line = new TLine(event_waveform_x_min, 1023, event_waveform_x_max, 1023);
                        _event_1023_line->SetLineColor(kPink + 1);
                        _event_1023_line->SetLineStyle(2);
                        _event_1023_line->SetLineWidth(1);
                        _event_1023_line->Draw("same");

                        auto _event_latex = new TLatex();
                        _event_latex->SetTextFont(42);
                        _event_latex->SetTextAlign(32);
                        _event_latex->SetTextSize(0.03);
                        _event_latex->SetTextColor(kBlue - 6);
                        _event_latex->SetNDC();
                        _event_latex->DrawLatex(0.88, 0.87, Form("FPGA %d Channel %d", _fpga_id, _unified_valid_channel_number));
                        _event_latex->DrawLatex(0.88, 0.83, Form("Event %d", _event_index));
                        _event_latex->DrawLatex(0.88, 0.79, Form("ToA: %.2f ns, Tot: %.2f", _event_toa_time, _event_tot_time));
                        _event_latex->DrawLatex(0.88, 0.75, Form("Fit A: %.2f, t0: %.2f ns", _event_fit_result_A, _event_fit_result_t0));
                        _event_latex->DrawLatex(0.88, 0.71, Form("Chi2: %.2f / %d", _event_fit_result_chi2, _event_fit_result_ndf));
                        _event_latex->DrawLatex(0.88, 0.67, Form("Optimal DAC: %d", int(_min_chi2_index * _dac_step_size)));

                        if (pdf_first_page) {
                            _event_canvas->Print((pdf_file_name + "(").c_str(), "pdf");
                            pdf_first_page = false;
                        } else {
                            _event_canvas->Print(pdf_file_name.c_str(), "pdf");
                        }
                        output_root->cd();
                        _event_canvas->Write(Form("FPGA_%d_Channel_%d_Event_%d", _fpga_id, _unified_valid_channel_number, _event_index));

                        _event_canvas->Close();
                        delete _event_canvas;
                        delete _event_graph_error;

                        _counter_event_printed++;
                    }


                } // end of event valid
                else 
                    continue;
            } // end of event loop

            // Draw amp* DAC -- ToT correlation
            auto _channel_fit_height_tot_correlation_th2d = new TH2D(Form("FPGA %d Channel %d ToT vs Amp*DAC", _fpga_id, _unified_valid_channel_number), Form("FPGA %d Channel %d ToT vs Amp*DAC", _fpga_id, _unified_valid_channel_number), 256, 0, 1024 * 2, 256, 0, 4096);
            for (int _event_index = 0; _event_index < _events_fit_amp_array.size(); _event_index++) {
                auto _event_fit_amp = _events_fit_amp_array[_event_index];
                auto _event_fit_dac = _events_fit_dac_array[_event_index];
                auto _event_fit_tot = _events_fit_tot_array[_event_index];
                if (_event_fit_tot > 0) {
                    _channel_fit_height_tot_correlation_th2d->Fill(_event_fit_amp * _event_fit_dac, _event_fit_tot);
                }
            } // end of event loop

            auto _channel_fit_height_tot_correlation_canvas = new TCanvas(Form("FPGA %d Channel %d ToT vs Amp*DAC", _fpga_id, _unified_valid_channel_number), Form("FPGA %d Channel %d ToT vs Amp*DAC", _fpga_id, _unified_valid_channel_number), 700, 600);
            _channel_fit_height_tot_correlation_canvas->SetLeftMargin(0.15);
            _channel_fit_height_tot_correlation_canvas->SetRightMargin(0.15);
            _channel_fit_height_tot_correlation_th2d->SetTitle("");
            _channel_fit_height_tot_correlation_th2d->GetXaxis()->SetTitle("Injection Fit Amplitude * DAC");
            _channel_fit_height_tot_correlation_th2d->GetYaxis()->SetTitle("ToT");
            _channel_fit_height_tot_correlation_th2d->GetXaxis()->SetRangeUser(0, 1024 * 2);
            _channel_fit_height_tot_correlation_th2d->GetYaxis()->SetRangeUser(0, 4096);

            _channel_fit_height_tot_correlation_th2d->SetStats(0);
            _channel_fit_height_tot_correlation_th2d->Draw("COLZ");

            _channel_fit_height_tot_correlation_canvas->Print(pdf_file_name.c_str(), "pdf");

            // Draw amp* DAC -- adc sum correlation
            auto _channel_fit_height_adc_sum_correlation_th2d = new TH2D(Form("FPGA %d Channel %d ADC Sum vs Amp*DAC", _fpga_id, _unified_valid_channel_number), Form("FPGA %d Channel %d ADC Sum vs Amp*DAC", _fpga_id, _unified_valid_channel_number), 256, 0, 1024 * 2, 256, 0, 1024 * 10);
            for (int _event_index = 0; _event_index < _events_fit_amp_array.size(); _event_index++) {
                auto _event_fit_amp = _events_fit_amp_array[_event_index];
                auto _event_fit_dac = _events_fit_dac_array[_event_index];
                auto _event_fit_adc_sum = _events_fit_adc_sum_array[_event_index];
                if (_event_fit_adc_sum > 0) {
                    _channel_fit_height_adc_sum_correlation_th2d->Fill(_event_fit_amp * _event_fit_dac, _event_fit_adc_sum);
                }
            } // end of event loop

            auto _channel_fit_height_adc_sum_correlation_canvas = new TCanvas(Form("FPGA %d Channel %d ADC Sum vs Amp*DAC", _fpga_id, _unified_valid_channel_number), Form("FPGA %d Channel %d ADC Sum vs Amp*DAC", _fpga_id, _unified_valid_channel_number), 700, 600);
            _channel_fit_height_adc_sum_correlation_canvas->SetLeftMargin(0.15);
            _channel_fit_height_adc_sum_correlation_canvas->SetRightMargin(0.15);
            _channel_fit_height_adc_sum_correlation_th2d->SetTitle("");
            _channel_fit_height_adc_sum_correlation_th2d->GetXaxis()->SetTitle("Injection Fit Amplitude * DAC");
            _channel_fit_height_adc_sum_correlation_th2d->GetYaxis()->SetTitle("ADC Sum");
            _channel_fit_height_adc_sum_correlation_th2d->GetXaxis()->SetRangeUser(0, 1024 * 2);
            _channel_fit_height_adc_sum_correlation_th2d->GetYaxis()->SetRangeUser(0, 1024 * 10);

            _channel_fit_height_adc_sum_correlation_th2d->SetStats(0);
            _channel_fit_height_adc_sum_correlation_th2d->Draw("COLZ");

            _channel_fit_height_adc_sum_correlation_canvas->Print(pdf_file_name.c_str(), "pdf");
        } // end of channel loop
    } // end of fpga loop

    // add dummy canvas to pdf
    auto dummy_canvas = new TCanvas("dummy_canvas", "dummy_canvas", 800, 600);
    dummy_canvas->Print((pdf_file_name + ")").c_str(), "pdf");
    output_root->Close();
    return 0;
}