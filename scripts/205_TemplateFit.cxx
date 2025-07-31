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

Double_t HalfGaussianFlat(Double_t *x, Double_t *par) {
    Double_t x0 = par[0];  // Central transition point
    Double_t x1 = par[1];  // Width of the flat region
    Double_t sigma = par[2];  // Standard deviation of the Gaussian
    Double_t amplitude = par[3];  // Amplitude of the Gaussian
    
    Double_t xval = x[0];

    if (xval < x0) {  // Left half-Gaussian
        return amplitude * TMath::Exp(-0.5 * TMath::Power((xval - x0) / sigma, 2));
    } 
    else if (xval >= x0 && xval <= x0 + x1) {  // Flat region
        return amplitude;
    } 
    else {  // Right half-Gaussian
        return amplitude * TMath::Exp(-0.5 * TMath::Power((xval - (x0 + x1)) / sigma, 2));
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
    ScriptOptions opts = parse_arguments_single_root_single_csv(argc, argv, "1.2");

    ROOT::EnableImplicitMT(30);
    gROOT->SetBatch(kTRUE);

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
    const int toa_time_hist_bins    = int((toa_time_hist_max - toa_time_hist_min)/phase_shift_time) * 2;

    const double toa_code_hist_min  = 0.0;
    const double toa_code_hist_max  = 1024.0;
    const int toa_code_hist_bins    = 256;

    const double time_peak_value    = 72.0;
    const double time_rising_ratio  = 1.60;
    const double time_falling_ratio = 1.75;

    const int example_fitting_channel = 145;
    const int module_to_check = 5;

    const double adc_sum_slice_min = 4000;
    const double adc_sum_slice_max = 4500;

    const double _text_line_height = 0.04;
    const double _text_line_start  = 0.85;
    const double _text_line_left   = 0.13;
    
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

    // for (int i = 0; i < template_DAC_values.size(); i++) {
    //     LOG(INFO) << "DAC " << template_DAC_values[i] << " has " << template_DAC_phase_results[i].size() << " entries";
    // }

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
    std::vector <std::vector<std::vector<double> >> channel_toa_time_toa_code_DNL_sums; // channel < toa code < toa time
    std::vector <std::vector<std::vector<double> >> channel_toa_time_toa_code_DNL_sums_count;
    TDirectory *channel_toa_time_toa_code_folder = output_root->mkdir("ChannelToaTimeToaCode");

    std::vector <TH2D*> channel_toa_time_adc_max_hist_list;
    std::vector <std::vector< std::vector<double> >> channel_toa_time_adc_max_DNL_sums; // channel < adc max < toa time
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
            std::vector <std::vector<double>> _toa_time_toa_code_DNL_sums_temp(toa_time_hist_bins, std::vector<double>(toa_code_hist_bins, 0));
            channel_toa_time_toa_code_DNL_sums.push_back(_toa_time_toa_code_DNL_sums_temp);
            std::vector <std::vector<double>> _toa_time_toa_code_DNL_sums_count_temp(toa_time_hist_bins, std::vector<double>(toa_code_hist_bins, 0));
            channel_toa_time_toa_code_DNL_sums_count.push_back(_toa_time_toa_code_DNL_sums_count_temp);
            _toa_time_toa_code_hist->SetDirectory(channel_toa_time_toa_code_folder);

            // -- Create the channel TOA time and ADC max histograms
            // ----------------------------------------------------------------
            auto *_toa_time_adc_max_hist = new TH2D(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa_time_adc_max").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA Time and ADC Max").c_str(), channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max, toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max);
            channel_toa_time_adc_max_hist_list.push_back(_toa_time_adc_max_hist);
            std::vector <std::vector<double>> _toa_time_adc_max_DNL_sums_temp(channel_adc_hist_bins, std::vector<double>(toa_time_hist_bins, 0));
            channel_toa_time_adc_max_DNL_sums.push_back(_toa_time_adc_max_DNL_sums_temp);
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

    std::vector <std::vector< std::vector <double> >> _input_adc_time_list_list;        // channel < event < sample
    std::vector <std::vector< std::vector <double> >> _input_adc_time_list_list_temp;   // channel < event < sample
    std::vector <std::vector< bool >> _input_adc_time_list_list_valid;                  // channel < event
    std::vector <std::vector< double >> _input_adc_time_toa_time_list;             // channel < event

    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
            auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_index);
            auto _channel_valid = get_valid_fpga_channel(_channel_index);
            
            std::vector <std::vector <double> > _input_adc_time_list_temp;
            std::vector <std::vector <double> > _input_adc_time_list_temp_temp;
            std::vector <bool> _input_adc_time_list_temp_valid;
            std::vector <double> _input_adc_time_toa_time_list_temp;
            
            _input_adc_time_list_temp.resize(entry_max);
            _input_adc_time_list_temp_temp.resize(entry_max);
            _input_adc_time_list_temp_valid.resize(entry_max);
            _input_adc_time_toa_time_list_temp.resize(entry_max);
            
            _input_adc_time_list_list.push_back(_input_adc_time_list_temp);
            _input_adc_time_list_list_temp.push_back(_input_adc_time_list_temp_temp);
            _input_adc_time_list_list_valid.push_back(_input_adc_time_list_temp_valid);
            _input_adc_time_toa_time_list.push_back(_input_adc_time_toa_time_list_temp);
        }
    }

    long counter_total_events = 0;
    long counter_hamming_code_error = 0;
    long counter_bad_daqh_start_end = 0;
    long counter_hamming_code_passed = 0;
    long counter_waveform_shape_passed = 0;
    long counter_adc_max_passed = 0;

    LOG(INFO) << "Start to read the events";
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
                
                if (_val2_max > 0 && _val0_max > 0) {
                    // LOG(DEBUG) << "Channel " << _unified_valid_channel_number << " has TOA value " << _val2_max;
                    auto _fine_DNL_index        = channel_to_fine_DNL_index[_unified_valid_channel_number];
                    // LOG(DEBUG) << "Channel " << _unified_valid_channel_number << " has fine DNL index " << _fine_DNL_index;
                    auto _fine_DNL_factor       = fine_DNL_factors[_fine_DNL_index];
                    auto _val2_max_fineTDC      = _val2_max % 8;

                    auto _fine_DNL_corase_index = channel_to_corase_DNL_index[_unified_valid_channel_number];
                    auto _fine_DNL_corase_factor = corase_DNL_factors[_fine_DNL_corase_index];
                    auto _val2_max_coraseTDC    = (_val2_max / 8) % 32;

                    auto _fine_DNL_factor_value = _fine_DNL_factor[_val2_max_fineTDC] * _fine_DNL_corase_factor[_val2_max_coraseTDC];
                    // LOG(DEBUG) << "Channel " << _unified_valid_channel_number << " has fine DNL factor " << _fine_DNL_factor_value;

                    auto _val2_time = decode_toa_value_ns(_val2_max);
                    double _toa_value_ns = _val2_max_index * sample_time + _val2_time;
                    // ! -- Manually fix the TOA value
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
                    int _hist_nx = channel_toa_time_toa_code_hist_list[_channel_hist_index]->GetXaxis()->FindBin(_toa_value_ns);
                    int _hist_ny = channel_toa_time_toa_code_hist_list[_channel_hist_index]->GetYaxis()->FindBin(_val2_max);
                    if (_hist_nx - 1 >= 0 && _hist_ny - 1 >= 0 && _hist_nx - 1 < toa_time_hist_bins && _hist_ny - 1 < toa_code_hist_bins) {
                        channel_toa_time_toa_code_DNL_sums[_channel_hist_index][_hist_nx - 1][_hist_ny - 1] += _fine_DNL_factor_value;
                        channel_toa_time_toa_code_DNL_sums_count[_channel_hist_index][_hist_nx - 1][_hist_ny - 1] += 1;
                    } else {
                        // LOG(WARNING) << "Channel " << _unified_valid_channel_number << " has TOA value " << _toa_value_ns << " TOA code " << _val2_max << " hist_nx " << _hist_nx << " hist_ny " << _hist_ny;
                    }
                    // channel_toa_time_toa_code_DNL_sums[_channel_hist_index][_hist_nx - 1][_hist_ny - 1] += _fine_DNL_factor_value;

                    channel_toa_time_adc_max_hist_list[_channel_hist_index]->Fill(_val0_max, _toa_value_ns);
                    int _hist_nx_adc_max = channel_toa_time_adc_max_hist_list[_channel_hist_index]->GetXaxis()->FindBin(_val0_max);
                    int _hist_ny_adc_max = channel_toa_time_adc_max_hist_list[_channel_hist_index]->GetYaxis()->FindBin(_toa_value_ns);
                    if (_hist_nx_adc_max - 1 >= 0 && _hist_ny_adc_max - 1 >= 0 && _hist_nx_adc_max - 1 < channel_adc_hist_bins && _hist_ny_adc_max - 1 < toa_time_hist_bins) {
                        channel_toa_time_adc_max_DNL_sums[_channel_hist_index][_hist_nx_adc_max - 1][_hist_ny_adc_max - 1] += _fine_DNL_factor_value;
                    } else {
                        // LOG(WARNING) << "Channel " << _unified_valid_channel_number << " has ADC max " << _val0_max << " TOA value " << _toa_value_ns << " hist_nx " << _hist_nx_adc_max << " hist_ny " << _hist_ny_adc_max;
                    }
                    // channel_toa_time_adc_max_DNL_sums[_channel_hist_index][_hist_nx_adc_max - 1][_hist_ny_adc_max - 1] += _fine_DNL_factor_value;

                    // get the adc samples
                    std::vector <double> _val0_samples;
                    std::vector <double> _val0_time;
                    double _val0_sum = 0;
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        double _sample_time = _sample_index * sample_time;
                        _sample_time -= (_toa_value_ns - 88.0);
                        _val0_samples.push_back(_val0);
                        _val0_time.push_back(_sample_time);
                        _val0_sum += _val0 - _val0_list[_channel_index];
                    } // end of the 2nd sample loop
                    auto _val0_min = std::min(_val0_samples[0], _val0_samples[1]);

                    counter_hamming_code_passed++;

                    bool _event_valid_for_adc_time_analysis = true;
                    // * filter on the ADC max value
                    if (_val0_sum < adc_sum_slice_min || _val0_sum > adc_sum_slice_max){
                        _event_valid_for_adc_time_analysis = false;
                    } else {
                        counter_adc_max_passed++;
                    }
                    // * filter on the waveform shape
                    if (_val0_samples[2] < _val0_min || _val0_samples[3] < _val0_min || _val0_samples[3] < _val0_samples[2] ||  _val0_samples[4] < _val0_samples[2] ||_val0_samples[4] < _val0_samples[5] || _val0_samples[5] < _val0_samples[6] || _val0_samples[6] < _val0_samples[7] || _val0_samples[7] < _val0_samples[8] || _val0_samples[8] < _val0_samples[9]) {
                        _event_valid_for_adc_time_analysis = false;
                    } else {
                        counter_waveform_shape_passed++;
                    }
                    _input_adc_time_list_list_valid[_channel_hist_index].push_back(_event_valid_for_adc_time_analysis);
                    _input_adc_time_toa_time_list[_channel_hist_index].push_back(_toa_value_ns);

                    auto _val0_time_temp = std::vector<double>(machine_gun_samples, 0);
                    auto _val0_samples_temp = std::vector<double>(machine_gun_samples, 0);

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        _val0_time_temp[_sample_index] = _val0_time[_sample_index];
                        _val0_samples_temp[_sample_index] = (_val0_samples[_sample_index] - _val0_min) + 60;
                    } // end of the 2nd sample loop

                    _input_adc_time_list_list[_channel_hist_index].push_back(_val0_time_temp);
                    _input_adc_time_list_list_temp[_channel_hist_index].push_back(_val0_samples_temp);
                } // end of if toa max > 0
            } // end of the 2nd channel loop
        } // end of fpga loop
    } // end of event loop

    LOG(INFO) << "Finish reading the events";
    LOG(INFO) << "---------------------------------------------------";

    LOG(INFO) << "Total events: " << counter_total_events;
    LOG(INFO) << "Hamming code error rate: " << double(counter_hamming_code_error) / double(counter_total_events) * 100.0 << "%";
    LOG(INFO) << "Bad daqh start end rate: " << double(counter_bad_daqh_start_end) / double(counter_total_events) * 100.0 << "%";

    LOG(INFO) << "Waveform filtering pass rate: " << double(counter_waveform_shape_passed) / double(counter_hamming_code_passed) * 100.0 << "%";
    LOG(INFO) << "ADC sum slice percentage    : " << double(counter_adc_max_passed) / double(counter_hamming_code_passed) * 100.0 << "%";

    LOG(INFO) << "---------------------------------------------------";

    std::vector <TCanvas*> canvas_channel_toa_time_toa_code_list;
    std::vector <TCanvas*> canvas_channel_toa_time_projection_list;

    std::vector <double> channel_toa_time_range_left;
    std::vector <double> channel_toa_time_range_right;
    

    TCanvas *canvas_example_toa_distribution = new TCanvas("ExampleToaDistribution", "Example TOA Distribution", 800, 600);
    for (int i = 0; i < channel_toa_time_toa_code_hist_list.size(); i++) {
        channel_toa_time_range_left.push_back(0);
        channel_toa_time_range_right.push_back(0);
        auto _canvas = new TCanvas(("ChannelTOATimeTOACode_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time and TOA Code " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        channel_toa_time_toa_code_hist_list[i]->SetStats(0);
        // normalize the histogram to 1
        // channel_toa_time_toa_code_hist_list[i]->Scale(1.0 / channel_toa_time_toa_code_hist_list[i]->GetMaximum());
        if (enable_fine_DNL_correction){
            int _nx = channel_toa_time_toa_code_hist_list[i]->GetNbinsX();
            int _ny = channel_toa_time_toa_code_hist_list[i]->GetNbinsY();
            for (int _bin_x = 1; _bin_x <= _nx; _bin_x++) {
                for (int _bin_y = 1; _bin_y <= _ny; _bin_y++) {
                    auto _fine_DNL_factor = channel_toa_time_toa_code_DNL_sums[i][_bin_x - 1][_bin_y - 1];
                    auto _fine_DNL_count = channel_toa_time_toa_code_DNL_sums_count[i][_bin_x - 1][_bin_y - 1];
                    if (_fine_DNL_factor > 0) {
                        auto _bin = channel_toa_time_toa_code_hist_list[i]->GetBin(_bin_x, _bin_y);
                        auto _bin_content = channel_toa_time_toa_code_hist_list[i]->GetBinContent(_bin);
                        auto _bin_corrected_content = _bin_content / _fine_DNL_factor * _fine_DNL_count;
                        channel_toa_time_toa_code_hist_list[i]->SetBinContent(_bin, _bin_corrected_content);
                    }
                }
            }
        }

        channel_toa_time_toa_code_hist_list[i]->Draw("COLZ");
        _canvas->Update();
        channel_toa_time_toa_code_folder->cd();
        channel_toa_time_toa_code_hist_list[i]->Write();
        
        canvas_channel_toa_time_toa_code_list.push_back(_canvas);

        auto _canvas_projection = new TCanvas(("ChannelTOATimeProjection_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time Projection " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);
        auto _channel_toa_time_toa_code_projection = channel_toa_time_toa_code_hist_list[i]->ProjectionX();
        _channel_toa_time_toa_code_projection->SetStats(0);
        _channel_toa_time_toa_code_projection->SetTitle("");
        _channel_toa_time_toa_code_projection->Draw("HIST");

        canvas_channel_toa_time_projection_list.push_back(_canvas_projection);


        // do the fitting
        if (true) {
            // LOG(INFO) << "Channel " << channel_index_to_unified_channel_number_map[i] << " has " << _channel_toa_time_toa_code_projection->GetEntries() << " entries";
            if  (_channel_toa_time_toa_code_projection->GetEntries() < 100) {
                continue;
            }

            // do gaussian fit over left edge
            double _fit_x_min = 0.0;
            double _fit_x_max = 200.0;

            // use the formula double HalfGaussianFlat
            TF1* fit_half_gaussian_flat = new TF1("fit_half_gaussian_flat", HalfGaussianFlat, _fit_x_min, _fit_x_max, 4);
            double _projection_x_mid = 0;
            double _projection_y_max = _channel_toa_time_toa_code_projection->GetMaximum();
            // get the average x by averaging the cross-half-points
            double _projection_x_cross_half_left = 0;
            double _projection_x_cross_half_right = 0;
            for (int j = 1; j < _channel_toa_time_toa_code_projection->GetNbinsX(); j++) {
                auto _bin_content = _channel_toa_time_toa_code_projection->GetBinContent(j);
                if (_bin_content > _projection_y_max / 4.0 && _projection_x_cross_half_left == 0) {
                    _projection_x_cross_half_left = _channel_toa_time_toa_code_projection->GetBinCenter(j);
                }
                if (_bin_content < _projection_y_max / 4.0 && _projection_x_cross_half_left != 0) {
                    _projection_x_cross_half_right = _channel_toa_time_toa_code_projection->GetBinCenter(j);
                    break;
                }
            }
            _projection_x_mid = (_projection_x_cross_half_left + _projection_x_cross_half_right) / 2.0;
            // draw lines
            auto _line = new TLine(_projection_x_mid, 0, _projection_x_mid, _projection_y_max);
            _line->SetLineColor(kRed - 9);
            _line->SetLineWidth(1);
            _line->SetLineStyle(3);
            _line->Draw("SAME");

            fit_half_gaussian_flat->SetParameters(_projection_x_mid-12.5, 25.0, 2.0, 0.8*_projection_y_max);
            fit_half_gaussian_flat->SetParLimits(0, _projection_x_mid - 42.5, _projection_x_mid + 17.5);
            fit_half_gaussian_flat->SetParLimits(1, 25.0, 30.0);
            fit_half_gaussian_flat->SetParLimits(2, 0.1, 10.0);
            fit_half_gaussian_flat->SetParLimits(3, 0.2*_projection_y_max, _projection_y_max);
            _channel_toa_time_toa_code_projection->Fit(fit_half_gaussian_flat, "RQN");

            double _fit_half_gaussian_flat_central = fit_half_gaussian_flat->GetParameter(0);
            double _fit_half_gaussian_flat_width = fit_half_gaussian_flat->GetParameter(1);
            double _fit_half_gaussian_flat_sigma = fit_half_gaussian_flat->GetParameter(2);
            double _fit_half_gaussian_flat_amplitude = fit_half_gaussian_flat->GetParameter(3);

            double _time_range_left  = _fit_half_gaussian_flat_central - 2.0 * _fit_half_gaussian_flat_sigma;
            double _time_range_right = _fit_half_gaussian_flat_central + 2.0 * _fit_half_gaussian_flat_sigma + _fit_half_gaussian_flat_width;

            channel_toa_time_range_left[i]  = _time_range_left;
            channel_toa_time_range_right[i] = _time_range_right;

            auto fit_half_gaussian_flat_latex = new TLatex();
            fit_half_gaussian_flat_latex->SetTextSize(0.04);
            fit_half_gaussian_flat_latex->SetTextColor(kBlack);
            fit_half_gaussian_flat_latex->SetTextAlign(12);
            fit_half_gaussian_flat_latex->SetNDC();

            fit_half_gaussian_flat_latex->DrawLatex(0.05, 0.88, Form("Channel %d", channel_index_to_unified_channel_number_map[i]));
            fit_half_gaussian_flat_latex->DrawLatex(0.05, 0.83, Form("Fit Central: %.2f ns", _fit_half_gaussian_flat_central));
            fit_half_gaussian_flat_latex->DrawLatex(0.05, 0.78, Form("Fit Width: %.2f ns", _fit_half_gaussian_flat_width));
            fit_half_gaussian_flat_latex->DrawLatex(0.05, 0.73, Form("Fit Sigma: %.2f ns", _fit_half_gaussian_flat_sigma));
            fit_half_gaussian_flat_latex->DrawLatex(0.05, 0.68, Form("Fit Amplitude: %.2f", _fit_half_gaussian_flat_amplitude));

            fit_half_gaussian_flat->SetLineColor(kRed);
            fit_half_gaussian_flat->SetLineWidth(1);
            fit_half_gaussian_flat->SetLineStyle(2);
            fit_half_gaussian_flat->Draw("SAME");

            if (channel_index_to_unified_channel_number_map[i] == example_fitting_channel) {
                canvas_example_toa_distribution->cd();
                auto example_toa_legend = new TLegend(0.65, 0.7, 0.89, 0.89);
                example_toa_legend->SetBorderSize(0);
                example_toa_legend->SetTextSize(0.03);
                example_toa_legend->SetTextFont(42);
                example_toa_legend->SetFillColor(0);
                example_toa_legend->SetFillStyle(0);

                TH1F* projection_clone = (TH1F*)_channel_toa_time_toa_code_projection->Clone("ExampleToaDistributionProjection");
                projection_clone->SetName("ExampleToaDistributionProjection");
                projection_clone->SetTitle("");
                projection_clone->GetXaxis()->SetTitle("TOA Time [ns]");
                projection_clone->GetYaxis()->SetTitle("Counts");
                projection_clone->SetStats(0);
    
                projection_clone->SetLineColor(kCyan + 3);
                projection_clone->SetLineWidth(2);
    
                projection_clone->Draw("HIST");
    
                TF1* fit_half_gaussian_flat_clone = (TF1*)fit_half_gaussian_flat->Clone("fit_half_gaussian_flat_clone");
                fit_half_gaussian_flat_clone->SetName("fit_half_gaussian_flat_clone");
                fit_half_gaussian_flat_clone->SetLineColor(kRed-4);
                fit_half_gaussian_flat_clone->SetLineWidth(2);

                // projection_clone->Fit(fit_half_gaussian_flat_clone, "RQN");

                fit_half_gaussian_flat_clone->SetLineWidth(2);
                fit_half_gaussian_flat_clone->SetLineStyle(2);
                fit_half_gaussian_flat_clone->Draw("SAME");
                example_toa_legend->AddEntry(projection_clone, "TOA Distribution", "l");
                

                auto fit_half_gaussian_result_central = fit_half_gaussian_flat_clone->GetParameter(0);
                auto fit_half_gaussian_result_width = fit_half_gaussian_flat_clone->GetParameter(1);
                auto fit_half_gaussian_result_sigma = fit_half_gaussian_flat_clone->GetParameter(2);
                auto fit_half_gaussian_result_amplitude = fit_half_gaussian_flat_clone->GetParameter(3);
                example_toa_legend->AddEntry(fit_half_gaussian_flat_clone, Form("Fit sigma: %.2f ns", fit_half_gaussian_result_sigma), "l");

                example_toa_legend->Draw("SAME");

                auto _toa_distribution_latex = new TLatex();
                _toa_distribution_latex->SetNDC();
                _toa_distribution_latex->SetTextSize(0.04);
                _toa_distribution_latex->SetTextFont(62);
                _toa_distribution_latex->DrawLatex(_text_line_left, _text_line_start, "FoCal-H Prototype II");
                _toa_distribution_latex->SetTextSize(0.03);
                _toa_distribution_latex->SetTextFont(42);
                _toa_distribution_latex->DrawLatex(_text_line_left, _text_line_start - _text_line_height, Form("Channel %d", channel_index_to_unified_channel_number_map[i]));
                _toa_distribution_latex->DrawLatex(_text_line_left, _text_line_start - 2 * _text_line_height, Form("Single Channel ToA Distribution"));
                _toa_distribution_latex->DrawLatex(_text_line_left, _text_line_start - 3 * _text_line_height, "September 2024");
                _toa_distribution_latex->SetTextFont(52);
                _toa_distribution_latex->SetTextColor(kGray + 3);
                _toa_distribution_latex->DrawLatex(_text_line_left, _text_line_start - 4* _text_line_height, "Work in Progress");

            } // end of if example channel
        } // end of if fit
        
        _canvas_projection->Update();
    }

    long counter_toa_time_passed = 0;
    long counter_toa_time_total = 0;

    // ! Fill the data to channel_adc_samples_hist_list
    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
            auto _channel_valid = get_valid_fpga_channel(_channel_index);
            auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
            auto _channel_hist_index = channel_unified_channel_number_to_index_map[_unified_valid_channel_number];
            // LOG(DEBUG) << "Channel " << _channel_hist_index;
            auto _toa_range_left = channel_toa_time_range_left[_channel_hist_index];
            auto _toa_range_right = channel_toa_time_range_right[_channel_hist_index];
            // LOG(DEBUG) << "TOA range left: " << _toa_range_left << " right: " << _toa_range_right;

            if (_channel_valid == -1){
                continue;
            }

            auto _val0_time = _input_adc_time_list_list[_channel_hist_index];
            auto _val0_samples = _input_adc_time_list_list_temp[_channel_hist_index];
            auto _val0_valid = _input_adc_time_list_list_valid[_channel_hist_index];
            auto _val0_toa_time = _input_adc_time_toa_time_list[_channel_hist_index];

            for (int _event_index = 0; _event_index < _val0_samples.size(); _event_index++) {
                auto _val0_samples_temp = _val0_samples[_event_index];
                auto _val0_time_temp = _val0_time[_event_index];
                auto _val0_toa_time_temp = _val0_toa_time[_event_index];
                auto _val0_valid_temp = _val0_valid[_event_index];

                if (_val0_valid_temp) {
                    counter_toa_time_total++;
                    if (_val0_toa_time_temp < _toa_range_left || _val0_toa_time_temp > _toa_range_right) {
                        continue;
                    } else {
                        counter_toa_time_passed++;
                    }
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        channel_adc_samples_hist_list[_channel_hist_index]->Fill(_val0_time_temp[_sample_index], _val0_samples_temp[_sample_index]);
                    }
                }
            }
        }
    }

    LOG(INFO) << "Finish data sample filling";
    LOG(INFO) << "---------------------------------------------------";

    LOG(INFO) << "ToA filtering pass rate: " << double(counter_toa_time_passed) / double(counter_toa_time_total) * 100.0 << "%";

    LOG(INFO) << "---------------------------------------------------";
        
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
    // Save only the central module histograms
    auto painter_canvas_toa_time_toa_code = global_painter->draw_module_channel_canvas(canvas_channel_toa_time_toa_code_list, channel_unified_channel_number_to_index_map, "ChannelTOATimeTOACode", "Channel TOA Time and TOA Code", module_to_check);
    if (painter_canvas_toa_time_toa_code == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_toa_time_toa_code->Write();
    painter_canvas_toa_time_toa_code->Print((pdf_file_name + "(").c_str());

    canvas_example_toa_distribution->Write();

    canvas_example_toa_distribution->Print(pdf_file_name.c_str());

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

    LOG(INFO) << "Hist Array Size: " << channel_toa_time_adc_max_hist_list.size();
    for (int i = 0; i < channel_toa_time_adc_max_hist_list.size(); i++) {
        auto _canvas = new TCanvas(("ChannelTOATimeADCMax_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel TOA Time and ADC Max " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        channel_toa_time_adc_max_hist_list[i]->SetStats(0);
        if (enable_fine_DNL_correction) {
            int _nx = channel_toa_time_adc_max_hist_list[i]->GetNbinsX();
            int _ny = channel_toa_time_adc_max_hist_list[i]->GetNbinsY();
            for (int _bin_x = 1; _bin_x <= _nx; _bin_x++) {
                for (int _bin_y = 1; _bin_y <= _ny; _bin_y++) {
                    auto _fine_DNL_factor = channel_toa_time_adc_max_DNL_sums[i][_bin_x - 1][_bin_y - 1];
                    if (_fine_DNL_factor > 0) {
                        auto _bin = channel_toa_time_adc_max_hist_list[i]->GetBin(_bin_x, _bin_y);
                        auto _bin_content = channel_toa_time_adc_max_hist_list[i]->GetBinContent(_bin);
                        auto _bin_corrected_content = _bin_content / _fine_DNL_factor;
                        channel_toa_time_adc_max_hist_list[i]->SetBinContent(_bin, _bin_corrected_content);
                    }
                }
            }
        }


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
                if (_xCentre < channel_adc_hist_min || _xCentre > channel_adc_hist_max) {
                    continue;
                }
                if (_bin_mean < 50 || _bin_mean > 170) {
                    continue;
                }
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
    _graph_threshold->Draw("AP same");
    _canvas_threshold->Update();

    auto _toa_theshold_latex = new TLatex();
    _toa_theshold_latex->SetNDC();
    _toa_theshold_latex->SetTextSize(0.04);
    _toa_theshold_latex->SetTextFont(62);
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
    std::vector <TF1*> template_DAC_fit_funcs;

    std::vector <double> template_minimal_index_list;
    std::vector <double> template_max_bin_list;

    std::vector <double> template_fit_result_A;
    std::vector <double> template_fit_result_x0;
    std::vector <int> template_fit_DAC_index;
    std::vector <int> template_fit_unified_channel_numbers;
    std::vector <double> template_fit_result_chi2;
    std::vector <double> template_fit_result_ndf;

    TGraphErrors *example_average_waveform = nullptr;
    TH2D *example_waveform_hist2d = nullptr;
    TF1 *example_fit_func = nullptr;

    for (int i = 0; i < channel_adc_samples_hist_list.size(); i++) {
        auto _canvas = new TCanvas(("ChannelADCSamples_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel ADC Samples " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);

        bool _is_example_channel = example_fitting_channel == channel_index_to_unified_channel_number_map[i];

        channel_adc_samples_hist_list[i]->SetStats(0);
        channel_adc_samples_hist_list[i]->Draw("COLZ");
        _canvas->Update();
        channel_adc_smaples_folder->cd();
        channel_adc_samples_hist_list[i]->Write();
        auto channel_event_number = channel_adc_samples_hist_list[i]->GetEntries();
        if (_is_example_channel) {
            example_waveform_hist2d = (TH2D*)channel_adc_samples_hist_list[i]->Clone();
        }
        
        canvas_channel_adc_samples_list.push_back(_canvas);

        auto _canvas_fit = new TCanvas(("ChannelADCSamplesFitted_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), ("Channel ADC Samples Fitted " + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(), 800, 600);
        auto _canvas_fit_legend = new TLegend(0.45, 0.6, 0.89, 0.89);
        _canvas_fit_legend->SetFillColor(kWhite);
        _canvas_fit_legend->SetLineColor(kWhite);
        _canvas_fit_legend->SetShadowColor(kWhite);
        _canvas_fit_legend->SetBorderSize(0);
        _canvas_fit_legend->SetTextSize(0.05);

        std::vector <double> _bin_mean_list;
        std::vector <double> _bin_error_y_list;
        std::vector <double> _bin_x_center_list;
        std::vector <double> _bin_error_x_list;

        double _max_y_value = -1.0;
        double _peak_position = -1.0;

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
                if (_bin_mean > _max_y_value) {
                    _max_y_value = _bin_mean;
                    _peak_position = _xCentre;
                }
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
        
        _canvas_fit->Update();

        // * Try doing the fitting here
        
        _canvas_fit->cd();
        // channel_adc_samples_hist_list[i]->Draw("scat");
        _graph->Draw("AP");
        _canvas_fit_legend->AddEntry(_graph, "Data", "p");
        // Draw the line
        TLine *_line = new TLine(0, _max_y_value, machine_gun_samples * sample_time, _max_y_value);
        _line->SetLineColor(kRed);
        _line->SetLineWidth(1);
        _line->SetLineStyle(2);
        _line->Draw("same");

        double _initial_x0 = _peak_position - (time_peak_value * time_rising_ratio);
        
        if (channel_event_number > 100) {
            std::vector<double> dac_fitting_chi2_list;
            for (int _dac_index = 0; _dac_index < template_DAC_phase_results.size(); _dac_index++) {
                auto _dac_sample = template_DAC_phase_results[_dac_index];
                // Capture by value rather than by reference:
                auto wrappedFitFunc = [=](double *x, double *par) {
                    return fitParameterizedFunction(x, par, template_time_column_values, _dac_sample);
                };
                TF1* _average_fit = new TF1(("fit_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(),
                                             wrappedFitFunc,
                                             template_time_column_values.front(),
                                             template_time_column_values.back(), 2);
                
                _average_fit->SetParameters(1.0, _initial_x0);
                // _average_fit->SetParLimits(0, 0.75, 1.25);
                _average_fit->FixParameter(0, 1.0);
                _average_fit->SetParLimits(1, -80.0, -20.0);
                // _average_fit->FixParameter(1, -44.0);
                
                _graph->Fit(_average_fit, "MULTITHREAD RNQ");

                template_DAC_fit_funcs.push_back(_average_fit);

                dac_fitting_chi2_list.push_back(_average_fit->GetChisquare());
                // Removed deletion code to avoid premature de-allocation:
                // _average_fit->SetParent(nullptr);
                _average_fit->Delete();
            }

            // LOG(DEBUG) << "Channel: " << channel_index_to_unified_channel_number_map[i] << " DAC Fitting Chi2 List: ";
            // for (auto _chi2 : dac_fitting_chi2_list) {
            //     LOG(DEBUG) << _chi2;
            // }

            // find the minimum chi2 value
            auto _min_chi2_index = std::min_element(dac_fitting_chi2_list.begin(), dac_fitting_chi2_list.end()) - dac_fitting_chi2_list.begin();
            auto _dac_sample = template_DAC_phase_results[_min_chi2_index];
            template_minimal_index_list.push_back(_min_chi2_index);
            template_max_bin_list.push_back(_max_y_value);
            // Capture by value:
            auto wrappedFitFunc = [=](double *x, double *par) {
                return fitParameterizedFunction(x, par, template_time_column_values, _dac_sample);
            };
            
            TF1* _average_fit_opt = new TF1(("fit_" + std::to_string(channel_index_to_unified_channel_number_map[i])).c_str(),
                                             wrappedFitFunc,
                                             template_time_column_values.front(),
                                             template_time_column_values.back(), 2);
            _average_fit_opt->SetParameters(1.0, -50.0);
            double amp_ratio_min = (double(_min_chi2_index) - 1.0) / double(_min_chi2_index);
            double amp_ratio_max = (double(_min_chi2_index) + 1.0) / double(_min_chi2_index);
            _average_fit_opt->SetParLimits(0, amp_ratio_min, amp_ratio_max);
            //  _average_fit_opt->SetParLimits(0, 0.75, 1.25);
            _average_fit_opt->SetParLimits(1, -80.0, -20.0);
            // _average_fit_opt->FixParameter(1, -44.0);
            _average_fit_opt->SetLineColor(kRed + 3);
            _average_fit_opt->SetLineWidth(2);
            _average_fit_opt->SetLineStyle(1);
            
            _graph->Fit(_average_fit_opt, "RNQ");
            _average_fit_opt->Draw("same");

            auto _fit_A = _average_fit_opt->GetParameter(0);
            auto _fit_x0 = _average_fit_opt->GetParameter(1);
            auto _fit_chi2 = _average_fit_opt->GetChisquare();
            auto _fit_ndf = _average_fit_opt->GetNDF();

            template_fit_result_A.push_back(_fit_A);
            template_fit_result_x0.push_back(_fit_x0);
            template_fit_DAC_index.push_back(_min_chi2_index);
            template_fit_unified_channel_numbers.push_back(channel_index_to_unified_channel_number_map[i]);
            template_fit_result_chi2.push_back(_fit_chi2);
            template_fit_result_ndf.push_back(_fit_ndf);

            _canvas_fit_legend->AddEntry(_average_fit_opt, ("Fit: A = " + std::to_string(_fit_A).substr(0, 5) + ", t_{offset} = " + std::to_string(_fit_x0).substr(0, 5)).c_str(), "l");
            _canvas_fit_legend->AddEntry(_average_fit_opt, ("DAC Sample: " + std::to_string(_min_chi2_index)).c_str(), "l");
            _canvas_fit_legend->AddEntry(_average_fit_opt, ("#chi^{2}/ndf = " + std::to_string(_fit_chi2).substr(0, 5) + "/" + std::to_string(_fit_ndf)).c_str(), "l");

            if (_is_example_channel) {
                example_fit_func = (TF1*)_average_fit_opt->Clone();
            }
            
        }
        
        _canvas_fit_legend->Draw();
        _canvas_fit->Update();

        channel_adc_smaples_folder->cd();
        _canvas_fit->Write();

        canvas_channel_adc_samples_fitted_list.push_back(_canvas_fit);
        if (_is_example_channel) {
            example_average_waveform = (TGraphErrors*)_graph->Clone();
        }
    }
    // draw the chi2/ndf distribution
    auto _canvas_chi2 = new TCanvas("TemplateDACChi2", "Template DAC Chi2", 800, 600);
    auto _graph_chi2 = new TH1D("TemplateDACChi2", "Template DAC Chi2", 50, 0, 100);
    for (int i = 0; i < template_fit_result_chi2.size(); i++) {
        _graph_chi2->Fill(template_fit_result_chi2[i] / template_fit_result_ndf[i]);
    }
    _graph_chi2->SetTitle("");
    _graph_chi2->SetStats(0);
    _graph_chi2->GetXaxis()->SetTitle("#chi^{2}/ndf");
    _graph_chi2->GetYaxis()->SetTitle("Entries");

    _graph_chi2->SetMarkerStyle(20);
    _graph_chi2->SetMarkerSize(0.2);
    _graph_chi2->SetMarkerColor(kBlack);
    _graph_chi2->SetLineWidth(2);
    _graph_chi2->SetLineColorAlpha(kBlack, 0.5);

    _graph_chi2->GetXaxis()->SetRangeUser(0, 100);
    _graph_chi2->GetYaxis()->SetRangeUser(0, _graph_chi2->GetMaximum() * 1.2);

    _graph_chi2->Draw("HIST");

    auto _chi2_latex = new TLatex();
    _chi2_latex->SetNDC();
    _chi2_latex->SetTextSize(0.04);
    _chi2_latex->SetTextFont(62);
    _chi2_latex->DrawLatex(_text_line_left, _text_line_start, "FoCal-H Prototype II");
    _chi2_latex->SetTextSize(0.03);
    _chi2_latex->SetTextFont(42);
    _chi2_latex->DrawLatex(_text_line_left, _text_line_start - _text_line_height, "Template DAC Chi2/ndf");
    _chi2_latex->DrawLatex(_text_line_left, _text_line_start - 2 * _text_line_height, "SPS H2 Beam Line");
    _chi2_latex->DrawLatex(_text_line_left, _text_line_start - 3 * _text_line_height, "September 2024");
    _chi2_latex->SetTextFont(52);
    _chi2_latex->SetTextColor(kGray + 3);
    _chi2_latex->DrawLatex(_text_line_left, _text_line_start - 54* _text_line_height, "Work in Progress");

    _canvas_chi2->Update();

    _canvas_chi2->Write();
    _canvas_chi2->Print(pdf_file_name.c_str());

    // Save only the central module histograms
    LOG(INFO) << "Canvas Channel ADC Samples List Size: " << canvas_channel_adc_samples_list.size();
    auto painter_canvas = global_painter->draw_module_channel_canvas(canvas_channel_adc_samples_list, channel_unified_channel_number_to_index_map, "ChannelADCSamples", "Channel ADC Samples", module_to_check);
    if (painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas->Write();
    painter_canvas->Print(pdf_file_name.c_str());

    LOG(DEBUG) << "Canvas Channel ADC Samples Fitted List Size: " << canvas_channel_adc_samples_fitted_list.size();
    auto painter_canvas_fitted = global_painter->draw_module_channel_canvas(canvas_channel_adc_samples_fitted_list, channel_unified_channel_number_to_index_map, "ChannelADCSamplesFitted", "Channel ADC Samples Fitted", module_to_check);
    if (painter_canvas_fitted == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }

    output_root->cd();
    painter_canvas_fitted->Write();
    painter_canvas_fitted->Print(pdf_file_name.c_str());

    // draw the template minimal index and maximum bin correlation
    auto _canvas_template = new TCanvas("TemplateMinimalIndex", "Template Minimal Index", 800, 600);
    auto _graph_template = new TGraphErrors(template_minimal_index_list.size(), &template_max_bin_list[0], &template_minimal_index_list[0], 0, 0);
    _graph_template->SetTitle("Template Minimal Index and Maximum Bin");
    _graph_template->SetMarkerStyle(20);
    _graph_template->SetMarkerSize(0.2);
    _graph_template->SetMarkerColor(kBlack);
    _graph_template->SetLineWidth(0);
    _graph_template->SetLineColorAlpha(kBlack, 0.5);
    _graph_template->GetXaxis()->SetRangeUser(0, 1024);
    _graph_template->GetYaxis()->SetRangeUser(0,  20);
    _graph_template->GetYaxis()->SetTitle("Minimal Index");
    _graph_template->GetXaxis()->SetTitle("Maximum Bin");
    _canvas_template->cd();
    _graph_template->Draw("AP same");
    _canvas_template->Update();
    output_root->cd();
    _canvas_template->Write();
    _canvas_template->Print(pdf_file_name.c_str());

    // draw the template DAC fit functions
    auto _canvas_template_data = new TCanvas("TemplateDACFit", "Template DAC Fit", 800, 600);
    auto _canvas_template_data_legend = new TLegend(0.15, 0.7, 0.89, 0.89);
    _canvas_template_data_legend->SetNColumns(2);        // two columns
    _canvas_template_data_legend->SetFillStyle(0);       // transparent (no fill color)
    _canvas_template_data_legend->SetBorderSize(0);      // no border
    _canvas_template_data_legend->SetTextSize(0.02);     // text size as desired

    std::vector <TGraphErrors*> template_DAC_graphs;
    for (int i = 0; i < template_DAC_phase_results.size(); i++) {
        auto _dac_sample = template_DAC_phase_results[i];
        auto _graph_template_data = new TGraphErrors(template_time_column_values.size(), &template_time_column_values[0], &_dac_sample[0], 0, 0);
        _graph_template_data->SetTitle(("Template DAC Fit " + std::to_string(i)).c_str());
        _graph_template_data->SetMarkerStyle(20);
        _graph_template_data->SetMarkerSize(0.2);
        _graph_template_data->SetMarkerColor(color_list[i % color_list.size()]);
        _graph_template_data->SetLineWidth(1);
        _graph_template_data->SetLineColorAlpha(color_list[i % color_list.size()], 1.0);
        _graph_template_data->GetXaxis()->SetRangeUser(template_time_column_values.front(), template_time_column_values.back());
        _graph_template_data->GetYaxis()->SetRangeUser(0, 1024 + 512);
        _graph_template_data->GetXaxis()->SetTitle("Time [ns]");
        _graph_template_data->GetYaxis()->SetTitle("ADC Value");
        _canvas_template_data->cd();
        if (i == 0) {
            _graph_template_data->Draw("APL");
        } else {
            _graph_template_data->Draw("PL same");
        }
        template_DAC_graphs.push_back(_graph_template_data);
        if (i < 10) {
            _canvas_template_data_legend->AddEntry(_graph_template_data, ("DAC Sample: " + std::to_string(template_DAC_values[i])).c_str(), "l");
        }
    }
    _canvas_template_data_legend->Draw();
    _canvas_template_data->Update();
    output_root->cd();
    _canvas_template_data->Write();
    _canvas_template_data->Print(pdf_file_name.c_str());

    // Draw the fitting of the example channel
    auto _canvas_example_channel = new TCanvas("ExampleChannelFit", "Example Channel Fit", 800, 600);
    auto _canvas_example_channel_legend = new TLegend(0.65, 0.7, 0.89, 0.89);
    _canvas_example_channel_legend->SetFillStyle(0);       // transparent (no fill color)
    _canvas_example_channel_legend->SetBorderSize(0);      // no border
    _canvas_example_channel_legend->SetTextSize(0.02);     // text size as desired

    // draw the 2D histogram
    _canvas_example_channel->cd();
    example_waveform_hist2d->SetTitle("Example Channel ADC Fitting");
    example_waveform_hist2d->SetStats(0);
    example_waveform_hist2d->GetXaxis()->SetTitle("Time [ns]");
    example_waveform_hist2d->GetYaxis()->SetTitle("ADC Value");
    example_waveform_hist2d->SetMarkerColorAlpha(kGray+1, 0.5);
    example_waveform_hist2d->Draw("scat");
    _canvas_example_channel->Update();

    // draw the average waveform
    example_average_waveform->SetMarkerStyle(20);
    example_average_waveform->SetMarkerSize(0.2);
    example_average_waveform->SetMarkerColor(kCyan - 2);
    example_average_waveform->SetLineWidth(1);
    example_average_waveform->SetLineColor(kCyan - 2);
    example_average_waveform->Draw("PE same");
    _canvas_example_channel_legend->AddEntry(example_average_waveform, "Average Waveform", "pe");

    // draw the fit function
    if (example_fit_func != nullptr) {
        example_fit_func->SetLineColor(kPink + 7);
        example_fit_func->SetLineWidth(2);
        example_fit_func->SetLineStyle(1);
        example_fit_func->Draw("same");
        _canvas_example_channel_legend->AddEntry(example_fit_func, ("Fit: A = " + std::to_string(example_fit_func->GetParameter(0)).substr(0, 5) + ", t_{offset} = " + std::to_string(example_fit_func->GetParameter(1)).substr(0, 5)).c_str(), "l");
        auto dummy_hist = new TH1D("dummy_hist", "dummy_hist", 100, 0, 100);
        dummy_hist->SetLineColor(kWhite);
        dummy_hist->SetLineWidth(0);

        _canvas_example_channel_legend->AddEntry(dummy_hist, ("Best injection DAC: " + std::to_string(template_DAC_values[template_minimal_index_list[0]])).c_str(), "l");
        _canvas_example_channel_legend->AddEntry(dummy_hist, ("#chi^{2}/ndf = " + std::to_string(example_fit_func->GetChisquare()).substr(0, 5) + "/" + std::to_string(example_fit_func->GetNDF())).c_str(), "l");
        _canvas_example_channel_legend->AddEntry(dummy_hist, ("Slice: " + std::to_string(int(adc_sum_slice_min)) + " - " + std::to_string(int(adc_sum_slice_max))).c_str(), "l");
    }

    _canvas_example_channel_legend->Draw();

    auto example_latex = new TLatex();
    example_latex->SetNDC();
    example_latex->SetTextSize(0.04);
    example_latex->SetTextFont(62);
    example_latex->DrawLatex(_text_line_left, _text_line_start, "FoCal-H Prototype II");
    example_latex->SetTextSize(0.03);
    example_latex->SetTextFont(42);
    example_latex->DrawLatex(_text_line_left, _text_line_start - _text_line_height, ("Channel " + std::to_string(example_fitting_channel) + " ADC Samples").c_str());
    example_latex->DrawLatex(_text_line_left, _text_line_start - 2 * _text_line_height, "SPS H2 Beam Line");
    example_latex->DrawLatex(_text_line_left, _text_line_start - 3 * _text_line_height, "September 2024");
    example_latex->SetTextFont(52);
    example_latex->SetTextColor(kGray + 3);
    example_latex->DrawLatex(_text_line_left, _text_line_start - 4 * _text_line_height, "Work in Progress");


    _canvas_example_channel->Update();
    output_root->cd();
    _canvas_example_channel->Write();
    _canvas_example_channel->Print(pdf_file_name.c_str());

    LOG(DEBUG) << "Template DAC Fit Funcs Size: " << template_DAC_fit_funcs.size();
    // Close the pdf file by saving a dummy canvas
    auto dummy_canvas = new TCanvas("dummy_canvas", "dummy_canvas", 800, 600);
    dummy_canvas->Print((pdf_file_name + ")").c_str());
    dummy_canvas->Close();

    output_root->Close();

    LOG(INFO) << "Output file: " << opts.output_file;
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

    json fit_result;
    fit_result["unified_channel_numbers"] = template_fit_unified_channel_numbers;
    fit_result["A"] = template_fit_result_A;
    fit_result["x0"] = template_fit_result_x0;
    fit_result["DAC_index"] = template_fit_DAC_index;
    fit_result["chi2"] = template_fit_result_chi2;
    fit_result["ndf"] = template_fit_result_ndf;

    std::string fit_result_file = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".json";
    std::ofstream fit_result_stream(fit_result_file);
    fit_result_stream << fit_result.dump(4);
    fit_result_stream.close();

    return 0;
}