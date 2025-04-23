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
    for (size_t i = 0; i < para_time.size() - 1; ++i) {
        if (_xx >= para_time[i] && _xx < para_time[i+1]) {
            double _slope = (para_values[i+1] - para_values[i]) / (para_time[i+1] - para_time[i]);
            double _intercept = para_values[i] + _slope * (_xx - para_time[i]);
            double _res = p[0] * _intercept;
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
    ScriptOptions opts = parse_arguments_single_json(argc, argv, "1.0");

    gROOT->SetBatch(kTRUE);
    ROOT::EnableImplicitMT(30);
    
    bool enable_focal_mapping = opts.focal;

    // * --- Read the configuration file ------------------------------------------------
    // * --------------------------------------------------------------------------------
    auto json_file_name = opts.input_file;
    json config_content;
    std::ifstream json_file(json_file_name);
    if (!json_file.is_open()) {
        LOG(ERROR) << "Failed to open json file " << json_file_name;
        return 1;
    }
    json_file >> config_content;
    json_file.close();

    std::vector <int> config_run_numbers            = config_content["Run Numbers"].get<std::vector<int>>();
    std::vector <std::string> config_run_files      = config_content["Run Files"].get<std::vector<std::string>>();
    std::vector <double> config_beam_energies       = config_content["Beam Energies"].get<std::vector<double>>();
    std::vector <std::string> config_beam_particles = config_content["Beam Particles"].get<std::vector<std::string>>();
    std::vector <int> config_target_machineguns     = config_content["MachineGuns"].get<std::vector<int>>();
  
    std::vector <std::string> config_plot_info      = config_content["Plot Info"].get<std::vector<std::string>>();
    std::string timewalk_json = config_content["TimeWalk JSON"].get<std::string>();
    std::string fitting_json = config_content["Fitting JSON"].get<std::string>();
    bool enable_working_in_progress = config_content["Working in progress"].get<bool>();

    int plot_sum_x_min = config_content["Plot Setting"]["Sum X Min"].get<int>();
    int plot_sum_x_max = config_content["Plot Setting"]["Sum X Max"].get<int>();
    int plot_sum_x_bin = config_content["Plot Setting"]["Sum X Bin"].get<int>();

    // * --- Read the timewalk correction file ------------------------------------------
    // * --------------------------------------------------------------------------------
    json timewalk_content;
    std::vector <int> timewalk_correction_uni_channel;
    std::vector <double> timewalk_correction_parA;
    std::vector <double> timewalk_correction_parB;
    std::vector <double> timewalk_correction_parC;
    std::vector <double> timewalk_correction_parx0;

    std::ifstream timewalk_file(timewalk_json);
    if (!timewalk_file.is_open()) {
        LOG(ERROR) << "Failed to open timewalk json file " << timewalk_json;
        return 1;
    }
    timewalk_file >> timewalk_content;
    timewalk_file.close();

    timewalk_correction_uni_channel = timewalk_content["unified_channel_numbers"].get<std::vector<int>>();
    timewalk_correction_parA        = timewalk_content["A"].get<std::vector<double>>();
    timewalk_correction_parB        = timewalk_content["B"].get<std::vector<double>>();
    timewalk_correction_parC        = timewalk_content["C"].get<std::vector<double>>();
    timewalk_correction_parx0       = timewalk_content["x0"].get<std::vector<double>>();

    // * --- Read the fitting parameters file -------------------------------------------
    // * --------------------------------------------------------------------------------
    json fitting_content;
    std::vector <int> fitting_unified_channel_numbers;
    std::vector <double> fitting_results_A;
    std::vector <double> fitting_results_x0;
    std::vector <double> fitting_results_chi2;
    std::vector <double> fitting_results_ndf;
    std::vector <int> fitting_results_DAC_index;
    
    std::ifstream fitting_file(fitting_json);
    if (!fitting_file.is_open()) {
        LOG(ERROR) << "Failed to open fitting json file " << fitting_json;
        return 1;
    }
    fitting_file >> fitting_content;
    fitting_file.close();

    fitting_unified_channel_numbers = fitting_content["unified_channel_numbers"].get<std::vector<int>>();
    fitting_results_A               = fitting_content["A"].get<std::vector<double>>();
    fitting_results_x0              = fitting_content["x0"].get<std::vector<double>>();
    fitting_results_chi2            = fitting_content["chi2"].get<std::vector<double>>();
    fitting_results_ndf             = fitting_content["ndf"].get<std::vector<double>>();
    fitting_results_DAC_index       = fitting_content["DAC_index"].get<std::vector<int>>();

    // * --- Read the parameterized function file ----------------------------------------
    // * --------------------------------------------------------------------------------
    std::string template_time_column_header = "Time[ns]";
    std::vector<double> template_time_column_values;
    std::vector<int> template_DAC_values;
    std::vector<std::vector<double>> template_DAC_phase_results;
    std::vector<std::vector<double>> template_DAC_phase_results_errors;

    const double time_peak_value    = 75.0;
    const double time_rising_ratio  = 1.8;
    const double time_falling_ratio = 1.8;

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
    
    // * --- Find the run with the lowest and highest beam energy -----------------------
    // * --------------------------------------------------------------------------------
    int run_index_lowest_energy = -1;
    int run_index_highest_energy = -1;
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _beam_energy = config_beam_energies[_run_index];
        if (run_index_lowest_energy == -1 || _beam_energy < config_beam_energies[run_index_lowest_energy]) {
            run_index_lowest_energy = _run_index;
        }
        if (run_index_highest_energy == -1 || _beam_energy > config_beam_energies[run_index_highest_energy]) {
            run_index_highest_energy = _run_index;
        }
    }

    const int channel_adc_hist_bins = 256;
    const double channel_adc_hist_min = 0;
    const double channel_adc_hist_max = 1024;

    const double toa_time_hist_min = 0.0;
    const double toa_time_hist_max = 200.0; // unit: ns
    const int toa_time_hist_bins = 200;

    const double toa_code_hist_min = 0.0;
    const double toa_code_hist_max = 1024.0;
    const int toa_code_hist_bins = 256;

    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 25.0 / 16.0; // unit: ns

    const double beam_energy_relative_error = 0.05;

    const double sample_time_errror = 0.025; // unit: ns
    const double sample_val0_error = 1.0 / sqrt(12.0); // unit: ADC

    const double adc_sum_hist_min = double(plot_sum_x_min);
    const double adc_sum_hist_max = double(plot_sum_x_max);
    const int adc_sum_hist_bins = plot_sum_x_bin;

    const int print_event_number = 10;

    GlobalChannelPainter *global_painter = nullptr;
    if (enable_focal_mapping){
        global_painter = new GlobalChannelPainter("data/SPS_2024/config/focalh_mapping.json", "data/SPS_2024/config/h2gcroc_mapping.json");
    } else {
        global_painter = new GlobalChannelPainter("data/DESY_2025/config/EEEMCal_Mapping_DESY_2025.json");
    }

    // * --- Read pedestal calibration profile ------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *pede_file = new TFile(opts.pedestal_file.c_str(), "READ");
    if (pede_file->IsZombie()) {
        LOG(ERROR) << "Failed to open pedestal file " << opts.pedestal_file;
        return 1;
    }

    TTree *pede_tree = (TTree*) pede_file->Get("ChannelPedestal");
    if (pede_tree == nullptr) {
        LOG(ERROR) << "Failed to get pedestal tree from pedestal file " << opts.pedestal_file;
        return 1;
    }

    std::vector <UInt_t> channel_pede_channel_number;
    std::vector <UInt_t> channel_pede_peak_counts;
    std::vector <Double_t> channel_pede_mean1;
    std::vector <Double_t> channel_pede_error1;
    std::vector <Double_t> channel_pede_mean2;
    std::vector <Double_t> channel_pede_error2;

    UInt_t _channel_pede_channel_number;
    UInt_t _channel_pede_peak_counts;
    Double_t _channel_pede_mean1;
    Double_t _channel_pede_error1;
    Double_t _channel_pede_mean2;
    Double_t _channel_pede_error2;

    pede_tree->SetBranchAddress("ChannelNumber", &_channel_pede_channel_number);
    pede_tree->SetBranchAddress("PeakCounts", &_channel_pede_peak_counts);
    pede_tree->SetBranchAddress("Mean1", &_channel_pede_mean1);
    pede_tree->SetBranchAddress("Error1", &_channel_pede_error1);
    pede_tree->SetBranchAddress("Mean2", &_channel_pede_mean2);
    pede_tree->SetBranchAddress("Error2", &_channel_pede_error2);

    LOG(INFO) << "Pedestal calibration profile loaded with " << pede_tree->GetEntries() << " entries";

    for (int _pede_entry = 0; _pede_entry < pede_tree->GetEntries(); _pede_entry++) {
        pede_tree->GetEntry(_pede_entry);
        channel_pede_channel_number.push_back(_channel_pede_channel_number);
        channel_pede_peak_counts.push_back(_channel_pede_peak_counts);
        channel_pede_mean1.push_back(_channel_pede_mean1);
        channel_pede_error1.push_back(_channel_pede_error1);
        channel_pede_mean2.push_back(_channel_pede_mean2);
        channel_pede_error2.push_back(_channel_pede_error2);
    }

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    std::vector <std::vector <TH2D*>> run_chn_hist2d; // [run_index][channel_index]
    TDirectory *run_chn_sample_hist2d_folder = output_root->mkdir("RunChnSampleHist2D");

    std::unordered_map <int, int> channel_unified_chn_num_to_index;
    std::unordered_map <int, int> channel_index_to_unified_chn_num;

    std::vector <std::vector <TH2D*>> run_chn_toatime_toacode_hist2d; // [run_index][channel_index]
    TDirectory *run_chn_toatime_toacode_hist2d_folder = output_root->mkdir("RunChnToaTimeToaCodeHist2D");

    std::vector <std::vector <TH2D*>> run_chn_toatime_adcmax_hist2d; // [run_index][channel_index]
    TDirectory *run_chn_toatime_adcmax_hist2d_folder = output_root->mkdir("RunChnToaTimeAdcMaxHist2D");

    std::vector <TH1D*> run_adc_sum_hist1d;
    TDirectory *run_adc_sum_hist1d_folder = output_root->mkdir("RunAdcSumHist1D");

    std::vector <TH1D*> run_adc_sum_hist1d_fitted;
    TDirectory *run_adc_sum_hist1d_fitted_folder = output_root->mkdir("RunAdcSumHist1DFitted");

    // * --- Go though all the run files ------------------------------------------------
    // * --------------------------------------------------------------------------------
    if (config_run_numbers.size() != config_run_files.size()) {
        LOG(ERROR) << "Run numbers and run files do not match!";
        return 1;
    }
    if (config_run_numbers.size() != config_beam_energies.size()) {
        LOG(ERROR) << "Run numbers and beam energies do not match!";
        return 1;
    }
    if (config_run_numbers.size() != config_beam_particles.size()) {
        LOG(ERROR) << "Run numbers and beam particles do not match!";
        return 1;
    }

    int fpga_count = -1;
    int machine_gun_samples = -1;
    int entry_max = -1;

    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _run_number    = config_run_numbers[_run_index];
        auto _run_file      = config_run_files[_run_index];
        auto _beam_energy   = config_beam_energies[_run_index];
        auto _beam_particle = config_beam_particles[_run_index];

        LOG(INFO) << "Run " << _run_number << " with beam energy " << _beam_energy << " GeV " << _beam_particle;

        TFile *_input_root = new TFile(_run_file.c_str(), "READ");
        if (_input_root->IsZombie()) {
            LOG(ERROR) << "Failed to open input file " << _run_file;
            return 1;
        }
        TTree *_input_tree = (TTree*) _input_root->Get("data_tree");
        if (_input_tree == nullptr) {
            LOG(ERROR) << "Failed to get data tree from input file " << _run_file;
            return 1;
        }

        int _entry_max = _input_tree->GetEntries();
        if (_entry_max == 0) {
            LOG(ERROR) << "No events in the input file!";
            return 1;
        }
        if (opts.n_events > 0 && opts.n_events < _entry_max) {
            _entry_max = opts.n_events;
        } else {
            if (opts.n_events > _entry_max) {
                LOG(WARNING) << "Requested number of events " << opts.n_events << " is larger than the number of events in the input file " << _entry_max;
            }
        }
        if (entry_max < 0 || entry_max > _entry_max) {
            entry_max = _entry_max;
        }

        TNamed *_legal_fpga_id_list_tnamed = (TNamed*) _input_root->Get("Rootifier_legal_fpga_id_list");
        if (_legal_fpga_id_list_tnamed == nullptr) {
            LOG(ERROR) << "Failed to get legal fpga id list from input file " << _run_file;
            return 1;
        }
        std::string _legal_fpga_id_list_str = _legal_fpga_id_list_tnamed->GetTitle();
        std::vector <UShort_t> _legal_fpga_id_list;
        std::istringstream _legal_fpga_id_list_stream(_legal_fpga_id_list_str);
        UShort_t _legal_fpga_id;
        while (_legal_fpga_id_list_stream >> _legal_fpga_id) {
            _legal_fpga_id_list.push_back(_legal_fpga_id);
        }
        if (fpga_count == -1) {
            fpga_count = _legal_fpga_id_list.size();
        } else {
            if (fpga_count != _legal_fpga_id_list.size()) {
                LOG(ERROR) << "Number of FPGAs in the input files do not match!";
                return 1;
            }
        }

        TNamed *_input_machine_gun_samples_tnamed = (TNamed*) _input_root->Get("EventRecon_machine_gun_samples");
        if (_input_machine_gun_samples_tnamed == nullptr) {
            LOG(ERROR) << "Failed to get machine gun samples from input file " << _run_file;
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

        if (fpga_count > 0 && machine_gun_samples > 0){
            double _max_time = machine_gun_samples * sample_time;
            int _time_bins = int (_max_time / phase_shift_time);
            // -- create 2D histograms array
            // ----------------------------------------------------------------
            std::vector <TH2D*> _chn_hist2d;
            std::vector <TH2D*> _chn_toatime_toacode_hist2d;
            std::vector <TH2D*> _chn_toatime_adcmax_hist2d;

            auto _adc_sum_hist1d = new TH1D(("Run_"+std::to_string(_run_number)+"_AdcSum").c_str(), ("Run " + std::to_string(_run_number) + " AdcSum").c_str(), adc_sum_hist_bins, adc_sum_hist_min, adc_sum_hist_max);
            _adc_sum_hist1d->SetDirectory(run_adc_sum_hist1d_folder);
            run_adc_sum_hist1d.push_back(_adc_sum_hist1d);
            

            auto _adc_sum_fitted_hist1d = new TH1D(("Run_"+std::to_string(_run_number)+"_AdcSumFitted").c_str(), ("Run " + std::to_string(_run_number) + " AdcSumFitted").c_str(), adc_sum_hist_bins, adc_sum_hist_min, adc_sum_hist_max);
            _adc_sum_fitted_hist1d->SetDirectory(run_adc_sum_hist1d_fitted_folder);
            run_adc_sum_hist1d_fitted.push_back(_adc_sum_fitted_hist1d);
            // ----------------------------------------------------------------
            if (run_chn_hist2d.size() == 0) {
                int _hist_index = 0;
                for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                    auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++) {
                        auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);
                        channel_unified_chn_num_to_index[_unified_valid_channel_number] = _hist_index;
                        channel_index_to_unified_chn_num[_hist_index] = _unified_valid_channel_number;

                        auto _hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)).c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number)).c_str(), _time_bins, 0, _max_time, channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max);
                        _chn_hist2d.push_back(_hist2d);
                        _hist2d->SetDirectory(run_chn_sample_hist2d_folder);

                        auto _toatime_toacode_hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)+"_ToaTimeToaCode").c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number) + " ToaTimeToaCode").c_str(), toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max, toa_code_hist_bins, toa_code_hist_min, toa_code_hist_max);
                        _chn_toatime_toacode_hist2d.push_back(_toatime_toacode_hist2d);
                        _toatime_toacode_hist2d->SetDirectory(run_chn_toatime_toacode_hist2d_folder);

                        auto _toatime_adcmax_hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)+"_ToaTimeAdcMax").c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number) + " ToaTimeAdcMax").c_str(), channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max, toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max);
                        _chn_toatime_adcmax_hist2d.push_back(_toatime_adcmax_hist2d);
                        _toatime_adcmax_hist2d->SetDirectory(run_chn_toatime_adcmax_hist2d_folder);
                        
                        _hist_index++;
                    }
                }
                run_chn_hist2d.push_back(_chn_hist2d);
                run_chn_toatime_toacode_hist2d.push_back(_chn_toatime_toacode_hist2d);
                run_chn_toatime_adcmax_hist2d.push_back(_chn_toatime_adcmax_hist2d);
            } else {
                for (int _hist_index = 0; _hist_index < run_chn_hist2d[0].size(); _hist_index++) {
                    auto _unified_valid_channel_number = channel_unified_chn_num_to_index[_hist_index];

                    auto _hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)).c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number)).c_str(), _time_bins, 0, _max_time, channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max);
                    _chn_hist2d.push_back(_hist2d);
                    _hist2d->SetDirectory(run_chn_sample_hist2d_folder);

                    auto _toatime_toacode_hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)+"_ToaTimeToaCode").c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number) + " ToaTimeToaCode").c_str(), toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max, toa_code_hist_bins, toa_code_hist_min, toa_code_hist_max);
                    _chn_toatime_toacode_hist2d.push_back(_toatime_toacode_hist2d);
                    _toatime_toacode_hist2d->SetDirectory(run_chn_toatime_toacode_hist2d_folder);

                    auto _toatime_adcmax_hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)+"_ToaTimeAdcMax").c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number) + " ToaTimeAdcMax").c_str(), channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max, toa_time_hist_bins, toa_time_hist_min, toa_time_hist_max);
                    _chn_toatime_adcmax_hist2d.push_back(_toatime_adcmax_hist2d);
                    _toatime_adcmax_hist2d->SetDirectory(run_chn_toatime_adcmax_hist2d_folder);
                }
                run_chn_hist2d.push_back(_chn_hist2d);
                run_chn_toatime_toacode_hist2d.push_back(_chn_toatime_toacode_hist2d);
                run_chn_toatime_adcmax_hist2d.push_back(_chn_toatime_adcmax_hist2d);
            }
        }
        _input_root->Close();
    }

    LOG(INFO) << "FPGA count: " << fpga_count;
    LOG(INFO) << "Machine gun samples: " << machine_gun_samples;
    LOG(INFO) << "Entry max: " << entry_max;

    // * --- Go though all the events ---------------------------------------------------
    // * --------------------------------------------------------------------------------
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _run_number    = config_run_numbers[_run_index];
        auto _run_file      = config_run_files[_run_index];
        auto _beam_energy   = config_beam_energies[_run_index];
        auto _beam_particle = config_beam_particles[_run_index];

        LOG(INFO) << "Run " << _run_number << " with beam energy " << _beam_energy << " GeV " << _beam_particle;

        TFile *_input_root = new TFile(_run_file.c_str(), "READ");
        if (_input_root->IsZombie()) {
            LOG(ERROR) << "Failed to open input file " << _run_file;
            return 1;
        }
        TTree *_input_tree = (TTree*) _input_root->Get("data_tree");
        if (_input_tree == nullptr) {
            LOG(ERROR) << "Failed to get data tree from input file " << _run_file;
            return 1;
        }

        TNamed *_legal_fpga_id_list_tnamed = (TNamed*) _input_root->Get("Rootifier_legal_fpga_id_list");
        if (_legal_fpga_id_list_tnamed == nullptr) {
            LOG(ERROR) << "Failed to get legal fpga id list from input file " << _run_file;
            return 1;
        }
        std::string _legal_fpga_id_list_str = _legal_fpga_id_list_tnamed->GetTitle();
        std::vector <UShort_t> _legal_fpga_id_list;
        std::istringstream _legal_fpga_id_list_stream(_legal_fpga_id_list_str);
        UShort_t _legal_fpga_id;
        while (_legal_fpga_id_list_stream >> _legal_fpga_id) {
            _legal_fpga_id_list.push_back(_legal_fpga_id);
        }
        if (fpga_count == -1) {
            fpga_count = _legal_fpga_id_list.size();
        } else {
            if (fpga_count != _legal_fpga_id_list.size()) {
                LOG(ERROR) << "Number of FPGAs in the input files do not match!";
                return 1;
            }
        }

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
            auto _fpga_id = _legal_fpga_id_list[_fpga_index];
            auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
            auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
            auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
            auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];
    
            _input_tree->SetBranchAddress(("timestamps_"    + std::to_string(_fpga_id)).c_str(), branch_timestamps);
            _input_tree->SetBranchAddress(("daqh_list_"     + std::to_string(_fpga_id)).c_str(), branch_daqh_list);
            _input_tree->SetBranchAddress(("tc_list_"       + std::to_string(_fpga_id)).c_str(), branch_tc_list);
            _input_tree->SetBranchAddress(("tp_list_"       + std::to_string(_fpga_id)).c_str(), branch_tp_list);
            _input_tree->SetBranchAddress(("val0_list_"     + std::to_string(_fpga_id)).c_str(), branch_val0_list);
            _input_tree->SetBranchAddress(("val1_list_"     + std::to_string(_fpga_id)).c_str(), branch_val1_list);
            _input_tree->SetBranchAddress(("val2_list_"     + std::to_string(_fpga_id)).c_str(), branch_val2_list);
            _input_tree->SetBranchAddress(("crc32_list_"    + std::to_string(_fpga_id)).c_str(), branch_crc32_list);
            _input_tree->SetBranchAddress(("last_heartbeat_"+ std::to_string(_fpga_id)).c_str(), branch_last_heartbeat);
    
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

        auto _run_adc_sum_hist1d = run_adc_sum_hist1d[_run_index];
        auto _run_adc_sum_hist1d_fitted = run_adc_sum_hist1d_fitted[_run_index];

        for (int _entry = 0; _entry < entry_max; _entry++) {
            _input_tree->GetEntry(_entry);

            double _run_adc_sum = 0;
            double _run_adc_sum_fitted = 0;

            for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                auto _fpga_id    = _legal_fpga_id_list[_fpga_index];
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
                    auto _hist_index = channel_unified_chn_num_to_index[_unified_valid_channel_number];
                    int _val0_max = -1;
                    int _val0_max_index = -1;
                    int _val1_max = -1;
                    int _val1_max_index = -1;
                    int _val2_max = -1;
                    int _val2_max_index = -1;

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val1 = _val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val2 = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        if (int(_val0) > _val0_max && int(_val0) < 1024) {
                            _val0_max = int(_val0);
                            if (_val0_max > 10000) {
                                LOG(ERROR) << "Run " << _run_number << " FPGA " << _fpga_id << " Channel " << _channel_valid << " Sample " << _sample_index << " Val0 " << _val0;
                            }
                            _val0_max_index = _sample_index;
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
                    }
                    _channel_val0_max_list.push_back(_val0_max);
                    _channel_val1_max_list.push_back(_val1_max);
                    _channel_val2_max_list.push_back(_val2_max);
                    _channel_val0_max_index_list.push_back(_val0_max_index);
                    _channel_val1_max_index_list.push_back(_val1_max_index);
                    _channel_val2_max_index_list.push_back(_val2_max_index);
                }

                for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                    auto _channel_valid = get_valid_fpga_channel(_channel_index);
                    if (_channel_valid == -1){
                        continue;
                    }
                    auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                    int _pedestal_peak_counts    = int(channel_pede_peak_counts[_unified_valid_channel_number]);
                    if (_pedestal_peak_counts <= 0) {
                        continue;
                    }

                    auto _asic_index = int(_unified_valid_channel_number / (FPGA_CHANNEL_NUMBER_VALID / 2));
                    auto _asic_half_index = int(_unified_valid_channel_number / (FPGA_CHANNEL_NUMBER_VALID / 4)) % 2;

                    auto _val0_max = _channel_val0_max_list[_channel_index];
                    auto _val1_max = _channel_val1_max_list[_channel_index];
                    auto _val2_max = _channel_val2_max_list[_channel_index];
                    auto _val0_max_index = _channel_val0_max_index_list[_channel_index];
                    auto _val1_max_index = _channel_val1_max_index_list[_channel_index];
                    auto _val2_max_index = _channel_val2_max_index_list[_channel_index];

                    if (_val2_max > 0) {
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
                        double _val0_max_double = double(_val0_max);
                        double _val0_max_double_copy = _val0_max_double;
                        auto _timewalk_correction_index = std::find(timewalk_correction_uni_channel.begin(), timewalk_correction_uni_channel.end(), _unified_valid_channel_number);
                        if (_timewalk_correction_index != timewalk_correction_uni_channel.end()) {
                            auto _timewalk_correction_index_int = std::distance(timewalk_correction_uni_channel.begin(), _timewalk_correction_index);
                            double _timewalk_correction_parA = timewalk_correction_parA[_timewalk_correction_index_int];
                            double _timewalk_correction_parB = timewalk_correction_parB[_timewalk_correction_index_int];
                            double _timewalk_correction_parC = timewalk_correction_parC[_timewalk_correction_index_int];
                            double _timewalk_correction_parx0 = timewalk_correction_parx0[_timewalk_correction_index_int];
                            double _param[4] = {_timewalk_correction_parA, _timewalk_correction_parx0, _timewalk_correction_parB, _timewalk_correction_parC};
                            double _timewalk_correction = ExpCorrection(&_val0_max_double_copy, _param);
                            _toa_value_ns += _timewalk_correction;
                        } // end of if (_timewalk_correction_index != timewalk_correction_uni_channel.end())
                    
                        auto _hist_index = channel_unified_chn_num_to_index[_unified_valid_channel_number];
                        auto _channel_sample_hist2d = run_chn_hist2d[_run_index][_hist_index];
                        auto _channel_toatime_toacode_hist2d = run_chn_toatime_toacode_hist2d[_run_index][_hist_index];
                        auto _channel_toatime_adcmax_hist2d = run_chn_toatime_adcmax_hist2d[_run_index][_hist_index];

                        int _pedestal_event = -1;

                        double _channel_adc_value = 0;

                        auto _event_waveform_graph = new TGraph();
                        _event_waveform_graph->SetNameTitle(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)+"_Event_"+std::to_string(_entry)).c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number) + " Event " + std::to_string(_entry)).c_str());

                        std::vector <double> _event_waveform_time;
                        std::vector <double> _event_waveform_val0;
                        // std::vector <double> _event_waveform_time_err;
                        // std::vector <double> _event_waveform_val0_err;
                        for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                            auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                            if (_sample_index == 0) {
                                _pedestal_event = int(_val0);
                            }
                            double _sample_time = _sample_index * sample_time;
                            _sample_time -= (_toa_value_ns - 88.0);
                            _event_waveform_time.push_back(_sample_time);
                            _event_waveform_val0.push_back(_val0);

                            // if (_hamming_code_pass_list[_sample_index] && _val2_max > 0) {
                            //     _channel_sample_hist2d->Fill(_sample_time, _val0);
                            // }

                            // _event_waveform_graph->SetPoint(_sample_index, _sample_time, _val0);
                            // _event_waveform_graph->SetPointError(_sample_index, sample_time_errror, sample_val0_error);
                        }

                        bool _event_valid = true;
                        if (_event_waveform_val0[3] < _event_waveform_val0[2]) {
                            _event_valid = false;
                        }
                        if (_event_waveform_val0[2] < _pedestal_event) {
                            _event_valid = false;
                        }
                        if (_event_valid) {
                            for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                                auto _val0 = _event_waveform_val0[_sample_index];
                                auto _sample_time = _event_waveform_time[_sample_index];
                                _channel_sample_hist2d->Fill(_sample_time, _val0);
                                _event_waveform_graph->SetPoint(_sample_index, _sample_time, _val0);
                            }

                            // ! -- Do the event fitting
                            // ! -- Fit the event waveform
                            std::vector<double> wave_fit_chi2_list;
                            std::vector<int> used_dac_index_list;
                            wave_fit_chi2_list.reserve(template_DAC_phase_results.size());
                            auto _channel_index_in_fit = std::find(fitting_unified_channel_numbers.begin(), fitting_unified_channel_numbers.end(), _unified_valid_channel_number);
                            int _dac_index_fitting_start = 0;
                            int _dac_index_fitting_end = template_DAC_phase_results.size();
                            if (_val0_max_double < 1000.0) {
                                _dac_index_fitting_start = 0;
                                _dac_index_fitting_end = 18;
                            } else {
                                _dac_index_fitting_start = 18;
                                _dac_index_fitting_end = template_DAC_phase_results.size();
                            }
                            for (int _dac_index = _dac_index_fitting_start; _dac_index < _dac_index_fitting_end; _dac_index++) {
                                auto _dac_sample = template_DAC_phase_results[_dac_index];
                                used_dac_index_list.push_back(_dac_index);
                                auto wrappedFitFunc = [=](double *x, double *par) {
                                    return fitParameterizedFunction(x, par, template_time_column_values, _dac_sample);
                                };

                                TF1* _wave_fit = new TF1(("fit_" + std::to_string(_unified_valid_channel_number)).c_str(),
                                                wrappedFitFunc,
                                                0, machine_gun_samples * sample_time, 2);
                                // _wave_fit->SetParameter(0, 1.0);
                                // _wave_fit->SetParLimits(0, 0.8, 1.2);
                                _wave_fit->FixParameter(0, 1.0);
                                // if the fittting results contains the channel number
                                
                                if (_channel_index_in_fit != fitting_unified_channel_numbers.end()) {
                                    auto _channel_index_in_fit_int = std::distance(fitting_unified_channel_numbers.begin(), _channel_index_in_fit);
                                    auto _fitting_x0 = fitting_results_x0[_channel_index_in_fit_int];
                                    _wave_fit->FixParameter(1, _fitting_x0);
                                } else {
                                    _wave_fit->SetParameter(1, (_val0_max_index * sample_time - (time_peak_value * time_rising_ratio)));
                                    _wave_fit->SetParLimits(1, -90.0, -10.0);
                                }

                                _event_waveform_graph->Fit(_wave_fit, "MULTITHREAD RNQ");

                                wave_fit_chi2_list.push_back(_wave_fit->GetChisquare());

                                _wave_fit->Delete();
                            }

                            auto _optimal_chi2_index = std::min_element(wave_fit_chi2_list.begin(), wave_fit_chi2_list.end()) - wave_fit_chi2_list.begin();
                            _optimal_chi2_index = used_dac_index_list[_optimal_chi2_index];
                            auto _optimal_dac_sample = template_DAC_phase_results[_optimal_chi2_index];
                            auto wrappedFitFunc = [=](double *x, double *par) {
                                return fitParameterizedFunction(x, par, template_time_column_values, _optimal_dac_sample);
                            };

                            TF1* _wave_fit_opt = new TF1(("fit_" + std::to_string(_unified_valid_channel_number)).c_str(),
                                            wrappedFitFunc,
                                            0, machine_gun_samples * sample_time, 2);
                            _wave_fit_opt->SetParameter(0, 1.0);
                            _wave_fit_opt->SetParLimits(0, 0.9, 1.1);
                            // if the fittting results contains the channel number
                            if (_channel_index_in_fit != fitting_unified_channel_numbers.end()) {
                                auto _channel_index_in_fit_int = std::distance(fitting_unified_channel_numbers.begin(), _channel_index_in_fit);
                                auto _fitting_x0 = fitting_results_x0[_channel_index_in_fit_int];
                                _wave_fit_opt->FixParameter(1, _fitting_x0);
                            } else {
                                _wave_fit_opt->SetParameter(1, (_val0_max_index * sample_time - (time_peak_value * time_rising_ratio)));
                                _wave_fit_opt->SetParLimits(1, -90.0, -10.0);
                            }

                            _event_waveform_graph->Fit(_wave_fit_opt, "MULTITHREAD RNQ");
                            
                            auto _wave_fit_res_A    = _wave_fit_opt->GetParameter(0);
                            auto _wave_fit_res_x0   = _wave_fit_opt->GetParameter(1);
                            auto _wave_fit_res_chi2 = _wave_fit_opt->GetChisquare();
                            auto _wave_fit_res_ndf  = _wave_fit_opt->GetNDF();

                            _run_adc_sum_fitted += _wave_fit_res_A * template_DAC_values[_optimal_chi2_index] * 3.0;

                            if (_entry < print_event_number) {
                                auto _event_waveform_canvas = new TCanvas(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)+"_Event_"+std::to_string(_entry)).c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number) + " Event " + std::to_string(_entry)).c_str(), 800, 600);
                                _event_waveform_canvas->cd();
                                auto _event_waveform_legend = new TLegend(0.6, 0.6, 0.89, 0.89);
                                _event_waveform_legend->AddEntry(_event_waveform_graph, "Event waveform", "l");
                                _event_waveform_legend->AddEntry(_wave_fit_opt, ("A = " + std::to_string(_wave_fit_res_A) + ", x0 = " + std::to_string(_wave_fit_res_x0) + ", #chi^{2}/NDF = " + std::to_string(_wave_fit_res_chi2) + "/" + std::to_string(_wave_fit_res_ndf)).c_str(), "l");
                                _event_waveform_legend->SetBorderSize(0);
                                _event_waveform_legend->SetFillStyle(0);

                                _event_waveform_graph->SetMarkerStyle(20);
                                _event_waveform_graph->SetMarkerSize(0.5);
                                _event_waveform_graph->SetMarkerColor(kBlack);
                                _event_waveform_graph->SetLineColor(kBlack);
                                _event_waveform_graph->SetLineWidth(1);
                                _event_waveform_graph->GetXaxis()->SetTitle("Time (ns)");
                                _event_waveform_graph->GetYaxis()->SetTitle("ADC value");
                                _event_waveform_graph->Draw("APL");
                                _wave_fit_opt->Draw("SAME");
                                _event_waveform_legend->Draw();
                                output_root->cd();
                                _event_waveform_canvas->Write();
                                _event_waveform_canvas->Close();
                            } // end of if (_entry < print_event_number)
                        } // end of if (_event_valid)
                        _channel_toatime_toacode_hist2d->Fill(_toa_value_ns, _val2_max);
                        _channel_toatime_adcmax_hist2d->Fill(_val0_max_double, _toa_value_ns);

                        // * --- Calculate pedestal subtracted ADC value -------------------------
                        // * --------------------------------------------------------------------
                        double _pedestal_mean1       = channel_pede_mean1[_unified_valid_channel_number];
                        double _pedestal_error1      = channel_pede_error1[_unified_valid_channel_number];
                        double _pedestal_mean2       = channel_pede_mean2[_unified_valid_channel_number]; 
                        double _pedestal_error2      = channel_pede_error2[_unified_valid_channel_number];
                        
                        double _pedestal_value = 0;
                        if (_pedestal_peak_counts == 1){
                            _pedestal_value = _pedestal_mean1;
                            _channel_adc_value = _val0_max_double - _pedestal_value;
                            if (_channel_adc_value < 0){
                                _channel_adc_value = 0;
                            }
                            _run_adc_sum += _channel_adc_value;
                        } else if (_pedestal_peak_counts == 2){
                            auto _pedestal1_dist = std::abs(_pedestal_event - _pedestal_mean1);
                            auto _pedestal2_dist = std::abs(_pedestal_event - _pedestal_mean2);
                            if (_pedestal1_dist < _pedestal2_dist){
                                _pedestal_value = _pedestal_mean1;
                            } else {
                                _pedestal_value = _pedestal_mean2;
                            }
                            _channel_adc_value = _val0_max_double - _pedestal_value;
                            if (_channel_adc_value < 0){
                                _channel_adc_value = 0;
                            }
                            _run_adc_sum += _channel_adc_value;
                        }
                    } // end of if (_val2_max > 0)
                }
            }
            LOG(INFO) << "Run " << _run_number << " event " << _entry << " with ADC sum " << _run_adc_sum << " and fitted ADC sum " << _run_adc_sum_fitted;
            _run_adc_sum_hist1d->Fill(_run_adc_sum);
            _run_adc_sum_hist1d_fitted->Fill(_run_adc_sum_fitted);
        }
        LOG(INFO) << "Run " << _run_number << " with " << _input_valid_tot_counter << " valid TOTs (" << (100.0 * _input_valid_tot_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%) and " << _input_double_tot_counter << " double TOTs (" << (100.0 * _input_double_tot_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%)";
        LOG(INFO) << "Run " << _run_number << " with " << _input_valid_toa_counter << " valid TOAs (" << (100.0 * _input_valid_toa_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%) and " << _input_double_toa_counter << " double TOAs (" << (100.0 * _input_double_toa_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%)";
        _input_root->Close();
    }

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
    std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".pdf";

    // * --- Save the histograms --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    bool is_first_pdf_page_written = false;
    LOG(INFO) << "Saving channel sample histograms ...";
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _channel_hist2d = run_chn_hist2d[_run_index];
        run_chn_sample_hist2d_folder->cd();
        for (int _channel_index = 0; _channel_index < _channel_hist2d.size(); _channel_index++) {
            auto _unified_valid_channel_number = channel_unified_chn_num_to_index[_channel_index];
            auto _hist2d = _channel_hist2d[_channel_index];
            _hist2d->SetStats(0);
            _hist2d->Write();
        }
        if (_run_index == run_index_highest_energy || _run_index == run_index_lowest_energy){
            global_painter->draw_global_channel_hists2D(_channel_hist2d, channel_unified_chn_num_to_index, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnSampleHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
            auto _canvas = global_painter->get_canvas();
            if (!is_first_pdf_page_written){
                _canvas->SaveAs((pdf_file_name + "(").c_str());
                is_first_pdf_page_written = true;
            } else {
                _canvas->SaveAs(pdf_file_name.c_str());
            }
            output_root->cd();
            _canvas->Write();
        }
    }

    LOG(INFO) << "Saving toa time vs toa code histograms ...";
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _channel_toatime_toacode_hist2d = run_chn_toatime_toacode_hist2d[_run_index];
        run_chn_toatime_toacode_hist2d_folder->cd();
        for (int _channel_index = 0; _channel_index < _channel_toatime_toacode_hist2d.size(); _channel_index++) {
            auto _unified_valid_channel_number = channel_unified_chn_num_to_index[_channel_index];
            auto _hist2d = _channel_toatime_toacode_hist2d[_channel_index];
            _hist2d->GetXaxis()->SetTitle("TOA Time [ns]");
            _hist2d->GetYaxis()->SetTitle("TOA Code");
            _hist2d->SetStats(0);
            _hist2d->Write();
        }
        if (_run_index == run_index_lowest_energy || _run_index == run_index_highest_energy){
            global_painter->draw_global_channel_hists2D(_channel_toatime_toacode_hist2d, channel_unified_chn_num_to_index, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnToaTimeToaCodeHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
            auto _canvas = global_painter->get_canvas();
            _canvas->SaveAs(pdf_file_name.c_str());
            output_root->cd();
            _canvas->Write();
        }
    }

    LOG(INFO) << "Saving toa time vs adc max histograms ...";
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _channel_toatime_adcmax_hist2d = run_chn_toatime_adcmax_hist2d[_run_index];
        run_chn_toatime_adcmax_hist2d_folder->cd();
        for (int _channel_index = 0; _channel_index < _channel_toatime_adcmax_hist2d.size(); _channel_index++) {
            auto _unified_valid_channel_number = channel_unified_chn_num_to_index[_channel_index];
            auto _hist2d = _channel_toatime_adcmax_hist2d[_channel_index];
            _hist2d->GetXaxis()->SetTitle("ADC Max");
            _hist2d->GetYaxis()->SetTitle("TOA Time [ns]");
            _hist2d->SetStats(0);
            _hist2d->Write();
        }
        if (_run_index == run_index_lowest_energy || _run_index == run_index_highest_energy){
            global_painter->draw_global_channel_hists2D(_channel_toatime_adcmax_hist2d, channel_unified_chn_num_to_index, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnToaTimeAdcMaxHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
            auto _canvas = global_painter->get_canvas();
            _canvas->Print(pdf_file_name.c_str());
            output_root->cd();
            _canvas->Write();
        }
    }

    // * --- Save the ADC sum histograms -------------------------------------------------
    // * --------------------------------------------------------------------------------
    int run_adc_sum_hist_max = 0;
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _adc_sum_hist1d = run_adc_sum_hist1d[_run_index];
        if (_adc_sum_hist1d->GetMaximum() > run_adc_sum_hist_max) {
            run_adc_sum_hist_max = _adc_sum_hist1d->GetMaximum();
        }
        run_adc_sum_hist1d_folder->cd();
        _adc_sum_hist1d->SetStats(0);
        _adc_sum_hist1d->Write();
    }
    auto adc_sum_hist1d_canvas = new TCanvas("adc_sum_hist1d_canvas", "adc_sum_hist1d_canvas", 800, 600);
    auto adc_sum_legend = new TLegend(0.42, 0.5, 0.89, 0.89);
    adc_sum_legend->SetBorderSize(0);
    adc_sum_legend->SetFillStyle(0);
    adc_sum_legend->SetTextSize(0.02);
    adc_sum_legend->SetTextFont(102);

    std::vector <double> adc_sum_fit_mu_list;
    std::vector <double> adc_sum_fit_sigma_list;
    std::vector <double> adc_sum_fit_mu_error_list;
    std::vector <double> adc_sum_fit_sigma_error_list;
    std::vector <double> adc_sum_beam_energy_list;
    std::vector <double> adc_sum_beam_energy_error_list;

    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _adc_sum_hist1d = run_adc_sum_hist1d[_run_index];
        _adc_sum_hist1d->SetMaximum(run_adc_sum_hist_max * 1.3);
        _adc_sum_hist1d->SetLineColor(color_list[_run_index % color_list.size()]);
        if (_run_index == 0) {
            _adc_sum_hist1d->SetTitle("ADC Sum for Different Beam Energies");
            _adc_sum_hist1d->GetXaxis()->SetTitle("ADC Sum");
            _adc_sum_hist1d->GetYaxis()->SetTitle("Counts");
            _adc_sum_hist1d->Draw();
        } else {
            _adc_sum_hist1d->Draw("SAME");
        }
        // * do the gaussian fit
        std::vector <double> _adc_sum_fit_range = {1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; // unit: sigma
        std::vector <double> _adc_sum_fit_offsets = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double _adc_sum_fit_range_min = 0;
        double _adc_sum_fit_range_max = 0;
        double _adc_sum_fit_offset_min = 0;
        double _adc_sum_fit_offset_max = 0;
        for (auto _range: _adc_sum_fit_range) {
            if (_range < _adc_sum_fit_range_min) {
                _adc_sum_fit_range_min = _range;
            }
            if (_range > _adc_sum_fit_range_max) {
                _adc_sum_fit_range_max = _range;
            }
        }
        for (auto _offset: _adc_sum_fit_offsets) {
            if (_offset < _adc_sum_fit_offset_min) {
                _adc_sum_fit_offset_min = _offset;
            }
            if (_offset > _adc_sum_fit_offset_max) {
                _adc_sum_fit_offset_max = _offset;
            }
        }

        // * find initial parameters
        double _adc_sum_fit_mean = _adc_sum_hist1d->GetBinCenter(_adc_sum_hist1d->GetMaximumBin());
        double _adc_sum_fit_rms = _adc_sum_hist1d->GetRMS();
        double _adc_sum_fit_initial_fit_min = _adc_sum_fit_mean - _adc_sum_fit_rms * 2;
        double _adc_sum_fit_initial_fit_max = _adc_sum_fit_mean + _adc_sum_fit_rms * 2;
        if (_adc_sum_fit_initial_fit_min < 0) {
            _adc_sum_fit_initial_fit_min = 0;
        }
        if (_adc_sum_fit_initial_fit_max > adc_sum_hist_max) {
            _adc_sum_fit_initial_fit_max = adc_sum_hist_max;
        }
        // * do the first round of fitting to get the mu and sigma
        TF1 *_adc_sum_fit_function = new TF1("adc_sum_fit_function", "gaus", _adc_sum_fit_initial_fit_min, _adc_sum_fit_initial_fit_max);
        _adc_sum_hist1d->Fit(_adc_sum_fit_function, "RQN");
        double _adc_sum_fit_mu = _adc_sum_fit_function->GetParameter(1);
        double _adc_sum_fit_sigma = _adc_sum_fit_function->GetParameter(2);

        // * do the fitting for pedestal peak
        TF1 *_adc_sum_pede_fit_function = new TF1("adc_sum_pede_fit_function", "gaus", adc_sum_hist_min, adc_sum_hist_min + (adc_sum_hist_max - adc_sum_hist_min) * 0.06);
        // set parameter range
        _adc_sum_pede_fit_function->SetParLimits(1, adc_sum_hist_min, adc_sum_hist_min + (adc_sum_hist_max - adc_sum_hist_min) * 0.1);
        _adc_sum_hist1d->Fit(_adc_sum_pede_fit_function, "RQN");
        double _adc_sum_pede_mu = _adc_sum_pede_fit_function->GetParameter(1);
        double _adc_sum_pede_sigma = _adc_sum_pede_fit_function->GetParameter(2);
        // draw the pede mu
        TLine *_adc_sum_pede_mu_line = new TLine(_adc_sum_pede_mu, 0, _adc_sum_pede_mu, run_adc_sum_hist_max * 1.3);
        _adc_sum_pede_mu_line->SetLineColor(color_list[_run_index % color_list.size()]);
        _adc_sum_pede_mu_line->SetLineStyle(2);
        _adc_sum_pede_mu_line->Draw("SAME");

        std::vector <double> _adc_sum_fit_res_mean;
        std::vector <double> _adc_sum_fit_res_sigma;
        std::vector <double> _adc_sum_fit_res_mean_error;
        std::vector <double> _adc_sum_fit_res_sigma_error;

        for (auto _fit_range: _adc_sum_fit_range) {
            for (auto _fit_offset: _adc_sum_fit_offsets) {
                double _adc_sum_fit_min = _adc_sum_fit_mu - _fit_range * _adc_sum_fit_sigma + _fit_offset * _adc_sum_fit_sigma;
                double _adc_sum_fit_max = _adc_sum_fit_mu + _fit_range * _adc_sum_fit_sigma + _fit_offset * _adc_sum_fit_sigma;
                if (_adc_sum_fit_min < 0) {
                    _adc_sum_fit_min = 0;
                }
                if (_adc_sum_fit_max > adc_sum_hist_max) {
                    _adc_sum_fit_max = adc_sum_hist_max;
                }
                auto _adc_sum_fit_function = new TF1(("adc_sum_fit_function_"+std::to_string(_run_index)+"_"+std::to_string(_fit_range)+"_"+std::to_string(_fit_offset)).c_str(), "gaus", _adc_sum_fit_min, _adc_sum_fit_max);
                _adc_sum_fit_function->SetParameters(_adc_sum_hist1d->GetMaximum(), _adc_sum_fit_mu, _adc_sum_fit_sigma);
                _adc_sum_hist1d->Fit(_adc_sum_fit_function, "RQN+");
                _adc_sum_fit_function->SetLineColorAlpha(color_list[_run_index % color_list.size()], 0.2);
                _adc_sum_fit_function->Draw("SAME");

                _adc_sum_fit_res_mean.push_back(_adc_sum_fit_function->GetParameter(1));
                _adc_sum_fit_res_sigma.push_back(_adc_sum_fit_function->GetParameter(2));
                _adc_sum_fit_res_mean_error.push_back(_adc_sum_fit_function->GetParError(1));
                _adc_sum_fit_res_sigma_error.push_back(_adc_sum_fit_function->GetParError(2));
            }
        }

        // * calculate the mean and sigma
        double _adc_sum_fit_res_mean_weighted  = 0;
        double _adc_sum_fit_res_sigma_weighted = 0;
        double _adc_sum_fit_res_mean_err_sys   = 0;
        double _adc_sum_fit_res_sigma_err_sys  = 0;
        double _adc_sum_fit_res_mean_err_stat  = 0;
        double _adc_sum_fit_res_sigma_err_stat = 0;

        // * calculate the weighted mean and sigma
        double _adc_sum_fit_res_mean_weighted_sum = 0;
        double _adc_sum_fit_res_sigma_weighted_sum = 0;
        double _adc_sum_fit_res_mean_weighted_sum_weight = 0;
        double _adc_sum_fit_res_sigma_weighted_sum_weight = 0;

        for (int _fit_index = 0; _fit_index < _adc_sum_fit_res_mean.size(); _fit_index++) {
            double _mean = _adc_sum_fit_res_mean[_fit_index];
            double _sigma = _adc_sum_fit_res_sigma[_fit_index];
            double _mean_error = _adc_sum_fit_res_mean_error[_fit_index];
            double _sigma_error = _adc_sum_fit_res_sigma_error[_fit_index];
            double _weight_mean = 1.0 / (_mean_error * _mean_error);
            double _weight_sigma = 1.0 / (_sigma_error * _sigma_error);
            _adc_sum_fit_res_mean_weighted_sum += _mean * _weight_mean;
            _adc_sum_fit_res_sigma_weighted_sum += _sigma * _weight_sigma;
            _adc_sum_fit_res_mean_weighted_sum_weight += _weight_mean;
            _adc_sum_fit_res_sigma_weighted_sum_weight += _weight_sigma;
        }

        _adc_sum_fit_res_mean_weighted = _adc_sum_fit_res_mean_weighted_sum / _adc_sum_fit_res_mean_weighted_sum_weight;
        _adc_sum_fit_res_sigma_weighted = _adc_sum_fit_res_sigma_weighted_sum / _adc_sum_fit_res_sigma_weighted_sum_weight;

        // * calculate the systematic error
        double _adc_sum_fit_res_mean_err_sys_sum = 0;
        double _adc_sum_fit_res_sigma_err_sys_sum = 0;
        for (int _fit_index = 0; _fit_index < _adc_sum_fit_res_mean.size(); _fit_index++) {
            double _mean = _adc_sum_fit_res_mean[_fit_index];
            double _sigma = _adc_sum_fit_res_sigma[_fit_index];
            double _mean_error = _adc_sum_fit_res_mean_error[_fit_index];
            double _sigma_error = _adc_sum_fit_res_sigma_error[_fit_index];
            _adc_sum_fit_res_mean_err_sys_sum += (_mean - _adc_sum_fit_res_mean_weighted) * (_mean - _adc_sum_fit_res_mean_weighted);
            _adc_sum_fit_res_sigma_err_sys_sum += (_sigma - _adc_sum_fit_res_sigma_weighted) * (_sigma - _adc_sum_fit_res_sigma_weighted);
        }
        _adc_sum_fit_res_mean_err_sys = std::sqrt(_adc_sum_fit_res_mean_err_sys_sum / (_adc_sum_fit_res_mean.size() - 1));
        _adc_sum_fit_res_sigma_err_sys = std::sqrt(_adc_sum_fit_res_sigma_err_sys_sum / (_adc_sum_fit_res_sigma.size() - 1));

        // subtract the pedestal
        // _adc_sum_fit_res_mean_weighted -= _adc_sum_pede_mu;

        // * calculate the statistical error
        double _adc_sum_fit_res_mean_err_stat_sum = 0;
        double _adc_sum_fit_res_sigma_err_stat_sum = 0;
        for (int _fit_index = 0; _fit_index < _adc_sum_fit_res_mean.size(); _fit_index++) {
            double _mean = _adc_sum_fit_res_mean[_fit_index];
            double _sigma = _adc_sum_fit_res_sigma[_fit_index];
            double _mean_error = _adc_sum_fit_res_mean_error[_fit_index];
            double _sigma_error = _adc_sum_fit_res_sigma_error[_fit_index];
            _adc_sum_fit_res_mean_err_stat_sum += 1.0 / (_mean_error * _mean_error);
            _adc_sum_fit_res_sigma_err_stat_sum += 1.0 / (_sigma_error * _sigma_error);
        }
        _adc_sum_fit_res_mean_err_stat = 1.0 / std::sqrt(_adc_sum_fit_res_mean_err_stat_sum);
        _adc_sum_fit_res_sigma_err_stat = 1.0 / std::sqrt(_adc_sum_fit_res_sigma_err_stat_sum);

        // fixed length string
        std::string _adc_sum_mean_str = std::to_string(_adc_sum_fit_res_mean_weighted).substr(0, 5);
        std::string _adc_sum_sigma_str = std::to_string(_adc_sum_fit_res_sigma_weighted).substr(0, 4);
        std::string _adc_sum_mean_err_stat_str = std::to_string(_adc_sum_fit_res_mean_err_stat).substr(0, 4);
        std::string _adc_sum_sigma_err_stat_str = std::to_string(_adc_sum_fit_res_sigma_err_stat).substr(0, 4);
        std::string _adc_sum_mean_err_sys_str = std::to_string(_adc_sum_fit_res_mean_err_sys).substr(0, 5);
        std::string _adc_sum_sigma_err_sys_str = std::to_string(_adc_sum_fit_res_sigma_err_sys).substr(0, 5);
        std::string _adc_sum_beam_energy_str;
        if (int(config_beam_energies[_run_index]) < 100) {
            _adc_sum_beam_energy_str = " " + std::to_string(int(config_beam_energies[_run_index]));
        } else {
            _adc_sum_beam_energy_str = std::to_string(int(config_beam_energies[_run_index]));
        }
        std::string _adc_sum_beam_energy_dummy_str = "     ";
        for (int _beam_energy_index = 0; _beam_energy_index < 3 - _adc_sum_beam_energy_str.size(); _beam_energy_index++) {
            _adc_sum_beam_energy_dummy_str += " ";
        }

        adc_sum_fit_mu_list.push_back(_adc_sum_fit_res_mean_weighted);
        adc_sum_fit_sigma_list.push_back(_adc_sum_fit_res_sigma_weighted);
        adc_sum_fit_mu_error_list.push_back(sqrt(_adc_sum_fit_res_mean_err_stat * _adc_sum_fit_res_mean_err_stat + _adc_sum_fit_res_mean_err_sys * _adc_sum_fit_res_mean_err_sys));
        adc_sum_fit_sigma_error_list.push_back(sqrt(_adc_sum_fit_res_sigma_err_stat * _adc_sum_fit_res_sigma_err_stat + _adc_sum_fit_res_sigma_err_sys * _adc_sum_fit_res_sigma_err_sys));
        adc_sum_beam_energy_list.push_back(config_beam_energies[_run_index]);
        adc_sum_beam_energy_error_list.push_back(config_beam_energies[_run_index] * beam_energy_relative_error);

        // * write the results to legend
        std::string _adc_sum_fit_res_mean_str = (_adc_sum_beam_energy_str + "GeV mu: " + _adc_sum_mean_str + " #pm " + _adc_sum_mean_err_stat_str + "(stat) #pm " + _adc_sum_mean_err_sys_str + "(sys)").c_str();
        std::string _adc_sum_fit_res_sigma_str = (_adc_sum_beam_energy_dummy_str + "sigma: " + _adc_sum_sigma_str + " #pm " + _adc_sum_sigma_err_stat_str + "(stat) #pm " + _adc_sum_sigma_err_sys_str + "(sys)").c_str();
        auto dummy_hist = new TH1D("dummy_hist", "dummy_hist", 1, 0, 1);
        // no line for the dummy hist
        dummy_hist->SetLineColor(kWhite);
        adc_sum_legend->AddEntry(_adc_sum_hist1d, _adc_sum_fit_res_mean_str.c_str(), "l");
        adc_sum_legend->AddEntry(dummy_hist, _adc_sum_fit_res_sigma_str.c_str(), "l");

        dummy_hist->Delete();
        
        // adc_sum_legend->AddEntry(_adc_sum_hist1d, (std::to_string(int(config_beam_energies[_run_index])) + " GeV " + config_beam_particles[_run_index]).c_str(), "l");
    }
    adc_sum_legend->Draw();
    auto adc_sum_latex = new TLatex();
    const double _text_line_height = 0.04;
    const double _text_line_start = 0.85;
    const double _text_line_left = 0.13;
    adc_sum_latex->SetNDC();
    adc_sum_latex->SetTextSize(0.04);
    adc_sum_latex->SetTextFont(62);
    adc_sum_latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    adc_sum_latex->SetTextSize(0.03);
    adc_sum_latex->SetTextFont(42);
    for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
        adc_sum_latex->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    }
    if (enable_working_in_progress){
        adc_sum_latex->SetTextFont(52);
        adc_sum_latex->SetTextColor(kGray+3);
        adc_sum_latex->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    }
    
    output_root->cd();
    adc_sum_hist1d_canvas->Print(pdf_file_name.c_str());
    adc_sum_hist1d_canvas->Write();

    // * --- Draw linearity plot --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    auto linearity_canvas = new TCanvas("linearity_canvas", "linearity_canvas", 800, 600);
    linearity_canvas->SetMargin(0.15, 0.1, 0.15, 0.1);
    auto linearity_dummy_graph = new TGraphErrors(1); // only for setting the axis
    linearity_dummy_graph->SetTitle("Linearity of ADC Sum");
    linearity_dummy_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    linearity_dummy_graph->GetYaxis()->SetTitle("ADC Sum Mean");
    linearity_dummy_graph->GetXaxis()->SetRangeUser(0, 400);
    linearity_dummy_graph->GetXaxis()->SetLimits(0, 400);
    linearity_dummy_graph->GetYaxis()->SetRangeUser(0, adc_sum_hist_max);
    linearity_dummy_graph->GetYaxis()->SetLimits(0, adc_sum_hist_max);  

    linearity_dummy_graph->Draw("AP");

    auto linearity_graph = new TGraphErrors(adc_sum_beam_energy_list.size(), &adc_sum_beam_energy_list[0], &adc_sum_fit_mu_list[0], &adc_sum_beam_energy_error_list[0], &adc_sum_fit_mu_error_list[0]);
    linearity_graph->SetMarkerStyle(20);
    linearity_graph->SetMarkerSize(0.5);
    linearity_graph->SetLineWidth(2);

    auto linearity_fit_function = new TF1("linearity_fit_function", "[0]*x + [1]", 0, 400);
    linearity_fit_function->SetParameter(0, 400.0/adc_sum_hist_max);
    linearity_fit_function->SetParameter(1, 0);
    linearity_graph->Fit(linearity_fit_function, "RQN");
    linearity_graph->Draw("PE SAME");

    linearity_fit_function->SetLineColor(kCyan+3);
    linearity_fit_function->SetLineStyle(2);
    linearity_fit_function->Draw("SAME");

    auto linearity_fit_confidence_band = new TH1F("linearity_fit_confidence_band", "linearity_fit_confidence_band", 100, 0, 400);
    TVirtualFitter::GetFitter()->GetConfidenceIntervals(linearity_fit_confidence_band, 0.68);
    linearity_fit_confidence_band->SetFillColorAlpha(kCyan+3, 0.5);
    linearity_fit_confidence_band->Draw("E3 SAME");

    auto linearity_legend = new TLegend(0.5, 0.7, 0.89, 0.89);
    linearity_legend->SetBorderSize(0);
    linearity_legend->SetFillStyle(0);
    linearity_legend->SetTextSize(0.02);
    linearity_legend->SetTextFont(102);

    double _linearity_fit_slope = linearity_fit_function->GetParameter(0);
    double _linearity_fit_slope_error = linearity_fit_function->GetParError(0);
    double _linearity_fit_intercept = linearity_fit_function->GetParameter(1);
    double _linearity_fit_intercept_error = linearity_fit_function->GetParError(1);
    double _linearity_fit_chi2 = linearity_fit_function->GetChisquare();
    double _linearity_fit_ndf = linearity_fit_function->GetNDF();

    auto linearity_dummy_hist = new TH1D("linearity_dummy_hist", "linearity_dummy_hist", 1, 0, 1);
    linearity_dummy_hist->SetLineColor(kWhite);
    linearity_legend->AddEntry(linearity_graph, "Gaussian Fit Mean", "epl");
    linearity_legend->AddEntry(linearity_fit_function, ("Fit: slope     = " + std::to_string(_linearity_fit_slope).substr(0,6) + " #pm " + std::to_string(_linearity_fit_slope_error).substr(0,3)).c_str(), "l");
    linearity_legend->AddEntry(linearity_dummy_hist, ("     intercept = " + std::to_string(_linearity_fit_intercept).substr(0,6) + " #pm " + std::to_string(_linearity_fit_intercept_error).substr(0,3)).c_str(), "l");
    linearity_legend->AddEntry(linearity_dummy_hist, ("     #chi^{2}/NDF     = " + std::to_string(_linearity_fit_chi2).substr(0,4) + "/" + std::to_string(int(_linearity_fit_ndf))).c_str(), "l");
    
    linearity_legend->Draw();

    auto linearity_latex = new TLatex();
    linearity_latex->SetNDC();
    linearity_latex->SetTextSize(0.04);
    linearity_latex->SetTextFont(62);
    linearity_latex->DrawLatex(_text_line_left+0.05, _text_line_start, (config_plot_info[0].c_str()));
    linearity_latex->SetTextSize(0.03);
    linearity_latex->SetTextFont(42);
    for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
        linearity_latex->DrawLatex(_text_line_left+0.05, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    }
    if (enable_working_in_progress){
        linearity_latex->SetTextFont(52);
        linearity_latex->SetTextColor(kGray+3);
        linearity_latex->DrawLatex(_text_line_left+0.05, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    }

    linearity_canvas->Print(pdf_file_name.c_str());
    linearity_canvas->Write();

    // * --- Draw resolution plot --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    auto resolution_canvas = new TCanvas("resolution_canvas", "resolution_canvas", 800, 600);
    auto resolution_dummy_graph = new TGraphErrors(1); // only for setting the axis
    resolution_dummy_graph->SetTitle("Resolution of ADC Sum");
    resolution_dummy_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    resolution_dummy_graph->GetYaxis()->SetTitle("#sigma_{E}/E");
    resolution_dummy_graph->GetXaxis()->SetRangeUser(0, 400);
    resolution_dummy_graph->GetYaxis()->SetRangeUser(0.1, 0.4);
    resolution_dummy_graph->GetYaxis()->SetLimits(0.1, 0.4);
    resolution_dummy_graph->GetXaxis()->SetLimits(0, 400);
    resolution_dummy_graph->Draw("AP");

    std::vector <double> adc_sum_fit_sigmaE_E_list;
    std::vector <double> adc_sum_fit_sigmaE_E_error_list;
    std::vector <double> reconstructed_beam_energy_list;
    std::vector <double> reconstructed_beam_energy_list_error;
    for (int _beam_energy_index = 0; _beam_energy_index < adc_sum_beam_energy_list.size(); _beam_energy_index++) {
        double _beam_energy = adc_sum_beam_energy_list[_beam_energy_index];
        double _adc_sum_sigma = adc_sum_fit_sigma_list[_beam_energy_index];
        double _adc_sum_sigma_error = adc_sum_fit_sigma_error_list[_beam_energy_index];
        double _adc_sum_mu = adc_sum_fit_mu_list[_beam_energy_index];
        double _adc_sum_mu_error = adc_sum_fit_mu_error_list[_beam_energy_index];
        double _adc_sum_sigmaE_E = _adc_sum_sigma / _adc_sum_mu;
        double _adc_sum_sigmaE_E_error = _adc_sum_sigmaE_E * std::sqrt((_adc_sum_sigma_error / _adc_sum_sigma) * (_adc_sum_sigma_error / _adc_sum_sigma) + (_adc_sum_mu_error / _adc_sum_mu) * (_adc_sum_mu_error / _adc_sum_mu));
        adc_sum_fit_sigmaE_E_list.push_back(_adc_sum_sigmaE_E);
        adc_sum_fit_sigmaE_E_error_list.push_back(_adc_sum_sigmaE_E_error);
        double _reconstructed_beam_energy = (_adc_sum_mu - _linearity_fit_intercept)/_linearity_fit_slope;
        double _reconstructed_beam_energy_error = _reconstructed_beam_energy * beam_energy_relative_error;
        reconstructed_beam_energy_list.push_back(_reconstructed_beam_energy);
        reconstructed_beam_energy_list_error.push_back(_reconstructed_beam_energy_error);
    }
    auto resolution_graph = new TGraphErrors(reconstructed_beam_energy_list.size(), &reconstructed_beam_energy_list[0], &adc_sum_fit_sigmaE_E_list[0], &reconstructed_beam_energy_list_error[0], &adc_sum_fit_sigmaE_E_error_list[0]);
    resolution_graph->SetTitle("Resolution of ADC Sum");
    resolution_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    resolution_graph->GetXaxis()->SetRangeUser(0, 400);
    resolution_graph->GetYaxis()->SetTitle("#sigma_{E}/E");
    resolution_graph->GetYaxis()->SetRangeUser(0.1, 0.4);
    resolution_graph->SetMarkerStyle(20);
    resolution_graph->SetMarkerSize(0.5);
    resolution_graph->SetLineWidth(2);
    
    auto resolution_fit_function = new TF1("resolution_fit_function", "[0]/sqrt(x) + [1]", 0, 400);
    resolution_fit_function->SetParameter(0, 0.1);
    resolution_fit_function->SetParameter(1, 0.01);
    resolution_graph->Fit(resolution_fit_function, "RQN");
    resolution_graph->Draw("PE SAME");

    resolution_fit_function->SetLineColor(kRed);
    resolution_fit_function->SetLineStyle(2);
    resolution_fit_function->Draw("SAME");

    auto resolution_fit_confidence_band = new TH1F("resolution_fit_confidence_band", "resolution_fit_confidence_band", 100, 0, 400);
    TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    if (fitter) {
        fitter->GetConfidenceIntervals(resolution_fit_confidence_band, 0.68);
    } else {
        LOG(WARNING) << "Fitter is not initialized. Confidence intervals will not be drawn.";
    }
    resolution_fit_confidence_band->SetFillColorAlpha(kRed, 0.3);
    resolution_fit_confidence_band->SetLineColor(kRed);
    resolution_fit_confidence_band->SetLineStyle(2);
    resolution_fit_confidence_band->Draw("e3 SAME");

    auto resolution_legend = new TLegend(0.5, 0.7, 0.89, 0.89);
    resolution_legend->SetBorderSize(0);
    resolution_legend->SetFillStyle(0);
    resolution_legend->SetTextSize(0.02);
    resolution_legend->SetTextFont(102);

    double _resolution_fit_a = resolution_fit_function->GetParameter(0);
    double _resolution_fit_a_error = resolution_fit_function->GetParError(0);
    double _resolution_fit_b = resolution_fit_function->GetParameter(1);
    double _resolution_fit_b_error = resolution_fit_function->GetParError(1);
    double _resolution_fit_chi2 = resolution_fit_function->GetChisquare();
    double _resolution_fit_ndf = resolution_fit_function->GetNDF();

    auto resolution_dummy_hist = new TH1D("resolution_dummy_hist", "resolution_dummy_hist", 1, 0, 1);
    resolution_dummy_hist->SetLineColor(kWhite);
    resolution_legend->AddEntry(resolution_graph, "Gaussian Fit Resolution", "epl");
    resolution_legend->AddEntry(resolution_fit_function, ("#sigma_{E}/E   = #frac{" + std::to_string(_resolution_fit_a).substr(0,4) + " #pm " + std::to_string(_resolution_fit_a_error).substr(0,3) + "}{#sqrt{E}} #oplus (" + std::to_string(_resolution_fit_b).substr(0,4) + " #pm " + std::to_string(_resolution_fit_b_error).substr(0,4) + ")").c_str(), "l");
    resolution_legend->AddEntry(resolution_dummy_hist, ("#chi^{2}/NDF = " + std::to_string(_resolution_fit_chi2).substr(0,4) + "/" + std::to_string(int(_resolution_fit_ndf))).c_str(), "l");
    resolution_legend->Draw();

    auto resolution_latex = new TLatex();
    resolution_latex->SetNDC();
    resolution_latex->SetTextSize(0.04);
    resolution_latex->SetTextFont(62);
    resolution_latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    resolution_latex->SetTextSize(0.03);
    resolution_latex->SetTextFont(42);
    for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
        resolution_latex->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    }

    if (enable_working_in_progress){
        resolution_latex->SetTextFont(52);
        resolution_latex->SetTextColor(kGray+3);
        resolution_latex->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    }

    resolution_canvas->Print(pdf_file_name.c_str());
    resolution_canvas->Write();

    // * --- Draw the Fitted ADC Sum 
    // * --------------------------------------------------------------------------------
    int run_adc_sum_hist_max_fitted = 0;
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _adc_sum_hist1d_fitted = run_adc_sum_hist1d_fitted[_run_index];
        if (_adc_sum_hist1d_fitted->GetMaximum() > run_adc_sum_hist_max_fitted) {
            run_adc_sum_hist_max_fitted = _adc_sum_hist1d_fitted->GetMaximum();
        }
        run_adc_sum_hist1d_fitted_folder->cd();
        _adc_sum_hist1d_fitted->SetStats(0);
        _adc_sum_hist1d_fitted->Write();
    }
    auto adc_sum_hist1d_fitted_canvas = new TCanvas("adc_sum_hist1d_fitted_canvas", "adc_sum_hist1d_fitted_canvas", 800, 600);
    auto adc_sum_legend_fitted = new TLegend(0.42, 0.5, 0.89, 0.89);
    adc_sum_legend_fitted->SetBorderSize(0);
    adc_sum_legend_fitted->SetFillStyle(0);
    adc_sum_legend_fitted->SetTextSize(0.02);
    adc_sum_legend_fitted->SetTextFont(102);

    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _adc_sum_hist1d_fitted = run_adc_sum_hist1d_fitted[_run_index];
        _adc_sum_hist1d_fitted->SetMaximum(run_adc_sum_hist_max_fitted * 1.3);
        _adc_sum_hist1d_fitted->SetLineColor(color_list[_run_index % color_list.size()]);
        if (_run_index == 0) {
            _adc_sum_hist1d_fitted->SetTitle("Fitted ADC Sum for Different Beam Energies");
            _adc_sum_hist1d_fitted->GetXaxis()->SetTitle("ADC Sum");
            _adc_sum_hist1d_fitted->GetYaxis()->SetTitle("Counts");
            _adc_sum_hist1d_fitted->Draw();
        } else {
            _adc_sum_hist1d_fitted->Draw("SAME");
        }
        // * write the results to legend
        // adc_sum_legend_fitted->AddEntry(_adc_sum_hist1d_fitted, (std::to_string(int(config_beam_energies[_run_index])) + " GeV " + config_beam_particles[_run_index]).c_str(), "l");
        // * do the gaussian fit
        std::vector <double> _adc_sum_fit_range = {1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; // unit: sigma
        std::vector <double> _adc_sum_fit_offsets = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double _adc_sum_fit_range_min = 0;
        double _adc_sum_fit_range_max = 0;
        double _adc_sum_fit_offset_min = 0;
        double _adc_sum_fit_offset_max = 0;
        for (auto _range: _adc_sum_fit_range) {
            if (_range < _adc_sum_fit_range_min) {
                _adc_sum_fit_range_min = _range;
            }
            if (_range > _adc_sum_fit_range_max) {
                _adc_sum_fit_range_max = _range;
            }
        }
        for (auto _offset: _adc_sum_fit_offsets) {
            if (_offset < _adc_sum_fit_offset_min) {
                _adc_sum_fit_offset_min = _offset;
            }
            if (_offset > _adc_sum_fit_offset_max) {
                _adc_sum_fit_offset_max = _offset;
            }
        }

        // * find initial parameters
        double _adc_sum_fit_mean = _adc_sum_hist1d_fitted->GetBinCenter(_adc_sum_hist1d_fitted->GetMaximumBin());
        double _adc_sum_fit_rms = _adc_sum_hist1d_fitted->GetRMS();
        double _adc_sum_fit_initial_fit_min = _adc_sum_fit_mean - _adc_sum_fit_rms * 2;
        double _adc_sum_fit_initial_fit_max = _adc_sum_fit_mean + _adc_sum_fit_rms * 2;
        if (_adc_sum_fit_initial_fit_min < 0) {
            _adc_sum_fit_initial_fit_min = 0;
        }
        if (_adc_sum_fit_initial_fit_max > adc_sum_hist_max) {
            _adc_sum_fit_initial_fit_max = adc_sum_hist_max;
        }
        // * do the first round of fitting to get the mu and sigma
        TF1 *_adc_sum_fit_function = new TF1("adc_sum_fit_function", "gaus", _adc_sum_fit_initial_fit_min, _adc_sum_fit_initial_fit_max);
        _adc_sum_hist1d_fitted->Fit(_adc_sum_fit_function, "RQN");
        double _adc_sum_fit_mu = _adc_sum_fit_function->GetParameter(1);
        double _adc_sum_fit_sigma = _adc_sum_fit_function->GetParameter(2);

        // * do the fitting for pedestal peak
        TF1 *_adc_sum_pede_fit_function = new TF1("adc_sum_pede_fit_function", "gaus", adc_sum_hist_min, adc_sum_hist_min + (adc_sum_hist_max - adc_sum_hist_min) * 0.06);
        // set parameter range
        _adc_sum_pede_fit_function->SetParLimits(1, adc_sum_hist_min, adc_sum_hist_min + (adc_sum_hist_max - adc_sum_hist_min) * 0.1);
        _adc_sum_hist1d_fitted->Fit(_adc_sum_pede_fit_function, "RQN");
        double _adc_sum_pede_mu = _adc_sum_pede_fit_function->GetParameter(1);
        double _adc_sum_pede_sigma = _adc_sum_pede_fit_function->GetParameter(2);
        // draw the pede mu
        TLine *_adc_sum_pede_mu_line = new TLine(_adc_sum_pede_mu, 0, _adc_sum_pede_mu, run_adc_sum_hist_max_fitted * 1.3);
        _adc_sum_pede_mu_line->SetLineColor(color_list[_run_index % color_list.size()]);
        _adc_sum_pede_mu_line->SetLineStyle(2);
        _adc_sum_pede_mu_line->Draw("SAME");

        std::vector <double> fit_adc_sum_fit_res_mean;
        std::vector <double> fit_adc_sum_fit_res_sigma;
        std::vector <double> fit_adc_sum_fit_res_mean_error;
        std::vector <double> fit_adc_sum_fit_res_sigma_error;

        for (auto _fit_range: _adc_sum_fit_range) {
            for (auto _fit_offset: _adc_sum_fit_offsets) {
                double _adc_sum_fit_min = _adc_sum_fit_mu - _fit_range * _adc_sum_fit_sigma + _fit_offset * _adc_sum_fit_sigma;
                double _adc_sum_fit_max = _adc_sum_fit_mu + _fit_range * _adc_sum_fit_sigma + _fit_offset * _adc_sum_fit_sigma;
                if (_adc_sum_fit_min < 0) {
                    _adc_sum_fit_min = 0;
                }
                if (_adc_sum_fit_max > adc_sum_hist_max) {
                    _adc_sum_fit_max = adc_sum_hist_max;
                }
                auto _adc_sum_fit_function = new TF1(("adc_sum_fit_function_"+std::to_string(_run_index)+"_"+std::to_string(_fit_range)+"_"+std::to_string(_fit_offset)).c_str(), "gaus", _adc_sum_fit_min, _adc_sum_fit_max);
                _adc_sum_fit_function->SetParameters(_adc_sum_hist1d_fitted->GetMaximum(), _adc_sum_fit_mu, _adc_sum_fit_sigma);
                _adc_sum_hist1d_fitted->Fit(_adc_sum_fit_function, "RQN+");
                _adc_sum_fit_function->SetLineColorAlpha(color_list[_run_index % color_list.size()], 0.2);
                _adc_sum_fit_function->Draw("SAME");

                fit_adc_sum_fit_res_mean.push_back(_adc_sum_fit_function->GetParameter(1));
                fit_adc_sum_fit_res_sigma.push_back(_adc_sum_fit_function->GetParameter(2));
                fit_adc_sum_fit_res_mean_error.push_back(_adc_sum_fit_function->GetParError(1));
                fit_adc_sum_fit_res_sigma_error.push_back(_adc_sum_fit_function->GetParError(2));
            }
        }

        // * calculate the mean and sigma
        double fit_adc_sum_fit_res_mean_weighted  = 0;
        double fit_adc_sum_fit_res_sigma_weighted = 0;
        double fit_adc_sum_fit_res_mean_err_sys   = 0;
        double fit_adc_sum_fit_res_sigma_err_sys  = 0;
        double fit_adc_sum_fit_res_mean_err_stat  = 0;
        double fit_adc_sum_fit_res_sigma_err_stat = 0;

        // * calculate the weighted mean and sigma
        double fit_adc_sum_fit_res_mean_weighted_sum = 0;
        double fit_adc_sum_fit_res_sigma_weighted_sum = 0;
        double fit_adc_sum_fit_res_mean_weighted_sum_weight = 0;
        double fit_adc_sum_fit_res_sigma_weighted_sum_weight = 0;

        for (int _fit_index = 0; _fit_index < fit_adc_sum_fit_res_mean.size(); _fit_index++) {
            double _mean = fit_adc_sum_fit_res_mean[_fit_index];
            double _sigma = fit_adc_sum_fit_res_sigma[_fit_index];
            double _mean_error = fit_adc_sum_fit_res_mean_error[_fit_index];
            double _sigma_error = fit_adc_sum_fit_res_sigma_error[_fit_index];
            double _weight_mean = 1.0 / (_mean_error * _mean_error);
            double _weight_sigma = 1.0 / (_sigma_error * _sigma_error);
            fit_adc_sum_fit_res_mean_weighted_sum += _mean * _weight_mean;
            fit_adc_sum_fit_res_sigma_weighted_sum += _sigma * _weight_sigma;
            fit_adc_sum_fit_res_mean_weighted_sum_weight += _weight_mean;
            fit_adc_sum_fit_res_sigma_weighted_sum_weight += _weight_sigma;
        }

        fit_adc_sum_fit_res_mean_weighted = fit_adc_sum_fit_res_mean_weighted_sum / fit_adc_sum_fit_res_mean_weighted_sum_weight;
        fit_adc_sum_fit_res_sigma_weighted = fit_adc_sum_fit_res_sigma_weighted_sum / fit_adc_sum_fit_res_sigma_weighted_sum_weight;

        // * calculate the systematic error
        double fit_adc_sum_fit_res_mean_err_sys_sum = 0;
        double fit_adc_sum_fit_res_sigma_err_sys_sum = 0;

        for (int _fit_index = 0; _fit_index < fit_adc_sum_fit_res_mean.size(); _fit_index++) {
            double _mean = fit_adc_sum_fit_res_mean[_fit_index];
            double _sigma = fit_adc_sum_fit_res_sigma[_fit_index];
            double _mean_error = fit_adc_sum_fit_res_mean_error[_fit_index];
            double _sigma_error = fit_adc_sum_fit_res_sigma_error[_fit_index];
            fit_adc_sum_fit_res_mean_err_sys_sum += (_mean - fit_adc_sum_fit_res_mean_weighted) * (_mean - fit_adc_sum_fit_res_mean_weighted);
            fit_adc_sum_fit_res_sigma_err_sys_sum += (_sigma - fit_adc_sum_fit_res_sigma_weighted) * (_sigma - fit_adc_sum_fit_res_sigma_weighted);
        }

        fit_adc_sum_fit_res_mean_err_sys = std::sqrt(fit_adc_sum_fit_res_mean_err_sys_sum / (fit_adc_sum_fit_res_mean.size() - 1));
        fit_adc_sum_fit_res_sigma_err_sys = std::sqrt(fit_adc_sum_fit_res_sigma_err_sys_sum / (fit_adc_sum_fit_res_sigma.size() - 1));

        // subtract the pedestal
        // fit_adc_sum_fit_res_mean_weighted -= _adc_sum_pede_mu;

        // * calculate the statistical error
        double fit_adc_sum_fit_res_mean_err_stat_sum = 0;
        double fit_adc_sum_fit_res_sigma_err_stat_sum = 0;

        for (int _fit_index = 0; _fit_index < fit_adc_sum_fit_res_mean.size(); _fit_index++) {
            double _mean = fit_adc_sum_fit_res_mean[_fit_index];
            double _sigma = fit_adc_sum_fit_res_sigma[_fit_index];
            double _mean_error = fit_adc_sum_fit_res_mean_error[_fit_index];
            double _sigma_error = fit_adc_sum_fit_res_sigma_error[_fit_index];
            fit_adc_sum_fit_res_mean_err_stat_sum += 1.0 / (_mean_error * _mean_error);
            fit_adc_sum_fit_res_sigma_err_stat_sum += 1.0 / (_sigma_error * _sigma_error);
        }

        fit_adc_sum_fit_res_mean_err_stat = 1.0 / std::sqrt(fit_adc_sum_fit_res_mean_err_stat_sum);
        fit_adc_sum_fit_res_sigma_err_stat = 1.0 / std::sqrt(fit_adc_sum_fit_res_sigma_err_stat_sum);

        // fixed length string
        std::string _adc_sum_mean_str = std::to_string(fit_adc_sum_fit_res_mean_weighted).substr(0, 5);
        std::string _adc_sum_sigma_str = std::to_string(fit_adc_sum_fit_res_sigma_weighted).substr(0, 4);
        std::string _adc_sum_mean_err_stat_str = std::to_string(fit_adc_sum_fit_res_mean_err_stat).substr(0, 4);
        std::string _adc_sum_sigma_err_stat_str = std::to_string(fit_adc_sum_fit_res_sigma_err_stat).substr(0, 4);
        std::string _adc_sum_mean_err_sys_str = std::to_string(fit_adc_sum_fit_res_mean_err_sys).substr(0, 5);
        std::string _adc_sum_sigma_err_sys_str = std::to_string(fit_adc_sum_fit_res_sigma_err_sys).substr(0, 5);
        std::string _adc_sum_beam_energy_str;

        if (int(config_beam_energies[_run_index]) < 100) {
            _adc_sum_beam_energy_str = " " + std::to_string(int(config_beam_energies[_run_index]));
        } else {
            _adc_sum_beam_energy_str = std::to_string(int(config_beam_energies[_run_index]));
        }
        std::string _adc_sum_beam_energy_dummy_str = "     ";
        for (int _beam_energy_index = 0; _beam_energy_index < 3 - _adc_sum_beam_energy_str.size(); _beam_energy_index++) {
            _adc_sum_beam_energy_dummy_str += " ";
        }

        adc_sum_fit_mu_list.push_back(fit_adc_sum_fit_res_mean_weighted);
        adc_sum_fit_sigma_list.push_back(fit_adc_sum_fit_res_sigma_weighted);
        adc_sum_fit_mu_error_list.push_back(sqrt(fit_adc_sum_fit_res_mean_err_stat * fit_adc_sum_fit_res_mean_err_stat + fit_adc_sum_fit_res_mean_err_sys * fit_adc_sum_fit_res_mean_err_sys));
        adc_sum_fit_sigma_error_list.push_back(sqrt(fit_adc_sum_fit_res_sigma_err_stat * fit_adc_sum_fit_res_sigma_err_stat + fit_adc_sum_fit_res_sigma_err_sys * fit_adc_sum_fit_res_sigma_err_sys));
        adc_sum_beam_energy_list.push_back(config_beam_energies[_run_index]);
        adc_sum_beam_energy_error_list.push_back(config_beam_energies[_run_index] * beam_energy_relative_error);

        // * write the results to legend
        std::string fit_adc_sum_fit_res_mean_str = (_adc_sum_beam_energy_str + "GeV mu: " + _adc_sum_mean_str + " #pm " + _adc_sum_mean_err_stat_str + "(stat) #pm " + _adc_sum_mean_err_sys_str + "(sys)").c_str();
        std::string fit_adc_sum_fit_res_sigma_str = (_adc_sum_beam_energy_dummy_str + "sigma: " + _adc_sum_sigma_str + " #pm " + _adc_sum_sigma_err_stat_str + "(stat) #pm " + _adc_sum_sigma_err_sys_str + "(sys)").c_str();
        auto dummy_hist = new TH1D("dummy_hist", "dummy_hist", 1, 0, 1);
        // no line for the dummy hist
        dummy_hist->SetLineColor(kWhite);
        adc_sum_legend_fitted->AddEntry(_adc_sum_hist1d_fitted, fit_adc_sum_fit_res_mean_str.c_str(), "l");
        adc_sum_legend_fitted->AddEntry(dummy_hist, fit_adc_sum_fit_res_sigma_str.c_str(), "l");
        
    }

    adc_sum_legend_fitted->Draw();
    auto adc_sum_latex_fitted = new TLatex();
    adc_sum_latex_fitted->SetNDC();
    adc_sum_latex_fitted->SetTextSize(0.04);
    adc_sum_latex_fitted->SetTextFont(62);
    adc_sum_latex_fitted->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    adc_sum_latex_fitted->SetTextSize(0.03);
    adc_sum_latex_fitted->SetTextFont(42);
    for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
        adc_sum_latex_fitted->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    }

    if (enable_working_in_progress){
        adc_sum_latex_fitted->SetTextFont(52);
        adc_sum_latex_fitted->SetTextColor(kGray+3);
        adc_sum_latex_fitted->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    }

    output_root->cd();
    adc_sum_hist1d_fitted_canvas->Print(pdf_file_name.c_str());
    adc_sum_hist1d_fitted_canvas->Write();

    // * --- Draw linearity plot --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    auto fit_linearity_fit_canvas = new TCanvas("fit_linearity_fit_canvas", "fit_linearity_fit_canvas", 800, 600);
    fit_linearity_fit_canvas->SetMargin(0.15, 0.1, 0.15, 0.1);
    auto fit_linearity_fit_dummy_graph = new TGraphErrors(1); // only for setting the axis
    fit_linearity_fit_dummy_graph->SetTitle("Linearity of Fitted ADC Sum");
    fit_linearity_fit_dummy_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    fit_linearity_fit_dummy_graph->GetYaxis()->SetTitle("ADC Sum Mean");
    fit_linearity_fit_dummy_graph->GetXaxis()->SetRangeUser(0, 400);
    fit_linearity_fit_dummy_graph->GetXaxis()->SetLimits(0, 400);
    fit_linearity_fit_dummy_graph->GetYaxis()->SetRangeUser(0, adc_sum_hist_max);
    fit_linearity_fit_dummy_graph->GetYaxis()->SetLimits(0, adc_sum_hist_max);
    fit_linearity_fit_dummy_graph->Draw("AP");

    auto fit_linearity_fit_graph = new TGraphErrors(adc_sum_beam_energy_list.size(), &adc_sum_beam_energy_list[0], &adc_sum_fit_mu_list[0], &adc_sum_beam_energy_error_list[0], &adc_sum_fit_mu_error_list[0]);
    fit_linearity_fit_graph->SetMarkerStyle(20);
    fit_linearity_fit_graph->SetMarkerSize(0.5);
    fit_linearity_fit_graph->SetLineWidth(2);

    auto fit_linearity_fit_function = new TF1("fit_linearity_fit_function", "[0]*x + [1]", 0, 400);
    fit_linearity_fit_function->SetParameter(0, 400.0/adc_sum_hist_max);
    fit_linearity_fit_function->SetParameter(1, 0);
    fit_linearity_fit_graph->Fit(fit_linearity_fit_function, "RQN");
    fit_linearity_fit_graph->Draw("PE SAME");

    fit_linearity_fit_function->SetLineColor(kCyan+3);
    fit_linearity_fit_function->SetLineStyle(2);
    fit_linearity_fit_function->Draw("SAME");

    auto fit_linearity_fit_confidence_band = new TH1F("fit_linearity_fit_confidence_band", "fit_linearity_fit_confidence_band", 100, 0, 400);
    TVirtualFitter::GetFitter()->GetConfidenceIntervals(fit_linearity_fit_confidence_band, 0.68);
    fit_linearity_fit_confidence_band->SetFillColorAlpha(kCyan+3, 0.5);
    fit_linearity_fit_confidence_band->SetMarkerColorAlpha(kCyan+3, 0.0);
    fit_linearity_fit_confidence_band->Draw("E3 SAME");

    auto fit_linearity_fit_legend = new TLegend(0.5, 0.7, 0.89, 0.89);
    fit_linearity_fit_legend->SetBorderSize(0);
    fit_linearity_fit_legend->SetFillStyle(0);
    fit_linearity_fit_legend->SetTextSize(0.02);
    fit_linearity_fit_legend->SetTextFont(102);

    double _fit_linearity_fit_slope_fitted = fit_linearity_fit_function->GetParameter(0);
    double _fit_linearity_fit_slope_error_fitted = fit_linearity_fit_function->GetParError(0);
    double _fit_linearity_fit_intercept_fitted = fit_linearity_fit_function->GetParameter(1);
    double _fit_linearity_fit_intercept_error_fitted = fit_linearity_fit_function->GetParError(1);
    double _fit_linearity_fit_chi2_fitted = fit_linearity_fit_function->GetChisquare();
    double _fit_linearity_fit_ndf_fitted = fit_linearity_fit_function->GetNDF();

    auto fit_linearity_dummy_hist = new TH1D("fit_linearity_dummy_hist", "fit_linearity_dummy_hist", 1, 0, 1);
    fit_linearity_fit_legend->AddEntry(fit_linearity_fit_graph, "Gaussian Fit Mean", "epl");
    fit_linearity_fit_legend->AddEntry(fit_linearity_fit_function, ("Fit: slope     = " + std::to_string(_fit_linearity_fit_slope_fitted).substr(0,6) + " #pm " + std::to_string(_fit_linearity_fit_slope_error_fitted).substr(0,3)).c_str(), "l");

    fit_linearity_fit_legend->AddEntry(fit_linearity_dummy_hist, ("     intercept = " + std::to_string(_fit_linearity_fit_intercept_fitted).substr(0,6) + " #pm " + std::to_string(_fit_linearity_fit_intercept_error_fitted).substr(0,3)).c_str(), "l");
    fit_linearity_fit_legend->AddEntry(fit_linearity_dummy_hist, ("     #chi^{2}/NDF     = " + std::to_string(_fit_linearity_fit_chi2_fitted).substr(0,4) + "/" + std::to_string(int(_fit_linearity_fit_ndf_fitted))).c_str(), "l");

    fit_linearity_fit_legend->Draw();

    auto fit_linearity_fit_latex = new TLatex();
    fit_linearity_fit_latex->SetNDC();
    fit_linearity_fit_latex->SetTextSize(0.04);
    fit_linearity_fit_latex->SetTextFont(62);
    fit_linearity_fit_latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    fit_linearity_fit_latex->SetTextSize(0.03);
    fit_linearity_fit_latex->SetTextFont(42);
    
    for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
        fit_linearity_fit_latex->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    }

    if (enable_working_in_progress){
        fit_linearity_fit_latex->SetTextFont(52);
        fit_linearity_fit_latex->SetTextColor(kGray+3);
        fit_linearity_fit_latex->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    }

    output_root->cd();
    fit_linearity_fit_canvas->Print(pdf_file_name.c_str());
    fit_linearity_fit_canvas->Write();

    // * --- Draw resolution plot --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    auto fit_resolution_fit_canvas = new TCanvas("fit_resolution_fit_canvas", "fit_resolution_fit_canvas", 800, 600);
    auto fit_resolution_fit_dummy_graph = new TGraphErrors(1); // only for setting the axis
    fit_resolution_fit_dummy_graph->SetTitle("Resolution of Fitted ADC Sum");
    fit_resolution_fit_dummy_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    fit_resolution_fit_dummy_graph->GetYaxis()->SetTitle("#sigma_{E}/E");
    fit_resolution_fit_dummy_graph->GetXaxis()->SetRangeUser(0, 400);
    fit_resolution_fit_dummy_graph->GetYaxis()->SetRangeUser(0.1, 0.4);
    fit_resolution_fit_dummy_graph->GetYaxis()->SetLimits(0.1, 0.4);
    fit_resolution_fit_dummy_graph->GetXaxis()->SetLimits(0, 400);
    fit_resolution_fit_dummy_graph->Draw("AP");

    std::vector <double> fit_adc_sum_fit_sigmaE_E_list;
    std::vector <double> fit_adc_sum_fit_sigmaE_E_error_list;
    std::vector <double> fit_reconstructed_beam_energy_list;
    std::vector <double> fit_reconstructed_beam_energy_list_error;

    for (int _beam_energy_index = 0; _beam_energy_index < adc_sum_beam_energy_list.size(); _beam_energy_index++) {
        double _beam_energy = adc_sum_beam_energy_list[_beam_energy_index];
        double _adc_sum_sigma = adc_sum_fit_sigma_list[_beam_energy_index];
        double _adc_sum_sigma_error = adc_sum_fit_sigma_error_list[_beam_energy_index];
        double _adc_sum_mu = adc_sum_fit_mu_list[_beam_energy_index];
        double _adc_sum_mu_error = adc_sum_fit_mu_error_list[_beam_energy_index];
        double _adc_sum_sigmaE_E = _adc_sum_sigma / _adc_sum_mu;
        double _adc_sum_sigmaE_E_error = _adc_sum_sigmaE_E * std::sqrt((_adc_sum_sigma_error / _adc_sum_sigma) * (_adc_sum_sigma_error / _adc_sum_sigma) + (_adc_sum_mu_error / _adc_sum_mu) * (_adc_sum_mu_error / _adc_sum_mu));
        fit_adc_sum_fit_sigmaE_E_list.push_back(_adc_sum_sigmaE_E);
        fit_adc_sum_fit_sigmaE_E_error_list.push_back(_adc_sum_sigmaE_E_error);
        double _fit_reconstructed_beam_energy = (_adc_sum_mu - _fit_linearity_fit_intercept_fitted)/_fit_linearity_fit_slope_fitted;
        double _fit_reconstructed_beam_energy_error = _fit_reconstructed_beam_energy * beam_energy_relative_error;
        fit_reconstructed_beam_energy_list.push_back(_fit_reconstructed_beam_energy);
        fit_reconstructed_beam_energy_list_error.push_back(_fit_reconstructed_beam_energy_error);
    }

    auto fit_resolution_fit_graph = new TGraphErrors(fit_reconstructed_beam_energy_list.size(), &fit_reconstructed_beam_energy_list[0], &fit_adc_sum_fit_sigmaE_E_list[0], &fit_reconstructed_beam_energy_list_error[0], &fit_adc_sum_fit_sigmaE_E_error_list[0]);

    fit_resolution_fit_graph->SetTitle("Resolution of Fitted ADC Sum");
    fit_resolution_fit_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    fit_resolution_fit_graph->GetXaxis()->SetRangeUser(0, 400);
    fit_resolution_fit_graph->GetYaxis()->SetTitle("#sigma_{E}/E");
    fit_resolution_fit_graph->GetYaxis()->SetRangeUser(0.1, 0.4);
    fit_resolution_fit_graph->SetMarkerStyle(20);
    fit_resolution_fit_graph->SetMarkerSize(0.5);
    fit_resolution_fit_graph->SetLineWidth(2);

    auto fit_resolution_fit_function = new TF1("fit_resolution_fit_function", "[0]/sqrt(x) + [1]", 0, 400);

    fit_resolution_fit_function->SetParameter(0, 0.1);
    fit_resolution_fit_function->SetParameter(1, 0.01);
    fit_resolution_fit_graph->Fit(fit_resolution_fit_function, "RQN");
    fit_resolution_fit_graph->Draw("PE SAME");

    fit_resolution_fit_function->SetLineColor(kRed);
    fit_resolution_fit_function->SetLineStyle(2);
    fit_resolution_fit_function->Draw("SAME");

    auto fit_resolution_fit_confidence_band = new TH1F("fit_resolution_fit_confidence_band", "fit_resolution_fit_confidence_band", 100, 0, 400);
    TVirtualFitter *fit_fitter = TVirtualFitter::GetFitter();
    if (fit_fitter) {
        fit_fitter->GetConfidenceIntervals(fit_resolution_fit_confidence_band, 0.68);
    } else {
        LOG(WARNING) << "Fitter is not initialized. Confidence intervals will not be drawn.";
    }

    fit_resolution_fit_confidence_band->SetFillColorAlpha(kRed, 0.3);
    fit_resolution_fit_confidence_band->SetLineColor(kRed);
    fit_resolution_fit_confidence_band->SetLineStyle(2);
    fit_resolution_fit_confidence_band->SetMarkerColorAlpha(kRed, 0.0);
    fit_resolution_fit_confidence_band->Draw("e3 SAME");

    auto fit_resolution_legend = new TLegend(0.5, 0.7, 0.89, 0.89);
    fit_resolution_legend->SetBorderSize(0);
    fit_resolution_legend->SetFillStyle(0);
    fit_resolution_legend->SetTextSize(0.02);
    fit_resolution_legend->SetTextFont(102);

    double _fit_resolution_fit_a = fit_resolution_fit_function->GetParameter(0);
    double _fit_resolution_fit_a_error = fit_resolution_fit_function->GetParError(0);
    double _fit_resolution_fit_b = fit_resolution_fit_function->GetParameter(1);
    double _fit_resolution_fit_b_error = fit_resolution_fit_function->GetParError(1);
    double _fit_resolution_fit_chi2 = fit_resolution_fit_function->GetChisquare();
    double _fit_resolution_fit_ndf = fit_resolution_fit_function->GetNDF();

    auto fit_resolution_dummy_hist = new TH1D("fit_resolution_dummy_hist", "fit_resolution_dummy_hist", 1, 0, 1);
    fit_resolution_dummy_hist->SetLineColor(kWhite);
    fit_resolution_legend->AddEntry(fit_resolution_fit_graph, "Gaussian Fit Resolution", "epl");
    fit_resolution_legend->AddEntry(fit_resolution_fit_function, ("#sigma_{E}/E   = #frac{" + std::to_string(_fit_resolution_fit_a).substr(0,4) + " #pm " + std::to_string(_fit_resolution_fit_a_error).substr(0,3) + "}{#sqrt{E}} #oplus (" + std::to_string(_fit_resolution_fit_b).substr(0,4) + " #pm " + std::to_string(_fit_resolution_fit_b_error).substr(0,4) + ")").c_str(), "l");

    fit_resolution_legend->AddEntry(fit_resolution_dummy_hist, ("#chi^{2}/NDF = " + std::to_string(_fit_resolution_fit_chi2).substr(0,4) + "/" + std::to_string(int(_fit_resolution_fit_ndf))).c_str(), "l");

    fit_resolution_legend->Draw();

    auto fit_resolution_latex = new TLatex();
    fit_resolution_latex->SetNDC();
    fit_resolution_latex->SetTextSize(0.04);
    fit_resolution_latex->SetTextFont(62);
    fit_resolution_latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    fit_resolution_latex->SetTextSize(0.03);
    fit_resolution_latex->SetTextFont(42);
    for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
        fit_resolution_latex->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    }

    if (enable_working_in_progress){
        fit_resolution_latex->SetTextFont(52);
        fit_resolution_latex->SetTextColor(kGray+3);
        fit_resolution_latex->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    }

    output_root->cd();

    fit_resolution_fit_canvas->Print(pdf_file_name.c_str());
    fit_resolution_fit_canvas->Write();




    // * --- Close the files ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TCanvas* dummy_canvas = new TCanvas("dummy_canvas", "dummy_canvas", 800, 600);
    dummy_canvas->SaveAs((pdf_file_name + ")").c_str());
    dummy_canvas->Close();

    output_root->Close();
    LOG(INFO) << "Output file " << opts.output_file << " has been saved.";
    LOG(INFO) << "PDF file " << pdf_file_name << " has been saved.";

    return 0;
}