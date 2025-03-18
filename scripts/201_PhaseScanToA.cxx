#include "H2GCROC_Common.hxx"

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

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_json(argc, argv, "1.0");

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

    std::vector <int> config_run_numbers    = config_content["Run Numbers"].get<std::vector<int>>();
    std::vector <std::string> config_run_files = config_content["Run Files"].get<std::vector<std::string>>();
    std::vector <int> config_phase_settings = config_content["Phase Settings"].get<std::vector<int>>();
    std::string config_beam_particle = config_content["Beam Particle"].get<std::string>();
    double config_beam_energy = config_content["Beam Energy"].get<double>();
    std::string config_calibration_profile = config_content["Calibration Profile"].get<std::string>();
    std::string config_info = config_content["Info"].get<std::string>();
    // std::vector <UInt_t> hist_cut_values = config_content["ADC Cut"].get<std::vector<UInt_t>>();
    std::vector <UInt_t> hist_cut_values = {};
    std::vector <double> toa_thresholds; // in adc
    bool enable_focal_mapping = opts.focal;

    LOG(INFO) << "Phase scan: " << config_info;

    const double phase_time  = 25.0 / 16.0;  // unit: ns
    const double sample_time = 25.0;         // unit: ns
    const int initial_phase = 7;
    const int hist_y_bins = 256;
    const double hist_x_bin_size = phase_time;
    const double hist_y_min = 0;
    const double hist_y_max = 1024;
    const double toa_ns_bin_size = phase_time;
    const double toa_ns_min = 25.0;
    const double toa_ns_max = 7*25.0;

    bool enable_timewalk_correction = false;
    std::string timewalk_correction_file = "dump/201_PhaseScanToA/FoCal/PhaseScan_01_fitting_results.json";
    std::ifstream timewalk_correction_json_file(timewalk_correction_file);
    if (!timewalk_correction_json_file.is_open()) {
        LOG(ERROR) << "Failed to open timewalk correction json file " << timewalk_correction_file;
        return 1;
    }
    json timewalk_correction_json;
    timewalk_correction_json_file >> timewalk_correction_json;
    timewalk_correction_json_file.close();

    std::vector <int> timewalk_correction_uni_channel = timewalk_correction_json["unified_channel_numbers"].get<std::vector<int>>();
    std::vector <double> timewalk_correction_parA = timewalk_correction_json["A"].get<std::vector<double>>();
    std::vector <double> timewalk_correction_parB = timewalk_correction_json["B"].get<std::vector<double>>();
    std::vector <double> timewalk_correction_parC = timewalk_correction_json["C"].get<std::vector<double>>();
    std::vector <double> timewalk_correction_parx0 = timewalk_correction_json["x0"].get<std::vector<double>>();

    std::vector <std::vector <TH2D*>> phase_scan_hists_cut;
    std::vector <TDirectory*> phase_scan_hists_cut_folders;

    std::vector <TH2D*> phase_scan_hists_toa_valid;
    std::vector <TH2D*> phase_scan_toa_hists;
    std::vector <TH2D*> phase_scan_hists_toa_shifted;
    std::vector <TH2D*> phase_scan_adc_toa_correlation;
    std::vector <std::vector <TH1D*>> phase_scan_toa_phase_hists;

    TDirectory *phase_scan_hists_toa_valid_folder;
    TDirectory *phase_scan_toa_hists_folder;
    TDirectory *phase_scan_hists_toa_shifted_folder;
    TDirectory *phase_scan_adc_toa_correlation_folder;
    TDirectory *phase_scan_toa_phase_hists_folder;

    std::unordered_map <int, int> phase_scan_hists_toa_valid_unified_channel_map;

    auto toa_toa_ns_2d_hist = new TH2D("toa_toa_ns_2d_hist", "TOA vs TOA_ns", 1000, 0, 500, 256, 0, 1024);
    toa_toa_ns_2d_hist->GetXaxis()->SetTitle("TOA [ns]");
    toa_toa_ns_2d_hist->GetYaxis()->SetTitle("TOA Value");
    toa_toa_ns_2d_hist->SetStats(0);

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    // * --- Go though all the run files ------------------------------------------------
    // * --------------------------------------------------------------------------------
    if (config_run_numbers.size() != config_phase_settings.size()) {
        LOG(ERROR) << "Run numbers and phase settings do not match!";
        return 1;
    }
    if (config_run_numbers.size() != config_run_files.size()) {
        LOG(ERROR) << "Run numbers and run files do not match!";
        return 1;
    }
    if (config_run_numbers.empty()) {
        LOG(ERROR) << "No run numbers in the configuration file!";
        return 1;
    }
    int fpga_count = -1;
    int machine_gun_samples = -1;

    double toa_target = 88.0;

    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _run_number = config_run_numbers[_run_index];
        auto _run_file = config_run_files[_run_index];
        auto _phase_setting = config_phase_settings[_run_index];
        auto _phase_shifted = _phase_setting - initial_phase;
        if (_phase_shifted < 0) {
            _phase_shifted += 16;
        }
        LOG(INFO) << "Run " << _run_number << " with phase setting " << _phase_setting;

        // * --- Read the input file ----------------------------------------------------
        // * ----------------------------------------------------------------------------
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

        if (fpga_count != 0 && machine_gun_samples != 0){
            if (phase_scan_hists_toa_valid.empty()){
                phase_scan_hists_toa_valid_folder  = output_root->mkdir(("PhaseScanHistsToAValid_Run" + std::to_string(_run_number)).c_str());
                phase_scan_toa_hists_folder     = output_root->mkdir(("PhaseScanTOAHists_Run" + std::to_string(_run_number)).c_str());
                phase_scan_hists_toa_shifted_folder = output_root->mkdir(("PhaseScanHistsTOAShifted_Run" + std::to_string(_run_number)).c_str());
                phase_scan_toa_phase_hists_folder = output_root->mkdir(("PhaseScanTOAPhaseHists_Run" + std::to_string(_run_number)).c_str());
                phase_scan_adc_toa_correlation_folder = output_root->mkdir(("PhaseScanADCTOACorrelation_Run" + std::to_string(_run_number)).c_str());
                int _hist_index = 0;
                for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                    auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++){
                        double _time_max = machine_gun_samples * sample_time;
                        int _hist_x_bins = _time_max / hist_x_bin_size;
                        auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);

                        auto *_hist = new TH2D(("phase_scan_hist_toa_valid_chn" + std::to_string(_unified_valid_channel_number)).c_str(), ("Phase Scan Channel # " + std::to_string(_unified_valid_channel_number)).c_str(), _hist_x_bins, 0, _time_max, hist_y_bins, hist_y_min, hist_y_max);
                        phase_scan_hists_toa_valid.push_back(_hist);
                        phase_scan_hists_toa_valid_unified_channel_map[_unified_valid_channel_number] = _hist_index;
                        _hist->SetDirectory(phase_scan_hists_toa_valid_folder);

                        auto *_hist_toa = new TH2D(("phase_scan_toa_hist_chn" + std::to_string(_unified_valid_channel_number)).c_str(), ("Phase Scan TOA Channel # " + std::to_string(_unified_valid_channel_number)).c_str(), _hist_x_bins, 0, _time_max, hist_y_bins, hist_y_min, hist_y_max);
                        phase_scan_toa_hists.push_back(_hist_toa);
                        _hist_toa->SetDirectory(phase_scan_toa_hists_folder);

                        auto *_hist_shifted = new TH2D(("phase_scan_hist_toa_shifted_chn" + std::to_string(_unified_valid_channel_number)).c_str(), ("Phase Scan TOA Shifted Channel # " + std::to_string(_unified_valid_channel_number)).c_str(), _hist_x_bins, 0, _time_max, hist_y_bins, hist_y_min, hist_y_max);
                        phase_scan_hists_toa_shifted.push_back(_hist_shifted);
                        _hist_shifted->SetDirectory(phase_scan_hists_toa_shifted_folder);

                        auto *_hist_adc_toa_correlation = new TH2D(("phase_scan_adc_toa_correlation_chn" + std::to_string(_unified_valid_channel_number)).c_str(), ("Phase Scan ADC TOA Correlation Channel # " + std::to_string(_unified_valid_channel_number)).c_str(), 256, 0, 1024, (toa_ns_max - toa_ns_min) / toa_ns_bin_size, toa_ns_min, toa_ns_max);
                        phase_scan_adc_toa_correlation.push_back(_hist_adc_toa_correlation);
                        _hist_adc_toa_correlation->SetDirectory(phase_scan_adc_toa_correlation_folder);

                        auto *_phase_hists_toa = new std::vector <TH1D*>();
                        int _hist_toa_time_bins = int((toa_ns_max - toa_ns_min) *2 / toa_ns_bin_size);
                        for (int _phase_index = 0; _phase_index < config_phase_settings.size(); _phase_index++) {
                            auto _phase = config_phase_settings[_phase_index];
                            auto *_hist_phase = new TH1D(("phase_scan_toa_phase_hist_chn" + std::to_string(_unified_valid_channel_number) + "_phase" + std::to_string(_phase)).c_str(), ("Phase Scan TOA Phase Channel # " + std::to_string(_unified_valid_channel_number) + " Phase " + std::to_string(_phase)).c_str(), _hist_toa_time_bins, toa_ns_min, toa_ns_max);
                            _phase_hists_toa->push_back(_hist_phase);
                            _hist_phase->SetDirectory(phase_scan_toa_phase_hists_folder);
                        }
                        phase_scan_toa_phase_hists.push_back(*_phase_hists_toa);

                        toa_thresholds.push_back(1024);

                        _hist_index++;
                    }
                }
                for (int _cut_index = 0; _cut_index < hist_cut_values.size(); _cut_index++) {
                    TDirectory *_phase_scan_hists_cut_folder = output_root->mkdir(("PhaseScanHistsCut_Run" + std::to_string(_run_number) + "_Cut" + std::to_string(hist_cut_values[_cut_index])).c_str());
                    phase_scan_hists_cut_folders.push_back(_phase_scan_hists_cut_folder);
                    phase_scan_hists_cut.push_back(std::vector <TH2D*>());
                    int _hist_index = 0;
                    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                        auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                        for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++){
                            double _time_max = machine_gun_samples * sample_time;
                            int _hist_x_bins = _time_max / hist_x_bin_size;
                            auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);

                            auto *_hist = new TH2D(("phase_scan_hist_cut_chn" + std::to_string(_unified_valid_channel_number)).c_str(), ("Phase Scan Channel # " + std::to_string(_unified_valid_channel_number)).c_str(), _hist_x_bins, 0, _time_max, hist_y_bins, hist_y_min, hist_y_max);
                            
                            phase_scan_hists_cut[_cut_index].push_back(_hist);
                            _hist->SetDirectory(_phase_scan_hists_cut_folder);
                            _hist_index++;
                        }
                    }
                }
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

        int _hamming_code_error_count = 0;
        int _hamming_code_error_count_total = 0;
        // * --- Read the events --------------------------------------------------------
        // * ----------------------------------------------------------------------------
        for (int _entry = 0; _entry < _entry_max; _entry++) {
            _input_tree->GetEntry(_entry);
           
            for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                auto _fpga_id    = _legal_fpga_id_list[_fpga_index];
                auto _timestamp  = branch_timestamps_list[_fpga_index][0];
                auto _daqh_list  = branch_daqh_list_list[_fpga_index];
                auto _tc_list    = branch_tc_list_list[_fpga_index];
                auto _tp_list    = branch_tp_list_list[_fpga_index];
                auto _val0_list  = branch_val0_list_list[_fpga_index];
                auto _val1_list  = branch_val1_list_list[_fpga_index];
                auto _val2_list  = branch_val2_list_list[_fpga_index];

                // * --- Process the event ------------------------------------------------
                // * ------------------------------------------------------------------------
                std::vector <bool> _hamming_code_pass_list;
                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    bool _hamming_code_pass = true;
                    for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
                        _hamming_code_error_count_total++;
                        auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
                        auto _h1h2h3 = (_daqh >> 4) & 0x7;
                        if (_h1h2h3 != 0x00){
                            _hamming_code_error_count++;
                            _hamming_code_pass = false;
                        }
                    }
                    _hamming_code_pass_list.push_back(_hamming_code_pass);
                }

                for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                    auto _channel_valid = get_valid_fpga_channel(_channel_index);
                    if (_channel_valid == -1){
                        continue;
                    }
                    auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                    auto _hist_index = phase_scan_hists_toa_valid_unified_channel_map[_unified_valid_channel_number];
                    UInt_t _val0_max = 0;
                    UInt_t _val2_max = 0;
                    bool _multiple_val2_max = false;
                    int _val2_max_machine_gun_sample = -1;
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _value = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        if (_value > _val0_max){
                            _val0_max = _value;
                        }
                        auto _toa_value = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        // if (_toa_value > _val2_max && _val2_max == 0){
                        if (_toa_value > _val2_max && _val2_max == 0 && _toa_value % 256 != 0){
                            if (_val2_max > 0){
                                _multiple_val2_max = true;
                            }
                            _val2_max = _toa_value;
                            _val2_max_machine_gun_sample = _sample_index;
                        }
                    }
                    if (_val2_max > 0) {
                        if (_val0_max < toa_thresholds[_unified_valid_channel_number]){
                            toa_thresholds[_unified_valid_channel_number] = _val0_max;
                        }
                    }
                    double _toa_value_ns = _val2_max_machine_gun_sample * sample_time + decode_toa_value_ns(_val2_max);
                    // if (_val2_max > 988){
                    //     _toa_value_ns -= 6.25;
                    // }
                    if (_val2_max > 730){
                        _toa_value_ns -= 25.0;
                    }
                    // if (_val2_max > 476){
                    //     _toa_value_ns -= 6.25;
                    // }
                    // if (_val2_max > 220){
                    //     _toa_value_ns -= 6.25;
                    // }
                    
                    // _toa_value_ns += _phase_setting * phase_time;
                    // if(_toa_value_ns > 115){
                    //     _toa_value_ns -= 25.0;
                    // }
                    // if(_toa_value_ns > 115){
                    //     _toa_value_ns -= 25.0;
                    // }
                    // for (int _subtract_time = 0; _subtract_time < 10; _subtract_time++){
                    //     if (_toa_value_ns > 115){
                    //         _toa_value_ns -= 25.0;
                    //     }
                    // }
                    _toa_value_ns += _phase_setting * phase_time;
                    if (_unified_valid_channel_number == 65){
                        toa_toa_ns_2d_hist->Fill(_toa_value_ns, _val2_max);
                    }
                    auto _hist_toa = phase_scan_toa_hists[_hist_index];
                    if (_val2_max > 0){
                        _hist_toa->Fill(_toa_value_ns, _val2_max);
                    }
                    if (enable_timewalk_correction){
                        auto _timewalk_correction_index = std::find(timewalk_correction_uni_channel.begin(), timewalk_correction_uni_channel.end(), _unified_valid_channel_number);
                        if (_timewalk_correction_index != timewalk_correction_uni_channel.end()){
                            auto _timewalk_correction_index_int = std::distance(timewalk_correction_uni_channel.begin(), _timewalk_correction_index);
                            auto _timewalk_correction_parA = timewalk_correction_parA[_timewalk_correction_index_int];
                            auto _timewalk_correction_parB = timewalk_correction_parB[_timewalk_correction_index_int];
                            auto _timewalk_correction_parC = timewalk_correction_parC[_timewalk_correction_index_int];
                            auto _timewalk_correction_parx0 = timewalk_correction_parx0[_timewalk_correction_index_int];
                            double _x_val = _val0_max;
                            auto _timewalk_correction = ExpShifted(&_x_val, new double[4]{_timewalk_correction_parA, _timewalk_correction_parx0, _timewalk_correction_parB, _timewalk_correction_parC});
                            _toa_value_ns -= _timewalk_correction;
                            _toa_value_ns += 90.0;
                        }
                    }
                    auto _hist_adc_toa_correlation = phase_scan_adc_toa_correlation[_hist_index];
                    _hist_adc_toa_correlation->Fill(_val0_max, _toa_value_ns);
                    // double _toa_value_ns = 2 * sample_time + decode_toa_value_ns(_val2_max) + _phase_setting * phase_time;
                    double _toa_target_shifted = toa_target;
                    double _toa_correction = _toa_target_shifted - _toa_value_ns;
                    if (_val2_max > 0){
                        phase_scan_toa_phase_hists[_hist_index][_phase_setting]->Fill(_toa_value_ns);
                    }
                    if (_multiple_val2_max){
                        LOG(WARNING) << "Multiple TOA max values found for channel " << _unified_valid_channel_number;
                        continue;   
                    }
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        if (!_hamming_code_pass_list[_sample_index]){
                            continue;
                        }
                        auto _time = _sample_index * sample_time + _phase_shifted * phase_time;
                        auto _time_shifted = _time + _toa_correction;
                        auto _value = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _toa_value = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _hist = phase_scan_hists_toa_valid[_hist_index];
                        auto _hist_shifted = phase_scan_hists_toa_shifted[_hist_index];
                        if (_val2_max > 0){
                            _hist->Fill(_time, _value);
                            _hist_shifted->Fill(_time_shifted, _value);
                        }
                        for (int _cut_index = 0; _cut_index < hist_cut_values.size(); _cut_index++) {
                            if (_val0_max > hist_cut_values[_cut_index]){
                                auto _hist = phase_scan_hists_cut[_cut_index][_hist_index];
                                _hist->Fill(_time, _value);
                            }
                        }
                    }
                }
            }
        }
        if (_hamming_code_error_count > 0){
            LOG(WARNING) << "Hamming code error count: " << _hamming_code_error_count << " (" << (double)_hamming_code_error_count / _hamming_code_error_count_total * 100 << "%)";
        }

        _input_root->Close();
    }

    // * --- Write the output file -------------------------------------------------------
    // * ------------------------------------------------------------------------------------
    phase_scan_hists_toa_valid_folder->cd();

    for (auto _hist : phase_scan_hists_toa_valid) {
        // format the histogram
        _hist->GetXaxis()->SetTitle("Time [ns]");
        _hist->GetYaxis()->SetTitle("ADC Value");
        _hist->SetStats(0);
        _hist->Write();
    }

    LOG(INFO) << "Drawing global phase scan canvas";
    GlobalChannelPainter *global_painter = nullptr;
    if (enable_focal_mapping){
        global_painter = new GlobalChannelPainter("data/SPS_2024/config/focalh_mapping.json", "data/SPS_2024/config/h2gcroc_mapping.json");
    } else {
        global_painter = new GlobalChannelPainter("data/DESY_2025/config/EEEMCal_Mapping_DESY_2025.json");
    }
    global_painter->draw_global_channel_hists2D(phase_scan_hists_toa_valid, phase_scan_hists_toa_valid_unified_channel_map, "PhaseScanHistsNoCut", "Phase Scan Channel # ");
    auto global_phase_scan_toa_valid_canvas = global_painter->get_canvas();

    output_root->cd();
    global_phase_scan_toa_valid_canvas->Write();

    phase_scan_toa_hists_folder->cd();

    for (auto _hist : phase_scan_toa_hists) {
        // format the histogram
        _hist->GetXaxis()->SetTitle("Time [ns]");
        _hist->GetYaxis()->SetTitle("TOA Value");
        _hist->SetStats(0);
        _hist->Write();
    }

    LOG(INFO) << "Drawing global toa phase scan canvas";
    global_painter->draw_global_channel_hists2D(phase_scan_toa_hists,phase_scan_hists_toa_valid_unified_channel_map, "PhaseScanTOAHists", "Phase Scan TOA Channel # ");
    auto global_phase_scan_toa_canvas = global_painter->get_canvas();

    output_root->cd();
    global_phase_scan_toa_canvas->Write();

    phase_scan_hists_toa_shifted_folder->cd();

    for (auto _hist : phase_scan_hists_toa_shifted) {
        // format the histogram
        _hist->GetXaxis()->SetTitle("Time [ns]");
        _hist->GetYaxis()->SetTitle("ADC Value");
        _hist->SetStats(0);
        _hist->Write();
    }

    LOG(INFO) << "Drawing global toa shifted phase scan canvas";
    global_painter->draw_global_channel_hists2D(phase_scan_hists_toa_shifted, phase_scan_hists_toa_valid_unified_channel_map, "PhaseScanHistsTOAShifted", "Phase Scan TOA Shifted Channel # ");
    auto global_phase_scan_toa_shifted_canvas = global_painter->get_canvas();

    output_root->cd();
    global_phase_scan_toa_shifted_canvas->Write();

    for (int _cut_index = 0; _cut_index < hist_cut_values.size(); _cut_index++) {
        phase_scan_hists_cut_folders[_cut_index]->cd();
        for (auto _hist : phase_scan_hists_cut[_cut_index]) {
            // format the histogram
            _hist->GetXaxis()->SetTitle("Time [ns]");
            _hist->GetYaxis()->SetTitle("ADC Value");
            _hist->SetStats(0);
            _hist->Write();
        }
    }

    LOG(INFO) << "Drawing global toa time phase scan canvas";
    std::vector <EColor> color_array = {kBlue, kRed, kGreen, kBlack, kMagenta, kCyan, kYellow, kOrange, kViolet, kTeal, kAzure, kGray, kSpring, kPink};
    std::vector <std::string> phase_labels;
    for (int i = 0; i < config_phase_settings.size(); i++){
        phase_labels.push_back("Phase " + std::to_string(config_phase_settings[i]));
    }
    double phase_scan_toa_hist_max = 0;
    for (auto _hist_vector : phase_scan_toa_phase_hists){
        for (auto _hist : _hist_vector){
            if (_hist->GetMaximum() > phase_scan_toa_hist_max){
                phase_scan_toa_hist_max = _hist->GetMaximum();
            }
        }
    }
    for (auto _hist_vector : phase_scan_toa_phase_hists){
        for (auto _hist : _hist_vector){
            _hist->GetYaxis()->SetRangeUser(0, phase_scan_toa_hist_max * 1.1);
        }
    }
    global_painter->draw_global_channel_hists1D_group(phase_scan_toa_phase_hists, phase_scan_hists_toa_valid_unified_channel_map, "PhaseScanTOAPhaseHists", "Phase Scan TOA Phase Channel # ", color_array, phase_labels);
    auto global_phase_scan_toa_phase_canvas = global_painter->get_canvas();

    output_root->cd();
    global_phase_scan_toa_phase_canvas->Write();

    phase_scan_toa_phase_hists_folder->cd();
    
    for (auto _hist_vector : phase_scan_toa_phase_hists){
        for (auto _hist : _hist_vector){
            // format the histogram
            _hist->GetXaxis()->SetTitle("Time [ns]");
            _hist->GetYaxis()->SetTitle("TOA Value");
            _hist->SetStats(0);
            _hist->Write();
        }
    }

    output_root->cd();
    toa_toa_ns_2d_hist->Write();

    std::vector <TCanvas*> phase_scan_adc_toa_correlation_canvas;

    for (int _hist_index = 0; _hist_index < phase_scan_adc_toa_correlation.size(); _hist_index++) {
        auto _canvas = new TCanvas(("phase_scan_adc_toa_correlation_chn" + std::to_string(_hist_index)).c_str(), ("Phase Scan ADC TOA Correlation Channel # " + std::to_string(_hist_index)).c_str(), 800, 600);
        phase_scan_adc_toa_correlation_canvas.push_back(_canvas);
    }

    phase_scan_adc_toa_correlation_folder->cd();
    std::vector <TF1*> phase_scan_adc_toa_correlation_fit;

    std::vector <int> fitting_unified_channel_numbers;
    std::vector <double> fitting_results_A;
    std::vector <double> fitting_results_x0;
    std::vector <double> fitting_results_B;
    std::vector <double> fitting_results_C;

    int fit_event_threshold = 100;

    for (int _hist_index = 0; _hist_index < phase_scan_adc_toa_correlation.size(); _hist_index++) {
        auto _hist = phase_scan_adc_toa_correlation[_hist_index];
        phase_scan_adc_toa_correlation_canvas[_hist_index]->cd();
        _hist->SetStats(0);

        auto _hist_profX = _hist->ProfileX();
        _hist_profX->SetMarkerStyle(20);
        _hist_profX->SetMarkerSize(0.5);
        _hist_profX->SetMarkerColor(kBlack);
        _hist_profX->SetLineColor(kBlack);
        _hist_profX->SetTitle("");
        _hist_profX->SetStats(0);
        _hist_profX->Draw("p");

        if (_hist_profX->GetEntries() > fit_event_threshold){
            auto _hist_min = _hist_profX->GetMinimum();
            TF1* _fit = new TF1(("fit_" + std::to_string(_hist_index)).c_str(), ExpShifted, toa_thresholds[_hist_index], 1024, 4);
            _fit->SetParameters(100, toa_thresholds[_hist_index], 0.2, _hist_min);
            _fit->SetParLimits(0, 0, 1000);
            _fit->SetParLimits(1, 0, 1024);
            _fit->SetParLimits(2, 0, 1);
            _fit->SetParLimits(3, 0, 1000);

            _fit->FixParameter(1, toa_thresholds[_hist_index]);
    
            _hist_profX->Fit(("fit_" + std::to_string(_hist_index)).c_str(), "RQN");
            _fit->Draw("same");
            phase_scan_adc_toa_correlation_fit.push_back(_fit);
    
            auto _fit_A = _fit->GetParameter(0);
            auto _fit_x0 = _fit->GetParameter(1);
            auto _fit_B = _fit->GetParameter(2);
            auto _fit_C = _fit->GetParameter(3);
    
            TLegend* _legend = new TLegend(0.65, 0.65, 0.89, 0.89);
            _legend->SetBorderSize(0);
            _legend->SetFillColor(0);
            _legend->SetFillStyle(0);

            // Add vertical line for the threshold
            TLine* _line = new TLine(toa_thresholds[_hist_index], 0, toa_thresholds[_hist_index], 1024);
            _line->SetLineColor(kPink);
            _line->SetLineStyle(2);
            _line->SetLineWidth(2);
            _line->Draw();

            TObject *dummy_obj = new TObject();
            _legend->AddEntry(_hist_profX, "Profile X", "p");
            _legend->AddEntry(_fit, ("A = " + std::to_string(_fit_A)).c_str());
            _legend->AddEntry(dummy_obj, ("x0 = " + std::to_string(_fit_x0)).c_str());
            _legend->AddEntry(dummy_obj, ("B = " + std::to_string(_fit_B)).c_str());
            _legend->AddEntry(dummy_obj, ("C = " + std::to_string(_fit_C)).c_str());

            _legend->Draw();

            fitting_unified_channel_numbers.push_back(_hist_index);
            fitting_results_A.push_back(_fit_A);
            fitting_results_x0.push_back(_fit_x0);
            fitting_results_B.push_back(_fit_B);
            fitting_results_C.push_back(_fit_C);
        }

        _hist_profX->SetDirectory(phase_scan_adc_toa_correlation_folder);

        phase_scan_adc_toa_correlation_canvas[_hist_index]->Update();
        phase_scan_adc_toa_correlation_canvas[_hist_index]->Write();
    }

    LOG(INFO) << "Drawing global adc toa correlation canvas";
    // global_painter->draw_global_channel_canvas(phase_scan_adc_toa_correlation_canvas, phase_scan_hists_toa_valid_unified_channel_map, "PhaseScanADCTOACorrelation", "Phase Scan ADC TOA Correlation Channel # ");
    // auto global_phase_scan_adc_toa_canvas = global_painter->get_canvas();

    // output_root->cd();
    // global_phase_scan_adc_toa_canvas->Write();

    output_root->Close();

    // save the fitting results to a json file
    if (!enable_timewalk_correction){
        json fitting_results;
        fitting_results["unified_channel_numbers"] = fitting_unified_channel_numbers;
        fitting_results["A"] = fitting_results_A;
        fitting_results["x0"] = fitting_results_x0;
        fitting_results["B"] = fitting_results_B;
        fitting_results["C"] = fitting_results_C;

        std::string fitting_results_file = opts.output_file.substr(0, opts.output_file.find_last_of(".")) + "_fitting_results.json";
        std::ofstream fitting_results_stream(fitting_results_file);
        fitting_results_stream << fitting_results.dump(4);
        fitting_results_stream.close();
    }
    
    return 0;
}