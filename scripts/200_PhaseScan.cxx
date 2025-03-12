#include "H2GCROC_Common.hxx"

INITIALIZE_EASYLOGGINGPP

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
    std::vector <UInt_t> hist_cut_values = config_content["ADC Cut"].get<std::vector<UInt_t>>();

    LOG(INFO) << "Phase scan: " << config_info;

    const double phase_time  = 25.0 / 16.0;  // unit: ns
    const double sample_time = 25.0;         // unit: ns
    const int initial_phase = 7;
    const int hist_y_bins = 256;
    const double hist_x_bin_size = phase_time;
    const double hist_y_min = 0;
    const double hist_y_max = 1024;

    std::vector <std::vector <TH2D*>> phase_scan_hists_cut;
    std::vector <TDirectory*> phase_scan_hists_cut_folders;

    std::vector <TH2D*> phase_scan_hists_no_cut;
    std::unordered_map <int, int> phase_scan_hists_no_cut_unified_channel_map;
    TDirectory *phase_scan_hists_no_cut_folder;

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
            if (phase_scan_hists_no_cut.empty()){
                phase_scan_hists_no_cut_folder = output_root->mkdir(("PhaseScanHistsNoCut_Run" + std::to_string(_run_number)).c_str());
                phase_scan_hists_no_cut_folder->cd();
                int _hist_index = 0;
                for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                    auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++){
                        double _time_max = machine_gun_samples * sample_time;
                        int _hist_x_bins = _time_max / hist_x_bin_size;
                        auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);

                        auto *_hist = new TH2D(("phase_scan_hist_no_cut_chn" + std::to_string(_unified_valid_channel_number)).c_str(), ("Phase Scan Channel # " + std::to_string(_unified_valid_channel_number)).c_str(), _hist_x_bins, 0, _time_max, hist_y_bins, hist_y_min, hist_y_max);
                        phase_scan_hists_no_cut.push_back(_hist);

                        phase_scan_hists_no_cut_unified_channel_map[_unified_valid_channel_number] = _hist_index;
                        _hist->SetDirectory(phase_scan_hists_no_cut_folder);
                        _hist_index++;
                    }
                }
                for (int _cut_index = 0; _cut_index < hist_cut_values.size(); _cut_index++) {
                    TDirectory *_phase_scan_hists_cut_folder = output_root->mkdir(("PhaseScanHistsCut_Run" + std::to_string(_run_number) + "_Cut" + std::to_string(hist_cut_values[_cut_index])).c_str());
                    _phase_scan_hists_cut_folder->cd();
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
                // for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                //     bool _hamming_code_pass = true;
                //     for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
                //         _hamming_code_error_count_total++;
                //         auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
                //         auto _h1h2h3 = (_daqh >> 4) & 0x7;
                //         if (_h1h2h3 != 0x00){
                //             _hamming_code_error_count++;
                //             _hamming_code_pass = false;
                //         }
                //     }
                //     if (!_hamming_code_pass){
                //         continue;
                //     }
                //     for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                //         auto _channel_valid = get_valid_fpga_channel(_channel_index);
                //         if (_channel_valid == -1){
                //             continue;
                //         }
                //         auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                //         auto _hist_index = phase_scan_hists_no_cut_unified_channel_map[_unified_valid_channel_number];
                //         auto _time = _sample_index * sample_time + _phase_shifted * phase_time;
                //         auto _value = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                //         auto _hist = phase_scan_hists_no_cut[_hist_index];
                //         _hist->Fill(_time, _value);
                //     }
                // }
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
                    auto _hist_index = phase_scan_hists_no_cut_unified_channel_map[_unified_valid_channel_number];
                    UInt_t _val0_max = 0;
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _value = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        if (_value > _val0_max){
                            _val0_max = _value;
                        }
                    }
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        if (!_hamming_code_pass_list[_sample_index]){
                            continue;
                        }
                        auto _time = _sample_index * sample_time + _phase_shifted * phase_time;
                        auto _value = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _hist = phase_scan_hists_no_cut[_hist_index];
                        _hist->Fill(_time, _value);
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
    phase_scan_hists_no_cut_folder->cd();

    for (auto _hist : phase_scan_hists_no_cut) {
        // format the histogram
        _hist->GetXaxis()->SetTitle("Time [ns]");
        _hist->GetYaxis()->SetTitle("ADC Value");
        _hist->SetStats(0);
        _hist->Write();
    }

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

    output_root->Close();
    return 0;
}