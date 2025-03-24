#include "H2GCROC_Common.hxx"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_json(argc, argv, "1.1");

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
    
    bool enable_working_in_progress = config_content["Working in progress"].get<bool>();
    bool enable_focal_mapping = opts.focal;

    int plot_sum_x_min = config_content["Plot Setting"]["Sum X Min"].get<int>();
    int plot_sum_x_max = config_content["Plot Setting"]["Sum X Max"].get<int>();
    int plot_sum_x_bin = config_content["Plot Setting"]["Sum X Bin"].get<int>();

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

    const double adc_sum_hist_min = double(plot_sum_x_min);
    const double adc_sum_hist_max = double(plot_sum_x_max);
    const int adc_sum_hist_bins = plot_sum_x_bin;

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

    std::unordered_map <int, int> channel_wise_hists_unified_channel_map;

    std::vector <std::vector <TH2D*>> run_chn_toatime_toacode_hist2d; // [run_index][channel_index]
    TDirectory *run_chn_toatime_toacode_hist2d_folder = output_root->mkdir("RunChnToaTimeToaCodeHist2D");

    std::vector <std::vector <TH2D*>> run_chn_toatime_adcmax_hist2d; // [run_index][channel_index]
    TDirectory *run_chn_toatime_adcmax_hist2d_folder = output_root->mkdir("RunChnToaTimeAdcMaxHist2D");

    std::vector <TH1D*> run_adc_sum_hist1d;
    TDirectory *run_adc_sum_hist1d_folder = output_root->mkdir("RunAdcSumHist1D");

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
            // ----------------------------------------------------------------
            if (run_chn_hist2d.size() == 0) {
                int _hist_index = 0;
                for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                    auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++) {
                        auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);
                        channel_wise_hists_unified_channel_map[_unified_valid_channel_number] = _hist_index;

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
                    auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_hist_index];

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

        for (int _entry = 0; _entry < entry_max; _entry++) {
            _input_tree->GetEntry(_entry);

            double _run_adc_sum = 0;

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

                std::vector <int> _channel_val1_max_list;
                std::vector <int> _channel_val2_max_list;
                std::vector <int> _channel_val1_max_index_list;
                std::vector <int> _channel_val2_max_index_list;

                for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                    auto _channel_valid = get_valid_fpga_channel(_channel_index);
                    if (_channel_valid == -1){
                        continue;
                    }
                    auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                    auto _hist_index = channel_wise_hists_unified_channel_map[_unified_valid_channel_number];
                    int _val1_max = -1;
                    int _val1_max_index = -1;
                    int _val2_max = -1;
                    int _val2_max_index = -1;

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val1 = _val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val2 = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
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
                    _channel_val1_max_list.push_back(_val1_max);
                    _channel_val2_max_list.push_back(_val2_max);
                    _channel_val1_max_index_list.push_back(_val1_max_index);
                    _channel_val2_max_index_list.push_back(_val2_max_index);
                }

                for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                    auto _channel_valid = get_valid_fpga_channel(_channel_index);
                    auto _val1_max = _channel_val1_max_list[_channel_index];
                    auto _val2_max = _channel_val2_max_list[_channel_index];
                    auto _val1_max_index = _channel_val1_max_index_list[_channel_index];
                    auto _val2_max_index = _channel_val2_max_index_list[_channel_index];

                    if (_channel_valid == -1){
                        continue;
                    }
                    auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                    auto _hist_index = channel_wise_hists_unified_channel_map[_unified_valid_channel_number];
                    auto _channel_sample_hist2d = run_chn_hist2d[_run_index][_hist_index];
                    auto _channel_toatime_toacode_hist2d = run_chn_toatime_toacode_hist2d[_run_index][_hist_index];
                    auto _channel_toatime_adcmax_hist2d = run_chn_toatime_adcmax_hist2d[_run_index][_hist_index];

                    int _val0_max = -1;
                    int _pedestal_event = -1;

                    double _channel_adc_value = 0;

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        if (_sample_index == 0) {
                            _pedestal_event = int(_val0);
                        }
                        if (int(_val0) > _val0_max) {
                            _val0_max = int(_val0);
                        }

                        if (_hamming_code_pass_list[_sample_index] && _val2_max > 0) {
                            _channel_sample_hist2d->Fill(_sample_index * sample_time, _val0);
                        }
                    }

                    // * --- Calculate pedestal subtracted ADC value -------------------------
                    // * --------------------------------------------------------------------
                    double _pedestal_mean1       = channel_pede_mean1[_unified_valid_channel_number];
                    double _pedestal_error1      = channel_pede_error1[_unified_valid_channel_number];
                    double _pedestal_mean2       = channel_pede_mean2[_unified_valid_channel_number]; 
                    double _pedestal_error2      = channel_pede_error2[_unified_valid_channel_number];
                    int _pedestal_peak_counts    = int(channel_pede_peak_counts[_unified_valid_channel_number]);

                    if (_val2_max > 0) {
                        double _toa_time = decode_toa_value_ns(_val2_max) + _val2_max_index * sample_time;
                        _channel_toatime_toacode_hist2d->Fill(_toa_time, _val2_max);
                        _channel_toatime_adcmax_hist2d->Fill(_val0_max, _toa_time);
                        // LOG(DEBUG) << "ADC Max: " << _val0_max << " TOA Time: " << _toa_time << " TOA Code: " << _val2_max;
                    }

                    if (_pedestal_peak_counts == 1){
                        auto _pedestal = _pedestal_mean1;
                        _channel_adc_value = _val0_max - _pedestal;
                        _run_adc_sum += _channel_adc_value;
                    } else if (_pedestal_peak_counts == 2){
                        auto _pedestal1_dist = std::abs(_pedestal_event - _pedestal_mean1);
                        auto _pedestal2_dist = std::abs(_pedestal_event - _pedestal_mean2);
                        auto _pedestal = 0;
                        if (_pedestal1_dist < _pedestal2_dist){
                            _pedestal = _pedestal_mean1;
                        } else {
                            _pedestal = _pedestal_mean2;
                        }
                        _channel_adc_value = _val0_max - _pedestal;
                        _run_adc_sum += _channel_adc_value;
                    }
                }
            }
            _run_adc_sum_hist1d->Fill(_run_adc_sum);
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
            auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
            auto _hist2d = _channel_hist2d[_channel_index];
            _hist2d->SetStats(0);
            _hist2d->Write();
        }
        if (_run_index == run_index_highest_energy || _run_index == run_index_lowest_energy){
            global_painter->draw_global_channel_hists2D(_channel_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnSampleHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
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
            auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
            auto _hist2d = _channel_toatime_toacode_hist2d[_channel_index];
            _hist2d->GetXaxis()->SetTitle("TOA Time [ns]");
            _hist2d->GetYaxis()->SetTitle("TOA Code");
            _hist2d->SetStats(0);
            _hist2d->Write();
        }
        if (_run_index == run_index_lowest_energy || _run_index == run_index_highest_energy){
            global_painter->draw_global_channel_hists2D(_channel_toatime_toacode_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnToaTimeToaCodeHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
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
            auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
            auto _hist2d = _channel_toatime_adcmax_hist2d[_channel_index];
            _hist2d->GetXaxis()->SetTitle("ADC Max");
            _hist2d->GetYaxis()->SetTitle("TOA Time [ns]");
            _hist2d->SetStats(0);
            _hist2d->Write();
        }
        if (_run_index == run_index_lowest_energy || _run_index == run_index_highest_energy){
            global_painter->draw_global_channel_hists2D(_channel_toatime_adcmax_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnToaTimeAdcMaxHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
            auto _canvas = global_painter->get_canvas();
            _canvas->Print(pdf_file_name.c_str());
            output_root->cd();
            _canvas->Write();
        }
    }

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
        std::vector <double> _adc_sum_fit_range = {1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; // unit: sigma
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