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

    std::vector <int> config_run_numbers            = config_content["Run Numbers"].get<std::vector<int>>();
    std::vector <std::string> config_run_files      = config_content["Run Files"].get<std::vector<std::string>>();
    std::vector <double> config_beam_energies       = config_content["Beam Energies"].get<std::vector<double>>();
    std::vector <std::string> config_beam_particles = config_content["Beam Particles"].get<std::vector<std::string>>();
    std::vector <int> config_target_machineguns     = config_content["MachineGuns"].get<std::vector<int>>();
    std::vector <std::string> config_plot_info      = config_content["Plot Info"].get<std::vector<std::string>>();
    bool enable_working_in_progress = config_content["Working in progress"].get<bool>();
    bool enable_focal_mapping = opts.focal;

    const int channel_adc_hist_bins = 256;
    const double channel_adc_hist_min = 0;
    const double channel_adc_hist_max = 1024;

    const double sample_time = 25.0; // unit: ns

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
            std::vector <TH2D*> _chn_hist2d;
            if (run_chn_hist2d.size() == 0) {
                int _hist_index = 0;
                for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                    auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++) {
                        auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);
                        auto *_hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)).c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number)).c_str(), machine_gun_samples, 0, _max_time, channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max);
                        _chn_hist2d.push_back(_hist2d);
                        _hist2d->SetDirectory(run_chn_sample_hist2d_folder);
                        channel_wise_hists_unified_channel_map[_unified_valid_channel_number] = _hist_index;
                        _hist_index++;
                    }
                }
                run_chn_hist2d.push_back(_chn_hist2d);
            } else {
                for (int _hist_index = 0; _hist_index < run_chn_hist2d[0].size(); _hist_index++) {
                    auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_hist_index];
                    auto *_hist2d = new TH2D(("Run_"+std::to_string(_run_number)+"_Chn_"+std::to_string(_unified_valid_channel_number)).c_str(), ("Run " + std::to_string(_run_number) + " Channel " + std::to_string(_unified_valid_channel_number)).c_str(), machine_gun_samples, 0, _max_time, channel_adc_hist_bins, channel_adc_hist_min, channel_adc_hist_max);
                    _chn_hist2d.push_back(_hist2d);
                    _hist2d->SetDirectory(run_chn_sample_hist2d_folder);
                }
                run_chn_hist2d.push_back(_chn_hist2d);
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

        for (int _entry = 0; _entry < entry_max; _entry++) {
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

                for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                    auto _channel_valid = get_valid_fpga_channel(_channel_index);
                    if (_channel_valid == -1){
                        continue;
                    }
                    auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                    auto _hist_index = channel_wise_hists_unified_channel_map[_unified_valid_channel_number];
                    auto _channel_sample_hist2d = run_chn_hist2d[_run_index][_hist_index];

                    UInt_t _val0_max = 0;
                    UInt_t _val1_max = 0;
                    UInt_t _val2_max = 0;

                    double _channel_adc_value = 0;

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val1 = _val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val2 = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        if (_val0 > _val0_max)  _val0_max = _val0;
                        if (_val1 > _val1_max)  _val1_max = _val1;
                        if (_val2 > _val2_max)  _val2_max = _val2;

                        if (_hamming_code_pass_list[_sample_index]){
                            _channel_sample_hist2d->Fill(_sample_index * sample_time, _val0);
                        }
                    }
                }
            }
        }
        _input_root->Close();
    }

    std::vector <EColor> color_list = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kTeal, kAzure, kGray, kPink, kSpring, kBrown};
    std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".pdf";

    // * --- Save the histograms --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _channel_hist2d = run_chn_hist2d[_run_index];
        run_chn_sample_hist2d_folder->cd();
        for (int _channel_index = 0; _channel_index < _channel_hist2d.size(); _channel_index++) {
            auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
            LOG(INFO) << "Saving channel " << _unified_valid_channel_number;
            auto _hist2d = _channel_hist2d[_channel_index];
            _hist2d->SetStats(0);
            // set logz
            _hist2d->SetContour(100);
            _hist2d->SetOption("colz");
            LOG(INFO) << "Entries: " << _hist2d->GetEntries();
            _hist2d->Write();
        }

        global_painter->draw_global_channel_hists2D(_channel_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnSampleHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
        auto _canvas = global_painter->get_canvas();
        if (_run_index == 0) {
            _canvas->SaveAs((pdf_file_name + "(").c_str());
        } else if (_run_index == config_run_numbers.size() - 1) {
            _canvas->SaveAs((pdf_file_name + ")").c_str());
        } else {
            _canvas->SaveAs(pdf_file_name.c_str());
        }
        output_root->cd();
        _canvas->Write();
    }

    output_root->Close();
    LOG(INFO) << "Output file " << opts.output_file << " has been saved.";
    LOG(INFO) << "PDF file " << pdf_file_name << " has been saved.";

    return 0;
}