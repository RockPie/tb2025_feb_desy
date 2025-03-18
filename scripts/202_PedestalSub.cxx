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

    const int channel_pedestal_adc_hist_bins = 256;
    const double channel_pedestal_adc_hist_min = 0;
    const double channel_pedestal_adc_hist_max = 256;

    const double gaussian_fit_max_sigma = 20.0;
    const double gaussian_fit_valid_min_mean = 0.5;

    const double good_fit_chi2ndf = 15.0;
    const double okay_fit_chi2ndf = 45.0;
    const double dual_peak_fit_chi2ndf_threshold = 5.0;

    // Set up the global channel painter mapping
    GlobalChannelPainter *global_painter = nullptr;
    if (enable_focal_mapping){
        global_painter = new GlobalChannelPainter("data/SPS_2024/config/focalh_mapping.json", "data/SPS_2024/config/h2gcroc_mapping.json");
    } else {
        global_painter = new GlobalChannelPainter("data/DESY_2025/config/EEEMCal_Mapping_DESY_2025.json");
    }
    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    TDirectory *channel_wise_hists_folder = output_root->mkdir("ChannelWiseHists");
    std::unordered_map <int, int> channel_wise_hists_unified_channel_map;
    std::vector <std::vector <TH1D*>> channel_wise_adc_hists;

    std::vector <int> channel_pedestal_peak_counts;

    std::vector <int> channel_single_peak_uni_index_list;
    std::vector <double> channel_single_peak_mean_list;
    std::vector <double> channel_single_peak_error_list;

    std::vector <int> channel_dual_peak_uni_index_list;
    std::vector <double> channel_dual_peak_mean1_list;
    std::vector <double> channel_dual_peak_mean2_list;
    std::vector <double> channel_dual_peak_error1_list;
    std::vector <double> channel_dual_peak_error2_list;

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
            std::vector <TH1D*> _channel_pede_hists;
            if (channel_wise_adc_hists.size() == 0) {
                int _hist_index = 0;
                for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
                    auto _fpga_id = _legal_fpga_id_list[_fpga_index];
                    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID; _valid_channel_index++) {
                        auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _valid_channel_index);
                        auto *_channel_pede_hist = new TH1D(("ChannelADC_" + std::to_string(_unified_valid_channel_number)).c_str(), ("Channel ADC " + std::to_string(_unified_valid_channel_number)).c_str(), channel_pedestal_adc_hist_bins, channel_pedestal_adc_hist_min, channel_pedestal_adc_hist_max);
                        _channel_pede_hists.push_back(_channel_pede_hist);
                        channel_wise_hists_unified_channel_map[_unified_valid_channel_number] = _hist_index;
                        channel_pedestal_peak_counts.push_back(-1);
                        _channel_pede_hist->SetDirectory(channel_wise_hists_folder);
                        _hist_index++;
                    }
                }
                channel_wise_adc_hists.push_back(_channel_pede_hists);
            } else {
                for (int _hist_index = 0; _hist_index < channel_wise_adc_hists[0].size(); _hist_index++) {
                    auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_hist_index];
                    auto *_channel_pede_hist = new TH1D(("ChannelADC_" + std::to_string(_unified_valid_channel_number)).c_str(), ("Channel ADC " + std::to_string(_unified_valid_channel_number)).c_str(), channel_pedestal_adc_hist_bins, channel_pedestal_adc_hist_min, channel_pedestal_adc_hist_max);
                    _channel_pede_hists.push_back(_channel_pede_hist);
                    channel_pedestal_peak_counts.push_back(-1);
                    _channel_pede_hist->SetDirectory(channel_wise_hists_folder);
                }
                channel_wise_adc_hists.push_back(_channel_pede_hists);
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
                    auto _channel_pede_hist = channel_wise_adc_hists[_run_index][_hist_index];

                    UInt_t _val0_max = 0;
                    UInt_t _val1_max = 0;
                    UInt_t _val2_max = 0;

                    double _channel_adc_value = 0;

                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val1 = _val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        auto _val2 = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
                        if (_val0 > _val0_max){
                            _val0_max = _val0;
                        }
                        if (_val1 > _val1_max){
                            _val1_max = _val1;
                        }
                        if (_val2 > _val2_max){
                            _val2_max = _val2;
                        }

                        // check if the sample index is in the machine gun list
                        auto _machine_gun_sample_index = std::find(config_target_machineguns.begin(), config_target_machineguns.end(), _sample_index);
                        if (_machine_gun_sample_index == config_target_machineguns.end()){
                            continue;
                        } else {
                            if (_hamming_code_pass_list[_sample_index]){
                                _channel_adc_value += _val0;
                            }
                        }

                    }
                    // add the val0 max if -1 is in the list
                    auto _machine_gun_sample_index = std::find(config_target_machineguns.begin(), config_target_machineguns.end(), -1);
                    if (_machine_gun_sample_index != config_target_machineguns.end()){
                        _channel_adc_value += _val0_max;
                    }

                    _channel_pede_hist->Fill(_channel_adc_value);
                }
            }
        }
        _input_root->Close();
    }

    std::vector <EColor> color_list = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kTeal, kAzure, kGray, kPink, kSpring, kBrown};
    std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".pdf";
    
    // * --- Save the histograms --------------------------------------------------------
    // * --------------------------------------------------------------------------------
    channel_wise_hists_folder->cd();
    std::vector <TCanvas*> channel_wise_adc_canvas;
    for (int _hist_index = 0; _hist_index < channel_wise_adc_hists[0].size(); _hist_index++) {
        auto _uni_channel_number = channel_wise_hists_unified_channel_map[_hist_index];
        auto _channel_adc_canvas = new TCanvas(("Pedestal Channel_" + std::to_string(_uni_channel_number)).c_str(), ("Channel " + std::to_string(_uni_channel_number) + " ADC").c_str(), 400, 400);
        channel_wise_adc_canvas.push_back(_channel_adc_canvas);

        _channel_adc_canvas->cd();
        auto _legend = new TLegend(0.55, 0.1, 0.89, 0.89);
        _legend->SetFillColor(0);
        _legend->SetLineColor(0);
        _legend->SetBorderSize(0);
        _legend->SetLineColorAlpha(0, 0);
        _legend->SetTextSize(0.035);
        int max_count = 0;
        for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
            if (channel_wise_adc_hists[_run_index][_hist_index]->GetMaximum() > max_count){
                max_count = channel_wise_adc_hists[_run_index][_hist_index]->GetMaximum();
            }
        }

        std::vector <double> single_gaussian_fit_mean;
        std::vector <double> single_gaussian_fit_sigma;
        std::vector <double> single_gaussian_fit_chi2;

        std::vector <double> dual_gaussian_fit_mean1;
        std::vector <double> dual_gaussian_fit_sigma1;
        std::vector <double> dual_gaussian_fit_mean2;
        std::vector <double> dual_gaussian_fit_sigma2;
        std::vector <double> dual_gaussian_fit_chi2;

        for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
            auto _run_number    = config_run_numbers[_run_index];
            auto _beam_energy   = config_beam_energies[_run_index];
            auto _beam_particle = config_beam_particles[_run_index];
            auto _channel_pede_hist = channel_wise_adc_hists[_run_index][_hist_index];
            
            _channel_pede_hist->SetStats(0);
            _channel_pede_hist->GetXaxis()->SetTitle("ADC Value");
            _channel_pede_hist->GetYaxis()->SetTitle("Counts");
            _channel_pede_hist->GetYaxis()->SetRangeUser(0, max_count * 1.1);
            _channel_pede_hist->SetLineColorAlpha(color_list[_run_index % color_list.size()], 0.5);
            if (_run_index == 0){
                _channel_pede_hist->Draw();
            } else {
                _channel_pede_hist->Draw("SAME");
            }
            
            _legend->AddEntry(_channel_pede_hist, ("Run " + std::to_string(_run_number) + " " + std::to_string(int(_beam_energy)) + " GeV").c_str(), "l");
            _channel_adc_canvas->Update();
            // do a gaussian fit
            TF1 *gaus_fit = new TF1("gaus_fit", "gaus", channel_pedestal_adc_hist_min, channel_pedestal_adc_hist_max);
            // set the initial parameters
            double _initial_mean = _channel_pede_hist->GetMean();
            double _initial_sigma = 5.0;
            double _initial_amplitude = _channel_pede_hist->GetMaximum();
            gaus_fit->SetParameter(0, _initial_amplitude);
            gaus_fit->SetParameter(1, _initial_mean);
            gaus_fit->SetParameter(2, _initial_sigma);
            

            gaus_fit->SetParLimits(1, channel_pedestal_adc_hist_min, channel_pedestal_adc_hist_max);
            gaus_fit->SetParLimits(2, 0.1, gaussian_fit_max_sigma);

            _channel_pede_hist->Fit("gaus_fit", "RQN+");
            gaus_fit->SetLineColorAlpha(color_list[_run_index % color_list.size()], 0.5);
            
            double fit_mean  = gaus_fit->GetParameter(1);
            double fit_sigma = gaus_fit->GetParameter(2);
            double fit_ndf = gaus_fit->GetNDF();
            double fit_chi2  = 1000.0;
            if (fit_ndf > 0){
                fit_chi2  = gaus_fit->GetChisquare() / fit_ndf;
                // normalize the chi2 to 4000 events
                fit_chi2 = fit_chi2 * 4000 / entry_max;
            }

            bool single_peak_fit_valid = (fit_chi2 >= 1.0 && fit_mean >= gaussian_fit_valid_min_mean);
            if (!single_peak_fit_valid && channel_pedestal_peak_counts[_hist_index] == -1){
                channel_pedestal_peak_counts[_hist_index] = 0;
            } else if (fit_chi2 < good_fit_chi2ndf && single_peak_fit_valid){
                gaus_fit->Draw("SAME");
                _legend->AddEntry(gaus_fit, ("#mu:" + std::to_string(fit_mean).substr(0, 5) + " #sigma:" + std::to_string(fit_sigma).substr(0, 5) + " #chi^{2}:" + std::to_string(fit_chi2).substr(0, 5)).c_str(), "l");
                if (channel_pedestal_peak_counts[_hist_index] <= 0) {
                    channel_pedestal_peak_counts[_hist_index] = 1;
                }
                single_gaussian_fit_mean.push_back(fit_mean);
                single_gaussian_fit_sigma.push_back(fit_sigma);
                single_gaussian_fit_chi2.push_back(fit_chi2);
                TLine* _line = new TLine(_initial_mean, 0, _initial_mean, max_count * 1.1);
                _line->SetLineColorAlpha(kPink, 0.5);
                _line->SetLineStyle(2);
                _line->SetLineWidth(1);
                _line->Draw();
            } else {
                // try to fit with dual gaussian
                TF1 *dual_gaus_fit = new TF1("dual_gaus_fit", "gaus(0) + gaus(3)", channel_pedestal_adc_hist_min, channel_pedestal_adc_hist_max);
                // set the initial parameters
                double _initial_mean1 = 10;
                double _initial_sigma1 = _initial_sigma;
                double _initial_amplitude1 = _initial_amplitude;
                double _initial_mean2 = 130;
                double _initial_sigma2 = _initial_sigma;
                double _initial_amplitude2 = _initial_amplitude / 3;
                dual_gaus_fit->SetParameter(0, _initial_amplitude1);
                dual_gaus_fit->SetParameter(1, _initial_mean1);
                dual_gaus_fit->SetParameter(2, _initial_sigma1);
                dual_gaus_fit->SetParameter(3, _initial_amplitude2);
                dual_gaus_fit->SetParameter(4, _initial_mean2);
                dual_gaus_fit->SetParameter(5, _initial_sigma2);

                dual_gaus_fit->SetParLimits(1, channel_pedestal_adc_hist_min, 180);
                dual_gaus_fit->SetParLimits(2, 0.1, gaussian_fit_max_sigma);
                dual_gaus_fit->SetParLimits(4, channel_pedestal_adc_hist_min, 200);
                dual_gaus_fit->SetParLimits(5, 0.1, gaussian_fit_max_sigma);
                
                _channel_pede_hist->Fit("dual_gaus_fit", "RQN+");
                dual_gaus_fit->SetLineColorAlpha(color_list[_run_index % color_list.size()], 0.5);

                double fit_mean1  = dual_gaus_fit->GetParameter(1);
                double fit_sigma1 = dual_gaus_fit->GetParameter(2);
                
                double fit_mean2  = dual_gaus_fit->GetParameter(4);
                double fit_sigma2 = dual_gaus_fit->GetParameter(5);

                double fit_ndf = dual_gaus_fit->GetNDF();
                double fit_chi2_dual = 1000.0;
                if (fit_ndf > 0){
                    fit_chi2_dual = dual_gaus_fit->GetChisquare() / fit_ndf;
                    // normalize the chi2 to 4000 events
                    fit_chi2_dual = fit_chi2_dual * 4000 / entry_max;
                }

                if (fit_chi2_dual < fit_chi2 - dual_peak_fit_chi2ndf_threshold && abs(fit_mean1 - fit_mean2) > 20 && fit_chi2_dual < okay_fit_chi2ndf){
                    dual_gaus_fit->Draw("SAME");
                    _legend->AddEntry(dual_gaus_fit, ("#mu1:" + std::to_string(fit_mean1).substr(0, 5) + " #mu2:" + std::to_string(fit_mean2).substr(0, 5) + " #chi2:" + std::to_string(fit_chi2_dual).substr(0, 5)).c_str(), "l");
                    if (channel_pedestal_peak_counts[_hist_index] != 2) {
                        channel_pedestal_peak_counts[_hist_index] = 2;
                    }
                    dual_gaussian_fit_mean1.push_back(fit_mean1);
                    dual_gaussian_fit_sigma1.push_back(fit_sigma1);
                    dual_gaussian_fit_mean2.push_back(fit_mean2);
                    dual_gaussian_fit_sigma2.push_back(fit_sigma2);
                    dual_gaussian_fit_chi2.push_back(fit_chi2_dual);

                    TLine* _line1 = new TLine(_initial_mean1, 0, _initial_mean1, max_count * 1.1);
                    _line1->SetLineColorAlpha(kOrange, 0.5);
                    _line1->SetLineStyle(2);
                    _line1->SetLineWidth(1);
                    _line1->Draw();

                    TLine* _line2 = new TLine(_initial_mean2, 0, _initial_mean2, max_count * 1.1);
                    _line2->SetLineColorAlpha(kBlue, 0.5);
                    _line2->SetLineStyle(2);
                    _line2->SetLineWidth(1);
                    _line2->Draw();
                } else {
                    if (fit_chi2 < okay_fit_chi2ndf && single_peak_fit_valid){
                        gaus_fit->Draw("SAME");
                        if (channel_pedestal_peak_counts[_hist_index] <= 0) {
                            channel_pedestal_peak_counts[_hist_index] = 1;
                        }
                        _legend->AddEntry(gaus_fit, ("#mu:" + std::to_string(fit_mean).substr(0, 5) + " #sigma:" + std::to_string(fit_sigma).substr(0, 5) + " #chi^{2}:" + std::to_string(fit_chi2).substr(0, 5)).c_str(), "l");
                        single_gaussian_fit_mean.push_back(fit_mean);
                        single_gaussian_fit_sigma.push_back(fit_sigma);
                        single_gaussian_fit_chi2.push_back(fit_chi2);

                        TLine* _line = new TLine(_initial_mean, 0, _initial_mean, max_count * 1.1);
                        _line->SetLineColorAlpha(kPink, 0.5);
                        _line->SetLineStyle(2);
                        _line->SetLineWidth(1);
                        _line->Draw();
                    }
                }
            }
        }
        _legend->Draw();

        auto _Latex = new TLatex();
        _Latex->SetNDC();
        _Latex->SetTextSize(0.04);
        _Latex->SetTextFont(52);
        if (channel_pedestal_peak_counts[_hist_index] == 1){
            _Latex->SetTextColor(kBlue+4);
            _Latex->DrawLatex(0.05, 0.95, "Single Peak");

            double weighted_sum = 0;
            double weight_sum = 0;

            for (int _run_index = 0; _run_index < single_gaussian_fit_mean.size(); _run_index++) {
                double _weight = 1.0 / (single_gaussian_fit_sigma[_run_index] * single_gaussian_fit_sigma[_run_index]);
                weighted_sum += single_gaussian_fit_mean[_run_index] * _weight;
                weight_sum += _weight;
            }

            double weighted_mean = weighted_sum / weight_sum;
            double weighted_error = 1.0 / sqrt(weight_sum);

            channel_single_peak_uni_index_list.push_back(_uni_channel_number);
            channel_single_peak_mean_list.push_back(weighted_mean);
            channel_single_peak_error_list.push_back(weighted_error);

            _Latex->DrawLatex(0.05, 0.90, ("Weighted Mean: " + std::to_string(weighted_mean).substr(0, 5) + " #pm " + std::to_string(weighted_error).substr(0, 5)).c_str());
            _Latex->DrawLatex(0.05, 0.85, ("Channel: " + std::to_string(_uni_channel_number)).c_str());

        } else if (channel_pedestal_peak_counts[_hist_index] == 0){
            _Latex->SetTextColor(kRed);
            _Latex->DrawLatex(0.05, 0.95, "No Pedestal");
        } else if (channel_pedestal_peak_counts[_hist_index] == -1){
            _Latex->SetTextColor(kOrange+2);
            _Latex->DrawLatex(0.05, 0.95, "Unprocessed");
        } else if (channel_pedestal_peak_counts[_hist_index] == 2){
            _Latex->SetTextColor(kTeal+3);
            _Latex->DrawLatex(0.05, 0.95, "Dual Peaks");

            double weighted_sum1 = 0;
            double weight_sum1 = 0;
            
            for (int _run_index = 0; _run_index < dual_gaussian_fit_mean1.size(); _run_index++) {
                double _weight = 1.0 / (dual_gaussian_fit_sigma1[_run_index] * dual_gaussian_fit_sigma1[_run_index]);
                weighted_sum1 += dual_gaussian_fit_mean1[_run_index] * _weight;
                weight_sum1 += _weight;
            }

            double weighted_mean1 = weighted_sum1 / weight_sum1;
            double weighted_error1 = 1.0 / sqrt(weight_sum1);

            double weighted_sum2 = 0;
            double weight_sum2 = 0;

            for (int _run_index = 0; _run_index < dual_gaussian_fit_mean2.size(); _run_index++) {
                double _weight = 1.0 / (dual_gaussian_fit_sigma2[_run_index] * dual_gaussian_fit_sigma2[_run_index]);
                weighted_sum2 += dual_gaussian_fit_mean2[_run_index] * _weight;
                weight_sum2 += _weight;
            }

            double weighted_mean2 = weighted_sum2 / weight_sum2;
            double weighted_error2 = 1.0 / sqrt(weight_sum2);

            channel_dual_peak_mean1_list.push_back(weighted_mean1);
            channel_dual_peak_error1_list.push_back(weighted_error1);
            channel_dual_peak_mean2_list.push_back(weighted_mean2);
            channel_dual_peak_error2_list.push_back(weighted_error2);
            channel_dual_peak_uni_index_list.push_back(_uni_channel_number);

            _Latex->DrawLatex(0.05, 0.90, ("Weighted Mean1: " + std::to_string(weighted_mean1).substr(0, 5) + " #pm " + std::to_string(weighted_error1).substr(0, 5)).c_str());
            _Latex->DrawLatex(0.05, 0.85, ("Weighted Mean2: " + std::to_string(weighted_mean2).substr(0, 5) + " #pm " + std::to_string(weighted_error2).substr(0, 5)).c_str());
            _Latex->DrawLatex(0.05, 0.80, ("Channel: " + std::to_string(_uni_channel_number)).c_str());
        }

        _channel_adc_canvas->Update();
        _channel_adc_canvas->Write();
    }

    std::vector <std::string> channel_wise_adc_legend_labels;
    for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
        auto _run_number    = config_run_numbers[_run_index];
        auto _beam_energy   = config_beam_energies[_run_index];
        auto _beam_particle = config_beam_particles[_run_index];
        channel_wise_adc_legend_labels.push_back("Run " + std::to_string(_run_number) + " " + std::to_string(int(_beam_energy)) + " GeV " + _beam_particle);
    } 
    global_painter->draw_global_channel_hists1D_run_group(channel_wise_adc_hists, channel_wise_hists_unified_channel_map, "ChannelWiseADC", "Channel Wise ADC Channel # ", color_list, channel_wise_adc_legend_labels);
    auto global_channel_adc_canvas = global_painter->get_canvas();
    output_root->cd();
    global_channel_adc_canvas->Write();
    global_channel_adc_canvas->SetTitle("Channel Wise ADC");
    global_channel_adc_canvas->SaveAs((pdf_file_name + "(").c_str());

    for (int _module_index = 1; _module_index < 10; _module_index++) {
        auto _module_canvas = global_painter->draw_module_channel_canvas(channel_wise_adc_canvas, channel_wise_hists_unified_channel_map, "ModuleCanvas" + std::to_string(_module_index), "Module Canvas " + std::to_string(_module_index), _module_index);
        _module_canvas->Write();

        _module_canvas->SaveAs(pdf_file_name.c_str());
        _module_canvas->Close();
    }

    // * --- Draw the distribution of the pedestal values --------------------------------
    // * --------------------------------------------------------------------------------
    TCanvas *channel_pede_distribution_canvas = new TCanvas("ChannelPedeDistribution", "Channel Pedestal Distribution", 1000, 600);
    channel_pede_distribution_canvas->cd();
    auto single_peak_channel_distribution_graph_error = new TGraphErrors(channel_single_peak_uni_index_list.size());
    for (int _index = 0; _index < channel_single_peak_uni_index_list.size(); _index++) {
        auto _uni_channel_number = channel_single_peak_uni_index_list[_index];
        auto _mean = channel_single_peak_mean_list[_index];
        auto _error = channel_single_peak_error_list[_index];
        single_peak_channel_distribution_graph_error->SetPoint(_index, _uni_channel_number, _mean);
        single_peak_channel_distribution_graph_error->SetPointError(_index, 0.0, _error);
    }
    auto dual_peak_channel_distribution_graph_error1 = new TGraphErrors(channel_dual_peak_uni_index_list.size());
    auto dual_peak_channel_distribution_graph_error2 = new TGraphErrors(channel_dual_peak_uni_index_list.size());
    for (int _index = 0; _index < channel_dual_peak_uni_index_list.size(); _index++) {
        auto _uni_channel_number = channel_dual_peak_uni_index_list[_index];
        auto _mean1 = channel_dual_peak_mean1_list[_index];
        auto _error1 = channel_dual_peak_error1_list[_index];
        auto _mean2 = channel_dual_peak_mean2_list[_index];
        auto _error2 = channel_dual_peak_error2_list[_index];
        if (_uni_channel_number > 250){
            LOG(DEBUG) << "Channel " << _uni_channel_number << " Mean1: " << _mean1 << " Error1: " << _error1 << " Mean2: " << _mean2 << " Error2: " << _error2;
        }
        dual_peak_channel_distribution_graph_error1->SetPoint(_index, _uni_channel_number, _mean1);
        dual_peak_channel_distribution_graph_error1->SetPointError(_index, 0.0, _error1);
        dual_peak_channel_distribution_graph_error2->SetPoint(_index, _uni_channel_number, _mean2);
        dual_peak_channel_distribution_graph_error2->SetPointError(_index, 0.0, _error2);
    }

    auto _legend = new TLegend(0.7, 0.7, 0.89, 0.89);
    _legend->SetFillColor(0);
    _legend->SetLineColor(0);
    _legend->SetBorderSize(0);
    _legend->SetLineColorAlpha(0, 0);
    _legend->SetTextSize(0.03);

    single_peak_channel_distribution_graph_error->SetMarkerStyle(20);
    single_peak_channel_distribution_graph_error->SetMarkerSize(0.2);
    single_peak_channel_distribution_graph_error->SetMarkerColor(kCyan+3);
    single_peak_channel_distribution_graph_error->SetLineColor(kCyan+3);
    single_peak_channel_distribution_graph_error->SetLineWidth(1);
    single_peak_channel_distribution_graph_error->SetTitle("");
    single_peak_channel_distribution_graph_error->GetXaxis()->SetTitle("Channel Number");
    single_peak_channel_distribution_graph_error->GetXaxis()->SetRangeUser(0, fpga_count * FPGA_CHANNEL_NUMBER_VALID + 1);
    single_peak_channel_distribution_graph_error->GetYaxis()->SetTitle("Pedestal Value [ADC]");
    single_peak_channel_distribution_graph_error->GetYaxis()->SetRangeUser(0, 256);
    single_peak_channel_distribution_graph_error->Draw("APE");
    _legend->AddEntry(single_peak_channel_distribution_graph_error, "Single Peak", "ple");
    
    dual_peak_channel_distribution_graph_error1->SetMarkerStyle(20);
    dual_peak_channel_distribution_graph_error1->SetMarkerSize(0.2);
    dual_peak_channel_distribution_graph_error1->SetMarkerColor(kOrange+2);
    dual_peak_channel_distribution_graph_error1->SetLineColor(kOrange+2);
    dual_peak_channel_distribution_graph_error1->SetLineWidth(1);
    dual_peak_channel_distribution_graph_error1->SetTitle("");
    dual_peak_channel_distribution_graph_error1->Draw("PE");
    _legend->AddEntry(dual_peak_channel_distribution_graph_error1, "Dual Peak #1", "ple");

    dual_peak_channel_distribution_graph_error2->SetMarkerStyle(20);
    dual_peak_channel_distribution_graph_error2->SetMarkerSize(0.2);
    dual_peak_channel_distribution_graph_error2->SetMarkerColor(kPink-5);
    dual_peak_channel_distribution_graph_error2->SetLineColor(kPink-5);
    dual_peak_channel_distribution_graph_error2->SetLineWidth(1);
    dual_peak_channel_distribution_graph_error2->SetTitle("");
    dual_peak_channel_distribution_graph_error2->Draw("PE");
    _legend->AddEntry(dual_peak_channel_distribution_graph_error2, "Dual Peak #2", "ple");

    _legend->Draw();

    auto _Latex = new TLatex();
    const double _text_line_height = 0.045;
    const double _text_line_start = 0.85;
    const double _text_line_left = 0.13;
    _Latex->SetNDC();
    _Latex->SetTextSize(0.04);
    _Latex->SetTextFont(62);
    _Latex->SetTextColor(kBlack);
    _Latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    _Latex->SetTextSize(0.03);
    _Latex->SetTextFont(42);
    
    for (int _line_index = 1; _line_index < config_plot_info.size(); _line_index++) {
        _Latex->DrawLatex(_text_line_left, _text_line_start - _text_line_height * _line_index, (config_plot_info[_line_index].c_str()));
    }
    if (enable_working_in_progress){
        _Latex->SetTextSize(0.04);
        _Latex->SetTextColor(kGray+3);
        _Latex->SetTextFont(52);
        _Latex->DrawLatex(_text_line_left, _text_line_start - _text_line_height * (config_plot_info.size()), "Work in Progress");
    }

    output_root->cd();
    channel_pede_distribution_canvas->Write();
    channel_pede_distribution_canvas->SaveAs((pdf_file_name + ")").c_str());

    // * --- Save a tree of channels to the output file ---------------------------------
    // * --------------------------------------------------------------------------------
    TTree *channel_pede_tree = new TTree("ChannelPedestal", "Channel Pedestal");
    std::vector <UInt_t> channel_pede_channel_number;
    std::vector <UInt_t> channel_pede_peak_counts;
    std::vector <Double_t> channel_pede_mean1;
    std::vector <Double_t> channel_pede_error1;
    std::vector <Double_t> channel_pede_mean2;
    std::vector <Double_t> channel_pede_error2;

    for (int _valid_channel_index = 0; _valid_channel_index < FPGA_CHANNEL_NUMBER_VALID * fpga_count; _valid_channel_index++) {
        channel_pede_channel_number.push_back(_valid_channel_index);
        auto _single_peak_index = std::find(channel_single_peak_uni_index_list.begin(), channel_single_peak_uni_index_list.end(), _valid_channel_index);
        if (_single_peak_index != channel_single_peak_uni_index_list.end()){
            auto _single_peak_index_int = std::distance(channel_single_peak_uni_index_list.begin(), _single_peak_index);
            channel_pede_peak_counts.push_back(1);
            channel_pede_mean1.push_back(channel_single_peak_mean_list[_single_peak_index_int]);
            channel_pede_error1.push_back(channel_single_peak_error_list[_single_peak_index_int]);
            channel_pede_mean2.push_back(0);
            channel_pede_error2.push_back(0);
        } else {
            auto _dual_peak_index = std::find(channel_dual_peak_uni_index_list.begin(), channel_dual_peak_uni_index_list.end(), _valid_channel_index);
            if (_dual_peak_index != channel_dual_peak_uni_index_list.end()){
                auto _dual_peak_index_int = std::distance(channel_dual_peak_uni_index_list.begin(), _dual_peak_index);
                channel_pede_peak_counts.push_back(2);
                channel_pede_mean1.push_back(channel_dual_peak_mean1_list[_dual_peak_index_int]);
                channel_pede_error1.push_back(channel_dual_peak_error1_list[_dual_peak_index_int]);
                channel_pede_mean2.push_back(channel_dual_peak_mean2_list[_dual_peak_index_int]);
                channel_pede_error2.push_back(channel_dual_peak_error2_list[_dual_peak_index_int]);
            } else {
                channel_pede_peak_counts.push_back(0);
                channel_pede_mean1.push_back(0);
                channel_pede_error1.push_back(0);
                channel_pede_mean2.push_back(0);
                channel_pede_error2.push_back(0);
            }
        }
    }

    UInt_t _channel_pede_channel_number;
    UInt_t _channel_pede_peak_counts;
    Double_t _channel_pede_mean1;
    Double_t _channel_pede_error1;
    Double_t _channel_pede_mean2;
    Double_t _channel_pede_error2;

    channel_pede_tree->Branch("ChannelNumber", &(_channel_pede_channel_number));
    channel_pede_tree->Branch("PeakCounts", &(_channel_pede_peak_counts));
    channel_pede_tree->Branch("Mean1", &(_channel_pede_mean1));
    channel_pede_tree->Branch("Error1", &(_channel_pede_error1));
    channel_pede_tree->Branch("Mean2", &(_channel_pede_mean2));
    channel_pede_tree->Branch("Error2", &(_channel_pede_error2));


    for (int _channel_index = 0; _channel_index < channel_pede_channel_number.size(); _channel_index++) {
        _channel_pede_channel_number = channel_pede_channel_number[_channel_index];
        _channel_pede_peak_counts = channel_pede_peak_counts[_channel_index];
        _channel_pede_mean1 = channel_pede_mean1[_channel_index];
        _channel_pede_error1 = channel_pede_error1[_channel_index];
        _channel_pede_mean2 = channel_pede_mean2[_channel_index];
        _channel_pede_error2 = channel_pede_error2[_channel_index];
        
        channel_pede_tree->Fill();
    }

    output_root->cd();
    channel_pede_tree->Write();

    output_root->Close();
    LOG(INFO) << "Output file " << opts.output_file << " has been saved.";
    LOG(INFO) << "PDF file " << pdf_file_name << " has been saved.";

    for (auto _channel_adc_canvas : channel_wise_adc_canvas) {
        if (_channel_adc_canvas != nullptr) {
            _channel_adc_canvas->Close();
        }
    }

    if (global_painter != nullptr) {
        delete global_painter;
    }
    return 0;
}