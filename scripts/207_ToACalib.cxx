#include "H2GCROC_Common.hxx"
#include "csv.hpp"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_root(argc, argv, "1.0");
    gROOT->SetBatch(kTRUE);

    bool enable_focal_mapping = opts.focal;
    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 25.0 / 16.0; // unit: ns

    const int code_hists_x_min  = -100;
    const int code_hists_x_max  = 1124;
    const int code_hists_x_bins = (code_hists_x_max - code_hists_x_min);

    const int code_grey_counter_hists_x_min  = 0;
    const int code_grey_counter_hists_x_max  = 4;
    const int code_grey_counter_hists_x_bins = (code_grey_counter_hists_x_max - code_grey_counter_hists_x_min);

    const int code_corase_tdc_hists_x_min  = 0;
    const int code_corase_tdc_hists_x_max  = 32;
    const int code_corase_tdc_hists_x_bins = (code_corase_tdc_hists_x_max - code_corase_tdc_hists_x_min);

    const int code_fine_tdc_hists_x_min  = 0;
    const int code_fine_tdc_hists_x_max  = 8;
    const int code_fine_tdc_hists_x_bins = (code_fine_tdc_hists_x_max - code_fine_tdc_hists_x_min);

    const int module_to_check = 5;

    // * --- Global channel painter -----------------------------------------------------
    // * --------------------------------------------------------------------------------
    GlobalChannelPainter *global_painter = nullptr;
    if (enable_focal_mapping){
        global_painter = new GlobalChannelPainter("data/SPS_2024/config/focalh_mapping.json", "data/SPS_2024/config/h2gcroc_mapping.json");
    } else {
        global_painter = new GlobalChannelPainter("data/DESY_2025/config/EEEMCal_Mapping_DESY_2025.json");
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

    std::vector <TH1I*> channel_adc_hist_list;
    TDirectory *channel_adc_folder = output_root->mkdir("ChannelADC");
    std::vector <TH1I*> channel_toa_hist_list;
    TDirectory *channel_toa_folder = output_root->mkdir("ChannelTOA");
    std::vector <TH1I*> channel_tot_hist_list;
    TDirectory *channel_tot_folder = output_root->mkdir("ChannelTOT");

    std::vector <TH1I*> channel_toa_grey_counter_hist_list;
    TDirectory *channel_toa_grey_counter_folder = output_root->mkdir("ChannelTOAGreyCounter");

    std::vector <TH1I*> channel_toa_corase_tdc_hist_list;
    TDirectory *channel_toa_corase_tdc_folder = output_root->mkdir("ChannelTOACoraseTDC");

    std::vector <TH1I*> channel_toa_fine_tdc_hist_list;
    TDirectory *channel_toa_fine_tdc_folder = output_root->mkdir("ChannelTOAFineTDC");

    std::vector <Long_t> channel_valid_toa_event_counters;
    std::vector <Long_t> channel_valid_tot_event_counters;
    std::vector <Long_t> channel_double_toa_event_counters;
    std::vector <Long_t> channel_double_tot_event_counters;

    std::vector <Long_t> channel_illegal_toa_event_counters;

    int channel_array_index = 0;
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

            // * --- ADC histogram ---------------------------------------------------------
            auto *_adc_hist = new TH1I(("channel_" + std::to_string(_unified_valid_channel_number) + "_adc").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " ADC").c_str(), code_hists_x_bins, code_hists_x_min, code_hists_x_max);
            channel_adc_hist_list.push_back(_adc_hist);
            _adc_hist->SetDirectory(channel_adc_folder);

            // * --- TOA histogram ---------------------------------------------------------
            auto *_toa_hist = new TH1I(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA").c_str(), code_hists_x_bins, code_hists_x_min, code_hists_x_max);
            channel_toa_hist_list.push_back(_toa_hist);
            _toa_hist->SetDirectory(channel_toa_folder);

            // * --- TOT histogram ---------------------------------------------------------
            auto *_tot_hist = new TH1I(("channel_" + std::to_string(_unified_valid_channel_number) + "_tot").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOT").c_str(), code_hists_x_bins, code_hists_x_min, code_hists_x_max);
            channel_tot_hist_list.push_back(_tot_hist);
            _tot_hist->SetDirectory(channel_tot_folder);

            // * --- TOA Grey Counter histogram ---------------------------------------------
            auto *_toa_grey_counter_hist = new TH1I(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa_grey_counter").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA Grey Counter").c_str(), code_grey_counter_hists_x_bins, code_grey_counter_hists_x_min, code_grey_counter_hists_x_max);
            channel_toa_grey_counter_hist_list.push_back(_toa_grey_counter_hist);
            _toa_grey_counter_hist->SetDirectory(channel_toa_grey_counter_folder);

            // * --- TOA Corase TDC histogram ----------------------------------------------
            auto *_toa_corase_tdc_hist = new TH1I(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa_corase_tdc").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA Corase TDC").c_str(), code_corase_tdc_hists_x_bins, code_corase_tdc_hists_x_min, code_corase_tdc_hists_x_max);
            channel_toa_corase_tdc_hist_list.push_back(_toa_corase_tdc_hist);
            _toa_corase_tdc_hist->SetDirectory(channel_toa_corase_tdc_folder);

            // * --- TOA Fine TDC histogram ------------------------------------------------
            auto *_toa_fine_tdc_hist = new TH1I(("channel_" + std::to_string(_unified_valid_channel_number) + "_toa_fine_tdc").c_str(), ("Channel " + std::to_string(_unified_valid_channel_number) + " TOA Fine TDC").c_str(), code_fine_tdc_hists_x_bins, code_fine_tdc_hists_x_min, code_fine_tdc_hists_x_max);
            channel_toa_fine_tdc_hist_list.push_back(_toa_fine_tdc_hist);
            _toa_fine_tdc_hist->SetDirectory(channel_toa_fine_tdc_folder);

            // * --- Event counters --------------------------------------------------------
            channel_valid_toa_event_counters.push_back(0);
            channel_valid_tot_event_counters.push_back(0);
            channel_double_toa_event_counters.push_back(0);
            channel_double_tot_event_counters.push_back(0);
            channel_illegal_toa_event_counters.push_back(0);

            channel_array_index++;
        } // end of channel loop
    } // end of fpga loop

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

            std::vector <bool> _hamming_code_pass_list; // for each machine gun sample
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
                auto _channel_array_index = channel_unified_channel_number_to_index_map[_unified_valid_channel_number];

                int _val0_max = -1;
                int _val0_max_index = -1;
                int _val1_max = -1;
                int _val1_max_index = -1;
                int _val2_max = -1;
                int _val2_max_index = -1;

                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    if (_hamming_code_pass_list[_sample_index] == false) {
                        continue;
                    }
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
                            channel_valid_tot_event_counters[_channel_array_index]++;
                        } else {
                            channel_double_tot_event_counters[_channel_array_index]++;
                        }
                    }
                    if (_val2 > 0) {
                        if (_val2_max == -1) {
                            _val2_max = _val2;
                            _val2_max_index = _sample_index;
                            channel_valid_toa_event_counters[_channel_array_index]++;
                        } else {
                            channel_double_toa_event_counters[_channel_array_index]++;
                        }
                    }
                } // end of sample loop

                if (_val0_max == -1) {
                    _val0_max = 0;
                }
                if (_val1_max == -1) {
                    _val1_max = 0;
                }
                if (_val2_max == -1) {
                    _val2_max = 0;
                }

                if (_val2_max % 250 < 16){
                    channel_illegal_toa_event_counters[_channel_array_index]++;
                }

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

                if (_channel_valid == -1){
                    continue;
                }

                auto _val0_max = _channel_val0_max_list[_channel_valid];
                auto _val1_max = _channel_val1_max_list[_channel_valid];
                auto _val2_max = _channel_val2_max_list[_channel_valid];

                int _toa_grey_counter = _val2_max / 256;
                int _toa_corase_tdc = _val2_max % 256 / 8;
                int _toa_fine_tdc = _val2_max % 8;

                if (_val0_max > 0){
                    channel_adc_hist_list[channel_unified_channel_number_to_index_map[_unified_valid_channel_number]]->Fill(_val0_max);
                }
                if (_val2_max > 0 && _val0_max > 400){
                    channel_toa_hist_list[channel_unified_channel_number_to_index_map[_unified_valid_channel_number]]->Fill(_val2_max);
                    channel_toa_grey_counter_hist_list[channel_unified_channel_number_to_index_map[_unified_valid_channel_number]]->Fill(_toa_grey_counter);
                    channel_toa_corase_tdc_hist_list[channel_unified_channel_number_to_index_map[_unified_valid_channel_number]]->Fill(_toa_corase_tdc);
                    channel_toa_fine_tdc_hist_list[channel_unified_channel_number_to_index_map[_unified_valid_channel_number]]->Fill(_toa_fine_tdc);
                }
                if (_val1_max > 0){
                    channel_tot_hist_list[channel_unified_channel_number_to_index_map[_unified_valid_channel_number]]->Fill(_val1_max);
                }
            } // end of 2nd channel loop
        } // end of fpga loop
    } // end of event loop

    input_root->Close();

    // * --- Write the histograms -------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of(".")) + ".pdf";

    LOG(INFO) << "Saving ADC code histograms...";
    channel_adc_folder->cd();
    std::vector <TCanvas*> channel_adc_canvas_list;
    for (int _channel_index = 0; _channel_index < channel_adc_hist_list.size(); _channel_index++) {
        auto _canvas = new TCanvas(("channel_" + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + "_adc_canvas").c_str(), ("Channel " + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + " ADC Canvas").c_str(), 800, 600);
        channel_adc_hist_list[_channel_index]->Draw();
        _canvas->Write();
        channel_adc_canvas_list.push_back(_canvas);
    }
    auto _painter_canvas = global_painter->draw_module_channel_canvas(channel_adc_canvas_list, channel_unified_channel_number_to_index_map, "ChannelADC", "Channel ADC", module_to_check);
    if (_painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }
    output_root->cd();
    _painter_canvas->Write();

    LOG(INFO) << "Saving TOA code histograms...";
    channel_toa_folder->cd();
    std::vector <TCanvas*> channel_toa_canvas_list;
    for (int _channel_index = 0; _channel_index < channel_toa_hist_list.size(); _channel_index++) {
        auto _canvas = new TCanvas(("channel_" + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + "_toa_canvas").c_str(), ("Channel " + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + " TOA Canvas").c_str(), 800, 600);
        channel_toa_hist_list[_channel_index]->SetStats(0);
        channel_toa_hist_list[_channel_index]->SetTitle("");
        channel_toa_hist_list[_channel_index]->Draw();
        
        auto _canvas_latex = new TLatex();
        _canvas_latex->SetNDC();
        _canvas_latex->SetTextSize(0.05);
        _canvas_latex->DrawLatex(0.1, 0.90, Form("Channel %d", channel_index_to_unified_channel_number_map[_channel_index]));

        auto _channel_toa_valid_event_counter = channel_valid_toa_event_counters[_channel_index];
        auto _channel_toa_double_event_counter = channel_double_toa_event_counters[_channel_index];
        auto _channel_toa_illegal_event_counter = channel_illegal_toa_event_counters[_channel_index];

        double _channel_toa_valid_event_ratio = double(_channel_toa_valid_event_counter) / double(entry_max) * 100.0;
        double _channel_toa_double_event_ratio = double(_channel_toa_double_event_counter) / double(entry_max) * 100.0;
        double _channel_toa_illegal_event_ratio = double(_channel_toa_illegal_event_counter) / double(entry_max) * 100.0;

        _canvas_latex->DrawLatex(0.1, 0.85, Form("Valid TOA Events: %ld (%.2f %%)", _channel_toa_valid_event_counter, _channel_toa_valid_event_ratio));
        _canvas_latex->DrawLatex(0.1, 0.80, Form("Double TOA Events: %ld (%.2f %%)", _channel_toa_double_event_counter, _channel_toa_double_event_ratio));
        _canvas_latex->DrawLatex(0.1, 0.75, Form("Illegal TOA Events: %ld (%.2f %%)", _channel_toa_illegal_event_counter, _channel_toa_illegal_event_ratio));
        _canvas_latex->Delete();

        _canvas->Write();
        channel_toa_canvas_list.push_back(_canvas);
    }
    _painter_canvas = global_painter->draw_module_channel_canvas(channel_toa_canvas_list, channel_unified_channel_number_to_index_map, "ChannelTOA", "Channel TOA", module_to_check);
    if (_painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }
    output_root->cd();
    _painter_canvas->Write();
    _painter_canvas->Print((pdf_file_name + "(").c_str());

    LOG(INFO) << "Saving TOT code histograms...";
    channel_tot_folder->cd();
    std::vector <TCanvas*> channel_tot_canvas_list;
    for (int _channel_index = 0; _channel_index < channel_tot_hist_list.size(); _channel_index++) {
        auto _canvas = new TCanvas(("channel_" + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + "_tot_canvas").c_str(), ("Channel " + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + " TOT Canvas").c_str(), 800, 600);
        channel_tot_hist_list[_channel_index]->SetStats(0);
        channel_tot_hist_list[_channel_index]->SetTitle("");
        channel_tot_hist_list[_channel_index]->Draw();

        auto _canvas_latex = new TLatex();
        _canvas_latex->SetNDC();
        _canvas_latex->SetTextSize(0.05);
        _canvas_latex->DrawLatex(0.1, 0.90, Form("Channel %d", channel_index_to_unified_channel_number_map[_channel_index]));
        auto _channel_tot_valid_event_counter = channel_valid_tot_event_counters[_channel_index];
        auto _channel_tot_double_event_counter = channel_double_tot_event_counters[_channel_index];
        double _channel_tot_valid_event_ratio = double(_channel_tot_valid_event_counter) / double(entry_max) * 100.0;
        double _channel_tot_double_event_ratio = double(_channel_tot_double_event_counter) / double(entry_max) * 100.0;
        _canvas_latex->DrawLatex(0.1, 0.85, Form("Valid TOT Events: %ld (%.2f %%)", _channel_tot_valid_event_counter, _channel_tot_valid_event_ratio));
        _canvas_latex->DrawLatex(0.1, 0.80, Form("Double TOT Events: %ld (%.2f %%)", _channel_tot_double_event_counter, _channel_tot_double_event_ratio));
        _canvas_latex->Delete();

        _canvas->Write();
        channel_tot_canvas_list.push_back(_canvas);
    }
    _painter_canvas = global_painter->draw_module_channel_canvas(channel_tot_canvas_list, channel_unified_channel_number_to_index_map, "ChannelTOT", "Channel TOT", module_to_check);
    if (_painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }
    output_root->cd();
    _painter_canvas->Write();
    _painter_canvas->Print(pdf_file_name.c_str());

    LOG(INFO) << "Saving TOA Grey Counter histograms...";
    channel_toa_grey_counter_folder->cd();
    std::vector <TCanvas*> channel_toa_grey_counter_canvas_list;
    for (int _channel_index = 0; _channel_index < channel_toa_grey_counter_hist_list.size(); _channel_index++) {
        auto _canvas = new TCanvas(("channel_" + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + "_toa_grey_counter_canvas").c_str(), ("Channel " + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + " TOA Grey Counter Canvas").c_str(), 800, 600);
        channel_toa_grey_counter_hist_list[_channel_index]->SetStats(0);
        channel_toa_grey_counter_hist_list[_channel_index]->Draw();
        _canvas->Write();
        channel_toa_grey_counter_canvas_list.push_back(_canvas);
    }
    _painter_canvas = global_painter->draw_module_channel_canvas(channel_toa_grey_counter_canvas_list, channel_unified_channel_number_to_index_map, "ChannelTOAGreyCounter", "Channel TOA Grey Counter", module_to_check);
    if (_painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }
    output_root->cd();
    _painter_canvas->Write();
    _painter_canvas->Print(pdf_file_name.c_str());

    LOG(INFO) << "Saving TOA Corase TDC histograms...";
    channel_toa_corase_tdc_folder->cd();
    std::vector <std::vector <double>> channel_toa_corase_tdc_DNL_correction_ratio_list;
    std::vector <TCanvas*> channel_toa_corase_tdc_canvas_list;
    for (int _channel_index = 0; _channel_index < channel_toa_corase_tdc_hist_list.size(); _channel_index++) {
        auto _canvas = new TCanvas(("channel_" + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + "_toa_corase_tdc_canvas").c_str(), ("Channel " + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + " TOA Corase TDC Canvas").c_str(), 800, 600);
        channel_toa_corase_tdc_hist_list[_channel_index]->SetStats(0);
        channel_toa_corase_tdc_hist_list[_channel_index]->Draw();
        // get the ratio
        std::vector <double> _channel_toa_corase_tdc_DNL_correction_ratio;
        for (int _bin_index = 1; _bin_index <= channel_toa_corase_tdc_hist_list[_channel_index]->GetNbinsX(); _bin_index++) {
            auto _bin_content = channel_toa_corase_tdc_hist_list[_channel_index]->GetBinContent(_bin_index);
            _channel_toa_corase_tdc_DNL_correction_ratio.push_back(1.0/double(_bin_content));
        }
        // normalize the ratio
        double _sum = 0;
        for (int _bin_index = 0; _bin_index < _channel_toa_corase_tdc_DNL_correction_ratio.size(); _bin_index++) {
            _sum += _channel_toa_corase_tdc_DNL_correction_ratio[_bin_index];
        }
        for (int _bin_index = 0; _bin_index < _channel_toa_corase_tdc_DNL_correction_ratio.size(); _bin_index++) {
            _channel_toa_corase_tdc_DNL_correction_ratio[_bin_index] /= (double)_sum / double(channel_toa_corase_tdc_hist_list[_channel_index]->GetNbinsX());
        }
        channel_toa_corase_tdc_DNL_correction_ratio_list.push_back(_channel_toa_corase_tdc_DNL_correction_ratio);
        // label the ratio for each bin
        auto _canvas_latex = new TLatex();
        _canvas_latex->SetNDC();
        _canvas_latex->SetTextSize(0.03);
        _canvas_latex->SetTextAlign(12);
        _canvas_latex->SetTextColor(kPink + 1);
        for (int _bin_index = 1; _bin_index <= channel_toa_corase_tdc_hist_list[_channel_index]->GetNbinsX(); _bin_index++) {
            auto _bin_center = channel_toa_corase_tdc_hist_list[_channel_index]->GetBinCenter(_bin_index);
            auto _bin_content = channel_toa_corase_tdc_hist_list[_channel_index]->GetBinContent(_bin_index);
            auto _bin_ratio = channel_toa_corase_tdc_DNL_correction_ratio_list[_channel_index][_bin_index - 1];
            _canvas_latex->DrawLatex(_bin_center/channel_toa_corase_tdc_hist_list[_channel_index]->GetNbinsX() - 0.02, 0.5+0.1*(_bin_index%2==0), Form("%.2f", _bin_ratio));
        }
        _canvas_latex->Delete();
        _canvas->Write();
        channel_toa_corase_tdc_canvas_list.push_back(_canvas);
    }
    _painter_canvas = global_painter->draw_module_channel_canvas(channel_toa_corase_tdc_canvas_list, channel_unified_channel_number_to_index_map, "ChannelTOACoraseTDC", "Channel TOA Corase TDC", module_to_check);
    if (_painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }
    output_root->cd();
    _painter_canvas->Write();
    _painter_canvas->Print(pdf_file_name.c_str());

    LOG(INFO) << "Saving TOA Fine TDC histograms...";
    channel_toa_fine_tdc_folder->cd();
    std::vector <std::vector <double>> channel_toa_fine_tdc_DNL_correction_ratio_list;
    std::vector <TCanvas*> channel_toa_fine_tdc_canvas_list;
    for (int _channel_index = 0; _channel_index < channel_toa_fine_tdc_hist_list.size(); _channel_index++) {
        auto _canvas = new TCanvas(("channel_" + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + "_toa_fine_tdc_canvas").c_str(), ("Channel " + std::to_string(channel_unified_channel_number_to_index_map[_channel_index]) + " TOA Fine TDC Canvas").c_str(), 800, 600);
        channel_toa_fine_tdc_hist_list[_channel_index]->SetStats(0);
        channel_toa_fine_tdc_hist_list[_channel_index]->Draw();
        // get the ratio of each bin
        std::vector <double> _channel_toa_fine_tdc_DNL_correction_ratio;
        for (int _bin_index = 1; _bin_index <= channel_toa_fine_tdc_hist_list[_channel_index]->GetNbinsX(); _bin_index++) {
            auto _bin_content = channel_toa_fine_tdc_hist_list[_channel_index]->GetBinContent(_bin_index);
            _channel_toa_fine_tdc_DNL_correction_ratio.push_back(1.0/double(_bin_content));
        }
        // normalize the ratio
        double _sum = 0;
        for (int _bin_index = 0; _bin_index < _channel_toa_fine_tdc_DNL_correction_ratio.size(); _bin_index++) {
            _sum += _channel_toa_fine_tdc_DNL_correction_ratio[_bin_index];
        }
        for (int _bin_index = 0; _bin_index < _channel_toa_fine_tdc_DNL_correction_ratio.size(); _bin_index++) {
            _channel_toa_fine_tdc_DNL_correction_ratio[_bin_index] /= (double)_sum / double(channel_toa_fine_tdc_hist_list[_channel_index]->GetNbinsX());
        }
        channel_toa_fine_tdc_DNL_correction_ratio_list.push_back(_channel_toa_fine_tdc_DNL_correction_ratio);
        // label the ratio for each bin
        auto _canvas_latex = new TLatex();
        _canvas_latex->SetNDC();
        // write the ratio
        _canvas_latex->SetTextSize(0.03);
        _canvas_latex->SetTextAlign(12);
        _canvas_latex->SetTextColor(kPink + 1);
        for (int _bin_index = 1; _bin_index <= channel_toa_fine_tdc_hist_list[_channel_index]->GetNbinsX(); _bin_index++) {
            auto _bin_center = channel_toa_fine_tdc_hist_list[_channel_index]->GetBinCenter(_bin_index);
            auto _bin_content = channel_toa_fine_tdc_hist_list[_channel_index]->GetBinContent(_bin_index);
            auto _bin_ratio = channel_toa_fine_tdc_DNL_correction_ratio_list[_channel_index][_bin_index - 1];
            _canvas_latex->DrawLatex(_bin_center/channel_toa_fine_tdc_hist_list[_channel_index]->GetNbinsX() - 0.02, 0.5, Form("%.2f", _bin_ratio));
        }
        _canvas_latex->Delete();
        _canvas->Write();
        channel_toa_fine_tdc_canvas_list.push_back(_canvas);
    }
    _painter_canvas = global_painter->draw_module_channel_canvas(channel_toa_fine_tdc_canvas_list, channel_unified_channel_number_to_index_map, "ChannelTOAFineTDC", "Channel TOA Fine TDC", module_to_check);
    if (_painter_canvas == nullptr) {
        LOG(ERROR) << "Failed to draw global channel canvas!";
        return 1;
    }
    output_root->cd();
    _painter_canvas->Write();
    _painter_canvas->Print(pdf_file_name.c_str());

    TCanvas dummy_canvas;
    dummy_canvas.Print((pdf_file_name + ")").c_str());

    // save the DNL correction ratio to a json file for each channel
    json channel_toa_fine_tdc_DNL_correction_ratio_json;
    for (int _channel_index = 0; _channel_index < channel_toa_fine_tdc_DNL_correction_ratio_list.size(); _channel_index++) {
        auto _unified_valid_channel_number = channel_unified_channel_number_to_index_map[_channel_index];
        auto _channel_toa_fine_tdc_DNL_correction_ratio = channel_toa_fine_tdc_DNL_correction_ratio_list[_channel_index];
        json _channel_toa_fine_tdc_DNL_correction_ratio_json;
        for (int _bin_index = 0; _bin_index < _channel_toa_fine_tdc_DNL_correction_ratio.size(); _bin_index++) {
            _channel_toa_fine_tdc_DNL_correction_ratio_json.push_back(_channel_toa_fine_tdc_DNL_correction_ratio[_bin_index]);
        }
        channel_toa_fine_tdc_DNL_correction_ratio_json[std::to_string(_unified_valid_channel_number)] = _channel_toa_fine_tdc_DNL_correction_ratio_json;
    }
    std::string json_file_name = opts.output_file.substr(0, opts.output_file.find_last_of(".")) + "_fineDNL.json";
    std::ofstream json_file(json_file_name);
    if (!json_file.is_open()) {
        LOG(ERROR) << "Failed to open json file " << json_file_name;
        return 1;
    }
    json_file << std::setw(4) << channel_toa_fine_tdc_DNL_correction_ratio_json << std::endl;
    json_file.close();
    LOG(INFO) << "Saved json file " << json_file_name;

    // save the DNL correction ratio to a json file for each channel
    json channel_toa_corase_tdc_DNL_correction_ratio_json;
    for (int _channel_index = 0; _channel_index < channel_toa_corase_tdc_DNL_correction_ratio_list.size(); _channel_index++) {
        auto _unified_valid_channel_number = channel_unified_channel_number_to_index_map[_channel_index];
        auto _channel_toa_corase_tdc_DNL_correction_ratio = channel_toa_corase_tdc_DNL_correction_ratio_list[_channel_index];
        json _channel_toa_corase_tdc_DNL_correction_ratio_json;
        for (int _bin_index = 0; _bin_index < _channel_toa_corase_tdc_DNL_correction_ratio.size(); _bin_index++) {
            _channel_toa_corase_tdc_DNL_correction_ratio_json.push_back(_channel_toa_corase_tdc_DNL_correction_ratio[_bin_index]);
        }
        channel_toa_corase_tdc_DNL_correction_ratio_json[std::to_string(_unified_valid_channel_number)] = _channel_toa_corase_tdc_DNL_correction_ratio_json;
    }
    std::string json_file_name2 = opts.output_file.substr(0, opts.output_file.find_last_of(".")) + "_coraseDNL.json";
    std::ofstream json_file2(json_file_name2);
    if (!json_file2.is_open()) {
        LOG(ERROR) << "Failed to open json file " << json_file_name2;
        return 1;
    }
    json_file2 << std::setw(4) << channel_toa_corase_tdc_DNL_correction_ratio_json << std::endl;
    json_file2.close();
    LOG(INFO) << "Saved json file " << json_file_name2;

    output_root->Close();
    return 0;
}