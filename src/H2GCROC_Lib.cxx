#include "H2GCROC_Lib.hxx"

bool readRootMetaData(const char* _root_file_path, TFile*& _input_root, TTree*& _input_tree, int& _fpga_count, int& _machine_gun_samples, int& _max_entry, std::vector<UShort_t>& _legal_fpga_id_list) {
    _input_root = new TFile(_root_file_path, "READ");
    if (_input_root->IsZombie()) {
        LOG(ERROR) << "Failed to open input file " << _root_file_path;
        return false;
    }
    _input_tree = (TTree*) _input_root->Get("data_tree");
    if (_input_tree == nullptr) {
        LOG(ERROR) << "Failed to get data tree from input file " << _root_file_path;
        return false;
    }
    _max_entry = _input_tree->GetEntries();
    if (_max_entry == 0) {
        LOG(ERROR) << "No events in the input file!";
        return false;
    }

    TNamed *legal_fpga_id_list_tnamed = (TNamed*) _input_root->Get("Rootifier_legal_fpga_id_list");
    if (legal_fpga_id_list_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get legal fpga id list from input file " << _root_file_path;
        return false;
    }
    std::string legal_fpga_id_list_str = legal_fpga_id_list_tnamed->GetTitle();
    std::istringstream legal_fpga_id_list_stream(legal_fpga_id_list_str);
    UShort_t legal_fpga_id;
    while (legal_fpga_id_list_stream >> legal_fpga_id) {
        _legal_fpga_id_list.push_back(legal_fpga_id);
    }
    _fpga_count = _legal_fpga_id_list.size();

    TNamed *_input_machine_gun_samples_tnamed = (TNamed*) _input_root->Get("EventRecon_machine_gun_samples");
    if (_input_machine_gun_samples_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get machine gun samples from input file " << _root_file_path;
        return false;
    }
    _machine_gun_samples = std::stoi(_input_machine_gun_samples_tnamed->GetTitle());

    return true;
}

bool drawTwoASIC(TCanvas *_canvas, TLegend *_legend, int _bin, std::vector<std::vector<double>> &_fpgas_data, std::vector<unsigned short> _fpga_id, std::vector<Color_t> &_fpga_colors, std::string _y_title) {
    _legend->SetBorderSize(0);
    _legend->SetFillColor(0);
    _legend->SetFillStyle(0);
    _canvas->cd();
    // LOG(INFO) << "Drawing two ASICs with bin size: " << _bin;
    auto _fpga_count  = _fpgas_data.size();
    auto _fpga_events = _fpgas_data[0].size();
    int _bin_size = _fpga_events / _bin;
    if (_bin_size <= 0) {
        LOG(ERROR) << "Bin size is zero or negative: " << _bin_size;
        return false;
    }
    if (_fpga_count % 2 != 0) {
        LOG(ERROR) << "FPGA data size is not even: " << _fpga_count;
        return false;
    }
    if (_fpga_count < 2) {
        LOG(ERROR) << "FPGA data size is less than 2: " << _fpga_count;
        return false;
    }
    _fpga_count /= 2; // we are drawing two ASICs, so we divide the count by 2

    double _fpgas_value_max = 0.0;
    std::vector<double> _plot_x_values; // binned x values
    std::vector<double> _plot_x_errors;
    bool _plot_x_set = false;

    std::vector<std::vector<double>> _fpgas_values_binned;
    std::vector<std::vector<double>> _fpgas_errors_binned;

    _fpgas_values_binned.resize(_fpga_count * 2);
    _fpgas_errors_binned.resize(_fpga_count * 2);

    for (int i = 0; i < _fpga_count; ++i) {
        // LOG(DEBUG) << "Processing FPGA " << i;
        // get the average and error for each bin
        auto &_fpga_data_a0 = _fpgas_data[i * 2];
        auto &_fpga_data_a1 = _fpgas_data[i * 2 + 1];

        for (int j = 0; j < _fpga_events; j+= _bin_size) {
            // LOG(DEBUG) << "Processing bin starting at index " << j;
            double _bin_sum_a0 = 0.0;
            double _bin_sum_a1 = 0.0;
            double _bin_sum_sq_a0 = 0.0;
            double _bin_sum_sq_a1 = 0.0;
            double _bin_count = 0.0;

            for (int k = 0; k < _bin_size && j + k < _fpga_events; k++) {
                if (j + k < _fpga_events) {
                    if (_fpga_data_a0[j + k] >= 0 && _fpga_data_a1[j + k] >= 0) {
                        _bin_sum_a0 += _fpga_data_a0[j + k];
                        _bin_sum_a1 += _fpga_data_a1[j + k];
                        _bin_sum_sq_a0 += _fpga_data_a0[j + k] * _fpga_data_a0[j + k];
                        _bin_sum_sq_a1 += _fpga_data_a1[j + k] * _fpga_data_a1[j + k];
                        _bin_count += 1.0;
                    }
                }
            }

            if (_bin_count > 0) {
                _fpgas_values_binned[i * 2].push_back(_bin_sum_a0 / _bin_count);
                _fpgas_values_binned[i * 2 + 1].push_back(_bin_sum_a1 / _bin_count);
                _fpgas_errors_binned[i * 2].push_back(sqrt(_bin_sum_sq_a0 / _bin_count - (_bin_sum_a0 / _bin_count) * (_bin_sum_a0 / _bin_count)) / sqrt(_bin_count));
                _fpgas_errors_binned[i * 2 + 1].push_back(sqrt(_bin_sum_sq_a1 / _bin_count - (_bin_sum_a1 / _bin_count) * (_bin_sum_a1 / _bin_count)) / sqrt(_bin_count));
                if (_fpgas_values_binned[i * 2].back() > _fpgas_value_max) {
                    _fpgas_value_max = _fpgas_values_binned[i * 2].back();
                }
                if (_fpgas_values_binned[i * 2 + 1].back() > _fpgas_value_max) {
                    _fpgas_value_max = _fpgas_values_binned[i * 2 + 1].back();
                }
            } else {
                _fpgas_values_binned[i * 2].push_back(0.0);
                _fpgas_values_binned[i * 2 + 1].push_back(0.0);
                _fpgas_errors_binned[i * 2].push_back(0.0);
                _fpgas_errors_binned[i * 2 + 1].push_back(0.0);
            }
            if (!_plot_x_set) {
                if (j + _bin_size < _fpga_events) {
                    _plot_x_values.push_back(j + _bin_size / 2.0); // center of the bin
                    _plot_x_errors.push_back(_bin_size / 2.0); // half of the bin size
                } else {
                    _plot_x_values.push_back(j + (_fpga_events - j) / 2.0); // center of the last bin
                    _plot_x_errors.push_back((_fpga_events - j) / 2.0); // half of the last bin size
                }
            }
        }
        _plot_x_set = true;
    }
    // LOG(INFO) << "Max FPGA value: " << _fpgas_value_max;

    // the actual drawing
    for (int i = 0; i < _fpga_count; ++i) {
        auto _id = _fpga_id[i];
        auto *_a0_graph_error = new TGraphErrors(_plot_x_values.size(), _plot_x_values.data(), _fpgas_values_binned[i * 2].data(), _plot_x_errors.data(), _fpgas_errors_binned[i * 2].data());
        auto *_a1_graph_error = new TGraphErrors(_plot_x_values.size(), _plot_x_values.data(), _fpgas_values_binned[i * 2 + 1].data(), _plot_x_errors.data(), _fpgas_errors_binned[i * 2 + 1].data());

        // Draw the average band
        double _a0_average = std::accumulate(_fpgas_values_binned[i * 2].begin(), _fpgas_values_binned[i * 2].end(), 0.0) / _fpgas_values_binned[i * 2].size();
        double _a1_average = std::accumulate(_fpgas_values_binned[i * 2 + 1].begin(), _fpgas_values_binned[i * 2 + 1].end(), 0.0) / _fpgas_values_binned[i * 2 + 1].size();
        double _a0_average_error = 0.0;
        double _a1_average_error = 0.0;
        if (_fpgas_errors_binned[i * 2].size() > 0) {
            _a0_average_error = std::accumulate(_fpgas_errors_binned[i * 2].begin(), _fpgas_errors_binned[i * 2].end(), 0.0) / _fpgas_errors_binned[i * 2].size();
        }
        if (_fpgas_errors_binned[i * 2 + 1].size() > 0) {
            _a1_average_error = std::accumulate(_fpgas_errors_binned[i * 2 + 1].begin(), _fpgas_errors_binned[i * 2 + 1].end(), 0.0) / _fpgas_errors_binned[i * 2 + 1].size();
        }

        std::string _legend_name_a0 = "FPGA " + std::to_string(_id) + "(A0), " + _y_title + ": " + std::to_string(_a0_average) + "#pm" + std::to_string(_a0_average_error);
        std::string _legend_name_a1 = "FPGA " + std::to_string(_id) + "(A1), " + _y_title + ": " + std::to_string(_a1_average) + "#pm" + std::to_string(_a1_average_error);

        _legend->AddEntry(_a0_graph_error, _legend_name_a0.c_str(), "lep");
        _legend->AddEntry(_a1_graph_error, _legend_name_a1.c_str(), "lep");

        _a0_graph_error->SetMarkerColor(_fpga_colors[_id % _fpga_colors.size()]);
        _a0_graph_error->SetLineColor(_fpga_colors[_id % _fpga_colors.size()]);
        _a0_graph_error->SetMarkerStyle(20);
        _a0_graph_error->SetLineWidth(2);
        _a0_graph_error->SetTitle("");
        _a0_graph_error->GetXaxis()->SetTitle("Event Index");
        _a0_graph_error->GetYaxis()->SetTitle(_y_title.c_str());
        _a0_graph_error->GetYaxis()->SetTitleOffset(1.2);
        _a0_graph_error->GetYaxis()->SetRangeUser(0, _fpgas_value_max * 1.2);
        _a0_graph_error->SetMarkerSize(0.5);

        _a1_graph_error->SetMarkerColor(_fpga_colors[_id % _fpga_colors.size()] + 1);
        _a1_graph_error->SetLineColor(_fpga_colors[_id % _fpga_colors.size()] + 1);
        _a1_graph_error->SetMarkerStyle(21);
        _a1_graph_error->SetLineWidth(2);
        _a1_graph_error->SetTitle("");
        _a1_graph_error->SetMarkerSize(0.5);

        if (i==0) {
            _a0_graph_error->Draw("AP");
        }
        else {
            _a0_graph_error->Draw("P");
        }
        _a1_graph_error->Draw("P same");

        
        TGraphErrors *average_band_a0 = new TGraphErrors(_plot_x_values.size());
        TGraphErrors *average_band_a1 = new TGraphErrors(_plot_x_values.size());
        for (size_t j = 0; j < _plot_x_values.size(); ++j) {
            average_band_a0->SetPoint(j, _plot_x_values[j], _a0_average);
            average_band_a0->SetPointError(j, _plot_x_errors[j], _a0_average_error);
            average_band_a1->SetPoint(j, _plot_x_values[j], _a1_average);
            average_band_a1->SetPointError(j, _plot_x_errors[j], _a1_average_error);
        }
        average_band_a0->SetFillColorAlpha(_fpga_colors[_id % _fpga_colors.size()], 0.3);
        average_band_a0->SetFillStyle(3001); // solid fill
        average_band_a1->SetFillColorAlpha(_fpga_colors[_id % _fpga_colors.size()] + 1, 0.3);
        average_band_a1->SetFillStyle(3001); // solid fill
        average_band_a0->Draw("3 same");
        average_band_a1->Draw("3 same");

    }

    _legend->Draw("same");

    return true;
}