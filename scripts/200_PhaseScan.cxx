#include "H2GCROC_Common.hxx"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments(argc, argv, "1.0");

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
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
    int fpga_number = legal_fpga_id_list.size();

    TNamed *input_machine_gun_samples_tnamed = (TNamed*) input_root->Get("EventRecon_machine_gun_samples");
    if (input_machine_gun_samples_tnamed == nullptr) {
        LOG(ERROR) << "Failed to get machine gun samples from input file " << opts.input_file;
        return 1;
    }
    int machine_gun_samples = std::stoi(input_machine_gun_samples_tnamed->GetTitle());

    LOG(INFO) << "Number of FPGAs: " << fpga_number;
    LOG(INFO) << "Machine gun samples: " << machine_gun_samples;

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
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

    for (int _fpga_index = 0; _fpga_index < fpga_number; _fpga_index++) {
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

        input_tree->SetBranchAddress(("timestamps_" + std::to_string(_fpga_id)).c_str(), branch_timestamps);
        input_tree->SetBranchAddress(("daqh_list_" + std::to_string(_fpga_id)).c_str(), branch_daqh_list);
        input_tree->SetBranchAddress(("tc_list_" + std::to_string(_fpga_id)).c_str(), branch_tc_list);
        input_tree->SetBranchAddress(("tp_list_" + std::to_string(_fpga_id)).c_str(), branch_tp_list);
        input_tree->SetBranchAddress(("val0_list_" + std::to_string(_fpga_id)).c_str(), branch_val0_list);
        input_tree->SetBranchAddress(("val1_list_" + std::to_string(_fpga_id)).c_str(), branch_val1_list);
        input_tree->SetBranchAddress(("val2_list_" + std::to_string(_fpga_id)).c_str(), branch_val2_list);
        input_tree->SetBranchAddress(("crc32_list_" + std::to_string(_fpga_id)).c_str(), branch_crc32_list);
        input_tree->SetBranchAddress(("last_heartbeat_" + std::to_string(_fpga_id)).c_str(), branch_last_heartbeat);

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

    int target_channel = 161;
    std::vector<double> sum_waveform_val0;
    sum_waveform_val0.resize(machine_gun_samples, 0);
    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        auto _channel_fpga = target_channel / FPGA_CHANNEL_NUMBER;
        auto _channel_fpga_index = -1;
        for (int i = 0; i < fpga_number; i++) {
            if (legal_fpga_id_list[i] == _channel_fpga) {
                _channel_fpga_index = i;
                break;
            }
        }
        auto _channel_index = target_channel % FPGA_CHANNEL_NUMBER;

        for (int _machinegun_index = 0; _machinegun_index < machine_gun_samples; _machinegun_index++) {
            auto _val0 = branch_val0_list_list[_channel_fpga_index][_machinegun_index * FPGA_CHANNEL_NUMBER + _channel_index];
            sum_waveform_val0[_machinegun_index] += _val0;
        }
    }

    TGraph *sum_waveform_val0_graph = new TGraph(machine_gun_samples);
    for (int _machinegun_index = 0; _machinegun_index < machine_gun_samples; _machinegun_index++) {
        sum_waveform_val0_graph->SetPoint(_machinegun_index, _machinegun_index, sum_waveform_val0[_machinegun_index]);
    }

    sum_waveform_val0_graph->SetTitle(("Sum waveform of channel " + std::to_string(target_channel)).c_str());
    sum_waveform_val0_graph->GetXaxis()->SetTitle("Time (ns)");
    sum_waveform_val0_graph->GetYaxis()->SetTitle("ADC counts");

    output_root->cd();
    sum_waveform_val0_graph->Write();

    input_root->Close();
    output_root->Close();
    return 0;
}