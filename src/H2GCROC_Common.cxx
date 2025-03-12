#include "H2GCROC_Common.hxx"

void set_easylogger(){
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s}[%levshort] (%fbase) %msg");
    defaultConf.set(el::Level::Info, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;34m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;33m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Error, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;31m%levshort\033[0m] (%fbase) %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

ScriptOptions parse_arguments_single_root(int argc, char **argv, const std::string& version) {
    ScriptOptions opts;

    opts.script_version = version;
    opts.script_name = argv[0];
    opts.script_name = opts.script_name.substr(opts.script_name.find_last_of("/\\") + 1);
    opts.script_name = opts.script_name.substr(0, opts.script_name.find_last_of("."));

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    argparse::ArgumentParser program(opts.script_name, opts.script_version);

    program.add_argument("-f", "--file").help("Input .root file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    program.add_argument("-v", "--verbose").help("Verbose mode").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
        opts.input_file = program.get<std::string>("--file");
        opts.output_file = program.get<std::string>("--output");
        auto n_events_str = program.get<std::string>("--events");
        opts.n_events = std::stoi(n_events_str);
        opts.verbose = program.get<bool>("--verbose");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        exit(1);
    }

    if (access(opts.input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << opts.input_file << " does not exist!";
        exit(1);
    }

    if (opts.input_file.substr(opts.input_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Input file " << opts.input_file << " should end with .root!";
        exit(1);
    }

    if (opts.output_file.substr(opts.output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << opts.output_file << " should end with .root!";
        exit(1);
    }

    opts.output_folder = opts.output_file.substr(0, opts.output_file.find_last_of("/\\"));
    if (opts.output_folder.empty()) {
        opts.output_folder = "./dump/" + opts.script_name;
    }

    if (access(opts.output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << opts.output_folder;
        if (mkdir(opts.output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << opts.output_folder;
            exit(1);
        }
    }

    if (access(opts.output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << opts.output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << opts.script_name;
    LOG(INFO) << "Input file: " << opts.input_file;
    LOG(INFO) << "Output file: " << opts.output_file << " in " << opts.output_folder;
    LOG(INFO) << "Number of events: " << opts.n_events;

    return opts;
}

ScriptOptions parse_arguments_single_json(int argc, char **argv, const std::string& version) {
    ScriptOptions opts;

    opts.script_version = version;
    opts.script_name = argv[0];
    opts.script_name = opts.script_name.substr(opts.script_name.find_last_of("/\\") + 1);
    opts.script_name = opts.script_name.substr(0, opts.script_name.find_last_of("."));

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    argparse::ArgumentParser program(opts.script_name, opts.script_version);

    program.add_argument("-f", "--file").help("Input .json file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    program.add_argument("-v", "--verbose").help("Verbose mode").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
        opts.input_file = program.get<std::string>("--file");
        opts.output_file = program.get<std::string>("--output");
        auto n_events_str = program.get<std::string>("--events");
        opts.n_events = std::stoi(n_events_str);
        opts.verbose = program.get<bool>("--verbose");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        exit(1);
    }

    if (access(opts.input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << opts.input_file << " does not exist!";
        exit(1);
    }

    if (opts.input_file.substr(opts.input_file.find_last_of(".") + 1) != "json") {
        LOG(ERROR) << "Input file " << opts.input_file << " should end with .json!";
        exit(1);
    }

    if (opts.output_file.substr(opts.output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << opts.output_file << " should end with .root!";
        exit(1);
    }

    opts.output_folder = opts.output_file.substr(0, opts.output_file.find_last_of("/\\"));
    if (opts.output_folder.empty()) {
        opts.output_folder = "./dump/" + opts.script_name;
    }

    if (access(opts.output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << opts.output_folder;
        if (mkdir(opts.output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << opts.output_folder;
            exit(1);
        }
    }

    if (access(opts.output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << opts.output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << opts.script_name;
    LOG(INFO) << "Input file: " << opts.input_file;
    LOG(INFO) << "Output file: " << opts.output_file << " in " << opts.output_folder;
    LOG(INFO) << "Number of events: " << opts.n_events;

    return opts;
}

int get_valid_fpga_channel(int fpga_channel){
    if (fpga_channel < 0 || fpga_channel >= FPGA_CHANNEL_NUMBER) {
        LOG(ERROR) << "Invalid FPGA channel " << fpga_channel;
        return -1;
    }
    std::vector<int> CM_channels = {0, 38, 76, 114};
    for (auto CM_channel : CM_channels) {
        if (fpga_channel == CM_channel) {
            return -1;
        }
    }
    std::vector<int> Calib_channels = {19, 57, 95, 133};
    for (auto Calib_channel : Calib_channels) {
        if (fpga_channel == Calib_channel) {
            return -1;
        }
    }
    if (fpga_channel < 19) {
        return fpga_channel - 1;
    }
    if (fpga_channel < 38) {
        return fpga_channel - 2;
    }
    if (fpga_channel < 57) {
        return fpga_channel - 3;
    }
    if (fpga_channel < 76) {
        return fpga_channel - 4;
    }
    if (fpga_channel < 95) {
        return fpga_channel - 5;
    }
    if (fpga_channel < 114) {
        return fpga_channel - 6;
    }
    if (fpga_channel < 133) {
        return fpga_channel - 7;
    }
    if (fpga_channel < 152) {
        return fpga_channel - 8;
    }
    return -1;
}

int get_total_fpga_channel(int fpga_channel){
    if (fpga_channel < 0 || fpga_channel >= FPGA_CHANNEL_NUMBER_VALID) {
        LOG(ERROR) << "Invalid FPGA channel " << fpga_channel;
        return -1;
    }
    if (fpga_channel <= 17) {
        return fpga_channel + 1;
    }
    if (fpga_channel <= 35) {
        return fpga_channel + 2;
    }
    if (fpga_channel <= 53) {
        return fpga_channel + 3;
    }
    if (fpga_channel <= 71) {
        return fpga_channel + 4;
    }
    if (fpga_channel <= 89) {
        return fpga_channel + 5;
    }
    if (fpga_channel <= 107) {
        return fpga_channel + 6;
    }
    if (fpga_channel <= 125) {
        return fpga_channel + 7;
    }
    if (fpga_channel <= 143) {
        return fpga_channel + 8;
    }
    return -1;
}