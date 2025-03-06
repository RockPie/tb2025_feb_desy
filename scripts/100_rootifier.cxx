#include <iostream>
#include <unistd.h>
#include "TCanvas.h" 
#include "TVectorD.h"
#include "TVector.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "easylogging++.h"
#include "argparse/argparse.hpp"

INITIALIZE_EASYLOGGINGPP

void set_easylogger();

int main(int argc, char **argv){
    std::string script_input_file, script_output_file;
    int script_n_events;
    std::string script_name = __FILE__;

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    argparse::ArgumentParser program(script_name);
    program.add_argument("-f", "--file").help("Input .h2g file").required();
    program.add_argument("-o", "--output").help("Output root file").required();
    program.add_argument("-n", "--events").help("Number of events").default_value(-1);

    try {
        program.parse_args(argc, argv);
        script_input_file  = program.get<std::string>("--file");
        script_output_file = program.get<std::string>("--output");
        script_n_events    = program.get<int>("--events");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        return 1;
    }

    LOG(INFO) << "Script name: " << script_name;
    LOG(INFO) << "Input file: " << script_input_file;
    LOG(INFO) << "Output file: " << script_output_file;
    LOG(INFO) << "Number of events: " << script_n_events;

    return 0;
}

void set_easylogger(){
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s}[%levshort] (%fbase) %msg");
    defaultConf.set(el::Level::Info,    el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;34m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;33m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Error,   el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;31m%levshort\033[0m] (%fbase) %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);
}