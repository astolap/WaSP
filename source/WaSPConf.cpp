/* WaSPConf.cpp */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#include "WaSPConf.hh"

WaSPConfig::WaSPConfig(int argc, char *argv[], const char *type) {

    if (!strcmp(type,"encoder") && !parseCommandLine_encoder(argc, argv)) {
        printf("\n HELP \n");
        exit(0);
    }

    if (!strcmp(type, "decoder") && !parseCommandLine_decoder(argc, argv)) {
        printf("\n HELP \n");
        exit(0);
    }

}

WaSPConfig::~WaSPConfig() {

}

bool WaSPConfig::parseCommandLine_decoder(int argc, char *argv[]) {

    if (argc < 7) {
        return false;
    }

    for (int32_t ii = 1; ii < argc - 1; ii += 2) {

        if (!strcmp(argv[ii], "-i")) {
            WaSP_encoder_setup.input_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--input")) {
            WaSP_encoder_setup.input_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-o")) {
            WaSP_encoder_setup.output_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--output")) {
            WaSP_encoder_setup.output_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-k")) {
            WaSP_encoder_setup.wasp_kakadu_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--kakadu")) {
            WaSP_encoder_setup.wasp_kakadu_directory = std::string(argv[ii + 1]);
        }

        else {
            return false;
        }

    }

    return true;

}

bool WaSPConfig::parseCommandLine_encoder(int argc, char *argv[]) {

    if (argc < 9) {
        return false;
    }

    for (int32_t ii = 1; ii < argc-1; ii+=2) {

        if (!strcmp(argv[ii], "-c")) {
            WaSP_encoder_setup.config_file = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--config")) {
            WaSP_encoder_setup.config_file = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-i")) {
            WaSP_encoder_setup.input_directory = std::string(argv[ii+1]);
        }

        else if (!strcmp(argv[ii], "--input")) {
            WaSP_encoder_setup.input_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-o")) {
            WaSP_encoder_setup.output_directory = std::string(argv[ii+1]);
        }

        else if (!strcmp(argv[ii], "--output")) {
            WaSP_encoder_setup.output_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-k")) {
            WaSP_encoder_setup.wasp_kakadu_directory = std::string(argv[ii+1]);
        }

        else if (!strcmp(argv[ii], "--kakadu")) {
            WaSP_encoder_setup.wasp_kakadu_directory = std::string(argv[ii + 1]);

        }

        else {
            return false;
        }

    }

    return true;

}