/* WaSPConf.cpp */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#include "WaSPConf.hh"

WaSPConfig::WaSPConfig(int argc, char *argv[], const char *type) {

    print_intro();

    if (!strcmp(type,"encoder") && !parseCommandLine_encoder(argc, argv)) {
        print_encoder_help();
        exit(0);
    }

    if (!strcmp(type, "decoder") && !parseCommandLine_decoder(argc, argv)) {
        print_decoder_help();
        exit(0);
    }

}

WaSPConfig::~WaSPConfig() {

}

bool WaSPConfig::parseCommandLine_decoder(int argc, char *argv[]) {

    if (argc < 7) {
        return false;
    }

    if (argc > 7) {
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

void WaSPConfig::print_encoder_help() {
    printf("\n\tUsage: wasp-encoder"
        "\n\t--input [INPUT DIRECTORY .PPM/.PGM]"
        "\n\t--output [OUTPUT DIRECTORY .PPM/.PGM]"
        "\n\t--config [JSON CONFIG]"
        "\n\t--kakadu [KAKADU BINARY DIRECTORY]\n");
    return;
}

void WaSPConfig::print_decoder_help() {
    printf("\n\tUsage: wasp-decoder"
        "\n\t--output [OUTPUT DIRECTORY .PPM/.PGM]"
        "\n\t--config [JSON CONFIG]"
        "\n\t--kakadu [KAKADU BINARY DIRECTORY]\n");
    return;
}

void WaSPConfig::print_intro() {
    printf("\n\t--------------------------------------------------------\n");
    printf(
        "\n\tWaSP - Warping and Sparse Prediction"
        "\n\tAuthor: Pekka Astola (pekka.astola@tuni.fi), 2018-2019"
        "\n\tWeb page: https://github.com/astolap/WaSP"
        "\n\tHome page of author: http://www.cs.tut.fi/~astolap/"
        "\n\tPublication: P. Astola and I. Tabus,"
        "\n\t\t\tWaSP: Hierarchical Warping, Merging, and Sparse Prediction for Light Field Image Compression,"
        "\n\t\t\t2018 7th European Workshop on Visual Information Processing (EUVIP), Tampere, 2018, pp. 1-6."
        "\n\n");
    printf("\n\t--------------------------------------------------------\n");
    return;
}

bool WaSPConfig::parseCommandLine_encoder(int argc, char *argv[]) {

    if (argc < 9) {
        return false;
    }

    if (argc > 9) {
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