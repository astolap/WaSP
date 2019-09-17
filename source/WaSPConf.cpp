/*BSD 2-Clause License
* Copyright(c) 2019, Pekka Astola
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
*
* 1. Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*     OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string.h>
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

    //if (argc < 7) {
    //    return false;
    //}

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

    if (WaSP_encoder_setup.input_directory.length() == 0) {
        printf("\n Input directory not set\n");
        return false;
    }

    if (WaSP_encoder_setup.output_directory.length() == 0) {
        printf("\n Output directory not set\n");
        return false;
    }

    if (WaSP_encoder_setup.wasp_kakadu_directory.length() == 0) {
        printf("\n Kakadu directory not set\n");
        return false;
    }

    return true;

}

void WaSPConfig::print_encoder_help() {
    printf("\n\tUsage: wasp-encoder"
        "\n\t--input [INPUT DIRECTORY .PPM/.PGM]"
        "\n\t--output [OUTPUT DIRECTORY .LF/.PPM/.PGM]"
        "\n\t--config [JSON CONFIG]"
        "\n\t--kakadu [KAKADU BINARY DIRECTORY]"
        "\n\t--sparse_subsampling [Subsampling factor when solving sparse filter."
        "\n\t\tneeds to be integer >0. Values 2 or 4 will increase encoder speed with some loss in PSNR.]\n\n");
    return;
}

void WaSPConfig::print_decoder_help() {
    printf("\n\tUsage: wasp-decoder"
        "\n\t--input [INPUT .LF]"
        "\n\t--output [OUTPUT DIRECTORY .PPM/.PGM]"
        "\n\t--kakadu [KAKADU BINARY DIRECTORY]\n\n");
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

    //if (argc < 9) {
    //    return false;
    //}

    //if (argc > 10) {
    //    return false;
    //}

    WaSP_encoder_setup.sparse_subsampling = 1;

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

        else if (!strcmp(argv[ii], "-s")) {
            WaSP_encoder_setup.sparse_subsampling = atoi(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--sparse_subsampling")) {
            WaSP_encoder_setup.sparse_subsampling = atoi(argv[ii + 1]);

        }

        else {
            return false;
        }

    }

    if (WaSP_encoder_setup.config_file.length() == 0) {
        printf("\n Config file (.json) not set\n");
        return false;
    }

    if (WaSP_encoder_setup.input_directory.length() == 0) {
        printf("\n Input directory not set\n");
        return false;
    }

    if (WaSP_encoder_setup.output_directory.length() == 0) {
        printf("\n Output directory not set\n");
        return false;
    }

    if (WaSP_encoder_setup.wasp_kakadu_directory.length() == 0) {
        printf("\n Kakadu directory not set\n");
        return false;
    }

    if (WaSP_encoder_setup.sparse_subsampling < 1) {
        printf("\n Sub sampling factor needs to be >= 1\n");
        return false;
    }

    return true;

}