/* wasp-encoder.cpp */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#include "WaSPConf.hh"
#include "WaSPEncoder.hh"

int main(int argc, char* argv[]) {

    WaSPConfig encoder_settings(argc, argv, "encoder");

    WaSPEncoder wasp_encoder(encoder_settings.WaSP_encoder_setup);

    wasp_encoder.encode();

    exit(0);
}