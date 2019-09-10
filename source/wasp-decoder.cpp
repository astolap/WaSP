/* wasp-decoder.cpp */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#include "WaSPDecoder.hh"

int main(int argc, char* argv[]) {

    WaSPConfig encoder_settings(argc, argv, "decoder");

    WaSPDecoder wasp_decoder(encoder_settings.WaSP_encoder_setup);

    wasp_decoder.decode();

    exit(0);
}