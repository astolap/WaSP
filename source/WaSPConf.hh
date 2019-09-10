/* WaSPConf.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#include <string>
#include <iostream>

using namespace std;

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#ifndef SOURCE_WASPCONF_HH_
#define SOURCE_WASPCONF_HH_

struct WaSPsetup {

    string input_directory;
    string output_directory;
    string config_file;
    string wasp_kakadu_directory;

};

class WaSPConfig {
 private:
 
  bool parseCommandLine_encoder(int argc, char *argv[]);
  bool parseCommandLine_decoder(int argc, char *argv[]);

  void print_encoder_help();
  void print_decoder_help();
  void print_intro();
 
  public:
  
  virtual ~WaSPConfig();

  WaSPConfig(int argc, char *argv[], const char *type);
  
  WaSPsetup WaSP_encoder_setup;
  
  
};
 
#endif