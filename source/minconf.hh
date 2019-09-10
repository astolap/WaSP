/* minconf.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef MINCONF_HH
#define MINCONF_HH

#include "view.hh"

struct minimal_config { /* this goes to bitstream */

  uint8_t r, c;  // SAI row,column subscript
  uint8_t mmode; /* view merging mode, 0=LS, 1=geometric weight, 2=median*/
  uint8_t level;
  uint16_t encoding_flags;

};

minimal_config makeMinimalConfig(view *view0);
void setup_form_minimal_config(minimal_config *mconf, view *view0);

#endif
