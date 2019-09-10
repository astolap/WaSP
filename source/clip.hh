/* clip.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef CLIP_HH
#define CLIP_HH

template<class T>
T clip(T in, const T min, const T max) {
  if (in > max)
    return max;
  if (in < min)
    return min;
  return in;
}

#endif
