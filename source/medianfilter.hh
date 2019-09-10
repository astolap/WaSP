/* medianfilter.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef MEDIANFILTER_HH
#define MEDIANFILTER_HH

#include <vector>
#include <algorithm>
#include <cmath>

template<class T>
T getMedian(std::vector<T> scores) {
  size_t s = scores.size();
  size_t n = s / 2;
  nth_element(scores.begin(), scores.begin() + n, scores.end());
  //sort(scores.begin(), scores.end());
  //double avg = (static_cast<double>(scores[n - 1]) + static_cast<double>(scores[n])) / 2.0;
  //T median = (s % 2 == 0) ? static_cast<T>(floor(avg+0.5)) : scores[n];
  T median = scores[n];
  //avg = floor(avg + 0.5);
  //T median = static_cast<T>(avg);
  return median;
}

template<class T>
void medfilt2D(
    T* input, 
    T* output, 
    const int32_t SZ, 
    const int32_t nr, 
    const int32_t nc) {

  int32_t dsz = (SZ / 2);
  std::vector<T> scores;

  for (int32_t y = 0; y < nr; y++) {
    for (int32_t x = 0; x < nc; x++) {
      scores.clear();
      for (int32_t dy = -dsz; dy <= dsz; dy++) {
        for (int32_t dx = -dsz; dx <= dsz; dx++) {
          if ((y + dy) >= 0 && (y + dy) < nr && (x + dx) >= 0 && (x + dx) < nc)
            scores.push_back(input[y + dy + (x + dx) * nr]);
        }
      }
      output[y + x * nr] = getMedian(scores);
    }
  }
}

template<class T>
T *medfilt2D(
    T* input, 
    const int32_t SZ,
    const int32_t nr,
    const int32_t nc) {

    T* output = new T[nr*nc]();

    int32_t dsz = (SZ / 2);
    std::vector<T> scores;

    for (int32_t y = 0; y < nr; y++) {
        for (int32_t x = 0; x < nc; x++) {
            scores.clear();
            for (int32_t dy = -dsz; dy <= dsz; dy++) {
                for (int32_t dx = -dsz; dx <= dsz; dx++) {
                    if ((y + dy) >= 0 && (y + dy) < nr && (x + dx) >= 0 && (x + dx) < nc)
                        scores.push_back(input[y + dy + (x + dx) * nr]);
                }
            }
            output[y + x * nr] = getMedian(scores);
        }
    }

    return output;
}

#endif
