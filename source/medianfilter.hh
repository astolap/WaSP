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
