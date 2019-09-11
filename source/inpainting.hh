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

#ifndef INPAINTING_HH
#define INPAINTING_HH

#include <cstdint>
#include <cstring>

#include "medianfilter.hh"

template<class T1, class T2>
uint32_t inpainting(
    T1 *pshort,
    const int32_t nr,
    const int32_t nc,
    const T2 maskval,
    T2 *mask_img) {

    std::vector<T1> neighbours;
    int32_t dsz = 1;

    //T2 *mask_im_work_copy = new T2[nr*nc]();
    //memcpy(mask_im_work_copy, mask_img, sizeof(T2)*nr*nc);

    uint32_t nholes_remaining = 0;

    for (int32_t ii = 0; ii < nr * nc; ii++) {
        int32_t y, x;
        y = ii % nr;
        x = (ii / nr);

        bool is_hole = 
            (static_cast<int32_t>(mask_img[ii]) == static_cast<int32_t>(maskval));

        if (is_hole) {
            neighbours.clear();
            for (int32_t dy = -dsz; dy <= dsz; dy++) {
                for (int32_t dx = -dsz; dx <= dsz; dx++) {
                    if (!(dy == 0 && dx == 0)) {

                        if (
                            (y + dy) >= 0 &&
                            (y + dy) < nr &&
                            (x + dx) >= 0 &&
                            (x + dx) < nc)
                        {

                            int32_t lin_ind = y + dy + (x + dx) * nr;

                            bool is_hole_again =
                                (static_cast<int32_t>(mask_img[lin_ind]) == static_cast<int32_t>(maskval));

                            if (!is_hole_again) 
                            {
                                neighbours.push_back(pshort[lin_ind]);
                            }

                        }
                    }
                }
            }
            if (neighbours.size() > 0) {
                pshort[ii] = getMedian(neighbours); // neighbours.at(2);//
                mask_img[ii] = mask_img[ii] + static_cast<T2>(1);
            }
            else {
                nholes_remaining++;
            }
        }
    }

    //delete[](mask_im_work_copy);

    //printf("holes remaining:\t%d\n", nholes_remaining);

    return nholes_remaining;
}

template<class T1, class T2>
uint32_t holefilling(
    T1 *pshort,
    const int32_t nr,
    const int32_t nc,
    const T2 maskval,
    const T2 *mask_img) {

    T2 *mask_img_work = new T2[nr*nc]();
    memcpy(mask_img_work, mask_img, sizeof(T2)*nr*nc);

    while (inpainting(
        pshort,
        nr,
        nc,
        maskval,
        mask_img_work) > 0);

    delete[](mask_img_work);

    return 0;
}

#endif
