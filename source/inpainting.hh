/* inpainting.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

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
