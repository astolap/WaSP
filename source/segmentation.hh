#ifndef SEGMENTATION_HH
#define SEGMENTATION_HH

#include <vector>

struct segmentation {

    std::vector<int32_t> seg;

    std::vector<int32_t> region_histogram;

    std::vector<std::pair<int32_t, int32_t>> region_sizes;

    int32_t number_of_regions;

};

segmentation normdispsegmentation(
    const std::vector<uint16_t> &img,
    const uint32_t iterations,
    const int32_t nr,
    const int32_t nc);

#endif