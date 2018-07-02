#ifndef MEDIANFILTER_HH
#define MEDIANFILTER_HH

#include <vector>
#include <algorithm>

template <class T>
T getMedian(std::vector<T> scores)
{
	T median;
	size_t size = scores.size();

	std::sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}

template <class T>
void medfilt2D(T* input, T* output, int SZ, int nr, int nc)
{
	int dsz = (SZ / 2);
	std::vector<T> scores;

	for (int y = 0; y < nr; y++) {
		for (int x = 0; x < nc; x++) {
			scores.clear();
			for (int dy = -dsz; dy < dsz; dy++) {
				for (int dx = -dsz; dx < dsz; dx++) {
					if ((y + dy) >= 0 && (y + dy) < nr
						&& (x + dx) >= 0 && (x + dx) < nc)
						scores.push_back(input[y + dy + (x + dx)*nr]);
				}
			}
			output[y + x*nr] = getMedian(scores);
		}
	}
}

#endif