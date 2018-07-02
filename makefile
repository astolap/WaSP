CC=g++
CFLAGS=-I. -std=c++11 -fopenmp -O3
HH_DEPS = source/bitdepth.hh source/clip.hh source/fastols.hh source/fileaux.hh source/inpainting.hh source/medianfilter.hh source/merging.hh source/minconf.hh source/ppm.hh source/predictdepth.hh source/psnr.hh source/residualjp2.hh source/sparsefilter.hh source/warping.hh source/view.hh source/ycbcr.hh
OBJ_ENCODER = source/wasp-encoder.o source/bitdepth.cpp source/clip.cpp source/fastols.cpp source/fileaux.cpp source/inpainting.cpp source/medianfilter.cpp source/merging.cpp source/minconf.cpp source/ppm.cpp source/predictdepth.cpp source/psnr.cpp source/residualjp2.cpp source/sparsefilter.cpp source/warping.cpp source/view.cpp source/ycbcr.cpp
OBJ_DECODER = source/wasp-decoder.o source/bitdepth.cpp source/clip.cpp source/fastols.cpp source/fileaux.cpp source/inpainting.cpp source/medianfilter.cpp source/merging.cpp source/minconf.cpp source/ppm.cpp source/predictdepth.cpp source/psnr.cpp source/residualjp2.cpp source/sparsefilter.cpp source/warping.cpp source/view.cpp source/ycbcr.cpp
 
%.o: %.cpp $(HH_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tut-hdca-encoder: $(OBJ_ENCODER)
	$(CC) -o $@ $^ $(CFLAGS)

tut-hdca-decoder: $(OBJ_DECODER)
	$(CC) -o $@ $^ $(CFLAGS)

all: tut-hdca-encoder tut-hdca-decoder
 
clean:
	rm TUT-HDCA-Decoder/tut-hdca-decoder.o TUT-HDCA-Encoder/tut-hdca-encoder.o
