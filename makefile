CC=g++
CFLAGS=-I. -std=c++11 -fopenmp -O2 -g -lstdc++fs -fPIC
HH_DEPS = source/bitdepth.hh source/clip.hh source/fastols.hh source/fileaux.hh source/inpainting.hh source/medianfilter.hh source/merging.hh source/minconf.hh source/ppm.hh source/predictdepth.hh source/psnr.hh source/residualjp2.hh source/sparsefilter.hh source/warping.hh source/view.hh source/ycbcr.hh source/WaSPDecoder.hh source/WaSPEncoder.hh source/WaSPConf.hh
OBJ_ENCODER = source/wasp-encoder.o source/WaSPEncoder.cpp source/WaSPConf.cpp source/fastols.cpp source/codestream.cpp source/fileaux.cpp source/merging.cpp source/minconf.cpp source/ppm.cpp source/predictdepth.cpp source/psnr.cpp source/residualjp2.cpp source/sparsefilter.cpp source/warping.cpp source/view.cpp source/ycbcr.cpp
OBJ_DECODER = source/wasp-decoder.o source/WaSPDecoder.cpp source/WaSPConf.cpp source/fastols.cpp source/codestream.cpp source/fileaux.cpp source/merging.cpp source/minconf.cpp source/ppm.cpp source/predictdepth.cpp source/psnr.cpp source/residualjp2.cpp source/sparsefilter.cpp source/warping.cpp source/view.cpp source/ycbcr.cpp
 
%.o: %.cpp $(HH_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

wasp-encoder-bin: $(OBJ_ENCODER)
	$(CC) -o $@ $^ $(CFLAGS)

wasp-decoder-bin: $(OBJ_DECODER)
	$(CC) -o $@ $^ $(CFLAGS)

all: wasp-encoder-bin wasp-decoder-bin
 
clean:
	rm source/*.o
