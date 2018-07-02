# WaSP light field compression

## Introduction

This is the WaSP (Warping and Sparse Prediction) light field compression software.

## Installation instructions

This software has been developed and tested on Windows 7,10 and has mostly been developed using Visual Studio. However, a makefile for compiling on Linux is provided and it has been tested on **Ubuntu 14.04 LTS**. Run 

>make all

to build both encoder and decoder.

### Kakadu installation on Linux

Some of the encoding is currently done by JPEG2000 and for that we use the Kakadu Software

[Download Kakadu for Linux] (http://kakadusoftware.com/downloads/) and see Kakadu's README.txt for instructions regarding **LD\_LIBRARY\_PATH**. 

If you encounter the error *"kakadu/kdu_compress: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.21'* not found" do the following,

>sudo add-apt-repository ppa:ubuntu-toolchain-r/test
 
>sudo apt-get update

>sudo apt-get install

>sudo apt-get install libstdc++6

## Demo

An encoding of the Bikes_01 at rate 0.75 bpp is provided as an example for the decoder. You can run the decoder with,

>./wasp-decoder /path/to/Bikes\_01_0.75.LF /path/to/output /path/to/Kakadu/

For encoding a sample configuration file is provided for the same Bikes_01 encoding. You can run that by,

>./wasp-encoder /path/to/original_images/ /path/to/output/directory /path/to/Kakadu/ /path/to/Bikes\_01\_0.75.conf

## Running the encoder

Encoder takes four arguments: path to input images, path to output directory, path to Kakadu directory and path to the configuration file. 

Encoder assumes that the images are named using col_row.ppm in %03d format, for example view row=7,col=3 should be in a file named **003\_007.ppm**. For the UNSW style inverse depths, the encoder assumes that depth files are named as **003\_007.pgm** and are found in the same folder as the color images.

Encoder also assumes that the directory pointing to Kakadu installation contains **kdu_compress** and **kdu_expand** executables.

### Configuration file

Currently, the configuration file is a collection of parameters (as 32-bit signed integers) generated with Matlab. The structure of the configuration file is the following (each line representing one int32), starting with the global (applied to each view) header information

1. total number of views to be encoded

2. YUV\_TRANSFORM, true (1) or false (0)

3. YUV\_RATIO_SEARCH, true (1) or false (0)

4. STD\_SEARCH, true (1) or false (0)

   continuing for each view to be encoded,

5. row index in light field, 7, (used to fetch .ppm files, the example here is view 003_007.ppm)

6. column index in light field, 3

7. camera position along the horizontal (e.g. from UNSW camera centers) multiplied with 100000

8. camera position along the vertical (e.g. from UNSW camera centers) multiplied with 100000

9. rate for the color component of the view multiplied with 100000

10. rate for the depth component of the view multiplied with 100000

11. Sparse filter order, for example 25

12. Sparse filter neighborhood, for example 3 would imply 7x7 neighborhood ( -3:3, -3:3 )

13. standard deviation used for geometry based view merging. Use this for only very low bit rates.

14. minimum inverse depth value to be subtracted from the inverse depth file. This is needed if the inverse depth file contains also negative values (as can be in the lenslet case). This value is represented in the integer range after multiplication by 2^14 and should always be positive.

15. number of reference views used in predicting the **color** component of this view, for example Nc,

   if Nc>0

   15+1:15+Nc. the indices of the reference views. Index for each view is one integer, and for example for the first view this integer == 0, the second view == 1 and so on. To generate a configuration file you need to keep track of the order.
	
16. (or 15+Nc+1) number of reference views used in predicting the **depth** component of this view, for example Nd,

   if Nd>0

   15+Nc+2:15+Nc+Nd+1. same instructions as for color views
	
17. (or 15+Nc+Nd+2) whether the view uses a segmentation, true (1) or false (0)

### Output of the encoder

The encoder outputs to the bitstream **output.LF** in the output directory. Additionally, it also decodes the views to the output directory.

## Running the decoder

The decoder takes three arguments: path to input file (for example path to the output.LF file), directory for outputting the decoded images (.ppm and .pgm) and path to the directory containing Kakadu executables.

## References

This work is based on academic research and any research based on this software should cite the following papers:

**P. Astola, I. Tabus, *Light Field Compression of HDCA Images Combining Linear Prediction and JPEG 2000*, EUSIPCO 2018**

**I. Tabus, P. Helin, P. Astola, *Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000*, ICIP 2017**