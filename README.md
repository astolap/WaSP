
# WaSP - light field compression
### Table of contents
 1. [Introduction](#introduction)
 2. [Installing](#installing-and-compiling)
    1. [Kakadu JPEG 2000 encoder](#kakadu-installation)
 3. [Running the software](#running-the-software)

## Introduction

This is the WaSP (Warping and Sparse Prediction) light field compression software. The program is intended for encoding light fields, such as the [JPEG Pleno datasets](https://jpeg.org/plenodb/lf/pleno_lf/). 

The program is developed and maintained by [Pekka Astola](http://www.cs.tut.fi/~astolap/).

If you plan to use this program for research, remember to cite the following publication,

```
@INPROCEEDINGS{8611756,  
author={P. {Astola} and I. {Tabus}},  
booktitle={2018 7th European Workshop on Visual Information Processing (EUVIP)},  
title={WaSP: Hierarchical Warping, Merging, and Sparse Prediction for Light Field Image Compression},  
year={2018},  
volume={},  
number={},  
pages={1-6},  
keywords={cameras;data compression;image coding;video codecs;warped versions;warped references;optimal LS merger;occluded pixels;merged image;original view;sparse predictor;plenoptic camera images;high density camera array data;JPEG Pleno test conditions;earlier scheme;JPEG Pleno Light Field coding standard;current view;reference views;improved view merging;generalized functionality;codec;earlier version;previous hierarchical levels;particular level;versatile light field compression scheme;Light Field image compression;sparse prediction;hierarchical warping;Encoding;Cameras;Image color analysis;Image coding;Transform coding;Codecs;Merging},  
doi={10.1109/EUVIP.2018.8611756},  
ISSN={2471-8963},  
month={Nov},}
```

## Installing and compiling

This software has been developed using Visual Studio on Windows 10. For Visual Studio a solution file is provided. 

However, a [makefile](#[https://github.com/astolap/WaSP/blob/master/makefile](https://github.com/astolap/WaSP/blob/master/makefile)) for compiling on Linux is provided, and it has been tested on **Ubuntu 14.04 LTS**. On Linux simply run, 

>make all

to build both encoder and decoder. 

### Kakadu installation

Normalized disparity and texture residual encoding is performed by JPEG 2000 using the [Kakadu](https://kakadusoftware.com/) encoder.

[Download Kakadu for Linux or Windows](http://kakadusoftware.com/downloads/) and see Kakadu's README.txt for instructions regarding **LD\_LIBRARY\_PATH**. 

If you encounter the error *"kakadu/kdu_compress: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.21'* not found" do the following,

>sudo add-apt-repository ppa:ubuntu-toolchain-r/test
 
>sudo apt-get update

>sudo apt-get install libstdc++6

## Running the software

Download the light field data sets from [JPEG Pleno database](https://jpeg.org/plenodb/lf/pleno_lf/), and use one of the [configuration files](https://github.com/astolap/WaSP/blob/master/configuration_files) provided. The syntax for the encoder is,
> wasp-encoder --input [INPUT DIRECTORY .PPM/.PGM --output [OUTPUT DIRECTORY .LF] --config [JSON CONFIG FILE] --kakadu [KAKADU BINARY DIRECTORY].

The syntax for the decoder is,
> wasp-decoder --input [INPUT .LF] --output [OUTPUT DIRECTORY .PPM/.PGM] --kakadu [KAKADU BINARY DIRECTORY].
