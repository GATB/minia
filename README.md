# Minia 

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

<!---
| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-macos-10.9.5-gcc-4.2.1/)
--->

# Before continuing..

If you are looking to do high-quality genome or metagenome assemblies, please go here: https://github.com/GATB/gatb-minia-pipeline This is a pipeline built on top of Minia that does a similar algorithm to metaSpades and MEGAHIT (multi-k assembly).

# Introduction

Minia is a short-read assembler based on a de Bruijn graph, capable of assembling a human genome on a desktop computer in a day. The output of Minia is a set of contigs. Back when it was released, Minia produced results of similar contiguity and accuracy to other de Bruijn assemblers (e.g. Velvet). Now (2015 onwards), genome assemblers have evolved and in order ot have high contiguity, see the previous section. 

# Getting the latest source code

## Instructions

It is recommended to use download the latest binary release (Linux or OSX) there: https://github.com/GATB/minia/releases

Otherwise, Minia may be compiled from sources as follows:

    # get a local copy of minia source code
    git clone --recursive https://github.com/GATB/minia.git
    
    # compile the code an run a simple test on your computer
    cd minia
    sh INSTALL

## Requirements

CMake 3.10+; see http://www.cmake.org/cmake/resources/software.html

C++11 compiler; (g++ version>=4.7 (Linux), clang version>=4.3 (Mac OSX))


# User manual	 

Type `minia` without any arguments for usage instructions.

A more complete manual is here: https://github.com/GATB/minia/raw/master/doc/manual.pdf

# What is new ? (2018)

Minia version 1 was implementing a rather unusual way to perform the assembly: traverse the graph and attempt to jump over errors and variants. This worked rather okay but not for e.g. repeated regions with many sequencing errors. Minia version 2 also followed the same philosophy, and had major improvements coming from the integration of the GATB library (mostly speed improvements) and cascading Bloom filter.  Minia version 3 uses newer techniques and has virtually nothing in common with Minia 1: there is no Bloom filter anymore (the data structure is based on unitigs produced by the BCALM software). The assembly is performed using graph simplifications that are heavily inspired by the SPAdes assembler.


# Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/
