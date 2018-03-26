# Minia 

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Minia/job/tool-minia-build-macos-10.9.5-gcc-4.2.1/)


# What is Minia ?

Minia is a short-read assembler based on a de Bruijn graph, capable of assembling a human genome on a desktop computer in a day. The output of Minia is a set of contigs. Minia produces results of similar contiguity and accuracy to other de Bruijn assemblers (e.g. Velvet).

# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

C++11 compiler; (g++ version>=4.7 (Linux), clang version>=4.3 (Mac OSX))

## Instructions

    # get a local copy of minia source code
    git clone --recursive https://github.com/GATB/minia.git
    
    # compile the code an run a simple test on your computer
    cd minia
    sh INSTALL

# User manual	 

Type `minia` without any arguments for usage instructions.

A more complete manual:

    cd doc 
    pdflatex manual.tex

If you cannot compile it: http://minia.genouest.org/files/minia.pdf

# What is new ? (2018)

Minia version 1 was implementing a rather unusual way to perform the assembly: traverse the graph and attempt to jump over errors and variants. This worked rather okay but not for e.g. repeated regions with many sequencing errors. Minia version 2 also followed the same philosophy, and had major improvements coming from the integration of the GATB library (mostly speed improvements) and cascading Bloom filter.  Minia version 3 uses newer techniques and has virtually nothing in common with Minia 1: there is no Bloom filter anymore (the data structure is based on unitigs produced by the BCALM software). The assembly is performed using graph simplifications that are heavily inspired by the SPAdes assembler.


# Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/
