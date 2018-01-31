# Minia 

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

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

#Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/
