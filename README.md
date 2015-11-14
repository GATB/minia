Type `bin/minia` for usage.

short manual: ./manual/manual.pdf

we have set up a Q&A website, see the Support section in http://minia.genouest.org

to contact an author directly: rayan.chikhi@ens-cachan.org


Compilation 
-----------

A C++11-compatible compiler is necesary. (e.g. GCC >= 4.7)

If cmake complains that pdflatex/bibtex are absent, compile with: 
    
    cmake -DSKIP_DOC=1


Disregard this
--------------

    cmake .. -DGFORGE_USER=chikhi -DMAJOR=2 -DMINOR=0 -DPATCH=1 && make -j 4
