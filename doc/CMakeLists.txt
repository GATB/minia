
################################################################################
#  DOC 
################################################################################

# http://public.kitware.com/Wiki/CMakeUserUseLATEX
INCLUDE (UseLATEX)

# We define our document. This will add a 'pdf' target in our makefile
ADD_LATEX_DOCUMENT (manual.tex)

# We tell what doc files to be installed
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf DESTINATION doc/minia/)
