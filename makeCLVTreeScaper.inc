# -*- Makefile -*-
####################################################################
#  TreeScaper make include file.                                   #
#  TreeScaper, Version 1.0.0                                       #
#  June 2020	                                                   #
####################################################################
#
#  The machine (platform) identifier to append to the CLAPACK library 
#  libraries. Make sure it is consistent with PLAT in CLAPACK/make.inc
#  and the postfix of the CLAPACK libraries lapack_XXX.a and bla
#
PLAT = _LINUX
#
#  CLAPACK path. It is strongly recommended to placed the CLPACK 
#  outside the GitHub repository to avoid recompling it.
#
CLAPPATH = ../CLAPACK-3.2.1
#
#  See more detailed instruction in /CLAPACK-3.2.1/README.install
#  for using different version of CLAPACK.
