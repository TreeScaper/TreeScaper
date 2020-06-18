# -*- Makefile -*-
####################################################################
#  TreeScaper make include file.                                   #
#  TreeScaper, Version 1.0.0                                       #
#  June 2020		                                           #
####################################################################
#
#  The machine (platform) identifier to append to the CLAPACK library 
#  libraries.
#
PLAT = _MAC
#
#  CLAPACK path. It is strongly recommended to placed the CLPACK 
#  outside the GitHub repository to avoid recompling it.
#
CLAPATH = ./CLAPACK-3.2.1
#
#  See more detailed instruction in /CLAPACK-3.2.1/README.install
#  for using different version of CLAPACK.