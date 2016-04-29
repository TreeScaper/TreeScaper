

//##########################################################################
//# This software is part of the Treescaper i
//# -- Version 0.1
//# Copyright (C) 2010 Wen Huang
//#
//# This program is free software; you can redistribute it and/or
//# modify it under the terms of the GNU General Public License
//# as published by the Free Software Foundation; either version 2
//# of the License, or (at your option) any later version.
//#
//# This program is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//# http://www.gnu.org/copyleft/gpl.html
//##########################################################################

This readme file describes how to compile the code on Windows, Linux(Ubontu) and 
MAC system.

Compiling:(Windows, Linux, Mac)
First, we need to change the file 'wdef.h'. If the platform is window, Linux, MAC respectively,
then _WINDOWS, _LINUX and _MAC need to be defined in 'wdef.h'.

Windows: 
1, download qt-opensource-windows-x86-mingw492-5.5.1.exe and install it. Put the bin 
directory of mingw in the system variable environment and restart computer.
2, download vtk-6.3.0.zip from http://www.vtk.org/download/
3, download cmake-3.5.0-rc3-win32-x86.zip from https://cmake.org/download/
4, download  mingw-get-setup.exe from http://www.mingw.org/wiki/Getting_Started and only install msys.
Note that do not install its mingw and use the mingw from QT.
 5, unzip vtk-6.3.0.zip and create ''bin-6.3.0'' and ''bin'' directory.
 6, open cmake-gui and choose vtk-6.3.0 as source code and choose bin-6.3.0 as build directory, then click configure.
 7, choose mingw makefile as generator. Do NOT uncheck BUILD_SHARED_LIBS. Modify CMAKE_INSTALL_PREFIX to be the directory of ''bin''.
 8, click generate. Open cmd and go to the directory of ''bin-6.3.0'' and use command mingw32-make to compile vtk.
  Note that you have have ''access denied'' error, which is caused by antivirus software. You can ignore it and rerun mingw32-make.
 9, Run mingw32-make install to install the library. All the head and library files will be in the ''bin'' directory, i.e., bin/include and
  bin/bin.
 10, download clapack-3.2.1-CMAKE.tgz from http://www.netlib.org/clapack/ and use the cmake and same way to generate library files.
 11, Open Windows_TreeScaper.pro using Qt. open the file wdef.h, make sure _WINDOWS is in the file.
 12, Copy all necessary vtk libraries and QT libraries to the directory of TreeScaper.exe

Linux:
Similar to MAC version
MAC:
1, download cmake, cmake-2.8.12-Darwin-universal.dmg, from: http://www.cmake.org/cmake/resources/software.html
2, download vtk from http://www.vtk.org/VTK/resources/software.html#latest2.
3, What I download are vtk-6.0.0.tar.gz and vtkdata-6.0.0.tar.gz.
4, untar vtk-6.0.0.tar.gz and vtkdata-6.0.0.tar.gz in a same directory. create a new fold 'bin' in the directory.
5, open cmake, choose vtk-6.0.0 as source code and choose bin as build directory. then click configure.
6, choose Unix-makefiles as generator. Uncheck BUILD_SHARED_LIBS and delete the value '-fobjc-gc' in VTK_REQUIRED_OBJCXX_FLAGS
7, now, we can click generate. After that, go to the bin directory that just is generated, use 'make' to build libs of vtk.
Modify code
8, open the file wdef.h, make sure _MAC is in the file.
9, open the file Mac_TreeScaper.pro. then change the directory for each library and head file.
package
10, use command "macdeployqt Treescaper.app/ -verbose=1 -dmg" to package. The macdeployqt is in the directory of QTSDK somewhere.
11, use "otool -L TreeScaper.app/Contents/MacOS/TreeScaper" to see the dependencies of the program.
​



