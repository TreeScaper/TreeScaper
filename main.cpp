
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

<<<<<<< HEAD
//#if VTK_MAJOR_VERSION <= 6 && VTK_MINOR_VERSION >= 3
//// for 6.0.0 VTK
=======
#if VTK_MAJOR_VERSION <= 6 && VTK_MINOR_VERSION >= 3
// for 6.0.0 VTK
>>>>>>> refs/remotes/origin/developer
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
//// for 6.3.0 VTK
//#else
//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//#endif

#ifdef _WIN32
//#include <windef.h>
#include <windows.h>
//#include <stddef.h>
#endif

#include "wImage.h"
#include "randgen.h"
#include <iostream>
#include <fstream>
#include "warray.h"
#include "wstring.h"
#include "wmatrix.cpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include "wfile.h"
#include "wNLDR.h"
#include "wDimEst.h"
//#include <QtGui/QApplication>
#undef max
#undef min
#include <QApplication>
#include <QDebug>
#include "treescaper.h"
#include "woutput.h"

int main(int argc, char *argv[])
{
//    std::cout << "h0" << std::endl;//--
//    return 0;//---
    QApplication a(argc, argv);
    TreeScaper w;

    w.show();

    return a.exec();
}
