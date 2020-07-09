
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

#ifndef WPLOTDIMTHREAD_CPP
#define WPLOTDIMTHREAD_CPP

#include "wdef.h"
#include "wImage.h"
#include "wPlotdimthread.h"
#include "treescaper.h"

Plotdimthread::Plotdimthread(QObject *parent)
    : QThread(parent)
{
};

Plotdimthread::~Plotdimthread()
{
};

void Plotdimthread::initialization(String fname)
{
    lpimage = new Image;
    lpimage->parameters.cluster_num = 1;
    File file(fname);
    int size = file.lines();
    int *indexs = new int[size];
    for(int i = 0; i < size; i++)
        indexs[i] = 0;
    lpimage->parameters.cluster_index = indexs;
    lpimage->parameters.sep_factor = 0;
    lpimage->parameters.index_size = size;
    lpimage->Initialize_Image(fname);
#ifdef _MAC
    lpimage->Load_Points();
    lpimage->Create_LinesActors();
    lpimage->Plot_Lines();
    this->exec();
    delete [] lpimage->parameters.cluster_index;
    delete lpimage;
#endif
};

void Plotdimthread::run()
{
#ifdef _LINUX
    lpimage->Load_Points();
    lpimage->Create_LinesActors();
    lpimage->Plot_Lines();
    delete [] lpimage->parameters.cluster_index;
    delete lpimage;
#endif

#ifdef _WINDOWS
    lpimage->Load_Points();
    lpimage->Create_LinesActors();
    lpimage->Plot_Lines();
    delete [] lpimage->parameters.cluster_index;
    delete lpimage;
#endif
}

#endif // WPLOTDIMTHREAD_CPP
