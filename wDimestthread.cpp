
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

#ifndef WDIMESTTHREAD_CPP
#define WDIMESTTHREAD_CPP

#include "treescaper.h"
#include "wDimestthread.h"

Dimestthread::Dimestthread(QObject *parent)
    : QThread(parent)
{
};

Dimestthread::~Dimestthread()
{
};

void Dimestthread::initialization(String fname, String estimator, String type, int num)
{
    dimest = new DimEst;
    cout << "initialzing" << endl;
    dimest->parameters.cor_n = num;
    dimest->parameters.mle_n = num;
    dimest->parameters.nn_n = num;
    dimest->parameters.distance_file_type = 1;
    dimest->initial_DimEst(fname, estimator, type, "");
};

void Dimestthread::run()
{
    cout << "computing" << endl;
    dimest->Compute_Dim();
    cout << "outputing to files" << endl;
    dimest->output_to_files();
    delete dimest;
}

#endif // WDIMESTTHREAD_CPP
