
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

#ifndef WNLDRTHREAD_CPP
#define WNLDRTHREAD_CPP

#include "treescaper.h"
#include "wNLDRthread.h"

NLDRthread::NLDRthread(QObject *parent)
    : QThread(parent)
{
};

NLDRthread::~NLDRthread()
{
};

void NLDRthread::initialization(String fname, double **dist, int size, String dim, String method, String algo, String init_meth, String flag, long seed, nldr_parameters nldrparas)
{
    filename = fname;
    dimension = dim;
    Cost = method;
    Algorithm = algo;
    init = init_meth;
    fileflag = flag;
    randseed = seed;
    lpnldr = new NLDR;
    memcpy(&(lpnldr->parameters), &(paras::nldrparas), sizeof(nldr_parameters));
    lpnldr->parameters.distance_file_type = 1;

    lpnldr->init_NLDR(filename, dist, size, dimension, Cost, Algorithm, init, fileflag, randseed, "");
};

void NLDRthread::run()
{
    lpnldr->Compute_NLDR();
    lpnldr->result_analysis();
    lpnldr->output_to_files();
    delete lpnldr;

    emit sendbuttonvisible(1);
}

#endif // WNLDRTHREAD_CPP
