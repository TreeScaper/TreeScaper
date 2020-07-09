
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

#ifndef WNLDRTHREAD_H
#define WNLDRTHREAD_H

#include <QThread>
#include "wstring.h"
#include "wNLDR.h"
#include <iostream>

class NLDRthread : public QThread
{
    Q_OBJECT
public:
    NLDRthread(QObject *parent = 0);
    ~NLDRthread();
    void initialization(String fname, double **dist, int size, String dim, String method, String algo, String init_meth, String flag, long seed, nldr_parameters nldrparas);

protected:
    void run();

signals:
    void sendbuttonvisible(int);

private:
    String filename;
    String dimension;
    String Cost;
    String Algorithm;
    String init;
    String fileflag;
    long randseed;
    NLDR *lpnldr;
};

#endif // WNLDRTHREAD_H
