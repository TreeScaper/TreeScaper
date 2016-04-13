
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

#ifndef WDIMESTTHREAD_H
#define WDIMESTTHREAD_H

#include <QThread>
#include "wstring.h"
#include "wDimEst.h"
#include <iostream>

class Dimestthread : public QThread
{
    Q_OBJECT
public:
    Dimestthread(QObject *parent = 0);
    ~Dimestthread();
    void initialization(String fname, double **dist, int sizeinput, int dim, String estimator, String type, int num);

protected:
    void run();

private:
    DimEst *dimest;
};

#endif // WDIMESTTHREAD_H
