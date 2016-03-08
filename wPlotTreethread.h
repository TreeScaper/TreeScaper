
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

#ifndef WPLOTTREETHREAD_H
#define WPLOTTREETHREAD_H

#include "wPlotTree.h"
#include <QThread>
#include "wstring.h"
#include <iostream>
#include "TreeOPE.h"

class PlotTreethread : public QThread
{
    Q_OBJECT
public:
    PlotTreethread(QObject *parent = 0);
    ~PlotTreethread();
    void initialization(image_parameters plotparas, const NEWICKTREE *tree, bool isrooted, bool isweighted, const LabelMap *leaveslm, String title);

protected:
    void run();

private:
    PlotTree *lpplottree;
    const NEWICKTREE *tree;
    bool isrooted;
    bool isweighted;
    const LabelMap *leaveslm;
};

#endif // WPLOTTHREAD_H
