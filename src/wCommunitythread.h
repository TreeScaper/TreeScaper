
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

#ifndef WCOMMUNITYTHREAD_H
#define WCOMMUNITYTHREAD_H

#include <QThread>
#include "wstring.h"
#include <iostream>
#include "Trees.h"

class Communitythread : public QThread
{
    Q_OBJECT
public:
    Communitythread(QObject *parent = 0);
    ~Communitythread();
    void initialization(Trees *TreesData, String in_str_matrix, int in_modelType, Array<double> in_param1, Array<double> in_param2,
                        string in_highfreq, string in_lowfreq, int in_buttonflag, bool in_isauto = true);

//    void stop();
protected:
    void run();

signals:
    void sendbuttonCommunityPlotenable(int);

private:
    Trees *lptrees;
    String str_matrix;
    int modelType;
    Array<double> param1;
    Array<double> param2;
    string highfreq;
    string lowfreq;
    int buttonflag;
    bool isauto;

//    volatile bool stopped;
};

#endif // WCOMMUNITYTHREAD_H
