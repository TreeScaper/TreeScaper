
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

#ifndef WPLOTTHREAD_H
#define WPLOTTHREAD_H

#include "wImage.h"
#include <QThread>
#include "wstring.h"
#include <iostream>

class Plotthread : public QThread
{
    Q_OBJECT
public:
    Plotthread(QObject *parent = 0);
    ~Plotthread();
    void initialization(String fname, image_parameters imageparas, Array<int> &selected_ids);

protected:
    void run();

private:
    String filename;
    Image *lpimage;
};

#endif // WPLOTTHREAD_H
