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

#ifndef NLDRSETTING_H
#define NLDRSETTING_H

#include <QDialog>
#include <QTableWidgetItem>

namespace Ui {
class NLDRsetting;
}

class NLDRsetting : public QDialog
{
    Q_OBJECT

public:
    explicit NLDRsetting(QWidget *parent = 0);
    ~NLDRsetting();

private:
    void loadnldrparas();

    Ui::NLDRsetting *ui;
    QTableWidgetItem **nldrvalues;
    int nnldrvalues;

private slots:
    void on_pushnldrreset_clicked();
    void on_pushnldrcancel_clicked();
    void on_pushnldrapply_clicked();
    void on_radionldrKruButton_clicked();
    void on_radionldrtopButton_clicked();
    void on_radionldrNorButton_clicked();
    void on_radionldrNLMButton_clicked();
    void on_radionldrCCAButton_clicked();
//    void on_pushnldrclose_clicked();
};

#endif // NLDRSETTING_H
