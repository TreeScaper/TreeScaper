
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

#ifndef SETTING_H
#define SETTING_H

#include <QDialog>
#include <QTableWidgetItem>

namespace Ui {
    class setting;
}

class setting : public QDialog
{
    Q_OBJECT

public:
    explicit setting(QWidget *parent = 0);
    ~setting();

private:
    void loadplotparas();

    Ui::setting *ui;
    QTableWidgetItem **paravalues;
    QTableWidgetItem **indexvalues;
    int nparavalues;
    int nindex;

private slots:
    void on_pushButton_clicked();
    void on_pushplotcancel_clicked();
    void on_checkBoxplotstepsize_stateChanged(int );
    void on_textplotstepsize_textChanged();
    void on_pushplotnum_clicked();
    void on_textplotnumber_textChanged();
    void on_pushplotapply_clicked();
    void on_pushplotclose_clicked();
};

#endif // SETTING_H
