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

#ifndef SETTING_CPP
#define SETTING_CPP

#include "treescaper.h"
#include "setting.h"
#include "ui_setting.h"
#include <QTableWidgetItem>
#include <iostream>
#include "warray.cpp"
#include <QMessageBox>
#include <QFile>

setting::setting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::setting)
{
    ui->setupUi(this);
    ui->tableWidgetplot->setColumnWidth(0, 60);
    ui->tableWidgetplot->setColumnWidth(1, 245);
    ui->tableWidgetplotindex->setColumnWidth(0, 70);
    ui->textplotstepsize->setReadOnly(true);
//    ui->pushplotapply->setEnabled(true);
    paravalues = NULL;
    indexvalues = NULL;
    loadplotparas();
}

setting::~setting()
{
    for(int i = 0; i < nparavalues; i++)
        delete paravalues[i];

    for(int i = 0; i < nindex; i++)
        delete indexvalues[i];

    delete [] paravalues;
    delete [] indexvalues;
    delete ui;
}

void setting::loadplotparas()
{
    QMessageBox msgBox;

    //================load plotparas==================
//    File fparas("plotparas");
//    if(!fparas.is_open())
//    {
//        msgBox.setText("error: can not open parameter files \"plotparas\"!");
//        msgBox.exec();
//        return;
//    }

//    Array<String> entries;
//    String element;
//    fparas >> element;
//    int countelements = 1;      // Skip newline at the end of the file
//	while(! fparas.is_end())
//    {
//        entries.add(element);
//        if(countelements > 1024)
//            break;
//        else
//            fparas >> element;
//        countelements += 1;
//    }
    ifstream plotfparas;
    plotfparas.open("plotparas");
    if(!plotfparas.is_open())
    {
        msgBox.setText("Error: Cannot open parameter file \"plotparas\"!");
        msgBox.exec();
        return;
    }
    string linefromfile;
    string elementfromline;
    Array<String> entries;
    while(getline(plotfparas, linefromfile))
    {
        std::istringstream iss(linefromfile);
        while(iss >> elementfromline)
        {
            String element(elementfromline.c_str());
            entries.add(element);
        }

    }

    if(entries[1] == (String) "Points")
        ui->comboBoxplottype->setCurrentIndex(0);
    else
        ui->comboBoxplottype->setCurrentIndex(1);

    if(entries[3] == (String) "1")
        ui->checkBoxplotmovie->setChecked(true);
    else
        ui->checkBoxplotmovie->setChecked(false);

    QString qtstr(entries[5]);
    ui->textplotnumber->setPlainText(qtstr);
    ui->pushplotnum->click();

    if(entries[7] == (String) "1")
        ui->checkBoxplotsort->setChecked(true);
    else
        ui->checkBoxplotsort->setChecked(false);

    if(entries[9] == (String) "1")
    {
        ui->radioplotgradual->setChecked(true);
        ui->radioplotmost->setChecked(false);
    } else
    {
        ui->radioplotgradual->setChecked(false);
        ui->radioplotmost->setChecked(true);
    }

    nparavalues = 8;
    paravalues = new QTableWidgetItem *[nparavalues];
    for(int i = 0; i < nparavalues; i++)
    {
        paravalues[i] = new QTableWidgetItem((char *) entries[i * 2 + 11]);
        ui->tableWidgetplot->setItem(i, 0, paravalues[i]);
    }

    nindex = qtstr.toInt();
    indexvalues = new QTableWidgetItem *[nindex];
    for(int i = 0; i < nindex; i++)
        indexvalues[i] = NULL;
    for(int i = 0; i < nindex; i++)
    {
        if(i * 2 + 27 >= entries.get_length())
        {
            msgBox.setText("Warning: The index is wrong! Please fix it.");
            msgBox.exec();
            break;
        }
        indexvalues[i] = new QTableWidgetItem((char *) entries[i * 2 + 27]);
        ui->tableWidgetplotindex->setItem(i, 0, indexvalues[i]);
    }
}

void setting::on_pushplotapply_clicked()
{
    File fparas("plotparas");
    QMessageBox msgBox;
    if(!fparas.is_open())
    {
        msgBox.setText("Error: Cannot open parameter files \"plotparas\"!");
        msgBox.exec();
        return;
    }

    fparas.clean();
    fparas.seek(0);

    QString qtplottype = ui->comboBoxplottype->currentText();
    string stdplottype =qtplottype.toStdString();
    String plottype(stdplottype.c_str());
    paras::plot_type = plottype;
    fparas << "plot_type\t" << plottype << endl;

    paras::plot_makemovie = ui->checkBoxplotmovie->isChecked();
    if(paras::plot_makemovie)
        fparas << "plot_makemovie\t" << 1 << endl;
    else
        fparas << "plot_makemovie\t" << 0 << endl;

    QString qtplotnumber = ui->textplotnumber->toPlainText();
    paras::imageparas.cluster_num = qtplotnumber.toInt();
    fparas << "cluster_num\t" << paras::imageparas.cluster_num << endl;
    if(! paras::imageparas.cluster_num)
    {
        msgBox.setText("warning: invalid cluster number, it must be greater than 0. It may cause errors. Please fix it and click apply again.");
        msgBox.exec();
        return;
    }

    paras::is_sort_volume = ui->checkBoxplotsort->isChecked();
    if(paras::is_sort_volume)
        fparas << "is_sort_volume\t" << 1 << endl;
    else
        fparas << "is_sort_volume\t" << 0 << endl;

    if(ui->radioplotgradual->isChecked())
        paras::imageparas.color_type = 1;
    else
        paras::imageparas.color_type = 0;
    fparas << "color_type\t" << paras::imageparas.color_type << endl;

    QString stritem;
    QTableWidgetItem *item;
    item = ui->tableWidgetplot->item(0, 0);
    stritem = item->text();
    paras::imageparas.outlier_e = stritem.toDouble();
    fparas << "outlier_e\t" << paras::imageparas.outlier_e << endl;

    item = ui->tableWidgetplot->item(1, 0);
    stritem = item->text();
    paras::imageparas.outlierpercent = stritem.toDouble();
    fparas << "outlierpercent\t" << paras::imageparas.outlierpercent << endl;

    item = ui->tableWidgetplot->item(2, 0);
    stritem = item->text();
    paras::imageparas.sep_factor = stritem.toDouble();
    fparas << "sep_factor\t" << paras::imageparas.sep_factor << endl;

    item = ui->tableWidgetplot->item(3, 0);
    stritem = item->text();
    paras::imageparas.point_size = stritem.toDouble();
    fparas << "point_size\t" << paras::imageparas.point_size << endl;

    item = ui->tableWidgetplot->item(4, 0);
    stritem = item->text();
    paras::imageparas.moviesizex = stritem.toInt();
    fparas << "moviesizex\t" << paras::imageparas.moviesizex << endl;

    item = ui->tableWidgetplot->item(5, 0);
    stritem = item->text();
    paras::imageparas.moviesizey = stritem.toInt();
    fparas << "moviesizey\t" << paras::imageparas.moviesizey << endl;

    item = ui->tableWidgetplot->item(6, 0);
    stritem = item->text();
    paras::imageparas.zoom = stritem.toDouble();
    fparas << "zoom\t" << paras::imageparas.zoom << endl;

    item = ui->tableWidgetplot->item(7, 0);
    stritem = item->text();
    paras::imageparas.frame_num = stritem.toInt();
    fparas << "frame_num\t" << paras::imageparas.frame_num << endl;

    int row = ui->tableWidgetplotindex->rowCount();

    int max_index = 0;
    QStringList indices, range;
    QString qtelement;
    for(int i = 0; i < row; i++)
    {
        item = ui->tableWidgetplotindex->item(i, 0);
        if(item == NULL)
        {
            msgBox.setText("Warning: Index is not allow to be empty! Please fill them in.");
            msgBox.exec();
            return;
        }

        qtelement = (QString) item->text();
        qtelement = qtelement.remove(QChar(' '));
        indices = qtelement.split(",", QString::SkipEmptyParts);
        for(int j = 0; j < indices.size(); j++)
        {
            range = indices[j].split("-", QString::SkipEmptyParts);
            for(int k = 0; k < range.size(); k++)
                if(range[k].toInt() > max_index)
                    max_index = range[k].toInt();
        }
    }

    delete [] paras::imageparas.cluster_index;
    paras::imageparas.cluster_index = new int [max_index + 1];
    paras::imageparas.index_size = max_index + 1;

    for(int i = 0; i < row; i++)
    {
        item = ui->tableWidgetplotindex->item(i, 0);
        qtelement = (QString) item->text();
        qtelement = qtelement.remove(QChar(' '));
        fparas << i << "\t" << qtelement.toStdString() << endl;

        indices = qtelement.split(",", QString::SkipEmptyParts);
        for(int j = 0; j < indices.size(); j++)
        {
            range = indices[j].split("-", QString::SkipEmptyParts);
            if(range.size() == 1)
                paras::imageparas.cluster_index[range[0].toInt()] = i;
            else
            {
                for(int k = 0; k < range.size(); k = k + 2)
                {
                    for(int h = range[k].toInt(); h <= range[k + 1].toInt(); h++)
                        paras::imageparas.cluster_index[h] = i;
                }
            }
        }
    }

//    ui->pushplotapply->setEnabled(false);
}

void setting::on_textplotnumber_textChanged()
{
}

void setting::on_pushplotnum_clicked()
{
    int row = ui->tableWidgetplotindex->rowCount();
    QString qtnum = ui->textplotnumber->toPlainText();
    int num = qtnum.toInt();

    if(num > row)
    {
        for(int i = 0; i < num - row; i++)
        {
            ui->tableWidgetplotindex->insertRow(i + row);
        }
    }
    if(num < row)
    {
        for(int i = 0; i < row - num; i++)
            ui->tableWidgetplotindex->removeRow(row - i - 1);
        for(int i = num; i < nindex; i++)
            indexvalues[i] = NULL;
    }
}

void setting::on_textplotstepsize_textChanged()
{
    if(ui->checkBoxplotstepsize->isChecked())
    {
        QString qtstepsize = ui->textplotstepsize->toPlainText();
        int stepsize = qtstepsize.toInt();
        int row = ui->tableWidgetplotindex->rowCount();
        QTableWidgetItem *element;
        QString qtstart, qtend, range;

        for(int i = 0; i < row; i++)
        {
            qtstart = QString::number(i * stepsize);
            qtend = QString::number((i + 1) * stepsize - 1);
            range = qtstart;
            range += "-";
            range += qtend;

            element = ui->tableWidgetplotindex->item(i, 0);
            if(element == NULL)
            {
                element = new  QTableWidgetItem(range);
                ui->tableWidgetplotindex->setItem(i, 0, element);
            } else
                element->setText(range);
        }
    }
}

void setting::on_checkBoxplotstepsize_stateChanged(int )
{
    if(ui->checkBoxplotstepsize->isChecked())
    {
        ui->textplotstepsize->setReadOnly(false);
    } else
    {
        ui->textplotstepsize->setPlainText("");
        ui->textplotstepsize->setReadOnly(true);
    }
}

void setting::on_pushplotcancel_clicked()
{
    loadplotparas();
//    ui->pushplotapply->setEnabled(true);
}

void setting::on_pushButton_clicked()
{
    QMessageBox MSmsgBox;
    MSmsgBox.setText("Restore all original setting for plot parameters?");
    MSmsgBox.addButton(QMessageBox::Yes);
    MSmsgBox.addButton(QMessageBox::No);
    int ret = MSmsgBox.exec();
    if(ret == QMessageBox::No)
        return;

    QFile defaultplotparas("plotparas_default");
    QMessageBox msgBox;

    if(!defaultplotparas.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        msgBox.setText("Error: Cannot open default parameter file \"plotparas_default\"!");
        msgBox.exec();
        return;
    }

    defaultplotparas.remove("plotparas");
    defaultplotparas.copy("plotparas");

    ui->pushplotcancel->click();
}

void setting::on_pushplotclose_clicked()
{
    QDialog::accept();
}

#endif
