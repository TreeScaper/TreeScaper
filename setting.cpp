
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
    ui->tablewidgetnldr->resizeColumnsToContents();
    ui->textplotstepsize->setReadOnly(true);
    paravalues = NULL;
    indexvalues = NULL;
    nldrvalues = NULL;
    loadnldrparas();
    loadplotparas();
}

setting::~setting()
{
    for(int i = 0; i < nparavalues; i++)
        delete paravalues[i];

    for(int i = 0; i < nindex; i++)
        delete indexvalues[i];

    for(int i = 0; i < nnldrvalues; i++)
        delete nldrvalues[i];

    delete [] paravalues;
    delete [] indexvalues;
    delete [] nldrvalues;
    delete ui;
}

void setting::loadnldrparas()
{
    QMessageBox msgBox;

//======================load nldrparas===================
    File nldrfparas("nldrparas");
    if(!nldrfparas.is_open())
    {
        msgBox.setText("error: can not open parameter files \"nldrparas\"!");
        msgBox.exec();
        return;
    }

    Array<String> nldrentries;
    String nldrelement;
    nldrfparas >> nldrelement;
	while(! nldrfparas.is_end())
    {
        nldrentries.add(nldrelement);
        nldrfparas >> nldrelement;
    }

    nnldrvalues = 51;
    nldrvalues = new QTableWidgetItem *[nnldrvalues];
    for(int i = 0; i < nnldrvalues; i++)
    {
        nldrvalues[i] = new QTableWidgetItem((char *) nldrentries[2 * i + 1]);
        ui->tablewidgetnldr->setItem(i, 1, nldrvalues[i]);
    }
}

void setting::loadplotparas()
{
    QMessageBox msgBox;

    //================load plotparas==================
    File fparas("plotparas");
    if(!fparas.is_open())
    {
        msgBox.setText("error: can not open parameter files \"plotparas\"!");
        msgBox.exec();
        return;
    }

    Array<String> entries;
    String element;
    fparas >> element;
	while(! fparas.is_end())
    {
        entries.add(element);
        fparas >> element;
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
            msgBox.setText("warning: the index is wrong. Please fix it.");
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
        msgBox.setText("error: can not open parameter files \"plotparas\"!");
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
            msgBox.setText("warning : index is not allow to be empty. Please fill them in.");
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
}

void setting::on_textplotnumber_textChanged()
{
}

void setting::on_pushnldrapply_clicked()
{
    File fparas("nldrparas");
    QMessageBox msgBox;
    if(!fparas.is_open())
    {
        msgBox.setText("error: can not open parameter files \"nldrparas\"!");
        msgBox.exec();
        return;
    }

    fparas.clean();
    fparas.seek(0);

    QString stritem;
    QTableWidgetItem *item;

    item = ui->tablewidgetnldr->item(0, 1);
    stritem = item->text();
    paras::nldrparas.random_start = stritem.toDouble();
    fparas << "random_start\t" << paras::nldrparas.random_start << endl;

    item = ui->tablewidgetnldr->item(1, 1);
    stritem = item->text();
    paras::nldrparas.random_end = stritem.toDouble();
    fparas << "random_end\t" << paras::nldrparas.random_end << endl;

    item = ui->tablewidgetnldr->item(2, 1);
    stritem = item->text();
    paras::nldrparas.length_tru = stritem.toInt();
    fparas << "length_tru\t" << paras::nldrparas.length_tru << endl;

    item = ui->tablewidgetnldr->item(3, 1);
    stritem = item->text();
    paras::nldrparas.interval_tru = stritem.toInt();
    fparas << "interval_tru\t" << paras::nldrparas.interval_tru << endl;

    item = ui->tablewidgetnldr->item(4, 1);
    stritem = item->text();
    paras::nldrparas.length_con = stritem.toInt();
    fparas << "length_con\t" << paras::nldrparas.length_con << endl;

    item = ui->tablewidgetnldr->item(5, 1);
    stritem = item->text();
    paras::nldrparas.iterval_con = stritem.toInt();
    fparas << "interval_con\t" << paras::nldrparas.iterval_con << endl;

    item = ui->tablewidgetnldr->item(6, 1);
    stritem = item->text();
    paras::nldrparas.KRU_LIN_e = stritem.toDouble();
    fparas << "KRU_LIN_e\t" << paras::nldrparas.KRU_LIN_e << endl;

    item = ui->tablewidgetnldr->item(7, 1);
    stritem = item->text();
    paras::nldrparas.KRU_LIN_max_iter = stritem.toDouble();
    fparas << "KRU_LIN_max_iter\t" << paras::nldrparas.KRU_LIN_max_iter << endl;

    item = ui->tablewidgetnldr->item(8, 1);
    stritem = item->text();
    paras::nldrparas.KRU_LIN_run_time = stritem.toInt();
    fparas << "KRU_LIN_run_time\t" << paras::nldrparas.KRU_LIN_run_time << endl;

    item = ui->tablewidgetnldr->item(9, 1);
    stritem = item->text();
    paras::nldrparas.KRU_MAJ_e = stritem.toDouble();
    fparas << "KRU_MAJ_e\t" << paras::nldrparas.KRU_MAJ_e << endl;

    item = ui->tablewidgetnldr->item(10, 1);
    stritem = item->text();
    paras::nldrparas.KRU_MAJ_max_iter = stritem.toInt();
    fparas << "KRU_MAJ_max_iter\t" << paras::nldrparas.KRU_MAJ_max_iter << endl;

    item = ui->tablewidgetnldr->item(11, 1);
    stritem = item->text();
    paras::nldrparas.KRU_MAJ_run_time = stritem.toInt();
    fparas << "KRU_MAJ_run_time\t" << paras::nldrparas.KRU_MAJ_run_time << endl;

    item = ui->tablewidgetnldr->item(12, 1);
    stritem = item->text();
    paras::nldrparas.KRU_GAU_e = stritem.toDouble();
    fparas << "KRU_GAU_e\t" << paras::nldrparas.KRU_GAU_e << endl;

    item = ui->tablewidgetnldr->item(13, 1);
    stritem = item->text();
    paras::nldrparas.KRU_GAU_max_iter = stritem.toInt();
    fparas << "KRU_GAU_max_iter\t" << paras::nldrparas.KRU_GAU_max_iter << endl;

    item = ui->tablewidgetnldr->item(14, 1);
    stritem = item->text();
    paras::nldrparas.KRU_GAU_run_time = stritem.toInt();
    fparas << "KRU_GAU_run_time\t" << paras::nldrparas.KRU_GAU_run_time << endl;

    item = ui->tablewidgetnldr->item(15, 1);
    stritem = item->text();
    paras::nldrparas.KRU_STO_epochs = stritem.toInt();
    fparas << "KRU_STO_epochs\t" << paras::nldrparas.KRU_STO_epochs << endl;

    item = ui->tablewidgetnldr->item(16, 1);
    stritem = item->text();
    paras::nldrparas.KRU_STO_alpha0 = stritem.toDouble();
    fparas << "KRU_STO_alpha0\t" << paras::nldrparas.KRU_STO_alpha0 << endl;

    item = ui->tablewidgetnldr->item(17, 1);
    stritem = item->text();
    paras::nldrparas.KRU_STO_alphan = stritem.toDouble();
    fparas << "KRU_STO_alphan\t" << paras::nldrparas.KRU_STO_alphan << endl;

    item = ui->tablewidgetnldr->item(18, 1);
    stritem = item->text();
    paras::nldrparas.NOR_MAJ_e = stritem.toDouble();
    fparas << "NOR_MAJ_e\t" << paras::nldrparas.NOR_MAJ_e << endl;

    item = ui->tablewidgetnldr->item(19, 1);
    stritem = item->text();
    paras::nldrparas.NOR_MAJ_max_iter = stritem.toInt();
    fparas << "NOR_MAJ_max_iter\t" << paras::nldrparas.NOR_MAJ_max_iter << endl;

    item = ui->tablewidgetnldr->item(20, 1);
    stritem = item->text();
    paras::nldrparas.NOR_MAJ_run_time = stritem.toInt();
    fparas << "NOR_MAJ_run_time\t" << paras::nldrparas.NOR_MAJ_run_time << endl;

    item = ui->tablewidgetnldr->item(21, 1);
    stritem = item->text();
    paras::nldrparas.NOR_GAU_e = stritem.toDouble();
    fparas << "NOR_GAU_e\t" << paras::nldrparas.NOR_GAU_e << endl;

    item = ui->tablewidgetnldr->item(22, 1);
    stritem = item->text();
    paras::nldrparas.NOR_GAU_max_iter = stritem.toInt();
    fparas << "NOR_GAU_max_iter\t" << paras::nldrparas.NOR_GAU_max_iter << endl;

    item = ui->tablewidgetnldr->item(23, 1);
    stritem = item->text();
    paras::nldrparas.NOR_GAU_run_time = stritem.toInt();
    fparas << "NOR_GAU_run_time\t" << paras::nldrparas.NOR_GAU_run_time << endl;

    item = ui->tablewidgetnldr->item(24, 1);
    stritem = item->text();
    paras::nldrparas.NOR_STO_epochs = stritem.toInt();
    fparas << "NOR_STO_epochs\t" << paras::nldrparas.NOR_STO_epochs << endl;

    item = ui->tablewidgetnldr->item(25, 1);
    stritem = item->text();
    paras::nldrparas.NOR_STO_alpha0 = stritem.toDouble();
    fparas << "NOR_STO_alpha0\t" << paras::nldrparas.NOR_STO_alpha0 << endl;

    item = ui->tablewidgetnldr->item(26, 1);
    stritem = item->text();
    paras::nldrparas.NOR_STO_alphan = stritem.toDouble();
    fparas << "NOR_STO_alphan\t" << paras::nldrparas.NOR_STO_alphan << endl;

    item = ui->tablewidgetnldr->item(27, 1);
    stritem = item->text();
    paras::nldrparas.NLM_MAJ_e = stritem.toDouble();
    fparas << "NLM_MAJ_e\t" << paras::nldrparas.NLM_MAJ_e << endl;

    item = ui->tablewidgetnldr->item(28, 1);
    stritem = item->text();
    paras::nldrparas.NLM_MAJ_max_iter = stritem.toInt();
    fparas << "NLM_MAJ_max_iter\t" << paras::nldrparas.NLM_MAJ_max_iter << endl;

    item = ui->tablewidgetnldr->item(29, 1);
    stritem = item->text();
    paras::nldrparas.NLM_MAJ_run_time = stritem.toInt();
    fparas << "NLM_MAJ_run_time\t" << paras::nldrparas.NLM_MAJ_run_time << endl;

    item = ui->tablewidgetnldr->item(30, 1);
    stritem = item->text();
    paras::nldrparas.NLM_GAU_e = stritem.toDouble();
    fparas << "NLM_GAU_e\t" << paras::nldrparas.NLM_GAU_e << endl;

    item = ui->tablewidgetnldr->item(31, 1);
    stritem = item->text();
    paras::nldrparas.NLM_GAU_max_iter = stritem.toInt();
    fparas << "NLM_GAU_max_iter\t" << paras::nldrparas.NLM_GAU_max_iter << endl;

    item = ui->tablewidgetnldr->item(32, 1);
    stritem = item->text();
    paras::nldrparas.NLM_GAU_run_time = stritem.toInt();
    fparas << "NLM_GAU_run_time\t" << paras::nldrparas.NLM_GAU_run_time << endl;

    item = ui->tablewidgetnldr->item(33, 1);
    stritem = item->text();
    paras::nldrparas.NLM_STO_epochs = stritem.toInt();
    fparas << "NLM_STO_epochs\t" << paras::nldrparas.NLM_STO_epochs << endl;

    item = ui->tablewidgetnldr->item(34, 1);
    stritem = item->text();
    paras::nldrparas.NLM_STO_alpha0 = stritem.toDouble();
    fparas << "NLM_STO_alpha0\t" << paras::nldrparas.NLM_STO_alpha0 << endl;

    item = ui->tablewidgetnldr->item(35, 1);
    stritem = item->text();
    paras::nldrparas.NLM_STO_alphan = stritem.toDouble();
    fparas << "NLM_STO_alphan\t" << paras::nldrparas.NLM_STO_alphan << endl;

    item = ui->tablewidgetnldr->item(36, 1);
    stritem = item->text();
    paras::nldrparas.CCA_MAJ_e = stritem.toDouble();
    fparas << "CCA_MAJ_e\t" << paras::nldrparas.CCA_MAJ_e << endl;

    item = ui->tablewidgetnldr->item(37, 1);
    stritem = item->text();
    paras::nldrparas.CCA_MAJ_max_iter = stritem.toInt();
    fparas << "CCA_MAJ_max_iter\t" << paras::nldrparas.CCA_MAJ_max_iter << endl;

    item = ui->tablewidgetnldr->item(38, 1);
    stritem = item->text();
    paras::nldrparas.CCA_MAJ_run_time = stritem.toInt();
    fparas << "CCA_MAJ_run_time\t" << paras::nldrparas.CCA_MAJ_run_time << endl;

    item = ui->tablewidgetnldr->item(39, 1);
    stritem = item->text();
    paras::nldrparas.CCA_MAJ_lambda0 = stritem.toDouble();
    fparas << "CCA_MAJ_lambda0\t" << paras::nldrparas.CCA_MAJ_lambda0 << endl;

    item = ui->tablewidgetnldr->item(40, 1);
    stritem = item->text();
    paras::nldrparas.CCA_MAJ_lambdan = stritem.toDouble();
    fparas << "CCA_MAJ_lambdan\t" << paras::nldrparas.CCA_MAJ_lambdan << endl;

    item = ui->tablewidgetnldr->item(41, 1);
    stritem = item->text();
    paras::nldrparas.CCA_GAU_e = stritem.toDouble();
    fparas << "CCA_GAU_e\t" << paras::nldrparas.CCA_GAU_e << endl;

    item = ui->tablewidgetnldr->item(42, 1);
    stritem = item->text();
    paras::nldrparas.CCA_GAU_max_iter = stritem.toInt();
    fparas << "CCA_GAU_max_iter\t" << paras::nldrparas.CCA_GAU_max_iter << endl;

    item = ui->tablewidgetnldr->item(43, 1);
    stritem = item->text();
    paras::nldrparas.CCA_GAU_run_time = stritem.toInt();
    fparas << "CCA_GAU_run_time\t" << paras::nldrparas.CCA_GAU_run_time << endl;

    item = ui->tablewidgetnldr->item(44, 1);
    stritem = item->text();
    paras::nldrparas.CCA_GAU_lambda0 = stritem.toDouble();
    fparas << "CCA_GAU_lambda0\t" << paras::nldrparas.CCA_GAU_lambda0 << endl;

    item = ui->tablewidgetnldr->item(45, 1);
    stritem = item->text();
    paras::nldrparas.CCA_GAU_lambdan = stritem.toDouble();
    fparas << "CCA_GAU_lambdan\t" << paras::nldrparas.CCA_GAU_lambdan << endl;

    item = ui->tablewidgetnldr->item(46, 1);
    stritem = item->text();
    paras::nldrparas.CCA_STO_epochs = stritem.toInt();
    fparas << "CCA_STO_epochs\t" << paras::nldrparas.CCA_STO_epochs << endl;

    item = ui->tablewidgetnldr->item(47, 1);
    stritem = item->text();
    paras::nldrparas.CCA_STO_lambda0 = stritem.toDouble();
    fparas << "CCA_STO_lambda0\t" << paras::nldrparas.CCA_STO_lambda0 << endl;

    item = ui->tablewidgetnldr->item(48, 1);
    stritem = item->text();
    paras::nldrparas.CCA_STO_lambdan = stritem.toDouble();
    fparas << "CCA_STO_lambdan\t" << paras::nldrparas.CCA_STO_lambdan << endl;

    item = ui->tablewidgetnldr->item(49, 1);
    stritem = item->text();
    paras::nldrparas.CCA_STO_alpha0 = stritem.toDouble();
    fparas << "CCA_STO_alpha0\t" << paras::nldrparas.CCA_STO_alpha0 << endl;

    item = ui->tablewidgetnldr->item(50, 1);
    stritem = item->text();
    paras::nldrparas.CCA_STO_alphan = stritem.toDouble();
    fparas << "CCA_STO_alphan\t" << paras::nldrparas.CCA_STO_alphan << endl;
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

void setting::on_pushnldrcancel_clicked()
{
    loadnldrparas();
}

void setting::on_pushplotcancel_clicked()
{
    loadplotparas();
}

void setting::on_pushnldrreset_clicked()
{
    QMessageBox MSmsgBox;
    MSmsgBox.setText("Restore all original setting for NLDR parameters?");
    MSmsgBox.addButton(QMessageBox::Yes);
    MSmsgBox.addButton(QMessageBox::No);
    int ret = MSmsgBox.exec();
    if(ret == QMessageBox::No)
        return;

    QFile defaultnldrparas("nldrparas_default");
    QMessageBox msgBox;

    if(!defaultnldrparas.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        msgBox.setText("error: can not open default parameter file \"nldrparas_default\"!");
        msgBox.exec();
        return;
    }

    defaultnldrparas.remove("nldrparas");
    defaultnldrparas.copy("nldrparas");

    ui->pushnldrcancel->click();
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
        msgBox.setText("error: can not open default parameter file \"plotparas_default\"!");
        msgBox.exec();
        return;
    }

    defaultplotparas.remove("plotparas");
    defaultplotparas.copy("plotparas");

    ui->pushplotcancel->click();
}

void setting::on_radionldrtopButton_clicked()
{
    ui->tablewidgetnldr->scrollToTop();
}

void setting::on_radionldrKruButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(6, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

void setting::on_radionldrNorButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(18, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

void setting::on_radionldrNLMButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(27, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

void setting::on_radionldrCCAButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(36, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

#endif
