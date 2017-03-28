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

#ifndef NLDRSETTING_CPP
#define NLDRSETTING_CPP

#include "treescaper.h"
#include "NLDRsetting.h"
#include "ui_NLDRsetting.h"
#include <QTableWidgetItem>
#include <iostream>
#include "warray.cpp"
#include <QMessageBox>
#include <QFile>

NLDRsetting::NLDRsetting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::NLDRsetting)
{
    ui->setupUi(this);
    ui->tablewidgetnldr->resizeColumnsToContents();
//    ui->pushnldrapply->setEnabled(true);
    nldrvalues = NULL;
    loadnldrparas();
}

NLDRsetting::~NLDRsetting()
{
    for(int i = 0; i < nnldrvalues; i++)
        delete nldrvalues[i];

    delete [] nldrvalues;
    delete ui;
}

void NLDRsetting::loadnldrparas()
{
    QMessageBox msgBox;

//======================load nldrparas===================
    File nldrfparas("nldrparas");
    if(!nldrfparas.is_open())
    {
        msgBox.setText("Error: Cannot open parameter files \"nldrparas\"!");
        msgBox.exec();
        return;
    }

    Array<String> nldrentries;
    String nldrelement;
    nldrfparas >> nldrelement;
    int countnldrelements = 1;      // Skip newline at the end of the file
    while(! nldrfparas.is_end())
    {
        nldrentries.add(nldrelement);
        if(countnldrelements > 101)
            break;
        else
            nldrfparas >> nldrelement;
        countnldrelements += 1;
    }

    nnldrvalues = 51;
    nldrvalues = new QTableWidgetItem *[nnldrvalues];
    for(int i = 0; i < nnldrvalues; i++)
    {
        nldrvalues[i] = new QTableWidgetItem((char *) nldrentries[2 * i + 1]);
        ui->tablewidgetnldr->setItem(i, 1, nldrvalues[i]);
    }
}

void NLDRsetting::on_pushnldrapply_clicked()
{
    File fparas("nldrparas");
    QMessageBox msgBox;
    if(!fparas.is_open())
    {
        msgBox.setText("Error: Cannot open parameter files \"nldrparas\"!");
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

//    ui->pushnldrapply->setEnabled(false);
    QDialog::accept();
}

void NLDRsetting::on_pushnldrcancel_clicked()
{
    loadnldrparas();
//    ui->pushnldrapply->setEnabled(true);
    QDialog::accept();
}

void NLDRsetting::on_pushnldrreset_clicked()
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
        msgBox.setText("Error: Cannot open default parameter file \"nldrparas_default\"!");
        msgBox.exec();
        return;
    }

    defaultnldrparas.remove("nldrparas");
    defaultnldrparas.copy("nldrparas");

    ui->pushnldrcancel->click();
}

void NLDRsetting::on_radionldrtopButton_clicked()
{
    ui->tablewidgetnldr->scrollToTop();
}

void NLDRsetting::on_radionldrKruButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(6, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

void NLDRsetting::on_radionldrNorButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(18, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

void NLDRsetting::on_radionldrNLMButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(27, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

void NLDRsetting::on_radionldrCCAButton_clicked()
{
    QTableWidgetItem *ptitem = ui->tablewidgetnldr->item(36, 0);
    ui->tablewidgetnldr->scrollToItem(ptitem, QAbstractItemView::PositionAtTop);
}

//void NLDRsetting::on_pushnldrclose_clicked()
//{
//    QDialog::accept();
//}

#endif
