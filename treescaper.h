
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

#ifndef TREESCAPER_H
#define TREESCAPER_H

#ifdef _WIN32
#undef LONGLONG __int64
//#include <windef.h>
#include <windows.h>
//#include <stddef.h>
#undef LONGLONG
#endif

#include "woutput.h"
#include "wImage.h"
#include <QMainWindow>
#include "wNLDRthread.h"
#include "wPlotthread.h"
#include "wDimestthread.h"
#include "wPlotTreethread.h"
#include "wPlotdimthread.h"
#include "wCommunitythread.h"
#include "wstring.h"
#include "wNLDR.h"
#include "Trees.h"
#include "plot2d.h"

namespace Ui {
    class TreeScaper;
}

namespace paras{
    extern nldr_parameters nldrparas;
    extern image_parameters imageparas;
    extern String plot_type;
    extern bool plot_makemovie;
    extern bool is_sort_volume;
};

class TreeScaper : public QMainWindow
{
    Q_OBJECT

public:
    explicit TreeScaper(QWidget *parent = 0, Qt::WindowFlags flags = 0);
    ~TreeScaper();

private:
    void initialize_paras(nldr_parameters &nldrparas, image_parameters &imageparas);

    Ui::TreeScaper *ui;
    NLDRthread nldrthd;
    Plotthread plotthd;
    Dimestthread dimestthd;
    Plotdimthread plotdimthd;
    PlotTreethread plottreethd;
    Communitythread communitythd;
    QTimer *timer;
    String CORplot;
    String DIMplot;
    Trees *TreesData;
    QDebugStream *qout;

signals:
    void sendplottree(const NEWICKTREE *tree, char *title, bool isrooted, bool isweighted, const LabelMap *lm);

private slots:
    void on_comboboxNLDRmethod_currentIndexChanged(QString );
    void on_pushdimplotfile_clicked();
    void on_pushNLDRplotfile_clicked();
    void on_pushDIMfile_clicked();
//    void on_pushNLDRfile_clicked();
    void on_pushdimdimest_clicked();
//    void on_pusdimclear_clicked();
    void on_pushdimplot_clicked();
    void on_pushdimrun_clicked();
    void on_pushNLDRsetting_clicked();
    void on_pushNLDRplot_clicked();
    void timerloadlog();
    void plottree(const NEWICKTREE *tree, char *title, bool isrooted, bool isweighted, const LabelMap *lm);
    void buttonvisible(int);
    void buttonCommunityPlotenable(int);
//    void HandleFileChange(QString);
    void on_pushNLDRclear_clicked();
    void on_pushNLDRrun_clicked();
    void on_comboBoxNLDRinit_currentIndexChanged(const QString &);
    void on_pushButton_clicked();
//    void on_pushDISTclear_clicked();
    void on_pushDISTrun_clicked();
    void on_pushDATAcov_clicked();
    void on_pushDATAButton_clicked();
    void on_pushDATAclear_clicked();
    void on_pushDATApartition_clicked();
    void on_pushDISTaffinity_clicked();
    void on_pushDATAloadtrees_clicked();
//    void on_pushDISToutput_clicked();
    void on_pushDISTloaddist_clicked();
//    void on_pushDISTcleardata_clicked();
    void on_pushDISTcomm_clicked();
//    void on_comboBoxDISTmemdata_currentIndexChanged(const QString &arg1);
    void on_pushNLDRplottrees_clicked();
    void on_pushButtonDATAcon_clicked();
    void on_pushButtonDATAconplot_clicked();
//    void on_pushCONVbipartfreqCumulative_clicked();
//    void on_pushCONVbipartfreqSlide_clicked();
    void on_pushButtonDATAloadlist_clicked();
    void on_radioButtonDATAindex_clicked();
    void on_radioButtonDATArange_clicked();
    void on_pushDATAButtonlistbrowse_clicked();
    void on_pushNLDRsaveindices_clicked();
    void on_pushDISTplotcom_clicked();
    void on_pushCOVAcomm_clicked();
//    void on_pushCONSclear_clicked();
//    void on_pushButtonCOVAlog_clicked();
    void on_pushCOVAplotcom_clicked();
    void on_comboBoxDISTmemdata_2_currentIndexChanged(const QString &arg1);
    void on_comboBoxCOVAmemdata_currentIndexChanged(const QString &arg1);
    void on_comboBoxNLDRdata_currentIndexChanged(const QString &arg1);
    void on_pushButtonCOVAbrowse_clicked();
    void on_pushCOVAloaddist_clicked();
    void Getnewlambda(double, Plot2D *, int);
#ifdef _MAC
    void closeEvent(QCloseEvent *bar){exit(0);};
#endif
    void on_pushCOVAclear_clicked();
    void on_textBrowserDATAlog_textChanged();
    void on_radioButtonCOVAauto_clicked();
    void on_radioButtonCOVAmanu_clicked();
    void on_radioButtonDISTauto_clicked();
    void on_radioButtonDISTmanu_clicked();
//    void buttonCONVcomboBox(int input);
//    void on_pushCONVclear_clicked();
    void comboBoxTreeData(int index);
    void comboBoxCD(int index);
    void comboBoxDimNLDR(int index);
    void on_comboBoxDISTmemdata_3_currentIndexChanged(const QString &arg1);
    void on_listDATAmem_itemActivated(QListWidgetItem *item);
    void on_listDATAmem_2_itemActivated(QListWidgetItem *item);
    void on_listDATAmem_3_itemActivated(QListWidgetItem *item);
};

#endif // TREESCAPER_H
