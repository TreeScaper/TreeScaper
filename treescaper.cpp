
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

//#include"DIST.h"
//#include "Bipartition.h"
#include <QtCore>
#include "wdef.h"
#include "wImage.h"
#include "treescaper.h"
#include "ui_treescaper.h"
#include <iostream>
#include <fstream>
#include "wstring.h"
#include "wNLDR.h"
#include <QFile>
#include <QTextStream>
#include "wNLDRthread.h"
#include "warray.cpp"
#include "setting.h"
#include "NLDRsetting.h"
#include "QMessageBox"
#include "sys/stat.h"
#include "QFileDialog"
#include <QStringList>
#include <QString>
#include <QChar>
#include "plot2d.h"
#include "plot2d_AWTY.h"
#include <QDebug>

/*
 "No data in memory": empty memory
 "Unweighted treeset" : Unweighted treeset in memory
 "Weighted treeset" : Weighted treeset in memory
 "Bipartition Matrix": Bipartitin matrix
 "Covariance Matrix": Covariance matrix
 "Unweighted RF-distance": Unweighted Robinson-Foulds distance
 "Weighted RF-distance": Weighted Robinson-Foulds distance
 "Matching-distance": Matching distance
 "SPR-distance": Subtree Prune and Regraft distance
 "Geodesic-distance": Geodesic distance (to do)
 "File-distance": distance from load file
 "File-coordinate": distance coordinates from load file
 "Affinity-URF": Affinity matrix for unweighted Robinson-Foulds distance
 "Affinity-RF": Affinity matrix for weighted Robinson-Foulds distance
 "Affinity-match": Affinity matrix for matching distance
 "Affinity-SPR": Affinity matrix for subtree prune and regraft distance
 "Affinity-geodesic": Affinity matrix for geodisic distance
 "Affinity-filedist": Affinity matrix for file distance
*/

namespace paras{
    nldr_parameters nldrparas;
    image_parameters imageparas;
    String plot_type;
    bool plot_makemovie;
    bool is_sort_volume;
};

TreeScaper::TreeScaper(QWidget *parent, Qt::WindowFlags flags) :
    QMainWindow(parent, flags),
    ui(new Ui::TreeScaper)
{
#ifdef _MAC
    QDir curDir("./");
    QString CurrentDirName = curDir.canonicalPath();
    QDir dir = QDir(QApplication::applicationDirPath());
    QString ModifiedDirName = dir.canonicalPath();
    QString logfileName;
    if(CurrentDirName == ModifiedDirName)
    {
        // Running App inside Qt Creator
        logfileName = CurrentDirName.append("/log.txt");
    }
    else
    {
        // Running App outside Qt Creator
        dir.cdUp();
        dir.cdUp();
        dir.cdUp();
        ModifiedDirName = dir.canonicalPath();
        dir.setCurrent(ModifiedDirName);
        logfileName = ModifiedDirName.append("/log.txt");
    }
    freopen(logfileName.toLocal8Bit().data(),"w",stdout);
    setvbuf(stdout, NULL, _IONBF, 0);
#elif defined(_WINDOWS)
    freopen("./log.txt","w",stdout);
#elif defined(_LINUX)
    freopen("./log.txt","w",stdout);
#endif


//    freopen("./log.txt","w",stdout);
//    setvbuf(stdout, NULL, _IONBF, 0);
//    QFileSystemWatcher * watcher = new QFileSystemWatcher(this);
//    watcher->addPath("./log.txt");
//    connect(watcher,SIGNAL(fileChanged(QString)),this,SLOT(HandleFileChange(QString)));
//    qDebug() << "hsdf\n";//---
//    setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);//-- | Qt::WindowMinMaxButtonsHint);
    ui->setupUi(this);
    timer = new QTimer(this);
    connect(timer,SIGNAL(timeout()),this,SLOT(timerloadlog()));
//    connect(this,SIGNAL(sendplottree(const NEWICKTREE *, String)),this,SLOT(plottree(const NEWICKTREE *, String)), Qt::QueuedConnection);
    connect(this,SIGNAL(sendplottree(const NEWICKTREE *, char *, bool , bool , const LabelMap *)),
            this,SLOT(plottree(const NEWICKTREE *, char *, bool , bool , const LabelMap *)), Qt::QueuedConnection);
    connect(&nldrthd, SIGNAL(sendbuttonvisible(int)), this, SLOT(buttonvisible(int)));
    connect(&communitythd, SIGNAL(sendbuttonCommunityPlotenable(int)), this, SLOT(buttonCommunityPlotenable(int)));
//    connect(ui->buttonCONVcomboBox, SIGNAL(activated(int)), this, SLOT(buttonCONVcomboBox(int)));
    connect(ui->comboBoxTreeData, SIGNAL(activated(int)), this, SLOT(comboBoxTreeData(int)));
    connect(ui->comboBoxCD, SIGNAL(activated(int)), this, SLOT(comboBoxCD(int)));
    connect(ui->comboBoxDimNLDR, SIGNAL(activated(int)), this, SLOT(comboBoxDimNLDR(int)));

    timer->start(500);
    TreesData = new Trees;
    initialize_paras(paras::nldrparas, paras::imageparas);
    ui->pushNLDRplot->setEnabled(false);
    ui->pushNLDRplottrees->setEnabled(false);
    ui->pushNLDRsaveindices->setEnabled(false);

    ui->comboBoxdimtype->setCurrentIndex(0);
    ui->comboboxNLDRmethod->setCurrentIndex(4);
    ui->comboBoxNLDRalgo->setCurrentIndex(2);
    ui->textdimstart->setReadOnly(true);
    ui->textdimend->setReadOnly(true);
    ui->textdimdimest->setReadOnly(true);
//    ui->pushDISToutput->setEnabled(false);
//    ui->pushDISTcleardata->setEnabled(false);
    ui->pushDISTcomm->setEnabled(false);
    ui->pushCOVAcomm->setEnabled(false);

    ui->labelDATAlistname->setEnabled(false);
    ui->textDATAlistfile->setEnabled(false);
    ui->pushDATAButtonlistbrowse->setEnabled(false);

    const char *cnull = "";
    QString qsnull(cnull);

    ui->textBrowserNLDRlog->setHtml(qsnull);
//    ui->textBrowserdimlog->setHtml(qsnull);
//    ui->textBrowserDISTlog->setHtml(qsnull);
    ui->textBrowserDATAlog->setHtml(qsnull);
//    ui->textBrowserCONVlog->setHtml(qsnull);

    ui->labelDISTaffinity->setEnabled(true);
    ui->comboBoxDISTaffinity->setEnabled(true);
    ui->pushDISTaffinity->setEnabled(false);
//    ui->labelDISTHighfreq->setEnabled(false);
//    ui->labelDISTLowfreq->setEnabled(false);
//    ui->textDISThighfreq->setEnabled(false);
//    ui->textDISTlowfreq->setEnabled(false);

    ui->pushCOVAplotcom->setEnabled(false);
    ui->pushDISTplotcom->setEnabled(false);
    this->resize(800, 500);

    qout = NULL;
//    qout = new QDebugStream (std::cout, ui->textBrowserDATAlog,
//                             ui->textBrowserCONSlog,
//                             ui->textBrowserCONVlog,
//                             ui->textBrowserCOVAlog,
//                             ui->textBrowserDISTlog,
//                             ui->textBrowserdimlog,
//                             ui->textBrowserNLDRlog,"./log.txt");

//    ui->textBrowserDATAlog->ensureCursorVisible();
//    ui->textBrowserCONSlog->ensureCursorVisible();
//    ui->textBrowserCONVlog->ensureCursorVisible();
//    ui->textBrowserCOVAlog->ensureCursorVisible();
//    ui->textBrowserDISTlog->ensureCursorVisible();
//    ui->textBrowserdimlog->ensureCursorVisible();
//    ui->textBrowserNLDRlog->ensureCursorVisible();

//    ui->labelCOVAlambda_2->setVisible(false);
//    ui->radioButtonCOVALp->setVisible(false);
//    ui->radioButtonCOVALn->setVisible(false);
//    ui->textCOVALfixed->setVisible(false);
//    ui->labelCOVAlambda->setVisible(false);
//    ui->labelCOVAfrom->setVisible(false);
//    ui->textCOVALstart->setVisible(false);
//    ui->labelCOVAto->setVisible(false);
//    ui->textCOVALend->setVisible(false);
//    ui->labelCOVAintv->setVisible(false);
//    ui->textCOVALinterval->setVisible(false);

    ui->labelCOVAlambda_2->setEnabled(false);
    ui->radioButtonCOVALp->setEnabled(false);
    ui->radioButtonCOVALn->setEnabled(false);
    ui->textCOVALfixed->setEnabled(false);
    ui->labelCOVAlambda->setEnabled(false);
    ui->labelCOVAfrom->setEnabled(false);
    ui->textCOVALstart->setEnabled(false);
    ui->labelCOVAto->setEnabled(false);
    ui->textCOVALend->setEnabled(false);
    ui->labelCOVAintv->setEnabled(false);
    ui->textCOVALinterval->setEnabled(false);

    ui->labelDISTlambda_2->setEnabled(false);
    ui->labelDISTfrom->setEnabled(false);
    ui->textDISTLstart->setEnabled(false);
    ui->labelDISTto->setEnabled(false);
    ui->textDISTLend->setEnabled(false);
    ui->labelDISTto_2->setEnabled(false);
    ui->textDISTLinterval->setEnabled(false);
}

void TreeScaper::buttonCommunityPlotenable(int input)
{
    if(input == 0)
    {
        ui->pushCOVAplotcom->setEnabled(false);
        ui->pushDISTplotcom->setEnabled(false);
    }
    if(input == 1)
    {
        ui->pushCOVAplotcom->setEnabled(false);
        ui->pushDISTplotcom->setEnabled(true);
    }
    if(input == 2)
    {
        ui->pushCOVAplotcom->setEnabled(true);
        ui->pushDISTplotcom->setEnabled(false);
    }
}

void TreeScaper::buttonvisible(int input)
{
    ui->pushNLDRplot->setEnabled(true);
    ui->pushNLDRplottrees->setEnabled(false);
    ui->pushNLDRsaveindices->setEnabled(true);
}



TreeScaper::~TreeScaper()
{
    if(qout != NULL)
        delete qout;
    delete timer;
    delete TreesData;
    delete [] paras::imageparas.cluster_index;
    delete ui;
}

void TreeScaper::initialize_paras(nldr_parameters &nldrparas, image_parameters &imageparas)
{
    QMessageBox msgBox;

//=====================load nldrparas==================
    File nldrfparas("nldrparas");
    if(!nldrfparas.is_open())
    {
        msgBox.setText("error: can not open parameter file \"nldrparas\"!");
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

    QString nldrqtstr[51];

    for(int i = 0; i < 51; i++)
        nldrqtstr[i] = (char *) nldrentries[2 * i + 1];

    nldrparas.random_start = nldrqtstr[0].toDouble();
    nldrparas.random_end = nldrqtstr[1].toDouble();
    nldrparas.length_tru = nldrqtstr[2].toInt();
    nldrparas.interval_tru = nldrqtstr[3].toInt();
    nldrparas.length_con = nldrqtstr[4].toInt();
    nldrparas.iterval_con = nldrqtstr[5].toInt();

    nldrparas.KRU_LIN_e = nldrqtstr[6].toDouble();
    nldrparas.KRU_LIN_max_iter = nldrqtstr[7].toInt();
    nldrparas.KRU_LIN_run_time = nldrqtstr[8].toInt();
    nldrparas.KRU_MAJ_e = nldrqtstr[9].toDouble();
    nldrparas.KRU_MAJ_max_iter = nldrqtstr[10].toInt();
    nldrparas.KRU_MAJ_run_time = nldrqtstr[11].toInt();
    nldrparas.KRU_GAU_e = nldrqtstr[12].toDouble();
    nldrparas.KRU_GAU_max_iter = nldrqtstr[13].toInt();
    nldrparas.KRU_GAU_run_time = nldrqtstr[14].toInt();
    nldrparas.KRU_STO_epochs = nldrqtstr[15].toInt();
    nldrparas.KRU_STO_alpha0 = nldrqtstr[16].toDouble();
    nldrparas.KRU_STO_alphan = nldrqtstr[17].toDouble();

    nldrparas.NOR_MAJ_e = nldrqtstr[18].toDouble();
    nldrparas.NOR_MAJ_max_iter = nldrqtstr[19].toInt();
    nldrparas.NOR_MAJ_run_time = nldrqtstr[20].toInt();
    nldrparas.NOR_GAU_e = nldrqtstr[21].toDouble();
    nldrparas.NOR_GAU_max_iter = nldrqtstr[22].toInt();
    nldrparas.NOR_GAU_run_time = nldrqtstr[23].toInt();
    nldrparas.NOR_STO_epochs = nldrqtstr[24].toInt();
    nldrparas.NOR_STO_alpha0 = nldrqtstr[25].toDouble();
    nldrparas.NOR_STO_alphan = nldrqtstr[26].toDouble();

    nldrparas.NLM_MAJ_e = nldrqtstr[27].toDouble();
    nldrparas.NLM_MAJ_max_iter = nldrqtstr[28].toInt();
    nldrparas.NLM_MAJ_run_time = nldrqtstr[29].toInt();
    nldrparas.NLM_GAU_e = nldrqtstr[30].toDouble();
    nldrparas.NLM_GAU_max_iter = nldrqtstr[31].toInt();
    nldrparas.NLM_GAU_run_time = nldrqtstr[32].toInt();
    nldrparas.NLM_STO_epochs = nldrqtstr[33].toInt();
    nldrparas.NLM_STO_alpha0 = nldrqtstr[34].toDouble();
    nldrparas.NLM_STO_alphan = nldrqtstr[35].toDouble();

    nldrparas.CCA_MAJ_e = nldrqtstr[36].toDouble();
    nldrparas.CCA_MAJ_max_iter = nldrqtstr[37].toInt();
    nldrparas.CCA_MAJ_run_time = nldrqtstr[38].toInt();
    nldrparas.CCA_MAJ_lambda0 = nldrqtstr[39].toDouble();
    nldrparas.CCA_MAJ_lambdan = nldrqtstr[40].toDouble();
    nldrparas.CCA_GAU_e = nldrqtstr[41].toDouble();
    nldrparas.CCA_GAU_max_iter = nldrqtstr[42].toInt();
    nldrparas.CCA_GAU_run_time = nldrqtstr[43].toInt();
    nldrparas.CCA_GAU_lambda0 = nldrqtstr[44].toDouble();
    nldrparas.CCA_GAU_lambdan = nldrqtstr[45].toDouble();
    nldrparas.CCA_STO_epochs = nldrqtstr[46].toInt();
    nldrparas.CCA_STO_lambda0 = nldrqtstr[47].toDouble();
    nldrparas.CCA_STO_lambdan = nldrqtstr[48].toDouble();
    nldrparas.CCA_STO_alpha0 = nldrqtstr[49].toDouble();
    nldrparas.CCA_STO_alphan = nldrqtstr[50].toDouble();

//================plotparas========================
    File plotfparas("plotparas");
    if(!plotfparas.is_open())
    {
        msgBox.setText("error: can not open parameter file \"plotparas\"!");
        msgBox.exec();
        return;
    }

    Array<String> entries;
    String element;
    plotfparas >> element;
	while(! plotfparas.is_end())
    {
        entries.add(element);
        plotfparas >> element;
    }

    paras::plot_type = entries[1];

    if(entries[3] == (String) "1")
        paras::plot_makemovie = true;
    else
        paras::plot_makemovie = false;

    QString qtstr(entries[5]);
    imageparas.cluster_num = qtstr.toInt();

    qtstr = (char *) entries[7];
    if(qtstr.toInt() == 1)
        paras::is_sort_volume = true;
    else
        paras::is_sort_volume = false;

    qtstr = (char *) entries[9];
    imageparas.color_type = qtstr.toInt();

    qtstr = (char *) entries[11];
    imageparas.outlier_e = qtstr.toDouble();

    qtstr = (char *) entries[13];
    imageparas.outlierpercent = qtstr.toDouble();

    qtstr = (char *) entries[15];
    imageparas.sep_factor = qtstr.toDouble();

    qtstr = (char *) entries[17];
    imageparas.point_size = qtstr.toDouble();

    qtstr = (char *) entries[19];
    imageparas.moviesizex = qtstr.toInt();

    qtstr = (char *) entries[21];
    imageparas.moviesizey = qtstr.toInt();

    qtstr = (char *) entries[23];
    imageparas.zoom = qtstr.toDouble();

    qtstr = (char *) entries[25];
    imageparas.frame_num = qtstr.toInt();

    int max_index = 0;
    QStringList indices, range;
    QString qtelement;
    for(int i = 0; i < imageparas.cluster_num; i++)
    {
        if(i * 2 + 27 >= entries.get_length())
        {
            msgBox.setText("error : the index setting of NLDR plot is wrong, please set again.");
            msgBox.exec();
            return;
        }
        qtelement = (char *) entries[i * 2 + 27];
        indices = qtelement.split(",", QString::SkipEmptyParts);
        for(int j = 0; j < indices.size(); j++)
        {
            range = indices[j].split("-", QString::SkipEmptyParts);
            for(int k = 0; k < range.size(); k++)
                if(range[k].toInt() > max_index)
                    max_index = range[k].toInt();
        }
    }

    imageparas.cluster_index = new int [max_index + 1];
	for(int i = 0; i <= max_index; i++)
		imageparas.cluster_index[i] = -1;
    imageparas.index_size = max_index + 1;

    for(int i = 0; i < imageparas.cluster_num; i++)
    {
        qtelement = (char *) entries[i * 2 + 27];
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
            /*for(int k = 0; k < range.size(); k = k + 2)
            {
                for(int h = range[k].toInt(); h <= range[k + 1].toInt(); h++)
                    imageparas.cluster_index[h] = i;
            }*/
        }
    }
}

void TreeScaper::on_pushNLDRrun_clicked()
{
    if(nldrthd.isRunning())
    {
        QMessageBox::warning(NULL, "warning", "NLDR is running!", 0, 0);
        return;
    }
    double **dist = NULL;
    int **distint = NULL;
    int size = 0;

    QString qtcost = ui->comboboxNLDRmethod->currentText();
    string stdcost = qtcost.toStdString();
    String cost(stdcost.c_str());
    if(cost == (String) "Classic Scaling")
        cost = "CLASSIC_MDS";
    else
    if(cost == (String) "Kruskal-1")
        cost = "KRUSKAL1";
    else
    if(cost == (String) "Normalized")
        cost = "NORMALIZED";
    else
    if(cost == (String) "NLM")
        cost = "SAMMON";

    QString qtalgo = ui->comboBoxNLDRalgo->currentText();
    string stdalgo = qtalgo.toStdString();
    String algo(stdalgo.c_str());
    if(algo == (String) "Linear Iteration MDS")
        algo = "LINEAR_ITERATION";
    else
    if(algo == (String) "Majorization")
        algo = "MAJORIZATION";
    else
    if(algo == (String) "Gauss Seidel")
        algo = "GAUSS_SEIDEL";
    else
    if(algo == (String) "Stochastic")
        algo = "STOCHASTIC";
    else
    if(algo == (String) "Metropolis")
        algo = "METROPOLIS";

    QString qtdim = ui->textNLDRdim->toPlainText();
    string stddim = qtdim.toStdString();
    String dim(stddim.c_str());

    QString qtflag = ui->textNLDRflag->toPlainText();
    string stdflag = qtflag.toStdString();
    String flag(stdflag.c_str());

    QString qtinit = ui->comboBoxNLDRinit->currentText();
    string stdinit = qtinit.toStdString();
    String init(stdinit.c_str());
    if(init == (String) "Random")
        init = "RAND";
    else
        init = "CLASSIC_MDS";

    QString qtseed = ui->textNLDRrandseed->toPlainText();
    int seed = qtseed.toInt();

    if(nldrthd.isRunning())
    {
        cout << "Warning : NLDR is running already!\n\n";
        return;
    }

    QString qtNLDRdata = ui->comboBoxNLDRdata->currentText();
    string stdNLDRdata = qtNLDRdata.toStdString();
    String NLDRdata(stdNLDRdata.c_str());

    QString qtfname = ui->textDATAfile->toPlainText();
    string stdfname = qtfname.toStdString();
    String fname(stdfname.c_str());
    File file(fname);
    CORplot = file.prefix_name();

    size = TreesData->Get_n_trees();
    if(NLDRdata == (String) "Unweighted RF-distance")
    {
        dist = TreesData->GetdistURF();
        CORplot += "_URF";
    } else if(NLDRdata == (String) "Weighted RF-distance")
    {
        dist = TreesData->GetdistRF();
        CORplot += "_RF";
    } else if(NLDRdata == (String) "Matching-distance")
    {
        distint = TreesData->GetdistMatch();
        CORplot += "_Match";
    } else if(NLDRdata == (String) "SPR-distance")
    {
        distint = TreesData->GetdistSPR();
        CORplot += "_SPR";
    }
    else if(NLDRdata == (String) "Geodesic-distance")
    {
        dist = TreesData->GetdistGeo();
        CORplot += "_Geo";
    } else if(NLDRdata == (String) "File-distance")
    {
        dist = TreesData->GetdistFile();
        size = TreesData->Get_filedistsize();

        QString qtdistfname = ui->textDISTfile_2->toPlainText();
        string stddistfname = qtdistfname.toStdString();
        String fdistname(stddistfname.c_str());
        fname = fdistname;
        File filedist(fdistname);
        CORplot = filedist.prefix_name();
    } else
    {
        cout << "Warning: There is no distance matrix in memory!\n\n";
        return;
    }

    if(NLDRdata == (String) "Matching-distance"
            || NLDRdata == (String) "SPR-distance")
    {
        dist = new double *[size];
        for(int i = 0; i < size; i++)
            dist[i] = new double [size];
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
            {
                dist[i][j] = (double) distint[i][j];
            }
    }

/*
        QString qttttfname = ui->textNLDRfile->toPlainText();
        string stdtttfname = qttttfname.toStdString();
        String ftttname(stdtttfname.c_str());
        File tttfile(ftttname);
        if(! tttfile.is_open())
        {
            cout << "error: can not open the data file!" << endl;
            return;
        }

        nldrthd.initialization(ftttname, (double **) NULL, 0, dim, cost, algo, init, flag, seed, paras::nldrparas);
*/
    nldrthd.initialization(CORplot, dist, size, dim, cost, algo, init, flag, seed, paras::nldrparas);
    nldrthd.start();
    CORplot += "_";
    CORplot += dim;
    CORplot += "D_";
    CORplot += cost;
    CORplot += "_COR_";
    CORplot += algo;
    CORplot += ".out";

    if(NLDRdata == (String) "Matching-distance"
       || NLDRdata == (String) "SPR-distance")
    {
        for(int i = 0; i < size; i++)
            delete [] dist[i];
        delete [] dist;
    }
}

void TreeScaper::on_pushNLDRclear_clicked()
{
    File fout("log.txt");
    fout.clean();
    fout.close();
    ui->textBrowserDATAlog->setText("");
//    ui->textBrowserCONSlog->setText("");
//    ui->textBrowserCONVlog->setText("");
    ui->textBrowserCOVAlog->setText("");
//    ui->textBrowserDISTlog->setText("");
//    ui->textBrowserdimlog->setText("");
    ui->textBrowserNLDRlog->setText("");
}

void TreeScaper::timerloadlog()
{
    char element;
    String line;
    Array<String> lines;
    int length = 0;
    String logtext;

    struct stat fstat;
    stat("log.txt",&fstat);
#ifdef _LINUX
    static long prensec = 0;
    static long presec = 0;
    if(prensec == fstat.st_mtim.tv_nsec && presec == fstat.st_mtim.tv_sec)
        return;
    prensec = fstat.st_mtim.tv_nsec;
    presec = fstat.st_mtim.tv_sec;
#endif
#ifdef _MAC
    static long prensec = 0;
    if(prensec == (long) fstat.st_size)
        return;
    prensec = (long) fstat.st_size;
#endif
#ifdef _WINDOWS
    static long prensec = 0;
    static long presec = 0;
    if(prensec == fstat.st_mtime && presec == fstat.st_size)
        return;
    prensec = fstat.st_mtime;
	presec = fstat.st_size;
#endif

    fstream file("log.txt", ios::out | ios::in);

    while(file.get(element))
    {
        if(element != '\n' && element != 13)
            line.add(element);

        if((element == '\n' || element == 13))
        {
            if(line != EMPTY_STRING)
            {
                lines.add(line);
                line = EMPTY_STRING;
            }
            else
                lines.add("\n");
        }
    }

    length = lines.get_length();

    for(int i = 0; i < length; i++)
    {
        logtext += lines[i];//--[length - 1 - i];
        logtext += "<br>";
    }

    const char *str = logtext;
    QString qtlogtext(str);

//    ui->textBrowserNLDRlog->setHtml(qtlogtext);
//    ui->textBrowserdimlog->setHtml(qtlogtext);
//    ui->textBrowserDISTlog->setHtml(qtlogtext);
//    ui->textBrowserDATAlog->setHtml(qtlogtext);

    ui->textBrowserNLDRlog->clear();
//    ui->textBrowserdimlog->clear();
//    ui->textBrowserDISTlog->clear();
    ui->textBrowserDATAlog->clear();
//    ui->textBrowserCONSlog->clear();
//    ui->textBrowserCONVlog->clear();
    ui->textBrowserCOVAlog->clear();

    ui->textBrowserNLDRlog->append(qtlogtext);
//    ui->textBrowserdimlog->append(qtlogtext);
//    ui->textBrowserDISTlog->append(qtlogtext);
    ui->textBrowserDATAlog->append(qtlogtext);
//    ui->textBrowserCONSlog->append(qtlogtext);
//    ui->textBrowserCONVlog->append(qtlogtext);
    ui->textBrowserCOVAlog->append(qtlogtext);
}
/*
void TreeScaper::HandleFileChange(QString qtlogtext)
{
    ui->textBrowserNLDRlog->setHtml(qtlogtext);
    ui->textBrowserdimlog->setHtml(qtlogtext);
    ui->textBrowserDISTlog->setHtml(qtlogtext);
    ui->textBrowserDATAlog->setHtml(qtlogtext);
}*/

void TreeScaper::on_pushNLDRplot_clicked()
{
    if(plotthd.isRunning())
    {
        cout << "Warning: The figure has already been plotted!\n\n";
        return;
    }

    String filename;
//    if(!ui->checkBoxNLDRplot->isChecked())
        filename = CORplot;
//    else
//    {
//        QString qtfname = ui->textNLDRplotfile->toPlainText();
//        string stdfname = qtfname.toStdString();
//       filename = stdfname.c_str();
//    }

    File file(filename);

    if(! file.is_open())
    {
        cout << "Error: The coordinates file \""<< filename << "\" can not be found! Run NLDR first!\n\n";//--- or select the coordinate file first by checking the box!" << endl;
        return;
    }

    if(file.lines() < paras::imageparas.index_size)
    {
        cout << "Error: The largest index should not be greater than the number of points!\n\n";
        return;
    }

    Array<int> selected_trees;
    if(TreesData->Get_n_trees() == 0)
    {
        selected_trees.resize(paras::imageparas.index_size + 1);
        selected_trees[0] = -1;
        plotthd.initialization(filename, paras::imageparas, selected_trees);
        plotthd.start();
        return;
    }

    ui->pushNLDRplottrees->setEnabled(true);
    cout << "\nTo plot the trees, press SHIFT + R and circle points from the figure. Then, click the 'Plot Trees' button.\n\n";

    plotthd.initialization(filename, paras::imageparas, TreesData->selected_trees);
    plotthd.start();

}

void TreeScaper::on_pushNLDRsetting_clicked()
{
    NLDRsetting nldrst(this);
    nldrst.show();
    nldrst.exec();
}

void TreeScaper::on_pushPlotsetting_clicked()
{
    setting st(this);
    st.show();
    st.exec();
}

void TreeScaper::on_pushdimrun_clicked()
{
    if(dimestthd.isRunning())
    {
        QMessageBox::warning(NULL, "Warning", "Dimension estimation is running!", 0, 0);
        return;
    }
    double **dist = NULL;
    int ** distint = NULL;
    int size = 0;
    int dim = 0;
    String type;

    QString qtest = ui->comboBoxdimestimator->currentText();
    string stdest = qtest.toStdString();
    String est(stdest.c_str());
    if(est == (String) "Correlation dimension")
        est = "CORR_DIM";
    else
    if(est == (String) "Maximum likelihood")
        est = "MLE_DIM";
    else
        est = "NN_DIM";

    QString qtnum = ui->textdimnumber->toPlainText();
    if(qtnum.toInt() < 2)
    {
        cout << "Warning: The number is at least 2!\n\n";
        return;
    }

    QString qtDIMdata = ui->comboBoxDIMdata->currentText();
    string stdDIMdata = qtDIMdata.toStdString();
    String DIMdata(stdDIMdata.c_str());

    QString qtfname = ui->textDATAfile->toPlainText();
    string stdfname = qtfname.toStdString();
    String fname(stdfname.c_str());
    File file(fname);
    DIMplot = file.prefix_name();
    if(DIMdata == (String) "Unweighted RF-distance")
    {
        dist = TreesData->GetdistURF();
        size = TreesData->Get_n_trees();
        DIMplot += "_URF";
        type = "DIS";
    } else if(DIMdata == (String) "Weighted RF-distance")
    {
        dist = TreesData->GetdistRF();
        size = TreesData->Get_n_trees();
        DIMplot += "_RF";
        type = "DIS";
    } else if(DIMdata == (String) "Matching-distance")
    {
        distint = TreesData->GetdistMatch();
        size = TreesData->Get_n_trees();
        DIMplot += "_Match";
        type = "DIS";
    } else if(DIMdata == (String) "SPR-distance")
    {
        distint = TreesData->GetdistSPR();
        size = TreesData->Get_n_trees();
        DIMplot += "_SPR";
        type = "DIS";
    }
    else if(DIMdata == (String) "Geodesic-distance")
    {
        dist = TreesData->GetdistGeo();
        size = TreesData->Get_n_trees();
        DIMplot += "_Geo";
        type = "DIS";
    }
    else if(DIMdata == (String) "File-distance")
    {
        dist = TreesData->GetdistFile();
        size = TreesData->Get_filedistsize();
        type = "DIS";

        QString qtdistfname = ui->textDISTfile_2->toPlainText();
        string stddistfname = qtdistfname.toStdString();
        String fdistname(stddistfname.c_str());
        File filedist(fdistname);
        DIMplot = filedist.prefix_name();
    }
    else if(DIMdata == (String) "File-coordinate")
    {
        cout << "Here1:" << endl; //-- MMtest
        dist = TreesData->GetcoordFile();
        size = TreesData->Get_filedcoordinatesize();
        dim = TreesData->Get_filedcoordinatedim();

        cout << "size: " << size << " dim: " << dim << endl;
        type = "COR";

        QString qtcoordfname = ui->textDISTfile_2->toPlainText();
        string stdcoordfname = qtcoordfname.toStdString();
        String fcoordname(stdcoordfname.c_str());
        File filecoord(fcoordname);
        DIMplot = filecoord.prefix_name();
    }
    else
    {
        cout << "Warning: There is no distance matrix or coordinate data in memory!\n\n";
        return;
    }

    if(DIMdata == (String) "Matching-distance" || DIMdata == (String) "SPR-distance")
    {
        dist = new double *[size];
        for(int i = 0; i < size; i++)
            dist[i] = new double [size];
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
            {
                dist[i][j] = (double) distint[i][j];
            }
    }


    dimestthd.initialization(DIMplot, dist, size, dim, est, type, qtnum.toInt());
    dimestthd.start();

//    DIMplot = file.prefix_name();
    DIMplot += "_";
    DIMplot += est;
    if(est == (String) "MLE_DIM")
        DIMplot += "_k.out";
    else
        DIMplot += "_logvslog.out";
    ui->textdimstart->setReadOnly(false);
    ui->textdimend->setReadOnly(false);
    if(ui->textdimstart->toPlainText() == "(disable)")
    {
        ui->textdimstart->setText("(input)");
        ui->textdimend->setText("(input)");
    }

    if(DIMdata == (String) "Matching-distance" || DIMdata == (String) "SPR-distance")
    {
        for(int i = 0; i < size; i++)
            delete [] dist[i];
        delete [] dist;
    }
}

//void TreeScaper::on_pushdimrun_clicked()
//{
//    QString qtfname = ui->textdimfile->toPlainText();
//    string stdfname = qtfname.toStdString();
//    String fname(stdfname.c_str());
//    File file(fname);
//    if(! file.is_open())
//    {
//        cout << "Error: Can not open the data file!\n\n";
//        return;
//    }

//    QString qttype = ui->comboBoxdimtype->currentText();
//    string stdtype = qttype.toStdString();
//    String type(stdtype.c_str());
//    if(type == (String) "Distance")
//        type = "DIS";
//    else
//        type = "COR";

//    QString qtest = ui->comboBoxdimestimator->currentText();
//    string stdest = qtest.toStdString();
//    String est(stdest.c_str());
//    if(est == (String) "Correlation dimension")
//        est = "CORR_DIM";
//    else
//    if(est == (String) "Maximum likelihood")
//        est = "MLE_DIM";
//    else
//        est = "NN_DIM";

//    QString qtnum = ui->textdimnumber->toPlainText();
//    if(qtnum.toInt() < 2)
//    {
//        cout << "Warning: The number is at least 2!\n\n";
//        return;
//    }

//    if(dimestthd.isRunning())
//    {
//        cout << "Warning: Estimating the dimension now!\n\n";
//        return;
//    }

//    dimestthd.initialization(fname, est, type, qtnum.toInt());
//    dimestthd.start();

//    DIMplot = file.prefix_name();
//    DIMplot += "_";
//    DIMplot += est;
//    if(est == (String) "MLE_DIM")
//        DIMplot += "_k.out";
//    else
//        DIMplot += "_logvslog.out";
//    ui->textdimstart->setReadOnly(false);
//    ui->textdimend->setReadOnly(false);
//    if(ui->textdimstart->toPlainText() == "(disable)")
//    {
//        ui->textdimstart->setText("(input)");
//        ui->textdimend->setText("(input)");
//    }
//}

void TreeScaper::on_pushdimplot_clicked()
{
    if(plotdimthd.isRunning())
    {
        cout << "Warning: The figure has already been plotted!\n\n";
        return;
    }

    String filename;
//    if(!ui->checkBoxdimplot->isChecked())
        filename = DIMplot;
//    else
//    {
//        QString qtfname = ui->textdimplotfile->toPlainText();
//        string stdfname = qtfname.toStdString();
//        filename = stdfname.c_str();
//    }

    File file(filename);
    if(! file.is_open())
    {
        cout << "Error: The dimension file was not specified!\n\n";
        return;
    }

    plotdimthd.initialization(filename);
    plotdimthd.start();
}

void TreeScaper::on_pushdimdimest_clicked()
{
    QString qtstart = ui->textdimstart->toPlainText();
    QString qtend = ui->textdimend->toPlainText();
    double start = qtstart.toDouble();
    double end = qtend.toDouble();
    if(start >= end)
    {
        cout << "Error: The initial value should be less than final value!\n\n";
        return;
    }


    String filename;
//    if(!ui->checkBoxdimplot->isChecked())
        filename = DIMplot;
//    else
//    {
//        QString qtfname = ui->textdimplotfile->toPlainText();
//        string stdfname = qtfname.toStdString();
//        filename = stdfname.c_str();
//    }

    File file(filename);
    if(! file.is_open())
    {
        cout << "Error: The dimension file was not specified!\n\n";
        return;
    }

    Array<String> entries;
    String element;
    file >> element;
    int num_element = 0;
	while(! file.is_end())
    {
        entries.add(element);
        num_element++;
        file >> element;
    }

    QString value, fvalue, slope;
    double x, y, fx, fy;
    double sum = 0;
    int num = 0;
    value = (char *) entries[0];
    fvalue = (char *) entries[1];
    x = value.toDouble();
    fx = fvalue.toDouble();
    value = (char *) entries[num_element - 2];
    fvalue = (char *) entries[num_element - 1];
    y = value.toDouble();
    fy = fvalue.toDouble();

    QString qtest = ui->comboBoxdimestimator->currentText();
    string stdest = qtest.toStdString();
    String est(stdest.c_str());

    if(est == (String) "Correlation dimension" || est == (String) "Nearest neighbor")
    {
        for(int i = 0; i < num_element / 2; i++)
        {
            if(start >= value.toDouble())
            {
                x = value.toDouble();
                fx = fvalue.toDouble();
            }
            value = (char *) entries[i * 2];
            fvalue = (char *) entries[i * 2 + 1];
            if(end >= value.toDouble())
            {
                y = value.toDouble();
                fy = fvalue.toDouble();
            }
        }
        slope = slope.number((fy - fx) / (y - x));
        ui->textdimdimest->setPlainText(slope);
    } else
    {
        for(int i = 0; i < num_element / 2; i++)
        {
            if(start <= value.toDouble())
            {
                sum += fvalue.toDouble();
                num++;
            }
            value = (char *) entries[i * 2];
            fvalue = (char *) entries[i * 2 + 1];
            if(end <= value.toDouble())
                break;
        }

        slope = slope.number(sum / num);
        ui->textdimdimest->setPlainText(slope);
    }
}
/*
void TreeScaper::on_checkBoxdimplot_stateChanged(int )
{
    if(ui->checkBoxdimplot->isChecked())
    {
        ui->labeldimplotfile->setVisible(true);
        ui->textdimplotfile->setVisible(true);
        ui->pushdimplotfile->setVisible(true);
        ui->textdimstart->setReadOnly(false);
        ui->textdimend->setReadOnly(false);
        if(ui->textdimstart->toPlainText() == "(disable)")
        {
            ui->textdimstart->setText("(input)");
            ui->textdimend->setText("(input)");
        }
    } else
    {
        ui->labeldimplotfile->setVisible(false);
        ui->textdimplotfile->setVisible(false);
        ui->pushdimplotfile->setVisible(false);
        if(DIMplot.is_empty())
        {
            ui->textdimstart->setReadOnly(true);
            ui->textdimend->setReadOnly(true);
            ui->textdimstart->setText("(disable)");
            ui->textdimend->setText("(disable)");
        }
    }
}
*/
/*
void TreeScaper::on_pushNLDRfile_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
//    fileDialog->setFilter(tr("Files(*)"));
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];
        ui->textNLDRfile->setPlainText(path);
    }
}*/

//void TreeScaper::on_pushDIMfile_clicked()
//{
//    QFileDialog *fileDialog = new QFileDialog(this);
//    fileDialog->setWindowTitle(tr("Open File"));
//    fileDialog->setDirectory(".");
////    fileDialog->setFilter(tr("Files(*)"));
//    if(fileDialog->exec() == QDialog::Accepted)
//    {
//        QString path = fileDialog->selectedFiles()[0];
//        ui->textdimfile->setPlainText(path);
//    }
//}

void TreeScaper::on_pushNLDRplotfile_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    QString path;
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
    
    Array<int> selected_trees(paras::imageparas.index_size + 1);
    selected_trees[0] = -1;
//    fileDialog->setFilter(tr("Files(*)"));
    if(fileDialog->exec() == QDialog::Accepted)
    {
        path = fileDialog->selectedFiles()[0];
        if(plotthd.isRunning())
        {
            cout << "Warning: The figure has already been plotted!\n\n";
            return;
        }

        String filename;
        string stdfname = path.toStdString();
        filename = stdfname.c_str();

        File file(filename);

        if(! file.is_open())
        {
            cout << "Error: The coordinates file cannot be found!\n\n";
            return;
        }

        if(file.lines() < paras::imageparas.index_size)
        {
            cout << "Error: The largest index should not be greater than the number of points!\n\n";
            return;
        }
        plotthd.initialization(filename, paras::imageparas, selected_trees);
        plotthd.start();
    }
//    delete [] selected_trees;
}

void TreeScaper::on_pushdimplotfile_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
//    fileDialog->setFilter(tr("Files(*)"));
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];

        if(plotdimthd.isRunning())
        {
            cout << "Warning: The figure has already been plotted!\n\n";
            return;
        }

        String filename;
        string stdfname = path.toStdString();
        filename = stdfname.c_str();
        DIMplot = filename;

        File file(filename);
        if(! file.is_open())
        {
            cout << "Error: The dimension file was not specified!\n\n";
            return;
        }

        cout << "Successfully read dimension file.\n\n";

        ui->textdimstart->setReadOnly(false);
        ui->textdimend->setReadOnly(false);
        if(ui->textdimstart->toPlainText() == "(disable)")
        {
            ui->textdimstart->setText("(input)");
            ui->textdimend->setText("(input)");
        }
    }
}

void TreeScaper::on_comboboxNLDRmethod_currentIndexChanged(QString )
{
    QString qtmethod = ui->comboboxNLDRmethod->currentText();
    string stdmethod = qtmethod.toStdString();
    String method(stdmethod.c_str());

    QString SVD = "SVD";
    QString Linear_Iteration = "Linear Iteration MDS";
    QString Majorization = "Majorization";
    QString Gauss_Seidel = "Gauss Seidel";
    QString Stochastic = "Stochastic";
    QString Metropolis = "Metropolis";


    ui->comboBoxNLDRalgo->clear();
    if(method == (String) "Classic Scaling")
    {
        ui->comboBoxNLDRalgo->insertItem(0, SVD);
    } else
    if(method == (String) "Kruskal-1")
    {
        ui->comboBoxNLDRalgo->insertItem(0, Linear_Iteration);
        ui->comboBoxNLDRalgo->insertItem(1, Majorization);
        ui->comboBoxNLDRalgo->insertItem(2, Gauss_Seidel);
        ui->comboBoxNLDRalgo->insertItem(3, Metropolis);
    } else
    if(method == (String) "Normalized")
    {
        ui->comboBoxNLDRalgo->insertItem(0, Majorization);
        ui->comboBoxNLDRalgo->insertItem(1, Gauss_Seidel);
        ui->comboBoxNLDRalgo->insertItem(2, Stochastic);
        ui->comboBoxNLDRalgo->insertItem(3, Metropolis);
    } else
    if(method == (String) "NLM")
    {
        ui->comboBoxNLDRalgo->insertItem(0, Majorization);
        ui->comboBoxNLDRalgo->insertItem(1, Gauss_Seidel);
        ui->comboBoxNLDRalgo->insertItem(2, Stochastic);
        ui->comboBoxNLDRalgo->insertItem(3, Metropolis);
    } else
    if(method == (String) "CCA")
    {
        ui->comboBoxNLDRalgo->insertItem(0, Majorization);
        ui->comboBoxNLDRalgo->insertItem(1, Gauss_Seidel);
        ui->comboBoxNLDRalgo->insertItem(2, Stochastic);
        ui->comboBoxNLDRalgo->insertItem(3, Metropolis);
    }
}

void TreeScaper::on_comboBoxNLDRinit_currentIndexChanged(const QString &qttext)
{
    string stdtext = qttext.toStdString();
    String text(stdtext.c_str());

    QString stdseed = ui->textNLDRrandseed->toPlainText();

    if(!(text == (String) "Random"))
    {
        ui->textNLDRrandseed->setText(stdseed + "(disable)");
        ui->textNLDRrandseed->setReadOnly(true);
    } else
    {
        ui->textNLDRrandseed->setReadOnly(false);
        stdseed.chop(9);
        ui->textNLDRrandseed->setText(stdseed);
    }
}

void TreeScaper::on_pushButton_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
//    fileDialog->setFilter(tr("Files(*)"));
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];
        ui->textDISTfile_2->setPlainText(path);
    }
}

void TreeScaper::on_pushDISTrun_clicked()
{
    QString qtmethod = ui->comboBoxDISTmethod->currentText();
    string stdmethod = qtmethod.toStdString();
    String method(stdmethod.c_str());
    bool dis;
    if(!TreesData->treesAreexisting())
    {
        cout << "Warning: There are no trees in the memory!\n\n";
        return;
    }
    else
    if(method == (String) "Unweighted Robinson Foulds")
    {
        if(!TreesData->bipartmatrixIsexisting())
        {
            cout << "Warning: There is no bipartition matrix in the memory! Please compute it first.\n\n";
            return;
        }
        dis = TreesData->Compute_RF_dist_by_hash(false);
    }
    else
    if(method == (String) "Weighted Robinson Foulds")
    {
        if(!TreesData->bipartmatrixIsexisting())
        {
            cout << "Warning: There is no bipartition matrix in the memory! Please compute it first.\n\n";
            return;
        }
        dis = TreesData->Compute_RF_dist_by_hash(true);
    }
    else
    if(method == (String) "Matching")
    {
        dis = TreesData->Compute_Matching_dist();
    }
    else
    if(method == (String) "Subtree Prune and Regraft")
    {
        dis = TreesData->Compute_SPR_dist();
    }

    if (dis)
    {
        cout << "Successfully computed " << method << " distance.\n\n";
        if(method == (String) "Unweighted Robinson Foulds")
        {
            QString qstr("Unweighted RF-distance");
            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxNLDRdata->findText("Unweighted RF-distance") == -1)
                ui->comboBoxNLDRdata->addItem(qstr);
            if (ui->comboBoxDIMdata->findText("Unweighted RF-distance") == -1)
                ui->comboBoxDIMdata->addItem(qstr);
            if (ui->comboBoxDISTmemdata_2->findText("Unweighted RF-distance") == -1)
                ui->comboBoxDISTmemdata_2->addItem(qstr);
        }
        if(method == (String) "Weighted Robinson Foulds")
        {
            QString qstr("Weighted RF-distance");
            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxNLDRdata->findText("Weighted RF-distance") == -1)
                ui->comboBoxNLDRdata->addItem(qstr);
            if (ui->comboBoxDIMdata->findText("Weighted RF-distance") == -1)
                ui->comboBoxDIMdata->addItem(qstr);
            if (ui->comboBoxDISTmemdata_2->findText("Weighted RF-distance") == -1)
                ui->comboBoxDISTmemdata_2->addItem(qstr);
        }
        else
        if(method == (String) "Matching")
        {
            QString qstr("Matching-distance");
            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxNLDRdata->findText("Matching-distance") == -1)
                ui->comboBoxNLDRdata->addItem(qstr);
            if (ui->comboBoxDIMdata->findText("Matching-distance") == -1)
                ui->comboBoxDIMdata->addItem(qstr);
            if (ui->comboBoxDISTmemdata_2->findText("Matching-distance") == -1)
                ui->comboBoxDISTmemdata_2->addItem(qstr);
        }
        else
        if(method == (String) "Subtree Prune and Regraft")
        {
            QString qstr("SPR-distance");
            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxNLDRdata->findText("SPR-distance") == -1)
                ui->comboBoxNLDRdata->addItem(qstr);
            if (ui->comboBoxDIMdata->findText("SPR-distance") == -1)
                ui->comboBoxDIMdata->addItem(qstr);
            if (ui->comboBoxDISTmemdata_2->findText("SPR-distance") == -1)
                ui->comboBoxDISTmemdata_2->addItem(qstr);
        }

        if(ui->comboBoxDISTmemdata_2->findText("No distance data in memory") != -1)
        {
            int idex = ui->comboBoxDISTmemdata_2->findText("No distance data in memory");
            ui->comboBoxDISTmemdata_2->removeItem(idex);
        }
        if(ui->comboBoxNLDRdata->findText("No distance data in memory") != -1)
        {
            int idex = ui->comboBoxNLDRdata->findText("No distance data in memory");
            ui->comboBoxNLDRdata->removeItem(idex);
        }
        if(ui->comboBoxDIMdata->findText("No distance/coordinate data in memory") != -1)
        {
            int idex = ui->comboBoxDIMdata->findText("No distance/coordinate data in memory");
            ui->comboBoxDIMdata->removeItem(idex);
        }
    }
    ui->pushDISTaffinity->setEnabled(true);

    cout << "Successfully computed distance matrix.\n\n";
}

void TreeScaper::on_pushDATAcov_clicked()
{
    if(!TreesData->bipartmatrixIsexisting())
    {
        cout << "Warning: There is no bipartition matrix in the memory! Please compute it first.\n\n";
        return;
    }
    TreesData->Compute_Bipart_Covariance();
    QString qstr("Covariance Matrix");
    QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item1.empty())
        ui->listDATAmem->addItem(qstr);

    QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item2.empty())
        ui->listDATAmem_2->addItem(qstr);

    QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item3.empty())
        ui->listDATAmem_3->addItem(qstr);

    if(ui->comboBoxCOVAmemdata->findText("Covariance Matrix") == -1)
        ui->comboBoxCOVAmemdata->addItem(qstr);

    if(ui->comboBoxCOVAmemdata->findText("No covariance data in memory") != -1)
    {
        int idex = ui->comboBoxCOVAmemdata->findText("No covariance data in memory");
        ui->comboBoxCOVAmemdata->removeItem(idex);
    }
    ui->pushCOVAcomm->setEnabled(true);
    cout << "Successfully computed covariance matrix of bipartitions.\n\n";
}

void TreeScaper::on_pushDATAButton_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
//    fileDialog->setFilter(tr("Files(*)"));
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];
        ui->textDATAfile->setPlainText(path);
    }

}

void TreeScaper::on_pushDATAclear_clicked()
{
    File fout("log.txt");
    fout.clean();
    fout.close();
    ui->textBrowserDATAlog->setText("");
//    ui->textBrowserCONSlog->setText("");
//    ui->textBrowserCONVlog->setText("");
    ui->textBrowserCOVAlog->setText("");
//    ui->textBrowserDISTlog->setText("");
//    ui->textBrowserdimlog->setText("");
    ui->textBrowserNLDRlog->setText("");
}

void TreeScaper::on_pushDATApartition_clicked()
{
    if(!TreesData->treesAreexisting())
    {
        cout << "There are no trees in the memory! Please load trees first.\n\n";
        return;
    }
    if(TreesData->bipartmatrixIsexisting())
    {
        TreesData->deleteBipartitionMatrix();
    }

    TreesData->Compute_Hash();
    TreesData->Compute_Bipart_Matrix();

    QString qstr("Bipartition Matrix");
    QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item1.empty())
        ui->listDATAmem->addItem(qstr);

    QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item2.empty())
        ui->listDATAmem_2->addItem(qstr);

    QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item3.empty())
        ui->listDATAmem_3->addItem(qstr);

    cout << "Successfully computed bipartition matrix.\n\n";
}

void TreeScaper::on_pushDATAloadtrees_clicked()
{
    TreesData->destructor();
    delete TreesData;
    TreesData = new Trees;

    QString qtfname = ui->textDATAfile->toPlainText();
    string stdfname = qtfname.toStdString();
    String fname(stdfname.c_str());
    File file(fname);
    if(!file.is_open())
    {
        std::cout << "Error: Cannot open the data file!\n\n";
        return;
    }
    TreesData->initialTrees(stdfname);

//    QString qtbipart = ui->comboBoxDATAinformat->currentText();
//    string stdbipart = qtbipart.toStdString();
//    String bipart_type(stdbipart.c_str());
//    if(bipart_type == (String) "unweighted")
//        TreesData->Settreeweighttype(false);
//    else if(bipart_type == (String) "weighted")
//        TreesData->Settreeweighttype(true);

    if(ui->checkBoxDATAweighted->isChecked())
        TreesData->Settreeweighttype(true);
    else
        TreesData->Settreeweighttype(false);

    if(ui->checkBoxDATArooted->isChecked())
        TreesData->Settreeroottype(true);
    else
        TreesData->Settreeroottype(false);

//    QString qttreetype = ui->comboBoxDATAtreetype->currentText();
//    string stdtreetype = qttreetype.toStdString();
//    String treetype(stdtreetype.c_str());
//    if(treetype == (String) "rooted")
//        TreesData->Settreeroottype(true);
//    else
//        TreesData->Settreeroottype(false);
    TreesData->ReadTrees();

//    int aaa[] = {350};
//    TreesData->Printf(aaa, 1);
    TreesData->compute_numofbipart();
    ui->listDATAmem->clear();
    ui->listDATAmem_2->clear();
    ui->listDATAmem_3->clear();
    ui->comboBoxCOVAmemdata->clear();
    ui->comboBoxDISTmemdata_2->clear();
    ui->comboBoxDISTmemdata_3->clear();
    ui->comboBoxNLDRdata->clear();
    ui->comboBoxDIMdata->clear();

    if(! ui->checkBoxDATAweighted->isChecked()) //(bipart_type == (String) "unweighted")
    {
        QString qstr("Unweighted treeset");
        QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item1.empty())
            ui->listDATAmem->addItem(qstr);

        QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item2.empty())
            ui->listDATAmem_2->addItem(qstr);

        QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item3.empty())
            ui->listDATAmem_3->addItem(qstr);

    } else if(ui->checkBoxDATAweighted->isChecked())
    {
        QString qstr("Weighted treeset");
        QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item1.empty())
            ui->listDATAmem->addItem(qstr);

        QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item2.empty())
            ui->listDATAmem_2->addItem(qstr);

        QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item3.empty())
            ui->listDATAmem_3->addItem(qstr);
    }

    Array<int> *treeidx = TreesData->getidxlist();
    treeidx->resize(TreesData->Get_n_trees());
    for(int i = 0; i < TreesData->Get_n_trees(); i++)
        (*treeidx)[i] = i;

    file.close();
//    ui->pushDISTcleardata->setEnabled(true);
//    ui->pushDISTaffinity->setEnabled(false);
//    ui->pushDISToutput->setEnabled(false);


    cout << "Successfully read " << TreesData->Get_n_trees() << " trees.\n\n";

//    Array<int> *tttidx = TreesData->getidxlist();
//    for(int i = 0; i < tttidx->get_length(); i++)
//        std::cout << "i:" << i << ", idx:" << (*tttidx)[i] << std::endl;//----

//    int arr[3] = {0,1,2};
//    TreesData->Printf(arr, 3);
}

void TreeScaper::on_pushDISTloaddist_clicked()
{
    int msg1 = QMessageBox::question(this, "Load File", "Are you sure the input file and format are correct?", "Yes", "No");
    if(msg1 == 1) // No
        return;

    QString qtfname = ui->textDISTfile_2->toPlainText();
    string stdfname = qtfname.toStdString();
    String fname(stdfname.c_str());
    File file(fname);
    if(!file.is_open())
    {
        cout << "Error: Cannot open the data file!\n\n";
        return;
    }

    QString qttype = ui->comboBoxdimtype->currentText();
    string stdtype = qttype.toStdString();
    String type(stdtype.c_str());

    if(type == (String) "Distance")
    {
        TreesData->load_distfile(stdfname);

        QString qstr("File-distance");
        QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item1.empty())
            ui->listDATAmem->addItem(qstr);

        QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item2.empty())
            ui->listDATAmem_2->addItem(qstr);

        QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item3.empty())
            ui->listDATAmem_3->addItem(qstr);

        QList<QListWidgetItem *> itemTmp1 = ui->listDATAmem->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
        for (int i = 0; i < itemTmp1.size(); i++)
            delete ui->listDATAmem->takeItem(ui->listDATAmem->row(itemTmp1[i]));

        QList<QListWidgetItem *> itemTmp2 = ui->listDATAmem_2->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
        for (int i = 0; i < itemTmp2.size(); i++)
            delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(itemTmp2[i]));

        QList<QListWidgetItem *> itemTmp3 = ui->listDATAmem_3->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
        for (int i = 0; i < itemTmp3.size(); i++)
            delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(itemTmp3[i]));

        if (ui->comboBoxNLDRdata->findText("File-distance") == -1)
            ui->comboBoxNLDRdata->addItem(qstr);

        int idex = ui->comboBoxNLDRdata->findText("No distance data in memory");
        ui->comboBoxNLDRdata->removeItem(idex);

        if (ui->comboBoxDIMdata->findText("File-distance") == -1)
            ui->comboBoxDIMdata->addItem(qstr);

        idex = ui->comboBoxDIMdata->findText("No distance/coordinate data in memory");
        ui->comboBoxDIMdata->removeItem(idex);

        if (ui->comboBoxDISTmemdata_2->findText("File-distance") == -1)
            ui->comboBoxDISTmemdata_2->addItem(qstr);

        idex = ui->comboBoxDISTmemdata_2->findText("No distance data in memory");
        ui->comboBoxDISTmemdata_2->removeItem(idex);

        ui->pushDISTaffinity->setEnabled(true);

        cout << "Successfully read distance file.\n\n";
    }
    else
    {
        TreesData->load_coordinatefile(stdfname);

        QString qstr("File-coordinate");
        QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item1.empty())
            ui->listDATAmem->addItem(qstr);

        QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item2.empty())
            ui->listDATAmem_2->addItem(qstr);

        QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item3.empty())
            ui->listDATAmem_3->addItem(qstr);

        QList<QListWidgetItem *> itemTmp1 = ui->listDATAmem->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
        for (int i = 0; i < itemTmp1.size(); i++)
            delete ui->listDATAmem->takeItem(ui->listDATAmem->row(itemTmp1[i]));

        QList<QListWidgetItem *> itemTmp2 = ui->listDATAmem_2->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
        for (int i = 0; i < itemTmp2.size(); i++)
            delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(itemTmp2[i]));

        QList<QListWidgetItem *> itemTmp3 = ui->listDATAmem_3->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
        for (int i = 0; i < itemTmp3.size(); i++)
            delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(itemTmp3[i]));

        if (ui->comboBoxDIMdata->findText("File-coordinate") == -1)
            ui->comboBoxDIMdata->addItem(qstr);

        int idex = ui->comboBoxDIMdata->findText("No distance/coordinate data in memory");
        ui->comboBoxDIMdata->removeItem(idex);

        cout << "Successfully read coordinate file.\n\n";
    }
    file.close();
}

//void TreeScaper::on_pushDISTloaddist_clicked()
//{
//    QString qtfname = ui->textDISTfile_2->toPlainText();
//    string stdfname = qtfname.toStdString();
//    String fname(stdfname.c_str());
//    File file(fname);
//    if(! file.is_open())
//    {
//        cout << "Error: Cannot open the data file!\n\n";
//        return;
//    }
//    TreesData->load_distfile(stdfname);

//    QString qstr("File-distance");
////    if (ui->comboBoxDISTmemdata->findText("File-distance") == -1)
////        ui->comboBoxDISTmemdata->addItem(qstr);

//    QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
//    if(item1.empty())
//        ui->listDATAmem->addItem(qstr);

//    QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
//    if(item2.empty())
//        ui->listDATAmem_2->addItem(qstr);

//    QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
//    if(item3.empty())
//        ui->listDATAmem_3->addItem(qstr);

////    int idex = ui->comboBoxDISTmemdata->findText("No data in memory");
////    ui->comboBoxDISTmemdata->removeItem(idex);

//    QList<QListWidgetItem *> itemTmp1 = ui->listDATAmem->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
//    for (int i = 0; i < itemTmp1.size(); i++)
//        delete ui->listDATAmem->takeItem(ui->listDATAmem->row(itemTmp1[i]));

//    QList<QListWidgetItem *> itemTmp2 = ui->listDATAmem_2->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
//    for (int i = 0; i < itemTmp2.size(); i++)
//        delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(itemTmp2[i]));

//    QList<QListWidgetItem *> itemTmp3 = ui->listDATAmem_3->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
//    for (int i = 0; i < itemTmp3.size(); i++)
//        delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(itemTmp3[i]));

//    if (ui->comboBoxNLDRdata->findText("File-distance") == -1)
//        ui->comboBoxNLDRdata->addItem(qstr);

//    int idex = ui->comboBoxNLDRdata->findText("No distance data in memory");
//    ui->comboBoxNLDRdata->removeItem(idex);

//    if (ui->comboBoxDISTmemdata_2->findText("File-distance") == -1)
//        ui->comboBoxDISTmemdata_2->addItem(qstr);

//    idex = ui->comboBoxDISTmemdata_2->findText("No distance data in memory");
//    ui->comboBoxDISTmemdata_2->removeItem(idex);

//    cout << "Successfully read distance file.\n\n";
//}

//void TreeScaper::on_pushDISTcleardata_clicked()
//{
//    QString qtmemorydata = ui->comboBoxDISTmemdata->currentText();
//    string stdmemorydata = qtmemorydata.toStdString();
//    String memorydata(stdmemorydata.c_str());

//    if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
//    {
//        int idex = ui->comboBoxDISTmemdata->findText(qtmemorydata);
//        ui->comboBoxDISTmemdata->removeItem(idex);
//        TreesData->deletetrees();
//    }
//    else if(memorydata == (String) "No data in memory")
//    {
//        cout << "warning: There is no data in memory." << endl;
//        return;
//    }
//    else
//    {
//        int idex = ui->comboBoxDISTmemdata->findText(qtmemorydata);
//        ui->comboBoxDISTmemdata->removeItem(idex);
//        idex = ui->comboBoxNLDRdata->findText(qtmemorydata);
//        ui->comboBoxNLDRdata->removeItem(idex);
//        idex = ui->comboBoxDISTmemdata_2->findText(qtmemorydata);
//        ui->comboBoxDISTmemdata_2->removeItem(idex);
//        idex = ui->comboBoxDISTmemdata_3->findText(qtmemorydata);
//        ui->comboBoxDISTmemdata_3->removeItem(idex);
//        idex = ui->comboBoxCOVAmemdata->findText(qtmemorydata);
//        ui->comboBoxCOVAmemdata->removeItem(idex);

//        TreesData->delete_matrix(memorydata);
//    }

//    if(ui->comboBoxDISTmemdata->count() == 0)
//    {
//        QString qstr("No data in memory");
//        ui->comboBoxDISTmemdata->addItem(qstr);
//        ui->listDATAmem->addItem(qstr);
//    }
//    if(ui->comboBoxNLDRdata->count() == 0)
//    {
//        QString qstr("No distance data in memory");
//        ui->comboBoxNLDRdata->addItem(qstr);
//    }
//}

void TreeScaper::on_pushDISTaffinity_clicked()
{
    QString qtmemorydata = ui->comboBoxDISTmemdata_2->currentText();
    string stdmemorydata = qtmemorydata.toStdString();
    String memorydata(stdmemorydata.c_str());

    QString qtmethod = ui->comboBoxDISTaffinity->currentText();
    string stdmethod = qtmethod.toStdString();
    String method(stdmethod.c_str());
    int type;
    if(method == (String) "Reciprocal")
        type = 1;
    else if(method == (String) "Exponential")
        type = 2;

    if(memorydata == (String) "No distance data in memory" && memorydata != (String) "Unweighted RF-distance"
            && memorydata != (String) "Weighted RF-distance" && memorydata != (String)"Matching-distance"
            && memorydata != (String) "SPR-distance" && memorydata != (String)"File-distance")
    {
        cout << "Warning: No Distance file in memory! Please compute first.\n\n";
        return;
    }
    else
    {
        TreesData->Compute_Affinity_dist(memorydata, type);
        if (memorydata == (String)"Unweighted RF-distance")
        {
            QString qstr("Affinity-URF");

            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxDISTmemdata_3->findText("Affinity-URF") == -1)
                ui->comboBoxDISTmemdata_3->addItem(qstr);
        }
        else if (memorydata == (String)"Weighted RF-distance")
        {
            QString qstr("Affinity-RF");

            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxDISTmemdata_3->findText("Affinity-RF") == -1)
                ui->comboBoxDISTmemdata_3->addItem(qstr);

        }
        else if (memorydata == (String)"Matching-distance")
        {
            QString qstr("Affinity-match");

            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxDISTmemdata_3->findText("Affinity-match") == -1)
                ui->comboBoxDISTmemdata_3->addItem(qstr);

        }
        else if (memorydata == (String)"SPR-distance")
        {
            QString qstr("Affinity-SPR");

            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxDISTmemdata_3->findText("Affinity-SPR") == -1)
                ui->comboBoxDISTmemdata_3->addItem(qstr);
        }
        else if (memorydata == (String)"File-distance")
        {
            QString qstr("Affinity-filedist");

            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item1.empty())
                ui->listDATAmem->addItem(qstr);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item2.empty())
                ui->listDATAmem_2->addItem(qstr);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
            if(item3.empty())
                ui->listDATAmem_3->addItem(qstr);

            if (ui->comboBoxDISTmemdata_3->findText("Affinity-filedist") == -1)
                ui->comboBoxDISTmemdata_3->addItem(qstr);
        }

        int idex = ui->comboBoxDISTmemdata_3->findText("No affinity data in memory");
        ui->comboBoxDISTmemdata_3->removeItem(idex);


        cout << "Successfully computed affinity matrix.\n\n";
    }
}

//void TreeScaper::on_pushDISToutput_clicked()
//{
//    QString qtmemorydata = ui->comboBoxDISTmemdata->currentText();
//    string stdmemorydata = qtmemorydata.toStdString();
//    String memorydata(stdmemorydata.c_str());
//    if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
//    {
//        QString qtfname = ui->textDATAfile->toPlainText();
//        string stdfname = qtfname.toStdString();
//        int reply = QMessageBox::question(this, "output trees", "Choose output tree format", "Newick", "Nexus");

//        QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
//        if(reply == 0) // newick
//            qtconvert = "Newick";
//        else if(reply == 1)
//            qtconvert = "Nexus";

//        string stdconvert = qtconvert.toStdString();
//        String convert(stdconvert.c_str());
//        string outName = Trees::WriteTreesFilename(stdfname, stdconvert);
//        if(!TreesData->treesAreexisting())
//        {
//            cout << "Warning: there is no tree in the memory!" << endl;
//            return;
//        }

//        if(convert == (String) "Newick")
//        {
//            TreesData->WriteTrees(outName, NEWICK);
//            cout << "successfully outputted Newick format trees to file: " << outName << endl;
//        }
//        else
//        if(convert == (String) "Nexus")
//        {
//            TreesData->WriteTrees(outName, NEXUS);
//            cout << "successfully outputted NEXUS format trees to file: " << outName << endl;
//        }
//    }
//    else if (memorydata == (String) "Bipartition Matrix")
//    {
//        QString qtfname = ui->textDATAfile->toPlainText();
//        string stdfname = qtfname.toStdString();

//        int reply = QMessageBox::question(this, "output bipartition matrix", "Sparse matrix: Choose output format", "List", "Matrix");

//        QString qtformat;
//        if(reply == 0) // List
//            qtformat = "List format";
//        else if(reply == 1) // Matrix
//            qtformat = "Matrix format";

//        string stdformat = qtformat.toStdString();
//        String format(stdformat.c_str());

//        string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname);

//        ofstream outBipartMatrix;
//        outBipartMatrix.open(namebipartmatrix.c_str());

//        if(format == (String) "List format")
//        {
//            cout << "Output Bipartition matrix in list format to " << namebipartmatrix << endl;
//            TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
//        } else
//        if(format == (String) "Matrix format")
//        {
//            cout << "Output Bipartition matrix in matrix format to " << namebipartmatrix << endl;
//            TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
//        }
//    }
//    else
//    {
//        TreesData->print_matrix(memorydata);
//        cout << "Successfully printed " << memorydata << " !" << endl;
//    }

//}


void TreeScaper::on_pushDISTcomm_clicked()
{
    if(communitythd.isRunning())
    {
        QMessageBox::warning(NULL, "warning", "Community detection is running!", 0, 0);
        return;
    }
    QString qtmemorydata = ui->comboBoxDISTmemdata_3->currentText();
    string stdmemorydata = qtmemorydata.toStdString();
    String memorydata(stdmemorydata.c_str());

//    if(memorydata == (String) "No affinity data in memory" || memorydata != (String) "Affinity-URF"
//            || memorydata != (String) "Affinity-RF" || memorydata != (String) "Affinity-match"
//            || memorydata != (String) "Affinity-SPR" || memorydata != (String) "Affinity-geodesic"
//            || memorydata != (String) "Affinity-filedist")
//    {
//        cout << "Warning: Cannot detect community!" << endl;
//        return;
//    }

    QString qtcdmodel = ui->comboBoxDISTcdmodel->currentText();
    string stdcdmodel = qtcdmodel.toStdString();
    String cdmodel(stdcdmodel.c_str());
    int modelType;
    if(cdmodel == (String) "Configuration Null Model")
        modelType = 3;
    else if(cdmodel == (String) "Constant Potts Model")
            modelType = 4;
    else if(cdmodel == (String) "Erdos-Renyi Null Model")
            modelType = 2;
    else if(cdmodel == (String) "No Null Model")
            modelType = 1;

    QString qtstartlambda = ui->textDISTLstart->toPlainText();
    string stdstartlambda = qtstartlambda.toStdString();
    double startlambda = atof(stdstartlambda.c_str());
    QString qtendlambda = ui->textDISTLend->toPlainText();
    string stdendlambda = qtendlambda.toStdString();
    double endlambda = atof(stdendlambda.c_str());
    QString qtintervallambda = ui->textDISTLinterval->toPlainText();
    string stdintervallambda = qtintervallambda.toStdString();
    double intervallambda = atof(stdintervallambda.c_str());
    Array<double> paramfixed(1);
    paramfixed[0] = 0;
    int lambdasize = (endlambda - startlambda) / intervallambda + 1;
    Array<double> paramnonfixed(lambdasize);
    for(int i = 0; i < lambdasize; i++)
        paramnonfixed[i] = startlambda + i * intervallambda;

    QString qtparam3 = ui->textDISThighfreq->toPlainText();
    string highfreq = qtparam3.toStdString();
    QString qtparam4 = ui->textDISTlowfreq->toPlainText();
    string lowfreq = qtparam4.toStdString();


    if(ui->radioButtonDISTauto->isChecked())
        communitythd.initialization(TreesData, memorydata, modelType, paramnonfixed, paramfixed, highfreq, lowfreq, 1);
    else
    if(ui->radioButtonDISTmanu->isChecked())
        communitythd.initialization(TreesData, memorydata, modelType, paramnonfixed, paramfixed, highfreq, lowfreq, 1, false);
    else
        cout << "Please specify a method to find the plateaus!\n\n";

    communitythd.start();
}

//void TreeScaper::on_comboBoxDISTmemdata_currentIndexChanged(const QString &textstr)
//{
//    string stdtext = textstr.toStdString();
//    String text(stdtext.c_str());

//    if (text == (String) "")
//        return;

////    if (text == (String) "No data in memory")
////    {
////        ui->pushDISTcleardata->setEnabled(false);
////        ui->pushDISToutput->setEnabled(false);
////    }
////    else if (text != (String) "No data in memory")
////    {
////        ui->pushDISTcleardata->setEnabled(true);
////        ui->pushDISToutput->setEnabled(true);
////    }
//}

void TreeScaper::on_pushNLDRplottrees_clicked()
{
    if(TreesData->Get_n_trees() == 0)
    {
        cout << "There are no trees in memory!\n\n";
        return;
    }

/*
    int idx = 1;
    stringstream s;
    s << idx;
    string str = s.str();
    String title(str.c_str());
    title += "-th tree";
    plottreethd.initialization(paras::imageparas,
                               TreesData->get_tree(idx),
                               TreesData->Get_isrooted(),
                               TreesData->Get_isweighted(),
                               TreesData->Get_labelmap(),
                               title);
    plottreethd.start();
    return;*/

    for(int i = 0; TreesData->selected_trees[i] != -1; i++)
    {
        cout << "id:" << TreesData->selected_trees[i] << endl;
        stringstream s;
        s << TreesData->selected_trees[i];
        string str = s.str();
        String Strtitle(str.c_str());
        Strtitle += "-th tree";
        char *title = new char [50];
        strcpy(title, Strtitle);
        emit sendplottree(TreesData->get_tree(TreesData->selected_trees[i]), title,
                          TreesData->Get_isrooted(), TreesData->Get_isweighted(), TreesData->Get_labelmap());//--TreesData->selected_trees[i]);
    }
}

void TreeScaper::plottree(const NEWICKTREE *tree, char *title, bool isrooted, bool isweighted, const LabelMap *lm)
{
    if(tree == NULL)
    {
        cout << "Error: Tree does not exist! Unable to plot it.\n\n";
        return;
    }
    plottreethd.initialization(paras::imageparas,
                               tree,
                               isrooted,
                               isweighted,
                               lm,
                               (String) title);
    delete [] title;
    plottreethd.start();
}

void TreeScaper::on_pushButtonDATAcon_clicked()
{
    if(!TreesData->bipartmatrixIsexisting())
    {
        cout << "Warning: There is no bipartition matrix in the memory! Please compute it first.\n\n";
        return;
    }

    if(TreesData->compute_consensus_tree(MAJORITYTREE, ""))
    {
        cout << "Successfully computed the majority consensus tree!\n\n";

        QString qstr("Consensus Tree");

        QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item1.empty())
            ui->listDATAmem->addItem(qstr);

        QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item2.empty())
            ui->listDATAmem_2->addItem(qstr);

        QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
        if(item3.empty())
            ui->listDATAmem_3->addItem(qstr);
    }
}

void TreeScaper::on_pushButtonDATAconplot_clicked()
{
    String Strtitle = "consensus tree";

    char *title = new char [50];
    strcpy(title, Strtitle);

    emit sendplottree(TreesData->get_contree(0), title,
                      TreesData->Get_isrooted(), true, TreesData->Get_labelmap());
}

//void TreeScaper::on_pushCONVbipartfreqCumulative_clicked()
//{
//    if(!TreesData->bipartmatrixIsexisting())
//    {
//        cout << "Warning: There is no bipartition matrix in the memory! Please compute it first." << endl;
//        return;
//    }

//    double** bipartFreq = NULL;
//    int* bipartFreqIdx = NULL;
//    int nt = TreesData->Get_n_trees();
//    int Unique_idx = TreesData->Get_treecov_size();
//    double xstart = round(nt*0.1);

//    QString qtnumbipart = ui->TextCONVnumbipartsCumulative->toPlainText();
//    QString qtnumpts = ui->TextCONVnumptsCumulative->toPlainText();
//    int splits = qtnumbipart.toInt();
//    int increments = qtnumpts.toInt();

//    if(splits > Unique_idx)
//    {
//        cout << "The number of main bipartitions cannot be larger than the number of unique bipartitions!" << endl;
//        return;
//    }

//    if(increments > nt)
//    {
//        cout << "The number of intervals cannot be larger than the number of trees!" << endl;
//        return;
//    }

//    if(increments > (nt-xstart))
//    {
//        cout << "Must account for 10% burn-in when choosing number of intervals." << endl;
//        return;
//    }

//    bipartFreq = new double*[splits];
//    for(int i = 0; i < splits; i++)
//        bipartFreq[i] = new double[increments];

//    for(int i = 0; i < splits; i++)
//        for(int j = 0; j < increments; j++)
//            bipartFreq[i][j] = 0;

//    bipartFreqIdx = new int[splits];
//    for(int i = 0; i < splits; i++)
//        bipartFreqIdx[i] = 0;

//    TreesData->Compute_Cumulative(bipartFreq, bipartFreqIdx, splits, increments);

//    std::cout << "Done!" << endl;
////    // output
////    cout << "Cumuative sum for max " << splits << " bipart freq: " << endl;
////    for(int i = 0; i < splits; i++)
////    {
////        for(int j = 0; j < increments; j++)
////            cout << bipartFreq[i][j] << " ";
////        cout << endl;
////    }
////    cout << endl;


//    double xend = nt;
//    double ymax = 0;
//    for(int i = 0; i < increments; i++)
//    {
//        if(ymax < bipartFreq[0][i])
//            ymax = bipartFreq[0][i];
//    }

//    Plot2D_AWTY p2d("Number of Trees", "Relative Frequency", xstart, xend, -0.05 * ymax, 1.05 * ymax, "Cumulative", this);
//    QVector<double> x(increments), y(increments);
//    for(int i = 0; i < increments; i++)
//        x[i] = xstart + round((nt-xstart)/increments)*(i+1);


//    for(int j = 0; j < splits; j++)
//    {
//        for(int i = 0; i < increments; i++)
//            y[i] = bipartFreq[j][i];

//        QString bipartLegend = QString("Bipartition %1").arg(bipartFreqIdx[j] + 1);
//        p2d.Plot2DCurves_AWTY(x, y, bipartLegend, QColor(qSin(j*0.3)*100+100, qSin(j*0.6+0.7)*100+100, qSin(j*0.4+0.6)*100+100));
//    }

//    p2d.show();
//    p2d.exec();

//    for(int i = 0; i < splits; i++)
//        delete [] bipartFreq[i];
//    delete [] bipartFreq;
//    bipartFreq = NULL;
//    delete [] bipartFreqIdx;
//    bipartFreqIdx = NULL;
//}

//void TreeScaper::on_pushCONVbipartfreqSlide_clicked()
//{
//    if(!TreesData->bipartmatrixIsexisting())
//    {
//        cout << "Warning: there is no bipartition matrix in the memory! Please compute it first." << endl;
//        return;
//    }

//    double** bipartFreq = NULL;
//    int* bipartFreqIdx = NULL;
//    int nt = TreesData->Get_n_trees();
//    int Unique_idx = TreesData->Get_treecov_size();
//    double xstart = round(nt*0.1);

//    QString qtnumbipart = ui->TextCONVnumbipartsSlide->toPlainText();
//    QString qtnumpts = ui->TextCONVnumptsSlide->toPlainText();
//    int splits = qtnumbipart.toInt();
//    int increments = qtnumpts.toInt();

//    if(splits > Unique_idx)
//    {
//        cout << "The number of main bipartitions cannot be larger than the number of unique bipartitions!" << endl;
//        return;
//    }

//    if(increments > nt)
//    {
//        cout << "The number of intervals cannot be larger than the number of trees!" << endl;
//        return;
//    }

//    if(increments > (nt-xstart))
//    {
//        cout << "Must account for 10% burn-in when choosing number of intervals." << endl;
//        return;
//    }

//    bipartFreq = new double*[splits];
//    for(int i = 0; i < splits; i++)
//        bipartFreq[i] = new double[increments];

//    for(int i = 0; i < splits; i++)
//        for(int j = 0; j < increments; j++)
//            bipartFreq[i][j] = 0;

//    bipartFreqIdx = new int[splits];
//    for(int i = 0; i < splits; i++)
//        bipartFreqIdx[i] = 0;

//    TreesData->Compute_Slide(bipartFreq, bipartFreqIdx, splits, increments);

//    std::cout << "Done!" << endl;
////    // output
////    cout << "Cumuative sum for max " << splits << " bipart freq: " << endl;
////    for(int i = 0; i < splits; i++)
////    {
////        for(int j = 0; j < increments; j++)
////            cout << bipartFreq[i][j] << " ";
////        cout << endl;
////    }
////    cout << endl;


//    double xend = nt;
//    double ymax = 0;
//    for(int i = 0; i < increments; i++)
//    {
//        if(ymax < bipartFreq[0][i])
//            ymax = bipartFreq[0][i];
//    }

//    Plot2D_AWTY p2d("Number of Trees", "Frequency", xstart, xend, -0.05 * ymax, 1.05 * ymax, "Slide", this);
//    QVector<double> x(increments), y(increments);
//    for(int i = 0; i < increments; i++)
//        x[i] = xstart + round((nt-xstart)/increments)*(i+1);


//    for(int j = 0; j < splits; j++)
//    {
//        for(int i = 0; i < increments; i++)
//            y[i] = bipartFreq[j][i];

//        QString bipartLegend = QString("Bipartition %1").arg(bipartFreqIdx[j] + 1);
//        p2d.Plot2DCurves_AWTY(x, y, bipartLegend, QColor(qSin(j*0.3)*100+100, qSin(j*0.6+0.7)*100+100, qSin(j*0.4+0.6)*100+100));
//    }

//    p2d.show();
//    p2d.exec();

//    for(int i = 0; i < splits; i++)
//        delete [] bipartFreq[i];
//    delete [] bipartFreq;
//    bipartFreq = NULL;
//    delete [] bipartFreqIdx;
//    bipartFreqIdx = NULL;
//}

void TreeScaper::on_pushButtonDATAloadlist_clicked()
{
    if(ui->radioButtonDATArange->isChecked())
    {
        if(TreesData->Get_n_trees() <= 0)
        {
            cout << "No trees in the memory!\n\n";
            return;
        }

        int s = 0;
        int t = TreesData->Get_n_trees() - 1;

        QString qtstart = ui->textDATAidx1->toPlainText();
        string stdstart = qtstart.toStdString();
        String start(stdstart.c_str());
        if(start != (String) "start")
            s = atoi(stdstart.c_str());

        QString qtend = ui->textDATAidx2->toPlainText();
        string stdend = qtend.toStdString();
        String end(stdend.c_str());
        if(end != (String) "end")
            t = atoi(stdend.c_str());
        if(t < s)
            cout << "Warning: Initial index should be smaller than or equal to the final index!\n\n";
        if(s < 0)
            cout << "Warning: Initial index should be at least 0!\n\n";
        if(t >= TreesData->Get_n_trees())
            cout << "Warning: The final index should be smaller than the number of trees - 1!\n\n";
        if(t < s || s < 0 || t >= TreesData->Get_n_trees())
        {
            cout << "Use all trees" << endl;
            s = 0;
            t = TreesData->Get_n_trees() - 1;
        }
        Array<int> *treeidx = TreesData->getidxlist();
        treeidx->resize(t - s + 1);
        for(int i = s; i < t + 1; i++)
            (*treeidx)[i - s] = i;
        cout << "Successfully loaded from " << s << " to " << t << ".\n\n";
    } else
    {
        QString qtfname = ui->textDATAlistfile->toPlainText();
        string stdfname = qtfname.toStdString();
        String fname(stdfname.c_str());
        File file(fname);
        if(! file.is_open())
        {
            cout << "Error: Cannot open the data file!\n\n";
            return;
        }
        file.seek(0);
        int num = file.lines();
        file.seek(0);

        Array<int> *treeidx = TreesData->getidxlist();
        treeidx->resize(num);

        for(int i = 0; i < num; i++)
        {
            file >> (*treeidx)[i];
        }
        cout << "Successfully loaded from file: " << fname << "; " << num << " trees.\n\n";
    }
}

void TreeScaper::on_radioButtonDATAindex_clicked()
{
    ui->labelDATAfrom->setEnabled(false);
    ui->labelDATAto->setEnabled(false);
    ui->textDATAidx1->setEnabled(false);
    ui->textDATAidx2->setEnabled(false);

    ui->labelDATAlistname->setEnabled(true);
    ui->textDATAlistfile->setEnabled(true);
    ui->pushDATAButtonlistbrowse->setEnabled(true);
}

void TreeScaper::on_radioButtonDATArange_clicked()
{
    ui->labelDATAfrom->setEnabled(true);
    ui->labelDATAto->setEnabled(true);
    ui->textDATAidx1->setEnabled(true);
    ui->textDATAidx2->setEnabled(true);

    ui->labelDATAlistname->setEnabled(false);
    ui->textDATAlistfile->setEnabled(false);
    ui->pushDATAButtonlistbrowse->setEnabled(false);
}

void TreeScaper::on_pushDATAButtonlistbrowse_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
//    fileDialog->setFilter(tr("Files(*)"));
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];
        ui->textDATAlistfile->setPlainText(path);
    }
}

void TreeScaper::on_pushNLDRsaveindices_clicked()
{
    cout << "Output selected indices to the file: " << TreesData->Print_selected_indices() << "\n\n";
}

void TreeScaper::on_pushDISTplotcom_clicked()
{
    double **com_info = NULL;
    int length = 0;
    TreesData->Get_community_info(com_info, length);
    length--;
    if(length == 0)
    {
        cout << "Warning: No communities in memory!\n\n";
        return;
    }

    double xstart = com_info[1][1] - 0.05 * (com_info[1][length] - com_info[1][1]);
    double xend = com_info[1][length] + 0.05 * (com_info[1][length] - com_info[1][1]);
    double ymax = 0;
    for(int i = 0; i < length; i++)
    {
        if(ymax < com_info[0][i + 1])
            ymax = com_info[0][i + 1];
        if(ymax < com_info[2][i + 1])
            ymax = com_info[2][i + 1];
        if(ymax < com_info[3][i + 1])
            ymax = com_info[3][i + 1];
    }

    Plot2D p2d("lambda", "values", xstart, xend, -0.05 * ymax, 1.05 * ymax, "Community information for affinity matrix", this);
    QVector<double> x(length), y(length);
    for(int i = 0; i < length; i++)
        x[i] = com_info[1][i + 1];

    for(int i = 0; i < length; i++)
        y[i] = com_info[0][i + 1];
    p2d.Plot2DCurves(x, y, "Labels for communities", QColor(255, 0, 0));
    for(int i = 0; i < length; i++)
        y[i] = com_info[2][i + 1];
    p2d.Plot2DCurves(x, y, "Number of communities", QColor(0, 255, 0));
    for(int i = 0; i < length; i++)
        y[i] = com_info[3][i + 1];
    p2d.Plot2DCurves(x, y, "Modularities", QColor(0, 0, 255));

    connect(&p2d, SIGNAL(sendnewlambda(double, Plot2D *, int)), this, SLOT(Getnewlambda(double, Plot2D *, int)));

    p2d.show();
    p2d.exec();
}

void TreeScaper::on_pushCOVAcomm_clicked()
{
    if(communitythd.isRunning())
    {
        QMessageBox::warning(NULL, "Warning", "Community detection is running!", 0, 0);
        return;
    }

    QString qtcdmodel = ui->comboBoxCOVAcdmodel->currentText();
    string stdcdmodel = qtcdmodel.toStdString();
    String cdmodel(stdcdmodel.c_str());
    int modelType;
    if(cdmodel == (String) "Configuration Null Model")
        modelType = 3;
    else if(cdmodel == (String) "Constant Potts Model")
            modelType = 4;
    else if(cdmodel == (String) "Erdos-Renyi Null Model")
            modelType = 2;
    else if(cdmodel == (String) "No Null Model")
            modelType = 1;

    QString qtdata = ui->comboBoxCOVAmemdata->currentText();
    string stddata = qtdata.toStdString();
    String datastr(stddata.c_str());

    QString qtfixedlambda = ui->textCOVALfixed->toPlainText();
    double fixedlambda = qtfixedlambda.toDouble();
    QString qtstartlambda = ui->textCOVALstart->toPlainText();
    double startlambda = qtstartlambda.toDouble();
    QString qtendlambda = ui->textCOVALend->toPlainText();
    double endlambda = qtendlambda.toDouble();
    QString qtintervallambda = ui->textCOVALinterval->toPlainText();
    double intervallambda = qtintervallambda.toDouble();
    Array<double> paramfixed(1);
    paramfixed[0] = fixedlambda;
    int lambdasize = (endlambda - startlambda) / intervallambda + 1;
    Array<double> paramnonfixed(lambdasize);
    for(int i = 0; i < lambdasize; i++)
        paramnonfixed[i] = startlambda + i * intervallambda;

    QString qtparam3 = ui->textCOVAhighfreq->toPlainText();
    string highfreq = qtparam3.toStdString();
    QString qtparam4 = ui->textCOVAlowfreq->toPlainText();
    string lowfreq = qtparam4.toStdString();

    if(ui->radioButtonCOVAauto->isChecked())
        communitythd.initialization(TreesData, datastr, modelType, paramfixed, paramnonfixed, highfreq, lowfreq, 2);
    else
    if(ui->radioButtonCOVAmanu->isChecked())
    {
        if(ui->radioButtonCOVALp->isChecked())
        {
            communitythd.initialization(TreesData, datastr, modelType, paramfixed, paramnonfixed, highfreq, lowfreq, 2, false);
        } else
        {
            communitythd.initialization(TreesData, datastr, modelType, paramnonfixed, paramfixed, highfreq, lowfreq, 2, false);
        }
    } else
    {
        cout << "Please specify a method to find the plateaus!\n\n";
    }

    communitythd.start();
}

//void TreeScaper::on_pushCONSclear_clicked()
//{
//    File fout("log.txt");
//    fout.clean();
//    fout.close();
//    ui->textBrowserDATAlog->setText("");
////    ui->textBrowserCONVlog->setText("");
//    ui->textBrowserCOVAlog->setText("");
//    ui->textBrowserdimlog->setText("");
//    ui->textBrowserNLDRlog->setText("");
//}

//void TreeScaper::on_pushCONVclear_clicked()
//{
//    File fout("log.txt");
//    fout.clean();
//    fout.close();
//    ui->textBrowserDATAlog->setText("");
//    ui->textBrowserCONSlog->setText("");
//    ui->textBrowserCOVAlog->setText("");
//    ui->textBrowserdimlog->setText("");
//    ui->textBrowserNLDRlog->setText("");
//}

//void TreeScaper::on_pushButtonCOVAlog_clicked()
//{
//    File fout("log.txt");
//    fout.clean();
//    fout.close();
//    char *text = "";
//    ui->textBrowserCOVAlog->setText(text);
//}

void TreeScaper::on_pushCOVAplotcom_clicked()
{
    double **com_info = NULL;
    int length = 0;
    TreesData->Get_community_info(com_info, length);
    length--;
    if(length == 0)
    {
        cout << "Warning: No communities in memory!" << std::endl << endl;
        return;
    }

    double xstart = com_info[1][1] - 0.05 * (com_info[1][length] - com_info[1][1]);
    double xend = com_info[1][length] + 0.05 * (com_info[1][length] - com_info[1][1]);
    double ymax = 0;
    for(int i = 0; i < length; i++)
    {
        if(ymax < com_info[0][i + 1])
            ymax = com_info[0][i + 1];
        if(ymax < com_info[2][i + 1])
            ymax = com_info[2][i + 1];
        if(ymax < com_info[3][i + 1])
            ymax = com_info[3][i + 1];
    }

    Plot2D p2d("lambda", "values", xstart, xend, -0.05 * ymax, 1.05 * ymax, "Community information for convariance matrix", this);
    QVector<double> x(length), y(length);
    for(int i = 0; i < length; i++)
        x[i] = com_info[1][i + 1];

    for(int i = 0; i < length; i++)
        y[i] = com_info[0][i + 1];
    p2d.Plot2DCurves(x, y, "Labels for communities", QColor(255, 0, 0));
    for(int i = 0; i < length; i++)
        y[i] = com_info[2][i + 1];
    p2d.Plot2DCurves(x, y, "Number of communities", QColor(0, 255, 0));
    for(int i = 0; i < length; i++)
        y[i] = com_info[3][i + 1];
    p2d.Plot2DCurves(x, y, "Modularities", QColor(0, 0, 255));

    connect(&p2d, SIGNAL(sendnewlambda(double, Plot2D *, int)), this, SLOT(Getnewlambda(double, Plot2D *, int)));

    p2d.show();
    p2d.exec();
}

void TreeScaper::Getnewlambda(double newlambda, Plot2D *p2d, int plottype)
{
    if(plottype == 1)
    {
        cout << "Reload the community information of the covariance matrix!\n\n";
        QString qtcdmodel = ui->comboBoxCOVAcdmodel->currentText();
        string stdcdmodel = qtcdmodel.toStdString();
        String cdmodel(stdcdmodel.c_str());
        int modelType;
        if(cdmodel == (String) "Configuration Null Model")
            modelType = 3;
        else if(cdmodel == (String) "Constant Potts Model")
                modelType = 4;
        else if(cdmodel == (String) "Erdos-Renyi Null Model")
                modelType = 2;
        else if(cdmodel == (String) "No Null Model")
                modelType = 1;

        QString qtdata = ui->comboBoxCOVAmemdata->currentText();
        string stddata = qtdata.toStdString();
        String datastr(stddata.c_str());

        QString qtparam3 = ui->textCOVAhighfreq->toPlainText();
        string highfreq = qtparam3.toStdString();
        QString qtparam4 = ui->textCOVAlowfreq->toPlainText();
        string lowfreq = qtparam4.toStdString();

        if(ui->radioButtonCOVAauto->isChecked())
            TreesData->compute_community_fixedlambda(datastr, modelType, newlambda, 0, highfreq, lowfreq);

        QString qtfixedlambda = ui->textCOVALfixed->toPlainText();
        string stdfixedlambda = qtfixedlambda.toStdString();
        double fixedlambda = atof(stdfixedlambda.c_str());

        if(ui->radioButtonCOVAmanu->isChecked())
        {
            if(ui->radioButtonCOVALp->isChecked())
                TreesData->compute_community_fixedlambda(datastr, modelType, fixedlambda, newlambda, highfreq, lowfreq);
            else
                TreesData->compute_community_fixedlambda(datastr, modelType, newlambda, fixedlambda, highfreq, lowfreq);
        }
    } else
    if(plottype == 2)
    {
        cout << "Reload the community information of the affinity matrix!\n\n";
        QString qtmemorydata = ui->comboBoxDISTmemdata_3->currentText();
        string stdmemorydata = qtmemorydata.toStdString();
        String memorydata(stdmemorydata.c_str());

        QString qtcdmodel = ui->comboBoxDISTcdmodel->currentText();
        string stdcdmodel = qtcdmodel.toStdString();
        String cdmodel(stdcdmodel.c_str());
        int modelType;
        if(cdmodel == (String) "Configuration Null Model")
            modelType = 3;
        else if(cdmodel == (String) "Constant Potts Model")
                modelType = 4;
        else if(cdmodel == (String) "Erdos-Renyi Null Model")
                modelType = 2;
        else if(cdmodel == (String) "No Null Model")
                modelType = 1;

        QString qtparam3 = ui->textDISThighfreq->toPlainText();
        string highfreq = qtparam3.toStdString();
        QString qtparam4 = ui->textDISTlowfreq->toPlainText();
        string lowfreq = qtparam4.toStdString();

        TreesData->compute_community_fixedlambda(memorydata, modelType, newlambda, 0, highfreq, lowfreq);
    } else
    {
        cout << "Warning: Incorrect type: " << plottype << "\n\n";
        return;
    }
    p2d->Clear2DCurves();

    double **com_info = NULL;
    int length = 0;
    TreesData->Get_community_info(com_info, length);
    length--;
    if(length == 0)
    {
        cout << "Warning: No communities in memory!\n\n";
        return;
    }

    QVector<double> x(length), y(length);
    for(int i = 0; i < length; i++)
        x[i] = com_info[1][i + 1];

    for(int i = 0; i < length; i++)
        y[i] = com_info[0][i + 1];
    p2d->Plot2DCurves(x, y, "Labels for communities", QColor(255, 0, 0));
    for(int i = 0; i < length; i++)
        y[i] = com_info[2][i + 1];
    p2d->Plot2DCurves(x, y, "Number of communities", QColor(0, 255, 0));
    for(int i = 0; i < length; i++)
        y[i] = com_info[3][i + 1];
    p2d->Plot2DCurves(x, y, "Modularities", QColor(0, 0, 255));
}

void TreeScaper::on_comboBoxDISTmemdata_2_currentIndexChanged(const QString &textstr)
{
    string stdtext = textstr.toStdString();
    String text(stdtext.c_str());

    if(ui->comboBoxDISTmemdata_2->count() == 0)
    {
        QString qstr("No distance data in memory");
        ui->comboBoxDISTmemdata_2->addItem(qstr);
    }

    if (text == (String) "")
        return;

    if (text == (String) "Unweighted RF-distance" || (text == (String) "Weighted RF-distance")
            || (text == (String) "Matching-distance") || (text == (String) "SPR-distance")
            || (text == (String) "Geodesic-distance") || (text == (String) "File-distance"))
    {
        ui->pushDISTaffinity->setEnabled(true);
    }
    else
        ui->pushDISTaffinity->setEnabled(false);
}

void TreeScaper::on_comboBoxCOVAmemdata_currentIndexChanged(const QString &textstr)
{
    if(ui->comboBoxCOVAmemdata->count() == 0)
    {
        QString qstr("No covariance data in memory");
        ui->comboBoxCOVAmemdata->addItem(qstr);
        ui->pushCOVAcomm->setEnabled(false);
    }
}

void TreeScaper::on_comboBoxNLDRdata_currentIndexChanged(const QString &arg1)
{
    if(ui->comboBoxNLDRdata->count() == 0)
    {
        QString qstr("No distance data in memory");
        ui->comboBoxNLDRdata->addItem(qstr);
    }
}

void TreeScaper::on_pushButtonCOVAbrowse_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];
        ui->textCOVAfile->setPlainText(path);
    }
}

void TreeScaper::on_pushCOVAloaddist_clicked()
{
    QString qtfname = ui->textCOVAfile->toPlainText();
    string stdfname = qtfname.toStdString();
    String fname(stdfname.c_str());
    File file(fname);
    if(! file.is_open())
    {
        cout << "Error: Cannot open the covariance file!\n\n";
        return;
    }
    TreesData->load_covariancefile(stdfname);

    QString qstr("File-covariance");
    if (ui->comboBoxCOVAmemdata->findText("File-covariance") == -1)
        ui->comboBoxCOVAmemdata->addItem(qstr);

    QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item1.empty())
        ui->listDATAmem->addItem(qstr);

    QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item2.empty())
        ui->listDATAmem_2->addItem(qstr);

    QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item3.empty())
        ui->listDATAmem_3->addItem(qstr);

    int idex = ui->comboBoxCOVAmemdata->findText("No covariance data in memory");
    ui->comboBoxCOVAmemdata->removeItem(idex);

    QList<QListWidgetItem *> itemTmp1 = ui->listDATAmem->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
    for (int i = 0; i < itemTmp1.size(); i++)
        delete ui->listDATAmem->takeItem(ui->listDATAmem->row(itemTmp1[i]));

    QList<QListWidgetItem *> itemTmp2 = ui->listDATAmem_2->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
    for (int i = 0; i < itemTmp2.size(); i++)
        delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(itemTmp2[i]));

    QList<QListWidgetItem *> itemTmp3 = ui->listDATAmem_3->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
    for (int i = 0; i < itemTmp3.size(); i++)
        delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(itemTmp3[i]));

    ui->pushCOVAcomm->setEnabled(true);

    cout << "Successfully read covariance file.\n\n";
}

void TreeScaper::on_pushCOVAclear_clicked()
{
    File fout("log.txt");
    fout.clean();
    fout.close();
    ui->textBrowserDATAlog->setText("");
//    ui->textBrowserCONSlog->setText("");
//    ui->textBrowserCONVlog->setText("");
    ui->textBrowserCOVAlog->setText("");
//    ui->textBrowserDISTlog->setText("");
//    ui->textBrowserdimlog->setText("");
    ui->textBrowserNLDRlog->setText("");
}

void TreeScaper::on_textBrowserDATAlog_textChanged()
{
    ui->textBrowserDATAlog->verticalScrollBar()->maximumHeight();
//    ui->textBrowserDATAlog->moveCursor(QTextCursor::End);
    /*
    QString text = ui->textBrowserDATAlog->toPlainText();
    ui->textBrowserCONSlog->setText(text);
    ui->textBrowserCONVlog->setText(text);
    ui->textBrowserCOVAlog->setText(text);
    ui->textBrowserDISTlog->setText(text);
    ui->textBrowserdimlog->setText(text);
    ui->textBrowserNLDRlog->setText(text);

    ui->textBrowserCONSlog->moveCursor(QTextCursor::End);
    ui->textBrowserCONVlog->setText(QTextCursor::End);
    ui->textBrowserCOVAlog->moveCursor(QTextCursor::End);
    ui->textBrowserDISTlog->moveCursor(QTextCursor::End);
    ui->textBrowserdimlog->moveCursor(QTextCursor::End);
    ui->textBrowserNLDRlog->moveCursor(QTextCursor::End);

    ui->textBrowserCONSlog->verticalScrollBar()->setValue(ui->textBrowserCONSlog->verticalScrollBar()->maximumHeight());
    ui->textBrowserCONVlog->verticalScrollBar()->setValue(ui->textBrowserCONVlog->verticalScrollBar()->maximumHeight());
    ui->textBrowserCOVAlog->verticalScrollBar()->setValue(ui->textBrowserCOVAlog->verticalScrollBar()->maximumHeight());
    ui->textBrowserDISTlog->verticalScrollBar()->setValue(ui->textBrowserDISTlog->verticalScrollBar()->maximumHeight());
    ui->textBrowserdimlog->verticalScrollBar()->setValue(ui->textBrowserdimlog->verticalScrollBar()->maximumHeight());
    ui->textBrowserNLDRlog->verticalScrollBar()->setValue(ui->textBrowserNLDRlog->verticalScrollBar()->maximumHeight());
*/
}

void TreeScaper::on_radioButtonCOVAauto_clicked()
{
    ui->labelCOVAlambda_2->setEnabled(false);
    ui->radioButtonCOVALp->setEnabled(false);
    ui->radioButtonCOVALn->setEnabled(false);
    ui->textCOVALfixed->setEnabled(false);
    ui->labelCOVAlambda->setEnabled(false);
    ui->labelCOVAfrom->setEnabled(false);
    ui->textCOVALstart->setEnabled(false);
    ui->labelCOVAto->setEnabled(false);
    ui->textCOVALend->setEnabled(false);
    ui->labelCOVAintv->setEnabled(false);
    ui->textCOVALinterval->setEnabled(false);
}

void TreeScaper::on_radioButtonCOVAmanu_clicked()
{
    ui->labelCOVAlambda_2->setEnabled(true);
    ui->radioButtonCOVALp->setEnabled(true);
    ui->radioButtonCOVALn->setEnabled(true);
    ui->textCOVALfixed->setEnabled(true);
    ui->labelCOVAlambda->setEnabled(true);
    ui->labelCOVAfrom->setEnabled(true);
    ui->textCOVALstart->setEnabled(true);
    ui->labelCOVAto->setEnabled(true);
    ui->textCOVALend->setEnabled(true);
    ui->labelCOVAintv->setEnabled(true);
    ui->textCOVALinterval->setEnabled(true);
}

void TreeScaper::on_radioButtonDISTauto_clicked()
{
    ui->labelDISTlambda_2->setEnabled(false);
    ui->labelDISTfrom->setEnabled(false);
    ui->textDISTLstart->setEnabled(false);
    ui->labelDISTto->setEnabled(false);
    ui->textDISTLend->setEnabled(false);
    ui->labelDISTto_2->setEnabled(false);
    ui->textDISTLinterval->setEnabled(false);
}

void TreeScaper::on_radioButtonDISTmanu_clicked()
{
    ui->labelDISTlambda_2->setEnabled(true);
    ui->labelDISTfrom->setEnabled(true);
    ui->textDISTLstart->setEnabled(true);
    ui->labelDISTto->setEnabled(true);
    ui->textDISTLend->setEnabled(true);
    ui->labelDISTto_2->setEnabled(true);
    ui->textDISTLinterval->setEnabled(true);
}

//void TreeScaper::buttonCONVcomboBox(int input)
//{
//    if(input == 0)
//    {
//        ui->stackedWidgetCONV->setCurrentIndex(0);
//    }
//    else if(input == 1)
//    {
//        ui->stackedWidgetCONV->setCurrentIndex(1);
//    }
//    else if(input == 2)
//    {
//        ui->stackedWidgetCONV->setCurrentIndex(2);
//    }
//}

void TreeScaper::comboBoxTreeData(int index)
{
    if(index == 0)
        ui->stackedWidgetTreeData->setCurrentIndex(0);
    else if(index == 1)
        ui->stackedWidgetTreeData->setCurrentIndex(1);
    else if(index == 2)
        ui->stackedWidgetTreeData->setCurrentIndex(2);
    else if(index == 3)
        ui->stackedWidgetTreeData->setCurrentIndex(3);
    else if(index == 4)
        ui->stackedWidgetTreeData->setCurrentIndex(4);
    else if(index == 5)
        ui->stackedWidgetTreeData->setCurrentIndex(5);
    else if(index == 6)
        ui->stackedWidgetTreeData->setCurrentIndex(6);
}

void TreeScaper::comboBoxDimNLDR(int index)
{
    if(index == 0)
        ui->stackedWidgetDimNLDR->setCurrentIndex(0);
    else if(index == 1)
        ui->stackedWidgetDimNLDR->setCurrentIndex(1);
}

void TreeScaper::comboBoxCD(int index)
{
    if(index == 0)
        ui->stackedWidgetCD->setCurrentIndex(0);
    else if(index == 1)
        ui->stackedWidgetCD->setCurrentIndex(1);
}

void TreeScaper::on_comboBoxDISTmemdata_3_currentIndexChanged(const QString &textstr)
{
    string stdtext = textstr.toStdString();
    String text(stdtext.c_str());

    if(ui->comboBoxDISTmemdata_3->count() == 0)
    {
        QString qstr("No affinity data in memory");
        ui->comboBoxDISTmemdata_3->addItem(qstr);
    }

    if (text == (String) "")
        return;

    if ((text == (String) "Affinity-URF")
            || (text == (String) "Affinity-RF") || (text == (String) "Affinity-match")
            || (text == (String) "Affinity-SPR") || (text == (String) "Affinity-geodesic")
            || (text == (String) "Affinity-filedist") || (text == (String) "File-affinity"))
    {
        ui->pushDISTcomm->setEnabled(true);
    }
    else
    {
        ui->pushDISTcomm->setEnabled(false);
    }
}

void TreeScaper::on_listDATAmem_itemActivated(QListWidgetItem *item)
{
    int msg1 = QMessageBox::question(this, "Data option", "Do you want to delete or output the data value?", "Output", "Delete", "Cancel");

    QString qtmemorydata = item->text();
    string stdmemorydata = qtmemorydata.toStdString();
    String memorydata(stdmemorydata.c_str());
    if(msg1 == 0) // output
    {
        if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "Output trees", "Choose output tree format", "Newick", "Nexus", "Cancel");

            QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
            if(reply == 0) // newick
                qtconvert = "Newick";
            else if(reply == 1)
                qtconvert = "Nexus";
            else if(reply == 2) // cancel
                return;

            string stdconvert = qtconvert.toStdString();
            String convert(stdconvert.c_str());
            string outName = TreesData->WriteTreesFilename(stdfname, stdconvert);
            if(!TreesData->treesAreexisting())
            {
                cout << "Warning: There are no trees in the memory!\n\n";
                return;
            }

            if(convert == (String) "Newick")
            {
                TreesData->WriteTrees(outName, NEWICK);
                cout << "Successfully outputted Newick format trees to file: " << outName << "\n\n";
            }
            else
            if(convert == (String) "Nexus")
            {
                TreesData->WriteTrees(outName, NEXUS);
                cout << "Successfully outputted NEXUS format trees to file: " << outName << "\n\n";
            }
        }
        else if (memorydata == (String) "Bipartition Matrix")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "Output bipartition matrix", "Sparse matrix: Choose output format", "List", "Matrix", "Cancel");

            QString qtformat;
            if(reply == 0) // List
                qtformat = "List format";
            else if(reply == 1) // Matrix
                qtformat = "Matrix format";
            else if(reply == 2) // cancel
                return;

            string stdformat = qtformat.toStdString();
            String format(stdformat.c_str());

            string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, format);

            ofstream outBipartMatrix;
            outBipartMatrix.open(namebipartmatrix.c_str());

            if(format == (String) "List format")
            {
                cout << "Outputted Bipartition matrix in list format to " << namebipartmatrix << "\n\n";
                TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
            } else
            if(format == (String) "Matrix format")
            {
                cout << "Outputted Bipartition matrix in matrix format to " << namebipartmatrix << "\n\n";
                TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
            }
        }
        else if (memorydata == (String) "Consensus Tree")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "Output Consensus tree", "Choose output consensus tree format", "Newick", "Nexus", "Cancel");

            QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
            if(reply == 0) // newick
                qtconvert = "Newick";
            else if(reply == 1)
                qtconvert = "Nexus";
            else if(reply == 2) // cancel
                return;

            string stdconvert = qtconvert.toStdString();
            String convert(stdconvert.c_str());
            string outName = TreesData->WriteConsensusTreeFilename(stdfname, stdconvert);
            if(!TreesData->consensusTreeIsexisting())
            {
                cout << "Warning: There is no majority consensus tree in the memory!\n\n";
                return;
            }

            if(convert == (String) "Newick")
            {
                TreesData->WriteConsensusTree(outName, NEWICK);
                cout << "Successfully outputted Newick format consensus trees to file: " << outName << "\n\n";
            }
            else
            if(convert == (String) "Nexus")
            {
                TreesData->WriteConsensusTree(outName, NEXUS);
                cout << "Successfully outputted NEXUS format consensus trees to file: " << outName << "\n\n";
            }
        }
        else
        {
            TreesData->print_matrix(memorydata);
            cout << "Successfully printed " << memorydata << " !\n\n";
        }
    }
    else if(msg1 == 1) // delete
    {

        if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
        {
            cout << "Deleting: Trees\n";

            int rowItm = ui->listDATAmem->row(item);
            delete ui->listDATAmem->takeItem(rowItm);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(item2[i]));

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(item3[i]));

            TreesData->deletetrees();
        }
        else if(memorydata == (String) "Consensus Tree")
        {
            cout << "Deleting: Consensus Tree\n";

            int rowItm = ui->listDATAmem->row(item);
            delete ui->listDATAmem->takeItem(rowItm);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(item2[i]));

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(item3[i]));

            TreesData->deleteConsensustree();
        }
        else if(memorydata == (String) "No data in memory")
        {
            cout << "Warning: There is no data in memory!\n\n";
            return;
        }
        else
        {
            int rowItm = ui->listDATAmem->row(item);
            delete ui->listDATAmem->takeItem(rowItm);

            QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(item2[i]));

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(item3[i]));

            int idex = ui->comboBoxNLDRdata->findText(qtmemorydata);
            ui->comboBoxNLDRdata->removeItem(idex);
            idex = ui->comboBoxDIMdata->findText(qtmemorydata);
            ui->comboBoxDIMdata->removeItem(idex);
            idex = ui->comboBoxDISTmemdata_2->findText(qtmemorydata);
            ui->comboBoxDISTmemdata_2->removeItem(idex);
            idex = ui->comboBoxDISTmemdata_3->findText(qtmemorydata);
            ui->comboBoxDISTmemdata_3->removeItem(idex);
            idex = ui->comboBoxCOVAmemdata->findText(qtmemorydata);
            ui->comboBoxCOVAmemdata->removeItem(idex);

            TreesData->delete_matrix(memorydata);
        }

        if(ui->listDATAmem->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem->addItem(qstr);
        }
        if(ui->listDATAmem_2->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem_2->addItem(qstr);
        }
        if(ui->listDATAmem_3->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem_3->addItem(qstr);
        }
        if(ui->comboBoxNLDRdata->count() == 0)
        {
            QString qstr("No distance data in memory");
            ui->comboBoxNLDRdata->addItem(qstr);
        }
        if(ui->comboBoxDIMdata->count() == 0)
        {
            QString qstr("No distance/coordinate data in memory");
            ui->comboBoxDIMdata->addItem(qstr);
        }
        if(ui->comboBoxCOVAmemdata->count() == 0)
        {
            QString qstr("No covariance data in memory");
            ui->comboBoxCOVAmemdata->addItem(qstr);
        }
        if(ui->comboBoxDISTmemdata_3->count() == 0)
        {
            QString qstr("No affinity data in memory");
            ui->comboBoxDISTmemdata_3->addItem(qstr);
        }
    }
    else if(msg1 == 2) // Cancel
    {
        return;
    }
}

void TreeScaper::on_listDATAmem_2_itemActivated(QListWidgetItem *item)
{
    int msg1 = QMessageBox::question(this, "Data option", "Do you want to delete or output the data value?", "Output", "Delete", "Cancel");

    QString qtmemorydata = item->text();
    string stdmemorydata = qtmemorydata.toStdString();
    String memorydata(stdmemorydata.c_str());
    if(msg1 == 0) // output
    {
        if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();
            int reply = QMessageBox::question(this, "output trees", "Choose output tree format", "Newick", "Nexus", "Cancel");

            QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
            if(reply == 0) // newick
                qtconvert = "Newick";
            else if(reply == 1)
                qtconvert = "Nexus";
            else if(reply == 2) // cancel
                return;

            string stdconvert = qtconvert.toStdString();
            String convert(stdconvert.c_str());
            string outName = TreesData->WriteTreesFilename(stdfname, stdconvert);
            if(!TreesData->treesAreexisting())
            {
                cout << "Warning: There are no trees in the memory!\n\n";
                return;
            }

            if(convert == (String) "Newick")
            {
                TreesData->WriteTrees(outName, NEWICK);
                cout << "Successfully outputted Newick format trees to file: " << outName << "\n\n";
            }
            else
            if(convert == (String) "Nexus")
            {
                TreesData->WriteTrees(outName, NEXUS);
                cout << "Successfully outputted NEXUS format trees to file: " << outName << "\n\n";
            }
        }
        else if (memorydata == (String) "Bipartition Matrix")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "output bipartition matrix", "Sparse matrix: Choose output format", "List", "Matrix", "Cancel");

            QString qtformat;
            if(reply == 0) // List
                qtformat = "List format";
            else if(reply == 1) // Matrix
                qtformat = "Matrix format";
            else if(reply == 2) // cancel
                return;

            string stdformat = qtformat.toStdString();
            String format(stdformat.c_str());

            string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, format);

            ofstream outBipartMatrix;
            outBipartMatrix.open(namebipartmatrix.c_str());

            if(format == (String) "List format")
            {
                cout << "Outputted Bipartition matrix in list format to " << namebipartmatrix << "\n\n";
                TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
            } else
            if(format == (String) "Matrix format")
            {
                cout << "Outputted Bipartition matrix in matrix format to " << namebipartmatrix << "\n\n";
                TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
            }
        }
        else if (memorydata == (String) "Consensus Tree")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "Output Consensus tree", "Choose output consensus tree format", "Newick", "Nexus", "Cancel");

            QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
            if(reply == 0) // newick
                qtconvert = "Newick";
            else if(reply == 1)
                qtconvert = "Nexus";
            else if(reply == 2) // cancel
                return;

            string stdconvert = qtconvert.toStdString();
            String convert(stdconvert.c_str());
            string outName = TreesData->WriteConsensusTreeFilename(stdfname, stdconvert);
            if(!TreesData->consensusTreeIsexisting())
            {
                cout << "Warning: There is no majority consensus tree in the memory!\n\n";
                return;
            }

            if(convert == (String) "Newick")
            {
                TreesData->WriteConsensusTree(outName, NEWICK);
                cout << "Successfully outputted Newick format consensus trees to file: " << outName << "\n\n";
            }
            else
            if(convert == (String) "Nexus")
            {
                TreesData->WriteConsensusTree(outName, NEXUS);
                cout << "Successfully outputted NEXUS format consensus trees to file: " << outName << "\n\n";
            }
        }
        else
        {
            TreesData->print_matrix(memorydata);
            cout << "Successfully printed " << memorydata << " !\n\n";
        }
    }
    else if(msg1 == 1) // delete
    {
        if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
        {
            cout << "Deleting: Trees\n";

            QList<QListWidgetItem *> item2 = ui->listDATAmem->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem->takeItem(ui->listDATAmem->row(item2[i]));

            int rowItm = ui->listDATAmem_2->row(item);
            delete ui->listDATAmem_2->takeItem(rowItm);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(item3[i]));

            TreesData->deletetrees();
        }
        else if(memorydata == (String) "Consensus Tree")
        {
            cout << "Deleting: Consensus Tree\n";

            QList<QListWidgetItem *> item2 = ui->listDATAmem->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem->takeItem(ui->listDATAmem->row(item2[i]));

            int rowItm = ui->listDATAmem_2->row(item);
            delete ui->listDATAmem_2->takeItem(rowItm);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(item3[i]));

            TreesData->deleteConsensustree();
        }
        else if(memorydata == (String) "No data in memory")
        {
            cout << "Warning: There is no data in memory!\n\n";
            return;
        }
        else
        {
            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item1.size(); i++)
                delete ui->listDATAmem->takeItem(ui->listDATAmem->row(item1[i]));

            int rowItm = ui->listDATAmem_2->row(item);
            delete ui->listDATAmem_2->takeItem(rowItm);

            QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(item3[i]));

            int idex = ui->comboBoxNLDRdata->findText(qtmemorydata);
            ui->comboBoxNLDRdata->removeItem(idex);
            idex = ui->comboBoxDIMdata->findText(qtmemorydata);
            ui->comboBoxDIMdata->removeItem(idex);
            idex = ui->comboBoxDISTmemdata_2->findText(qtmemorydata);
            ui->comboBoxDISTmemdata_2->removeItem(idex);
            idex = ui->comboBoxDISTmemdata_3->findText(qtmemorydata);
            ui->comboBoxDISTmemdata_3->removeItem(idex);
            idex = ui->comboBoxCOVAmemdata->findText(qtmemorydata);
            ui->comboBoxCOVAmemdata->removeItem(idex);

            TreesData->delete_matrix(memorydata);
        }

        if(ui->listDATAmem->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem->addItem(qstr);
        }
        if(ui->listDATAmem_2->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem_2->addItem(qstr);
        }
        if(ui->listDATAmem_3->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem_3->addItem(qstr);
        }
        if(ui->comboBoxNLDRdata->count() == 0)
        {
            QString qstr("No distance data in memory");
            ui->comboBoxNLDRdata->addItem(qstr);
        }
        if(ui->comboBoxDIMdata->count() == 0)
        {
            QString qstr("No distance/coordinate data in memory");
            ui->comboBoxDIMdata->addItem(qstr);
        }
        if(ui->comboBoxCOVAmemdata->count() == 0)
        {
            QString qstr("No covariance data in memory");
            ui->comboBoxCOVAmemdata->addItem(qstr);
        }
        if(ui->comboBoxDISTmemdata_3->count() == 0)
        {
            QString qstr("No affinity data in memory");
            ui->comboBoxDISTmemdata_3->addItem(qstr);
        }
    }
    else if(msg1 == 2) // Cancel
    {
        return;
    }
}

void TreeScaper::on_listDATAmem_3_itemActivated(QListWidgetItem *item)
{
    int msg1 = QMessageBox::question(this, "Data option", "Do you want to delete or output the data value?", "Output", "Delete", "Cancel");

    QString qtmemorydata = item->text();
    string stdmemorydata = qtmemorydata.toStdString();
    String memorydata(stdmemorydata.c_str());
    if(msg1 == 0) // output
    {
        if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();
            int reply = QMessageBox::question(this, "output trees", "Choose output tree format", "Newick", "Nexus", "Cancel");

            QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
            if(reply == 0) // newick
                qtconvert = "Newick";
            else if(reply == 1)
                qtconvert = "Nexus";
            else if(reply == 2) // cancel
                return;

            string stdconvert = qtconvert.toStdString();
            String convert(stdconvert.c_str());
            string outName = TreesData->WriteTreesFilename(stdfname, stdconvert);
            if(!TreesData->treesAreexisting())
            {
                cout << "Warning: There are no trees in the memory!\n\n";
                return;
            }

            if(convert == (String) "Newick")
            {
                TreesData->WriteTrees(outName, NEWICK);
                cout << "Successfully outputted Newick format trees to file: " << outName << "\n\n";
            }
            else
            if(convert == (String) "Nexus")
            {
                TreesData->WriteTrees(outName, NEXUS);
                cout << "Successfully outputted NEXUS format trees to file: " << outName << "\n\n";
            }
        }
        else if (memorydata == (String) "Bipartition Matrix")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "output bipartition matrix", "Sparse matrix: Choose output format", "List", "Matrix", "Cancel");

            QString qtformat;
            if(reply == 0) // List
                qtformat = "List format";
            else if(reply == 1) // Matrix
                qtformat = "Matrix format";
            else if(reply == 2) // cancel
                return;

            string stdformat = qtformat.toStdString();
            String format(stdformat.c_str());

            string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, format);

            ofstream outBipartMatrix;
            outBipartMatrix.open(namebipartmatrix.c_str());

            if(format == (String) "List format")
            {
                cout << "Outputted Bipartition matrix in list format to " << namebipartmatrix << "\n\n";
                TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
            } else
            if(format == (String) "Matrix format")
            {
                cout << "Outputted Bipartition matrix in matrix format to " << namebipartmatrix << "\n\n";
                TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
            }
        }
        else if (memorydata == (String) "Consensus Tree")
        {
            QString qtfname = ui->textDATAfile->toPlainText();
            string stdfname = qtfname.toStdString();

            int reply = QMessageBox::question(this, "Output Consensus tree", "Choose output consensus tree format", "Newick", "Nexus", "Cancel");

            QString qtconvert;//--- = ui->comboBoxDATAoutformat->currentText();
            if(reply == 0) // newick
                qtconvert = "Newick";
            else if(reply == 1)
                qtconvert = "Nexus";
            else if(reply == 2) // cancel
                return;

            string stdconvert = qtconvert.toStdString();
            String convert(stdconvert.c_str());
            string outName = TreesData->WriteConsensusTreeFilename(stdfname, stdconvert);
            if(!TreesData->consensusTreeIsexisting())
            {
                cout << "Warning: There is no majority consensus tree in the memory!\n\n";
                return;
            }

            if(convert == (String) "Newick")
            {
                TreesData->WriteConsensusTree(outName, NEWICK);
                cout << "Successfully outputted Newick format consensus trees to file: " << outName << "\n\n";
            }
            else
            if(convert == (String) "Nexus")
            {
                TreesData->WriteConsensusTree(outName, NEXUS);
                cout << "Successfully outputted NEXUS format consensus trees to file: " << outName << "\n\n";
            }
        }
        else
        {
            TreesData->print_matrix(memorydata);
            cout << "Successfully printed " << memorydata << " !\n\n";
        }
    }
    else if(msg1 == 1) // delete
    {
        if(memorydata == (String) "Unweighted treeset" || memorydata == (String) "Weighted treeset")
        {
            cout << "Deleting: Trees\n";

            QList<QListWidgetItem *> item2 = ui->listDATAmem->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem->takeItem(ui->listDATAmem->row(item2[i]));

            QList<QListWidgetItem *> item3 = ui->listDATAmem_2->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(item3[i]));

            int rowItm = ui->listDATAmem_3->row(item);
            delete ui->listDATAmem_3->takeItem(rowItm);

            TreesData->deletetrees();
        }
        else if(memorydata == (String) "Consensus Tree")
        {
            cout << "Deleting: Consensus Tree\n";

            QList<QListWidgetItem *> item2 = ui->listDATAmem->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item2.size(); i++)
                delete ui->listDATAmem->takeItem(ui->listDATAmem->row(item2[i]));

            QList<QListWidgetItem *> item3 = ui->listDATAmem_2->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(item3[i]));

            int rowItm = ui->listDATAmem_3->row(item);
            delete ui->listDATAmem_3->takeItem(rowItm);

            TreesData->deleteConsensustree();
        }
        else if(memorydata == (String) "No data in memory")
        {
            cout << "Warning: There is no data in memory!\n\n";
            return;
        }
        else
        {
            QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item1.size(); i++)
                delete ui->listDATAmem->takeItem(ui->listDATAmem->row(item1[i]));

            QList<QListWidgetItem *> item3 = ui->listDATAmem_2->findItems(qtmemorydata, Qt::MatchCaseSensitive | Qt::MatchExactly);
            for (int i = 0; i < item3.size(); i++)
                delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(item3[i]));

            int rowItm = ui->listDATAmem_3->row(item);
            delete ui->listDATAmem_3->takeItem(rowItm);

            int idex = ui->comboBoxNLDRdata->findText(qtmemorydata);
            ui->comboBoxNLDRdata->removeItem(idex);
            idex = ui->comboBoxDIMdata->findText(qtmemorydata);
            ui->comboBoxDIMdata->removeItem(idex);
            idex = ui->comboBoxDISTmemdata_2->findText(qtmemorydata);
            ui->comboBoxDISTmemdata_2->removeItem(idex);
            idex = ui->comboBoxDISTmemdata_3->findText(qtmemorydata);
            ui->comboBoxDISTmemdata_3->removeItem(idex);
            idex = ui->comboBoxCOVAmemdata->findText(qtmemorydata);
            ui->comboBoxCOVAmemdata->removeItem(idex);

            TreesData->delete_matrix(memorydata);
        }

        if(ui->listDATAmem->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem->addItem(qstr);
        }
        if(ui->listDATAmem_2->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem_2->addItem(qstr);
        }
        if(ui->listDATAmem_3->count() == 0)
        {
            QString qstr("No data in memory");
            ui->listDATAmem_3->addItem(qstr);
        }
        if(ui->comboBoxNLDRdata->count() == 0)
        {
            QString qstr("No distance data in memory");
            ui->comboBoxNLDRdata->addItem(qstr);
        }
        if(ui->comboBoxDIMdata->count() == 0)
        {
            QString qstr("No distance/coordinate data in memory");
            ui->comboBoxDIMdata->addItem(qstr);
        }
        if(ui->comboBoxCOVAmemdata->count() == 0)
        {
            QString qstr("No covariance data in memory");
            ui->comboBoxCOVAmemdata->addItem(qstr);
        }
        if(ui->comboBoxDISTmemdata_3->count() == 0)
        {
            QString qstr("No affinity data in memory");
            ui->comboBoxDISTmemdata_3->addItem(qstr);
        }
    }
    else if(msg1 == 2) // Cancel
    {
        return;
    }
}

void TreeScaper::on_comboBoxDIMdata_currentIndexChanged(const QString &arg1)
{
    if(ui->comboBoxDIMdata->count() == 0)
    {
        QString qstr("No distance/coordinate data in memory");
        ui->comboBoxDIMdata->addItem(qstr);
    }
}

void TreeScaper::on_pushButtonAffinitybrowse_clicked()
{
    QFileDialog *fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setDirectory(".");
    if(fileDialog->exec() == QDialog::Accepted)
    {
        QString path = fileDialog->selectedFiles()[0];
        ui->textAfinityfile->setPlainText(path);
    }
}

void TreeScaper::on_pushAffinityloaddist_clicked()
{
    QString qtfname = ui->textAfinityfile->toPlainText();
    string stdfname = qtfname.toStdString();
    String fname(stdfname.c_str());
    File file(fname);
    if(! file.is_open())
    {
        cout << "Error: Cannot open the affinity file!\n\n";
        return;
    }
    TreesData->load_affinityfile(stdfname);

    QString qstr("File-affinity");
    if (ui->comboBoxDISTmemdata_3->findText("File-affinity") == -1)
        ui->comboBoxDISTmemdata_3->addItem(qstr);

    QList<QListWidgetItem *> item1 = ui->listDATAmem->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item1.empty())
        ui->listDATAmem->addItem(qstr);

    QList<QListWidgetItem *> item2 = ui->listDATAmem_2->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item2.empty())
        ui->listDATAmem_2->addItem(qstr);

    QList<QListWidgetItem *> item3 = ui->listDATAmem_3->findItems(qstr, Qt::MatchCaseSensitive | Qt::MatchExactly);
    if(item3.empty())
        ui->listDATAmem_3->addItem(qstr);

    int idex = ui->comboBoxDISTmemdata_3->findText("No affinity data in memory");
    ui->comboBoxDISTmemdata_3->removeItem(idex);

    QList<QListWidgetItem *> itemTmp1 = ui->listDATAmem->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
    for (int i = 0; i < itemTmp1.size(); i++)
        delete ui->listDATAmem->takeItem(ui->listDATAmem->row(itemTmp1[i]));

    QList<QListWidgetItem *> itemTmp2 = ui->listDATAmem_2->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
    for (int i = 0; i < itemTmp2.size(); i++)
        delete ui->listDATAmem_2->takeItem(ui->listDATAmem_2->row(itemTmp2[i]));

    QList<QListWidgetItem *> itemTmp3 = ui->listDATAmem_3->findItems(QString("No data in memory"), Qt::MatchCaseSensitive | Qt::MatchExactly);
    for (int i = 0; i < itemTmp3.size(); i++)
        delete ui->listDATAmem_3->takeItem(ui->listDATAmem_3->row(itemTmp3[i]));

    ui->pushDISTcomm->setEnabled(true);

    cout << "Successfully read affinity file.\n\n";
}
