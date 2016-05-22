
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

#ifndef WPLOTTREETHREAD_CPP
#define WPLOTTREETHREAD_CPP

#include "wdef.h"
#include "wImage.h"
#include "wPlotTreethread.h"
#include "treescaper.h"

PlotTreethread::PlotTreethread(QObject *parent)
    : QThread(parent)
{
};

PlotTreethread::~PlotTreethread()
{
};

void PlotTreethread::initialization(image_parameters plotparas, const NEWICKTREE *intree, bool inisrooted, bool inisweighted, const LabelMap *inleaveslm, String title)
{
    lpplottree = new PlotTree;
    memcpy(&(lpplottree->parameters), &(paras::imageparas), sizeof(image_parameters));
    tree = intree;
    isrooted = inisrooted;
    isweighted = inisweighted;
    leaveslm = inleaveslm;
#if defined(_MAC) || defined(_LINUX)
    lpplottree->Initialize_PlotTree(tree, isrooted, isweighted, leaveslm, title);
    this->exec();
    delete lpplottree;
#endif
};

void PlotTreethread::run()
{
    /*
#ifdef _LINUX
    if(paras::plot_type == (String) "Points")
    {
        lpimage->Load_points_with_selection();
        if(paras::is_sort_volume)
            lpimage->sort_volume();
        lpimage->Create_PointsActors();
        lpimage->Plot_Points();
    } else
    {
        lpimage->Load_points_with_selection();
        if(paras::is_sort_volume)
            lpimage->sort_volume();
        lpimage->Create_HullsActors();
        lpimage->Plot_Hulls();
    }

    if(paras::plot_makemovie)
        lpimage->make_movie();

    delete lpimage;
#endif

#ifdef _WINDOWS
    if(paras::plot_type == (String) "Points")
    {
        lpimage->Load_points_with_selection();
        if(paras::is_sort_volume)
            lpimage->sort_volume();
        lpimage->Create_PointsActors();
        lpimage->Plot_Points();
    } else
    {
        lpimage->Load_points_with_selection();
        if(paras::is_sort_volume)
            lpimage->sort_volume();
        lpimage->Create_HullsActors();
        lpimage->Plot_Hulls();
    }

    if(paras::plot_makemovie)
        lpimage->make_movie();

    delete lpimage;
#endif*/
}

#endif // WPlotTreethread_CPP
