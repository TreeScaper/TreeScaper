
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

#ifndef WPLOTTHREAD_CPP
#define WPLOTTHREAD_CPP

#include "wdef.h"
#include "wImage.h"
#include "wPlotthread.h"
#include "treescaper.h"

Plotthread::Plotthread(QObject *parent)
    : QThread(parent)
{
};

Plotthread::~Plotthread()
{
};

void Plotthread::initialization(String fname, image_parameters imageparas, Array<int> &selected_ids)
{
    filename = fname;
    lpimage = new Image;
    memcpy(&(lpimage->parameters), &(paras::imageparas), sizeof(image_parameters));
    lpimage->Initialize_Image(filename);

    lpimage->Set_SelectedIDs(selected_ids);
#ifdef _MAC
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
    this->exec();
//    delete lpimage;
#endif
};

void Plotthread::run()
{
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
#endif
}

#endif // WPLOTTHREAD_CPP
