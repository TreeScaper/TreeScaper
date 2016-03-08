
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

#ifndef WCOMMUNITYTHREAD_CPP
#define WCOMMUNITYTHREAD_CPP

#include "treescaper.h"
#include "wCommunitythread.h"

Communitythread::Communitythread(QObject *parent)
    : QThread(parent)
{
};

Communitythread::~Communitythread()
{
};

void Communitythread::initialization(Trees *TreesData, String in_str_matrix, int in_modelType, Array<double> in_param1, Array<double> in_param2,
                                     string in_highfreq, string in_lowfreq, int in_buttonflag, bool in_isauto)
{
    lptrees = TreesData;
    str_matrix = in_str_matrix;
    modelType = in_modelType;
    param1 = in_param1;
    param2 = in_param2;
    highfreq = in_highfreq;
    lowfreq = in_lowfreq;
    buttonflag = in_buttonflag;
    isauto = in_isauto;
};

void Communitythread::run()
{
    emit sendbuttonCommunityPlotenable(0);
    if(isauto)
        if(lptrees->compute_community_automatically(str_matrix, modelType, highfreq, lowfreq))
            std::cout << "Successfully detected communities of " << str_matrix << " by " << modelType << "!" << std::endl;

    if(!isauto)
        if(lptrees->compute_community_manually(str_matrix, modelType, param1, param2, highfreq, lowfreq))
            std::cout << "Successfully detected communities of " << str_matrix << " by " << modelType << "!" << std::endl;

    emit sendbuttonCommunityPlotenable(buttonflag);
}

#endif // WNLDRTHREAD_CPP
