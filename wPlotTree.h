
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

//wPlotTree.h
// using vtk to plot points
//                by whuang 15/May/2010

#ifndef WPLOTTREE_H
#define WPLOTTREE_H

#include "wfile.h"
#include "wstring.h"
#include "warray.cpp"
#include "wImage.h"
#include "TreeOPE.h"
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPoints.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointsProjectedHull.h>
#include <vtkPolyLine.h>
#include "vtkGlyph3D.h"
#include "vtkSphereSource.h"
#include "vtkMassProperties.h"
#include "vtkPlaneSource.h"
#include <vtkTriangleFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>

//#include "vtkMPEG2Writer.h"
#include <vtkWindowToImageFilter.h>
#include <vtkAVIWriter.h>
//#include <vtkFFMPEGWriter.h>
#include <vtkCamera.h>
#include <vtkWindow.h>
#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkImageMandelbrotSource.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkLineSource.h"
#include <vtkLegendBoxActor.h>

#include <vtkCircularLayoutStrategy.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkGraphLayoutView.h>
#include <vtkIntArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTree.h>


#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkDataSetAttributes.h>
#include <vtkGraphLayoutView.h>
#include <vtkIntArray.h>
#include <vtkLookupTable.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkTree.h>
#include <vtkTreeDFSIterator.h>
#include <vtkUnsignedCharArray.h>
#include <vtkViewTheme.h>

#include <vtkGraphLayoutView.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTree.h>
#include <vtkStringArray.h>

class PlotTree{
public:
    PlotTree(){};

    ~PlotTree(){destructor();};

    void destructor(){};

    PlotTree(const NEWICKTREE *tree, bool inisrooted, bool inisweighted, const LabelMap *inleaveslm, String title)
	{
        Initialize_PlotTree(tree, inisrooted, inisweighted, inleaveslm, title);
	};

    void maketree(const NEWICKTREE *tree)
    {
        vtkgraph = vtkSmartPointer<vtkMutableDirectedGraph>::New();

        if(isweighted)
        {
            edgeweights = vtkSmartPointer<vtkDoubleArray>::New();
            edgeweights->SetNumberOfComponents(1);
            edgeweights->SetName("Weights");
        }

        vertexIDs = vtkSmartPointer<vtkIntArray>::New();
        vertexIDs->SetNumberOfComponents(1);
        vertexIDs->SetName("VertexIDs");
        vertexNames = vtkSmartPointer<vtkStringArray>::New();
        vertexNames->SetNumberOfComponents(1);
        vertexNames->SetName("vertexNames");

        int depth = 0;
        makesubtree(tree->root, depth);
    };

    vtkIdType makesubtree(const NEWICKNODE *node, int depth)
    {
        vtkIdType cidx;
        vtkIdType vidx = vtkgraph->AddVertex();
        if(depth == 0)
        {
            if(isrooted)
                vertexNames->InsertNextValue("root");
            else
                vertexNames->InsertNextValue("psuedo-root");
            vertexIDs->InsertNextValue(-2);
        } else if(node->Nchildren != 0)
        {
            vertexIDs->InsertNextValue(-1);
            vertexNames->InsertNextValue("");
        } else
        {
            vertexIDs->InsertNextValue((int) atoi(node->label));
            vertexNames->InsertNextValue(leaveslm->name((int) atoi(node->label) - 1));
//            vertexNames->InsertNextValue(node->label);
        }

        depth++;
        for(int i = 0; i < node->Nchildren; i++)
        {
            cidx = makesubtree(node->child[i], depth);
            vtkgraph->AddEdge(vidx, cidx);
            if(isweighted)
            {
                edgeweights->InsertNextValue(node->child[i]->weight);
            }
        }
        return vidx;
    };

    void Initialize_PlotTree(const NEWICKTREE *tree, bool inisrooted, bool inisweighted, const LabelMap *inleaveslm, String title)
    {
        isrooted = inisrooted;
        isweighted = inisweighted;
        leaveslm = inleaveslm;

        maketree(tree);

        // Add the edge weight array to the graph
        vtkgraph->GetEdgeData()->AddArray(edgeweights);
        vtkgraph->GetVertexData()->AddArray(vertexIDs);
        vtkgraph->GetVertexData()->AddArray(vertexNames);

        vtkSmartPointer<vtkTree> vtktree = vtkSmartPointer<vtkTree>::New();
        if(! vtktree->CheckedShallowCopy(vtkgraph))
        {
            std::cout << "no a tree" << std::endl;
            return;
        }

        vtkSmartPointer<vtkGraphLayoutView> graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
        graphLayoutView->AddRepresentationFromInput(vtktree);
//        graphLayoutView->AddRepresentationFromInput(vtkgraph);

//        vtkSmartPointer<vtkCircularLayoutStrategy> circularLayoutStrategy = vtkSmartPointer<vtkCircularLayoutStrategy>::New();
//        graphLayoutView->SetLayoutStrategy(circularLayoutStrategy);
//        graphLayoutView->SetLayoutStrategyToCone();

//        graphLayoutView->SetLayoutStrategyToFast2D();
//        graphLayoutView->SetLayoutStrategyToCommunity2D();

        graphLayoutView->SetLayoutStrategyToTree();
//        graphLayoutView->SetLayoutStrategyToSpanTree();
        graphLayoutView->SetVertexLabelVisibility(true);
//        graphLayoutView->SetHideVertexLabelsOnInteraction(true);
        graphLayoutView->SetEdgeLabelVisibility(true);
        graphLayoutView->SetEdgeLabelArrayName("Weights"); //default is "labels"
//        graphLayoutView->SetVertexLabelArrayName("VertexIDs"); //default is "labels"
        graphLayoutView->SetVertexLabelArrayName("vertexNames"); //default is "labels"
        graphLayoutView->ResetCamera();
        graphLayoutView->Render();
//        graphLayoutView->GetRenderer()->SetBackground(1, 1, 1);
        graphLayoutView->GetRenderWindow()->SetWindowName((char *) title);
        graphLayoutView->GetRenderWindow()->SetSize(400, 400);
        graphLayoutView->GetInteractor()->Start();

        return;
    }

    image_parameters parameters;

private:

    bool isrooted;
    bool isweighted;
    const LabelMap *leaveslm;
    vtkSmartPointer<vtkMutableDirectedGraph> vtkgraph;
    vtkSmartPointer<vtkDoubleArray> edgeweights;
    vtkSmartPointer<vtkIntArray> vertexIDs;
    vtkSmartPointer<vtkStringArray> vertexNames;
};

#endif
