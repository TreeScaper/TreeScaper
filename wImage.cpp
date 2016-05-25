
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

//wplot.cpp
// using vtk to plot points
//                by whuang 15/May/2010

#ifndef WIMAGE_CPP
#define WIMAGE_CPP

#undef max
#undef min

#define mymin(a,b) ((a)<(b)?(a):(b))
#define mymax(a,b) ((a)>(b)?(a):(b))

#include "wfile.h"
#include "wImage.h"
#include <iostream>
#include "vtkLineSource.h"
#include "treescaper.h"

#include "PointsSource.h"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkIdTypeArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkRendererCollection.h>
#include <vtkProperty.h>
#include <vtkPlanes.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPointSource.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkAreaPicker.h>
#include <vtkExtractGeometry.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkIdFilter.h>
#include "vtkCallbackCommand.h"

// Define interaction style
class InteractorStyle : public vtkInteractorStyleRubberBandPick
{
  public:
    static InteractorStyle* New();
    vtkTypeMacro(InteractorStyle,vtkInteractorStyleRubberBandPick);

    InteractorStyle()
    {
        this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
        this->SelectedActor = vtkSmartPointer<vtkActor>::New();
        this->SelectedActor->SetMapper(SelectedMapper);
        this->PreSelectedPoints = new vtkSmartPointer<vtkPoints> [1];
        this->PreSelectedPoints[0] = vtkSmartPointer<vtkPoints>::New();
        iskeydown = false;
        preMtime = 0;
        idx = 0;
    }

    void SetSelectedPointsSize(int size)
    {
        this->pointsize = size;
    }

    virtual void OnKeyDown()
    {
        iskeydown = true;
    }

    virtual void OnKeyUp()
    {
        iskeydown = false;
    }

    virtual void OnLeftButtonUp()
    {
      // Forward events
      vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

      vtkPlanes* frustum = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetFrustum();
      if(preMtime == frustum->GetMTime())
          return;

      preMtime = frustum->GetMTime();

      vtkSmartPointer<vtkExtractGeometry> extractGeometry =
        vtkSmartPointer<vtkExtractGeometry>::New();
      extractGeometry->SetImplicitFunction(frustum);
#if VTK_MAJOR_VERSION <= 5
      extractGeometry->SetInput(this->Points);
#else
      extractGeometry->SetInputData(this->Points);
#endif
      extractGeometry->Update();

      vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
        vtkSmartPointer<vtkVertexGlyphFilter>::New();
      glyphFilter->SetInputConnection(extractGeometry->GetOutputPort());
      glyphFilter->Update();

      vtkPolyData* selected = glyphFilter->GetOutput();

      if(!iskeydown)
      {
          PreSelectedPoints[0]->SetNumberOfPoints(0);
          idx = 0;
      }

      vtkIdTypeArray* ids = vtkIdTypeArray::SafeDownCast(selected->GetPointData()->GetArray("OriginalIds"));
      int new_pro, j;
      cout << "new selected pts:" << endl;
      for(vtkIdType i = 0; i < ids->GetNumberOfTuples(); i++)
      {
          new_pro = IDMaps[ids->GetValue(i)];
          for(j = 0; j < idx; j++)
          {
              if((*selected_pts)[j] == new_pro)
                  break;
          }
          if(j == idx)
          {
              (*selected_pts)[idx] = new_pro;
              idx++;
          }
          cout << "i:" << IDMaps[ids->GetValue(i)] << ", ";
      }
      cout << endl;
      (*selected_pts)[idx] = -1;

      cout << "all selected pts:" << endl;
      for(int i = 0; i < idx; i++)
          cout << "i:" << (*selected_pts)[i] << ", ";
      cout << endl;

      vtkPoints *selectedpts = selected->GetPoints();
      for(vtkIdType i = 0; i < selectedpts->GetNumberOfPoints(); i++)
      {
          PreSelectedPoints[0]->InsertNextPoint(selectedpts->GetPoint(i));
      }

      vtkSmartPointer<PointSource> PS = vtkSmartPointer<PointSource>::New();
      vtkSmartPointer<vtkIdFilter> IDF = vtkSmartPointer<vtkIdFilter>::New();
      vtkSmartPointer<vtkDataSetSurfaceFilter> SF = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
      PS->SetPoints(PreSelectedPoints, 1);
      PS->Update();
      IDF->SetInputConnection(PS->GetOutputPort());
      IDF->Update();
      SF->SetInputConnection(IDF->GetOutputPort());
      SF->Update();
      selected = SF->GetOutput();

#if VTK_MAJOR_VERSION <= 5
      this->SelectedMapper->SetInput(selected);
#else
      this->SelectedMapper->SetInputData(selected);
#endif

      this->SelectedMapper->ScalarVisibilityOff();
      this->SelectedActor->GetProperty()->SetColor(0.0, 0.0, 0.0); //(R,G,B)
      this->SelectedActor->GetProperty()->SetPointSize(pointsize);

      this->CurrentRenderer->AddActor(SelectedActor);
      this->GetInteractor()->GetRenderWindow()->Render();
      this->HighlightProp(NULL);
    }

    void SetPoints(vtkSmartPointer<vtkPolyData> points) {this->Points = points;}
    void SetIDMaps(int *maps){IDMaps = maps;}
    void SetOutputPts(Array<int> *Pts){selected_pts = Pts;}
  private:
    vtkSmartPointer<vtkPolyData> Points;
    vtkSmartPointer<vtkPoints> *PreSelectedPoints;
    vtkSmartPointer<vtkActor> SelectedActor;
    vtkSmartPointer<vtkDataSetMapper> SelectedMapper;
    int pointsize;
    unsigned long preMtime;
    bool iskeydown;
    Array<int> * selected_pts;
    int idx;
    int *IDMaps;
};

vtkStandardNewMacro(InteractorStyle);

//------------

#include "vtkCommand.h"
#include "vtkRandomGraphSource.h"
#include "vtkGraph.h"
#include "PlotGraph.h"
#include "vtkVariant.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkContextView.h"
#include "vtkContextScene.h"

#include "vtkObjectFactory.h"

//#include "vtkRegressionTestImage.h" //---for windows---
//#include "newick.h"
#include "label-map.hh"

class GraphAnimate : public vtkCommand
{
public:
  static GraphAnimate *New() { return new GraphAnimate(); }
  vtkTypeMacro(GraphAnimate, vtkCommand);
  virtual void Execute(vtkObject *, unsigned long, void *)
    {
    this->GraphItem->UpdatePositions();
    this->View->Render();
    this->View->GetRenderWindow()->GetInteractor()->CreateOneShotTimer(10);
    }
  vtkGraphItem* GraphItem;
  vtkContextView* View;
};

int testfunction()
{

    //******************** Read a tree and parsing the tree for taxon label collection **************//
/*
    NEWICKTREE *newickTree;
    int err;
    FILE *fp;
    fp = fopen("triangle.txt", "r");

    if (!fp) {cout << "File open error\n"; exit(0);}

    newickTree = loadnewicktree2(fp, &err);
    if (!newickTree)
    {
        switch (err)
        {
            case -1:
                printf("Out of memory \n");
                break;
            case -2:
                printf("Pase error \n");
                break;
            case -3:
                printf("Cannot load file \n");
                break;
            default:
                printf("Error %d\n", err);
        }
    }

    LabelMap lm;

    try {
        GetTaxaLabels(newickTree->root, lm);
    }
    catch (LabelMap::AlreadyPushedEx ex) {
        cerr << "Error : The label ' " << ex.label << "' appeared twice" << endl;
        exit(2);
    }

    NUM_Taxa = lm.size();
    cout << "Number of taxa in the input file: " << NUM_Taxa << endl;
    killnewicktree(newickTree);
    fclose(fp);

    bool bUbid = false;// for counting the number pf unique bipartitions
*/
// ************* Collect bipartitions ****************************//
/*
    fp = fopen(fname.c_str(), "r");
    HashRFMap vec_hashrf; // Class HashRFMap

    unsigned long long M1 = 0;
    unsigned long long M2 = 0;

    if (NEWSEED != 1000)
    {
        vec_hashrf.uhashfunc_init(Num_tree, NUM_Taxa, C, NEWSEED);
    }
    else
        vec_hashrf.uhashfunc_init(Num_tree, NUM_Taxa, C);

    M1 = vec_hashrf._HF.getM1();
    M2 = vec_hashrf._HF.getM2();
    vec_hashrf._hashtab2.resize(M1);

    if (!fp) {cout << "File open error\n"; exit(0);}

    // fixed a leaf to be one side of the root
    std::string leafroot = lm.name(0); // normalized unrooted tree by WH
    NEWICKNODE *lrpt = NULL; // normalized unrooted tree by WH
    int indexchild = -1; // normalized unrooted tree by WH

    for (unsigned int treeIdx = 0; treeIdx < Num_tree; ++treeIdx)
    {
        newickTree = loadnewicktree2(fp, &err);

        // find the leaf
        if(!isrooted)
        {
            lrpt = findleaf(leafroot, newickTree->root, NULL, &indexchild); // normalized unrooted tree by WH
    //        std::cout << "pretree" << std::endl;//----
    //        printnewicktree(newickTree);

            // build the normalized tree
            normailzedTree(lrpt, newickTree, indexchild); // normalized unrooted tree by WH
    //        std::cout << "posttree" << std::endl;//----
    //        printnewicktree(newickTree);
        }

        if (!newickTree)
        {
            switch (err)
            {
                case -1:
                    printf("Out of Memory\n");
                    break;
                case -2:
                    printf("Parse Error\n");
                    break;
                case -3:
                    printf("Cann't load file\n");
                    break;
                default:
                    printf("Error %d\n", err);
            }
        }
        else {
            unsigned int numBitstr = 0;
            dfs_compute_hash(newickTree->root, lm, vec_hashrf, treeIdx, numBitstr, M1, M2, WEIGHTED,NUM_Taxa);
            killnewicktree(newickTree);
        }
    }
    treefile.close();



*/

  // Set up a 2D context view, context test object and add it to the scene

  vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
  view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
  view->GetRenderWindow()->SetSize(800, 600);

  vtkSmartPointer<vtkRandomGraphSource> source = vtkSmartPointer<vtkRandomGraphSource>::New();
  source->SetNumberOfVertices(100);
  source->SetNumberOfEdges(0);
  source->StartWithTreeOn();
  source->Update();
  vtkSmartPointer<vtkGraphItem> item = vtkSmartPointer<vtkGraphItem>::New();
  item->SetGraph(source->GetOutput());
  view->GetScene()->AddItem(item);

  vtkSmartPointer<GraphAnimate> anim = vtkSmartPointer<GraphAnimate>::New();
  anim->View = view;
  anim->GraphItem = item;
  view->GetRenderWindow()->GetInteractor()->Initialize();
  view->GetRenderWindow()->GetInteractor()->CreateOneShotTimer(10);
  view->GetRenderWindow()->GetInteractor()->AddObserver(vtkCommand::TimerEvent, anim);

  view->GetRenderWindow()->GetInteractor()->Start();
  return 0;
}



//-----------

void Image::Initialize_Image(String filename)
{
//---    testfunction();
    fname = filename;
    File fcor(filename);
    if(!fcor.is_open())
    {
        cout << "error: file \"" << filename << "\" can not be open." << endl;
        exit(0);
    }
    fcor.seek(0);
    size = fcor.lines();
    fcor.seek(0);
    dim = fcor.cols();
    fcor.seek(0);

    centers = new double[parameters.cluster_num * dim];
    COR = new double [size * dim];
    selectpts = new bool[size];
    for(int i = 0; i < size; i++)
        selectpts[i] = true;

    points = new vtkSmartPointer<vtkPoints> [parameters.cluster_num];
    Hull_points = new vtkSmartPointer<vtkPoints> [parameters.cluster_num];
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        points[i] = vtkSmartPointer<vtkPoints>::New();
        Hull_points[i] = vtkSmartPointer<vtkPoints>::New();
    }

    ranks = new int [parameters.cluster_num];
    volumes = new double [parameters.cluster_num];
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        ranks[i] = i;
        volumes[i] = i;
    }

//    PointsActors = new vtkSmartPointer<vtkActor> [parameters.cluster_num];
    PointsActors = new vtkSmartPointer<vtkActor> [parameters.cluster_num];
    sphere = new vtkSphereSource *[parameters.cluster_num];
    model = new vtkPolyData *[parameters.cluster_num];
//    glyph = new vtkGlyph3D *[parameters.cluster_num];
    pointMapper = new vtkPolyDataMapper *[parameters.cluster_num];
    VP = new vtkProperty *[parameters.cluster_num];

    HullsActors = new vtkSmartPointer<vtkActor> [parameters.cluster_num];
    HullMapper = new vtkPolyDataMapper *[parameters.cluster_num];
    polydata = new vtkSmartPointer<vtkPolyData> [parameters.cluster_num];
    delaunay = new vtkSmartPointer<vtkDelaunay3D> [parameters.cluster_num];
    surfaceFilter = new vtkSmartPointer<vtkDataSetSurfaceFilter> [parameters.cluster_num];

    for(int i = 0; i < parameters.cluster_num; i++)
    {
        sphere[i] = vtkSphereSource::New();
        sphere[i]->SetThetaResolution(7);
        sphere[i]->SetPhiResolution(7);
        sphere[i]->Update();
        model[i] = vtkPolyData::New();
//		glyph[i] = vtkGlyph3D::New();
        pointMapper[i] = vtkPolyDataMapper::New();
        VP[i] = vtkProperty::New();

        HullMapper[i] = vtkPolyDataMapper::New();
        polydata[i] = vtkSmartPointer<vtkPolyData>::New();
        delaunay[i] = vtkSmartPointer<vtkDelaunay3D>::New();
        surfaceFilter[i] = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    }

    legend = vtkSmartPointer<vtkLegendBoxActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    triangleFilter = new vtkTriangleFilter *[parameters.cluster_num];
    massproperties = new vtkMassProperties *[parameters.cluster_num];
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        triangleFilter[i] = vtkTriangleFilter::New();
        massproperties[i] = vtkMassProperties::New();
    }

    lines = new vtkLineSource *[size - 1];
    linemappers = new vtkPolyDataMapper *[size - 1];
    lineActors = new vtkActor *[size - 1];
    for(int i = 0; i < size - 1; i++)
    {
        lines[i] = vtkLineSource::New();
        linemappers[i] = vtkPolyDataMapper::New();
        lineActors[i] = vtkActor::New();
    }

    Colors = new double *[parameters.cluster_num];
    for(int i = 0; i < parameters.cluster_num; i++)
        Colors[i] = new double [3];
    Generate_Colors(parameters.cluster_num, Colors);

    IDMaps = new int [parameters.index_size];
    Load_COR(filename);
};

void Image::Create_Legend()
{
    char i_str[10];
    legend->SetHeight(1);
    legend->SetDisplayPosition(1, 1);
    legend->SetNumberOfEntries(parameters.cluster_num);
    vtkSmartPointer<vtkSphereSource> legendSphereSource = vtkSmartPointer<vtkSphereSource>::New();

    vtkSmartPointer<vtkPolyData> legendSphere = legendSphereSource->GetOutput();
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        int_to_string(i + 1, i_str, 10);
        legend->SetEntry(i, legendSphere, i_str, Colors[i]);
    }
};

bool Image::Load_points_with_selection()
{
    String filename = fname;
    select_Points(filename);
    Load_Points();
    return true;
}

void Image::sort_volume()
{
    Create_HullsActors();
    sort_HullsActors();
    Load_Points();
}

void Image::Generate_Colors(int num, double **Colors)
{
    int type = parameters.color_type;
    if(type == 0)
    {
        double interval = (double) 6 / num;
        double x = 0;
        for(int i = 0; i < num; i++)
        {
            if(x < 1)
            {
                Colors[i][0] = 1;
                Colors[i][1] = x;
                Colors[i][2] = 0;
            } else
            if(x < 2)
            {
                Colors[i][0] = 2 - x;
                Colors[i][1] = 1;
                Colors[i][2] = 0;
            } else
            if(x < 3)
            {
                Colors[i][0] = 0;
                Colors[i][1] = 1;
                Colors[i][2] = x - 2;
            } else
            if(x < 4)
            {
                Colors[i][0] = 0;
                Colors[i][1] = 4 - x;
                Colors[i][2] = 1;
            } else
            if(x < 5)
            {
                Colors[i][0] = x - 4;
                Colors[i][1] = 0;
                Colors[i][2] = 1;
            } else
            if(x < 6)
            {
                Colors[i][0] = 1;
                Colors[i][1] = 0;
                Colors[i][2] = 6 - x;
            }
            x +=interval;
        }
    } else
    {
        double interval = (double) 2 / (num - 1);
        double x = 0;
        for(int i = 0; i < num; i++)
        {
            if(x < 1)
            {
                Colors[i][0] = 1 - x;
                Colors[i][1] = x;
                Colors[i][2] = 0;
            } else
            {
                Colors[i][0] = 0;
                Colors[i][1] = 2 - x;
                Colors[i][2] = x - 1;
            }
            x += interval;
        }
    }
/*
    Colors[0][0] = 0.5;
    Colors[0][1] = 0;
    Colors[0][2] = 0.5;

    Colors[1][0] = 0;
    Colors[1][1] = 1;
    Colors[1][2] = 0;

    Colors[2][0] = 0;
    Colors[2][1] = 0;
    Colors[2][2] = 0;

    Colors[3][0] = 0.5;
    Colors[3][1] = 0.5;
    Colors[3][2] = 0;

    Colors[4][0] = 1;
    Colors[4][1] = 0;
    Colors[4][2] = 0;

    Colors[5][0] = 0;
    Colors[5][1] = 0.5;
    Colors[5][2] = 0.5;

    Colors[6][0] = 1;
    Colors[6][1] = 0;
    Colors[6][2] = 1;

    Colors[7][0] = 0.5;
    Colors[7][1] = 0;
    Colors[7][2] = 0;

    Colors[8][0] = 0;
    Colors[8][1] = 0;
    Colors[8][2] = 1;

    Colors[9][0] = 0;
    Colors[9][1] = 1;
    Colors[9][2] = 1;

    Colors[10][0] = 0;
    Colors[10][1] = 0;
    Colors[10][2] = 0.5;

    Colors[11][0] = 0;
    Colors[11][1] = 0.5;
    Colors[11][2] = 0;

    Colors[12][0] = 1;
    Colors[12][1] = 0.55;
    Colors[12][2] = 0;

    Colors[13][0] = 1;
    Colors[13][1] = 0.85;
    Colors[13][2] = 0;

    Colors[14][0] = 0.5;
    Colors[14][1] = 0.5;
    Colors[14][2] = 0.5;*/

};

bool Image::Load_COR(String filename)
{
    File fcor(filename);
    fcor.seek(0);
    size = fcor.lines();
    fcor.seek(0);
    dim = fcor.cols();
    fcor.seek(0);
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < dim; j++)
        {
//            COR[dim * i + j] = (double) i;//-------
            fcor >> COR[dim * i + j];
        }
    }

    return true;
};

bool Image::select_Points(String filename)
{
    double *mean = new double [dim];
    double var1 = 0;
    double var2 = 0;
    int max_n1 = 0;
    int max_n2 = 0;
    int exclude_num = 0;
    int *index = NULL;

    index = parameters.cluster_index;

    int *numbers_in_cluster = new int[parameters.cluster_num];

    for(int i = 0; i < parameters.cluster_num; i++)
        numbers_in_cluster[i] = 0;
    for(int i = 0; i < parameters.index_size; i++)
    {
        if((int) index[i] < parameters.cluster_num && (int) index[i] >= 0)
            numbers_in_cluster[index[i]]++;
    }

    for(int i = 0; i < parameters.cluster_num; i++)
    {
        compute_mean(mean, index, i);
        for(int j = 0; j < dim; j++)
            centers[i * dim + j] = mean[j];
        compute_var(mean, index, i, var1, max_n1);
        selectpts[max_n1] = false;
        compute_var(mean, index, i, var2, max_n2);
        selectpts[max_n1] = true;
        exclude_num = 0;
        while(fabs((var2 - var1) / var1) > parameters.outlier_e && ((double) exclude_num / numbers_in_cluster[i]) < parameters.outlierpercent)
        {
            selectpts[max_n1] = false;
            var1 = var2;
            max_n1 = max_n2;
            selectpts[max_n1] = false;
            compute_var(mean, index, i, var2, max_n2);
            selectpts[max_n1] = true;
            exclude_num++;
        }
    }
    delete [] mean;
    return true;
};

void Image::compute_mean(double *mean, int *index, int i)
{
    int select_n = 0;

    for(int j = 0; j < dim; j++)
        mean[j] = 0;

    select_n = 0;
    for(int j = 0; j < mymin(parameters.index_size, size); j++)
        if(selectpts[j] && (int) index[j] == i)
        {
            for(int k = 0; k < dim; k++)
                mean[k] += COR[dim * j + k];
            select_n++;
        }

    for(int j = 0; j < dim; j++)
        mean[j] /= select_n;
};

void Image::compute_var(double *mean, int *index, int i, double &var, int &max_n)
{
    double dis = 0;
    double max = 0;
    int select_n = 0;
    var = 0;
    for(int j = 0; j < mymin(parameters.index_size, size); j++)
    {
        if(selectpts[j] && index[j] == i)
        {
            dis = 0;
            for(int k = 0; k < dim; k++)
                dis += (mean[k] - COR[dim * j + k]) * (mean[k] - COR[dim * j + k]);
            if(dis > max)
            {
                max = dis;
                max_n = j;
            }
            var += dis;
            select_n++;
        }
    }
    var /= select_n;
};

bool Image::Load_Points()
{
    String filename = fname;
    File fcor(filename);
    fcor.seek(0);
    size = fcor.lines();
    fcor.seek(0);
    dim = fcor.cols();
    fcor.seek(0);

    int *index = NULL;

    index = parameters.cluster_index;

    for(int i = 0; i < parameters.cluster_num; i++)
    {
        points[i] = vtkSmartPointer<vtkPoints>::New();
        Hull_points[i] = vtkSmartPointer<vtkPoints>::New();
    }

    double cor[3] = {0};
    double away = parameters.sep_factor;
    int idx = 0;
    for(int i = 0; i < mymin(size, parameters.index_size); i++)
    {
        if(index[i] >= 0 && selectpts[i])
        {
            IDMaps[idx] = i;
            idx++;
            for(int j = 0; j < dim; j++)
                cor[j] = COR[i * dim + j];
            if(dim == 2)
            {
                points[index[i]]->InsertNextPoint(cor[0] + centers[index[i] * dim + 0] * away, cor[1] + centers[index[i] * dim + 0] * away, cor[2] + ranks[index[i]] / parameters.cluster_num);
                Hull_points[index[i]]->InsertNextPoint(cor[0] + centers[index[i] * dim + 0] * away, cor[1] + centers[index[i] * dim + 0] * away, cor[2] + ranks[index[i]] / parameters.cluster_num);
                Hull_points[index[i]]->InsertNextPoint(cor[0] + centers[index[i] * dim + 0] * away, cor[1] + centers[index[i] * dim + 0] * away, cor[2] + ranks[index[i]] / parameters.cluster_num + 1.0 / parameters.cluster_num);
            } else
            {
                points[index[i]]->InsertNextPoint(cor[0] + centers[index[i] * dim + 0] * away, cor[1] + centers[index[i] * dim + 1] * away, cor[2] + centers[index[i] * dim + 2] * away);
                Hull_points[index[i]]->InsertNextPoint(cor[0] + centers[index[i] * dim + 0] * away, cor[1] + centers[index[i] * dim + 1] * away, cor[2] + centers[index[i] * dim + 2] * away);
            }
        }
    }
    return true;
};

void Image::Create_PointsActors()
{
    vtkSmartPointer<PointSource> *PSs = new vtkSmartPointer<PointSource> [parameters.cluster_num];
    vtkSmartPointer<vtkIdFilter> *IDFs = new vtkSmartPointer<vtkIdFilter> [parameters.cluster_num];
    vtkSmartPointer<vtkDataSetSurfaceFilter> *SFs = new vtkSmartPointer<vtkDataSetSurfaceFilter> [parameters.cluster_num];

    for(int i = 0; i < parameters.cluster_num; i++)
    {
        PSs[i] = vtkSmartPointer<PointSource>::New();
        PSs[i]->SetPoints(points + i, 1);
        PSs[i]->Update();

        IDFs[i] = vtkSmartPointer<vtkIdFilter>::New();
        IDFs[i]->SetInputConnection(PSs[i]->GetOutputPort());
        IDFs[i]->SetIdsArrayName("OriginalIds");
        IDFs[i]->Update();

        SFs[i] = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        SFs[i]->SetInputConnection(IDFs[i]->GetOutputPort());
        SFs[i]->Update();

        //model[i] = SFs[i]->GetOutput();

      #if VTK_MAJOR_VERSION <= 5
        pointMapper[i]->SetInputConnection(model[i]->GetProducerPort());
      #else
        //pointMapper[i]->SetInputData(model[i]);
        pointMapper[i]->SetInputData(SFs[i]->GetOutput());
      #endif
        pointMapper[i]->ScalarVisibilityOff();
        pointMapper[i]->Update();

        PointsActors[i] = vtkSmartPointer<vtkActor>::New();
        PointsActors[i]->SetMapper(pointMapper[i]);
        VP[i]->SetColor(Colors[i]);
        VP[i]->SetPointSize(parameters.point_size);
        PointsActors[i]->SetProperty(VP[i]);
/*
        PointsActors[i] = vtkSmartPointer<vtkActor>::New();
        model[i]->SetPoints(points[i]);
        glyph[i]->SetInputData(model[i]);
        glyph[i]->SetSourceData(sphere[i]->GetOutput());
        glyph[i]->SetVectorModeToUseNormal();
        glyph[i]->SetScaleModeToScaleByVector();
        glyph[i]->SetScaleFactor(parameters.point_size);
        glyph[i]->Update();
        pointMapper[i]->SetInputData(glyph[i]->GetOutput());
        pointMapper[i]->Update();
        PointsActors[i]->SetMapper(pointMapper[i]);
        VP[i]->SetColor(Colors[i]);
        PointsActors[i]->SetProperty(VP[i]);
        */
    }

    delete [] PSs;
    delete [] IDFs;
    delete [] SFs;
};

void Image::Create_HullsActors()
{
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        HullsActors[i] = vtkSmartPointer<vtkActor>::New();
        polydata[i]->SetPoints(Hull_points[i]);
        delaunay[i]->SetInputData(polydata[i]);
        delaunay[i]->Update();

        surfaceFilter[i]->SetInputConnection(delaunay[i]->GetOutputPort());
        surfaceFilter[i]->Update();
        HullMapper[i]->SetInputData(surfaceFilter[i]->GetOutput());
        HullsActors[i]->SetMapper(HullMapper[i]);
        VP[i]->SetColor(Colors[i]);
        HullsActors[i]->SetProperty(VP[i]);

        triangleFilter[i]->SetInputConnection(surfaceFilter[i]->GetOutputPort());
        massproperties[i]->SetInputConnection(triangleFilter[i]->GetOutputPort());
        volumes[i] = massproperties[i]->GetVolume();
    }
};

void Image::Create_LinesActors()
{
    for(int i = 0; i < size - 1; i++)
    {
        lines[i]->SetPoint1(points[0]->GetPoint(i));
        lines[i]->SetPoint2(points[0]->GetPoint(i + 1));
        lines[i]->SetResolution(20);
        lines[i]->Update();
        linemappers[i]->SetInputData(lines[i]->GetOutput());
        linemappers[i]->Update();
        lineActors[i]->SetMapper(linemappers[i]);
        lineActors[i]->GetProperty()->SetDiffuseColor(Colors[0]);
    }
};

void Image::sort_HullsActors()
{
    double c = 0;
    int *order = new int[parameters.cluster_num];

    for(int i = 0; i < parameters.cluster_num; i++)
        order[i] = i;

    for(int i = 0; i < parameters.cluster_num; i++)
    {
        for(int j = i; j < parameters.cluster_num; j++)
        {
            if(volumes[i] < volumes[j])
            {
                c = volumes[i];
                volumes[i] = volumes[j];
                volumes[j] = c;
                c = order[i];
                order[i] = order[j];
                order[j] = c;
            }
        }
    }

    for(int i = 0; i < parameters.cluster_num; i++)
        ranks[order[i]] = i;
    delete [] order;
};

void Image::Plot_Points()
{

    vtkProperty *allVP = vtkProperty::New();

    double Colors[] = {1,1,1};
    allVP->SetColor(Colors);
    allVP->SetPointSize(1);

    vtkSmartPointer<PointSource> allpointSource =
      vtkSmartPointer<PointSource>::New();
    allpointSource->SetPoints(points, parameters.cluster_num);
    allpointSource->Update();

    vtkSmartPointer<vtkIdFilter> allidFilter =
      vtkSmartPointer<vtkIdFilter>::New();
    allidFilter->SetInputConnection(allpointSource->GetOutputPort());
    allidFilter->SetIdsArrayName("OriginalIds");
    allidFilter->Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> allsurfaceFilter =
      vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    allsurfaceFilter->SetInputConnection(allidFilter->GetOutputPort());
    allsurfaceFilter->Update();

    vtkPolyData* allinput = allsurfaceFilter->GetOutput();

    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> allmapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();

  #if VTK_MAJOR_VERSION <= 5
    allmapper->SetInputConnection(allinput->GetProducerPort());
  #else
    allmapper->SetInputData(allinput);
  #endif
    allmapper->ScalarVisibilityOff();

    vtkSmartPointer<vtkActor> allactor =
      vtkSmartPointer<vtkActor>::New();
    allactor->SetMapper(allmapper);
    allactor->SetProperty(allVP);

    renderWindow->AddRenderer(renderer);

    renderWindow->SetSize(500, 500);

    renderWindow->SetPosition(100, 100);

    vtkSmartPointer<vtkAreaPicker> areaPicker =
      vtkSmartPointer<vtkAreaPicker>::New();
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetPicker(areaPicker);

    renderWindowInteractor->SetRenderWindow(renderWindow);

    Create_Legend();
    renderer->AddActor(allactor);
    renderer->AddActor(legend);
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        renderer->AddActor(PointsActors[i]);
    }

    renderer->SetBackground(1,1,1);

    renderWindow->Render();

    vtkSmartPointer<InteractorStyle> allstyle = vtkSmartPointer<InteractorStyle>::New();
    allstyle->SetPoints(allinput);
    allstyle->SetIDMaps(IDMaps);
    allstyle->SetOutputPts(SelectedIDs);
    allstyle->SetSelectedPointsSize(parameters.point_size + 1);
    renderWindowInteractor->SetInteractorStyle( allstyle );

    renderer->GetActiveCamera()->Zoom(parameters.zoom);

//    if(dim == 2)
//        renderWindowInteractor->EnableRenderOff();

    if(! paras::plot_makemovie)
        renderWindowInteractor->Start();
};

void Image::Plot_Hulls()
{
    renderWindow->AddRenderer(renderer);

    renderWindow->SetSize(500, 500);

    renderWindow->SetPosition(100, 100);

    renderWindowInteractor->SetRenderWindow(renderWindow);

    Create_Legend();
    renderer->AddActor(legend);
    for(int i = 0; i < parameters.cluster_num; i++)
        renderer->AddActor(HullsActors[i]);

    renderer->SetBackground(1,1,1);

    renderWindow->Render();

    renderer->GetActiveCamera()->Zoom(parameters.zoom);

    if(! paras::plot_makemovie)
    renderWindowInteractor->Start();
};

void Image::Plot_Lines()
{
    renderWindow->AddRenderer(renderer);

    renderWindow->SetSize(499, 499);

    renderWindow->SetPosition(100, 100);

    renderWindowInteractor->SetRenderWindow(renderWindow);

    for(int i = 0; i < size - 1; i++)
        renderer->AddActor(lineActors[i]);

    renderer->SetBackground(0,0,0);

    renderWindow->Render();

    vtkCubeAxesActor2D *axes = vtkCubeAxesActor2D::New();

    vtkLineSource *rangesquare = vtkLineSource::New();
    double *arr0 = points[0]->GetPoint(0);
    double xstart = arr0[0], ystart = arr0[1], zstart = arr0[2];
    double *arr1 = points[0]->GetPoint(size - 1);
    double xend = arr1[0], yend = arr1[1], zend = arr1[2];
    double xmove = (xend - xstart) * 0.05;
    double ymove = (yend - ystart) * 0.05;

    rangesquare->SetPoint1(xstart - xmove, ystart - ymove, zstart);
    rangesquare->SetPoint2(xend + xmove, yend + ymove, zend);
    rangesquare->SetResolution(20);
    rangesquare->Update();

    axes ->SetInputData(rangesquare->GetOutput());
    axes->SetCamera(renderer->GetActiveCamera());
    axes->SetLabelFormat("%6.2f");
    axes->SetFlyModeToOuterEdges();
    axes->SetFontFactor(1.2);
    axes->SetZLabel(" ");
    axes->SetNumberOfLabels(5);
    axes->ScalingOn();
    renderer->AddViewProp(axes);
    renderer->GetActiveCamera()->Zoom(0.85);;

    renderWindowInteractor->UpdateSize(500, 500);
    renderWindowInteractor->Render();
    renderWindowInteractor->EnableRenderOff();
    renderWindowInteractor->Start();
    axes->Delete();
    rangesquare->Delete();
};

void Image::make_movie()
{
    double el = 0;
    double sign = 1;
    double ratio = 0.5;
    renderer->ResetCamera();
    renderWindow->SetSize(parameters.moviesizex, parameters.moviesizey);
    renderer->GetActiveCamera()->Zoom(parameters.zoom);

    vtkWindowToImageFilter *filter = vtkWindowToImageFilter::New();

    filter->SetInput(renderer->GetVTKWindow());
/*
    vtkSmartPointer<vtkMPEG2Writer> writer = vtkSmartPointer<vtkMPEG2Writer>::New();
//    vtkSmartPointer<vtkFFMPEGWriter> writer = vtkSmartPointer<vtkFFMPEGWriter>::New();
    writer->SetInputConnection(filter->GetOutputPort());
    writer->SetFileName("test.avi");

    writer->Start();
    for(int i = 0; i < parameters.frame_num; i++)
    {
        renderWindow->Render();
        renderer->GetActiveCamera()->Azimuth(3 * ratio);
        renderer->GetActiveCamera()->Roll(ratio);
        el = el + sign;
        if(el > 89.5 || el < -89.5)
        {
            el = el - 2 * sign * ratio;
            sign = -sign;
        }
        renderer->GetActiveCamera()->Elevation(sign * ratio);

        filter->Modified();
        writer->Write();
    }
    writer->End();
    */
};

void Image::destructor()
{
//    std::cout << "h1" << std::endl;//---
    //delete [] PointsActors;
//    std::cout << "h2" << std::endl;//---
    for(int i = 0; i < parameters.cluster_num; i++)
    {
        sphere[i]->Delete();
        model[i]->Delete();

//        glyph[i]->Delete();
        pointMapper[i]->Delete();
        VP[i]->Delete();

        HullMapper[i]->Delete();
        delete [] Colors[i];
        triangleFilter[i]->Delete();
        massproperties[i]->Delete();
    }
//    std::cout << "h3" << std::endl;//---
    for(int i = 0; i < size - 1; i++)
    {
        lines[i]->Delete();
        linemappers[i]->Delete();
        lineActors[i]->Delete();
    }
//    std::cout << "h4" << std::endl;//---
    delete [] lines;
    delete [] linemappers;
    delete [] lineActors;

//    std::cout << "h5" << std::endl;//---
    delete [] PointsActors; // in VTK 6.3.0 version, deleting this and model[i] together first causes garbage collector errors.
    delete [] points;
    delete [] Hull_points;

//    std::cout << "h6" << std::endl;//---
    delete [] sphere;
    delete [] model;
//	delete [] glyph;
    delete [] pointMapper;
    delete [] VP;
    delete [] HullsActors;
    delete [] HullMapper;
    delete [] polydata;
    delete [] delaunay;
    delete [] surfaceFilter;
//    std::cout << "h7" << std::endl;//---
    delete [] Colors;
    delete [] volumes;
    delete [] ranks;
    delete [] triangleFilter;
    delete [] massproperties;
    delete [] COR;
    delete [] selectpts;
    delete [] centers;
    delete [] IDMaps;



//    renderer->Delete();//---
//    renderWindow->Delete();
//    renderWindowInteractor->Delete();//---
//    renderWindowInteractor->UnRegister(0);//----
//    std::cout << "h8" << std::endl;//---
};

void Image::int_to_string(int n, char *str, int length)
{
    int *num = new int[length - 1];
    int digits = 0;
    for(int i = 0; i < length - 1; i++)
        num[i] = 0;

    for(int i = 0; i < length - 1; i++)
    {
        num[i] = n % 10;
        n = n / 10;
        if(n == 0)
        {
            digits = i;
            break;
        }
    }

    for(int i = 0; i < digits + 1; i++)
        str[i] = num[digits - i] + 48;
    str[digits + 1] = 0;

    delete [] num;
};

#endif
