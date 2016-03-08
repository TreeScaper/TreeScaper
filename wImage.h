
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

//wplot.h
// using vtk to plot points
//                by whuang 15/May/2010

#ifndef WIMAGE_H
#define WIMAGE_H

#include "wfile.h"
#include "wstring.h"
#include "warray.cpp"
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

struct image_parameters{
    int cluster_num; // 15
    int *cluster_index; //
    int index_size; //
    double outlier_e; // 0.01
    double outlierpercent; // 0.5
    double sep_factor; // 0
    double point_size; // 5
    int moviesizex; // 700
    int moviesizey; // 700
    double zoom; // 1.35
    int color_type; // 0
    int frame_num; // 100
};

class Image{
public:
	Image(){};

    ~Image(){destructor();};

	Image(String filename)
	{
		Initialize_Image(filename);
	};

    void Initialize_Image(String filename);

    void Set_SelectedIDs(Array<int> &selected_ids){SelectedIDs = &selected_ids;}

    void Generate_Colors(int num, double **Colors);

	bool Load_COR(String filename);

    bool Load_Points();

    bool Load_points_with_selection();

    void sort_volume();

    void Create_Legend();

	void Create_PointsActors();

	void Create_HullsActors();

    void Create_LinesActors();

	void sort_HullsActors();

	bool select_Points(String filename);

    void compute_mean(double *mean, int *index, int i);

    void compute_var(double *mean, int *index, int i, double &var, int &max_n);

    void Plot_Points();

	void Plot_Hulls();

    void Plot_Lines();

	void make_movie();

	void destructor();

    image_parameters parameters;

private:
    void int_to_string(int n, char *str, int length);

    String fname;
	double **Colors;
	double *COR;
	bool *selectpts;
	double *centers;
    Array<int> *SelectedIDs;
    int *IDMaps;
	int dim;
	int size;
	vtkSmartPointer<vtkPoints> *points;
	vtkSmartPointer<vtkPoints> *Hull_points;

	vtkSmartPointer<vtkActor> *PointsActors;
	vtkSphereSource **sphere;
	vtkPolyData **model;
//	vtkGlyph3D **glyph;
	vtkPolyDataMapper **pointMapper;
	vtkProperty **VP;

    vtkSmartPointer<vtkLegendBoxActor> legend;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

	vtkSmartPointer<vtkActor> *HullsActors;
	vtkPolyDataMapper **HullMapper;
	vtkSmartPointer<vtkPolyData> *polydata;
	vtkSmartPointer<vtkDelaunay3D> *delaunay;
	vtkSmartPointer<vtkDataSetSurfaceFilter> *surfaceFilter;
	vtkTriangleFilter **triangleFilter;
	vtkMassProperties **massproperties;
	int *ranks;
	double *volumes;

    vtkLineSource **lines;
    vtkPolyDataMapper **linemappers;
    vtkActor **lineActors;
};

#endif
