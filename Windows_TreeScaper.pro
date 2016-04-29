#-------------------------------------------------
#
# Project created by QtCreator 2010-07-11T13:11:32
#
#-------------------------------------------------

CONFIG += qt
QT += core gui widgets printsupport
QT += webkitwidgets

TARGET = TreeScaper
TEMPLATE = app
#CONFIG += staticlib

#CONFIG += debug
#CONFIG -= app_bundle
#CONFIG += console
#QT -= gui

SOURCES += main.cpp\
    randgen.cpp \
    treescaper.cpp \
    wstring.cpp \
    warray.cpp \
    wmapping.cpp \
    wmix.cpp \
    wfile.cpp \
    wimport_form.cpp \
    wDimEst.cpp \
    wmatrix.cpp \
    wNLDR.cpp \
    wImage.cpp \
    setting.cpp \
    wNLDRthread.cpp \
    wPlotthread.cpp \
    wPlotTreethread.cpp \
    wPlotTree.cpp \
    wDimestthread.cpp \
    wPlotdimthread.cpp \
    hash.cc \
    hashfunc.cc \
    label-map.cc \
#    newick.c \
    hungarian.c \
#    DIST.cpp \
#    Vconverter.cpp \
#    Bipartition.cpp \
    Sparse_matrix.cpp \
#    Affinity.cpp \
    PointsSource.cpp \
    community.cpp \
    graph.cpp \
    greedy_louvain.cpp \
    info.cpp \
    MersenneTwister.cpp \
    slicer.cpp \
    PlotGraph.cpp \
    Trees.cpp \
    TreeOPE.cpp \
    plot2d.cpp \
    qcustomplot.cpp \
    wCommunitythread.cpp \
    woutput.cpp \
    plot2d_awty.cpp \
    ClusterForest.cpp \
    ClusterInstance.cpp \
    Forest.cpp \
    rspr.cpp \
    SiblingPair.cpp \
    SPRLCA.cpp \
    SPRNode.cpp \
    UndoMachine.cpp

HEADERS  += treescaper.h \
    randgen.h \
    warray.h \
    wDimEst.h \
    wfile.h \
    wimport_form.h \
    wmapping.h \
    wmatrix.h \
    wmix.h \
    wNLDR.h \
    wstring.h \
    wImage.h \
    setting.h \
    wNLDRthread.h \
    wPlotthread.h \
    wPlotTreethread.h \
    wPlotTree.h \
    wDimestthread.h \
    wPlotdimthread.h \
    wdef.h \
    hash.hh \
    hashfunc.hh \
    label-map.hh \
#    newick.h \
    randomc.h \
    hungarian.h \
#    DIST.h \
#    Vconverter.h \
#    Bipartition.h \
    Sparse_matrix.h \
#    Affinity.h \
    PointsSource.h \
    community.h \
    graph.h \
    greedy_louvain.h \
    info.h \
    MersenneTwister.h \
    slicer.h \
    PlotGraph.h \
    Trees.h \
    TreeOPE.h \
    plot2d.h \
    qcustomplot.h \
    wCommunitythread.h \
    woutput.h \
    plot2d_awty.h \
    ClusterForest.hpp \
    ClusterInstance.hpp \
    Forest.hpp \
    rspr.hpp \
    SiblingPair.hpp \
    SPRLCA.hpp \
    SPRNode.cpp \
    UndoMachine.hpp

FORMS    += treescaper.ui \
    setting.ui \
    plot2d.ui

# for WH Windows system --WH
VTKPATH = D:\VTK\BIN

INCLUDEPATH += D:\VTK\clapack-3.2.1-CMAKE\INCLUDE

INCLUDEPATH += $${VTKPATH}\include\vtk-6.3

LIBS += D:\VTK\clapack\SRC\liblapack.a
LIBS += D:\VTK\clapack\BLAS\SRC\libblas.a
LIBS += D:\VTK\clapack\F2CLIBS\libf2c\libf2c.a

LIBS += $${VTKPATH}\bin\libvtkalglib-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkChartsCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonColor-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonComputationalGeometry-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonDataModel-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonExecutionModel-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonMath-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonMisc-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonSystem-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkCommonTransforms-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkDICOMParser-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkDomainsChemistry-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkexoIIc-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkexpat-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersAMR-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersExtraction-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersFlowPaths-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersGeneral-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersGeneric-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersGeometry-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersHybrid-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersHyperTree-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersImaging-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersModeling-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersParallel-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersParallelImaging-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersProgrammable-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersSelection-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersSMP-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersSources-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersStatistics-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersTexture-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkFiltersVerdict-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkfreetype-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkftgl-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkGeovisCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkgl2ps-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingColor-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingFourier-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingGeneral-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingHybrid-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingMath-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingMorphological-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingSources-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingStatistics-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkImagingStencil-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkInfovisCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkInfovisLayout-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkInteractionImage-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkInteractionStyle-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkInteractionWidgets-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOAMR-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOEnSight-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOExodus-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOExport-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOGeometry-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOImage-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOImport-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOInfovis-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOLegacy-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOLSDyna-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOMINC-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOMovie-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIONetCDF-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOParallel-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOParallelXML-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOPLY-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOSQL-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOVideo-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOXML-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkIOXMLParser-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkjpeg-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkjsoncpp-6.3.dll
LIBS += $${VTKPATH}\bin\libvtklibxml2-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkmetaio-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkNetCDF-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkNetCDF_cxx-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkoggtheora-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkParallelCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkpng-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkproj4-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingAnnotation-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingContext2D-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingContextOpenGL-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingFreeType-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingGL2PS-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingImage-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingLabel-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingLIC-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingLOD-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingOpenGL-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingVolume-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkRenderingVolumeOpenGL-6.3.dll
LIBS += $${VTKPATH}\bin\libvtksys-6.3.dll
LIBS += $${VTKPATH}\bin\libvtktiff-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkverdict-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkViewsContext2D-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkViewsCore-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkViewsInfovis-6.3.dll
LIBS += $${VTKPATH}\bin\libvtkzlib-6.3.dll
