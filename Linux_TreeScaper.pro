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

#CONFIG += debug
#CONFIG -= app_bundle

#CONFIG += console

#QT -= gui

SOURCES += main.cpp\
    randgen.cpp \
    treescaper.cpp \
    wstring.cpp \
    warray.cpp \
#    wmapping.cpp \
#    wmix.cpp \
    wfile.cpp \
#    wimport_form.cpp \
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
    UndoMachine.cpp \
    NLDRsetting.cpp

HEADERS  += treescaper.h \
    randgen.h \
    warray.h \
    wDimEst.h \
    wfile.h \
#    wimport_form.h \
#    wmapping.h \
    wmatrix.h \
#    wmix.h \
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
    UndoMachine.hpp \
    NLDRsetting.h

FORMS    += treescaper.ui \
    setting.ui \
    plot2d.ui \
    NLDRsetting.ui

#VTKBINPATH = /Users/GZhou/Desktop/New/bin
#VTKVTKPATH = /Users/GZhou/Desktop/New/VTK6.0.0
#CLAPACKPATH = /Users/GZhou/Desktop/TreeScaper/clapack-3.2.1-CMAKE

#VTKBINPATH = /Users/zhouguifang/Desktop/TreeScaperEvri/bin
#VTKVTKPATH = /Users/zhouguifang/Desktop/TreeScaperEvri/VTK6.0.0
#CLAPACKPATH = /Users/zhouguifang/Desktop/TreeScaperEvri/CLAPACK-3.2.1

#VTKBINPATH = /Users/whuang/TreeScaper_dev/bin6.0.0
#VTKVTKPATH = /Users/whuang/TreeScaper_dev/VTK6.0.0
#CLAPACKPATH = /Users/whuang/TreeScaper_dev/clapack-3.2.1-CMAKE


# for WH linux system --WH
VTKBINPATH = /home/whuang/VTK/Bin-6.3.0
VTKVTKPATH = /home/whuang/VTK/VTK-6.3.0
CLAPACKPATH = /home/whuang/CLAPACK-3.2.1

## for Guifang linux system --WH
#VTKBINPATH = /home/zhouzhou/vtk6.0.0/bin6.0.0
#VTKVTKPATH = /home/zhouzhou/vtk6.0.0/VTK6.0.0
#CLAPACKPATH = /home/zhouzhou/Documents/c++programs/CLAPACK-3.2.1

INCLUDEPATH += $${CLAPACKPATH}/INCLUDE

#INCLUDEPATH += /home/zhouzhou/vtk6.0.0/VTK6.0.0/Common/Core
#INCLUDEPATH += /home/zhouzhou/vtk6.0.0/bin6.0.0/Common/Core

INCLUDEPATH += $${VTKVTKPATH}/Common/Math
INCLUDEPATH += $${VTKBINPATH}/Common/Math
INCLUDEPATH += $${VTKVTKPATH}/Testing/Rendering
INCLUDEPATH += $${VTKBINPATH}/Testing/Rendering
INCLUDEPATH += $${VTKBINPATH}/Rendering/Volume
INCLUDEPATH += $${VTKBINPATH}/Filters/Hybrid
INCLUDEPATH += $${VTKBINPATH}/Interaction/Widgets
INCLUDEPATH += $${VTKVTKPATH}/Views/Context2D
INCLUDEPATH += $${VTKBINPATH}/Views/Context2D
INCLUDEPATH += $${VTKVTKPATH}/Views/Core
INCLUDEPATH += $${VTKBINPATH}/Views/Core
INCLUDEPATH += $${VTKVTKPATH}/Views/Infovis
INCLUDEPATH += $${VTKBINPATH}/Views/Infovis
INCLUDEPATH += $${VTKVTKPATH}/Infovis/Core
INCLUDEPATH += $${VTKVTKPATH}/Infovis/Layout
INCLUDEPATH += $${VTKVTKPATH}/Infovis/BoostGraphAlgorithms
INCLUDEPATH += $${VTKBINPATH}/Infovis/Core
INCLUDEPATH += $${VTKBINPATH}/Infovis/Layout
INCLUDEPATH += $${VTKBINPATH}/Rendering/OpenGL
INCLUDEPATH += $${VTKBINPATH}/Rendering/Context2D
INCLUDEPATH += $${VTKVTKPATH}/Rendering/Context2D
INCLUDEPATH += $${VTKVTKPATH}/Common/Transforms
INCLUDEPATH += $${VTKBINPATH}/Common/Transforms
INCLUDEPATH += $${VTKVTKPATH}/Filters
INCLUDEPATH += $${VTKBINPATH}
INCLUDEPATH += $${VTKBINPATH}/Common
INCLUDEPATH += $${VTKBINPATH}/Common/Core
INCLUDEPATH += $${VTKVTKPATH}/Common/Core
INCLUDEPATH += $${VTKBINPATH}/Common/DataModel
INCLUDEPATH += $${VTKVTKPATH}/Common/DataModel
INCLUDEPATH += $${VTKVTKPATH}/Common/ExecutionModel
INCLUDEPATH += $${VTKBINPATH}/Common/ExecutionModel
INCLUDEPATH += $${VTKBINPATH}/Utilities
INCLUDEPATH += $${VTKBINPATH}/Rendering
INCLUDEPATH += $${VTKBINPATH}/Rendering/VolumeOpenGL
INCLUDEPATH += $${VTKBINPATH}/Charts/Core
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/alglib
INCLUDEPATH += $${VTKVTKPATH}/IO
INCLUDEPATH += $${VTKVTKPATH}/Rendering/Annotation
INCLUDEPATH += $${VTKVTKPATH}/Geovis/Core
INCLUDEPATH += $${VTKVTKPATH}/Views
INCLUDEPATH += $${VTKVTKPATH}/IO/XML
INCLUDEPATH += $${VTKBINPATH}/IO/XML
INCLUDEPATH += $${VTKBINPATH}/IO/Geometry
INCLUDEPATH += $${VTKVTKPATH}/Filters/Geometry
INCLUDEPATH += $${VTKBINPATH}/Common/ExecutionModel
INCLUDEPATH += $${VTKBINPATH}/Filters/Geometry
INCLUDEPATH += $${VTKBINPATH}/Rendering/Core
INCLUDEPATH += $${VTKBINPATH}/Filters/Extraction
INCLUDEPATH += $${VTKBINPATH}/Filters/Statistics
INCLUDEPATH += $${VTKBINPATH}/IO/Image
INCLUDEPATH += $${VTKVTKPATH}/Filters/Sources
INCLUDEPATH += $${VTKBINPATH}/Filters/Sources
INCLUDEPATH += $${VTKBINPATH}/Filters/Imaging
INCLUDEPATH += $${VTKVTKPATH}/Imaging/Hybrid
INCLUDEPATH += $${VTKBINPATH}/Imaging/Hybrid
INCLUDEPATH += $${VTKVTKPATH}/Common/Misc
INCLUDEPATH += $${VTKBINPATH}/Common/Misc
INCLUDEPATH += $${VTKVTKPATH}/IO/Movie
INCLUDEPATH += $${VTKBINPATH}/IO/Movie
INCLUDEPATH += $${VTKVTKPATH}/IO/FFMPEG
INCLUDEPATH += $${VTKVTKPATH}/Imaging/Core
INCLUDEPATH += $${VTKBINPATH}/Imaging/Core
INCLUDEPATH += $${VTKVTKPATH}/Imaging/Sources
INCLUDEPATH += $${VTKBINPATH}/Imaging/Sources
INCLUDEPATH += $${VTKBINPATH}/Rendering/Annotation
INCLUDEPATH += $${VTKBINPATH}/IO/Export
INCLUDEPATH += $${VTKBINPATH}/Rendering/FreeType
INCLUDEPATH += $${VTKBINPATH}/Rendering/Label
INCLUDEPATH += $${VTKVTKPATH}/Interaction/Widgets
INCLUDEPATH += $${VTKVTKPATH}/Rendering/Core
INCLUDEPATH += $${VTKVTKPATH}/Charts/Core
INCLUDEPATH += $${VTKVTKPATH}/Rendering/OpenGL/Testing/Cxx
INCLUDEPATH += $${VTKVTKPATH}/IO/Core
INCLUDEPATH += $${VTKVTKPATH}/IO/Image
INCLUDEPATH += $${VTKVTKPATH}/Filters/General
INCLUDEPATH += $${VTKVTKPATH}/Filters/Modeling
INCLUDEPATH += $${VTKBINPATH}/Filters/Modeling
INCLUDEPATH += $${VTKVTKPATH}/Filters/Core
INCLUDEPATH += $${VTKBINPATH}/Filters/Core
INCLUDEPATH += $${VTKVTKPATH}/Filters/Generic
INCLUDEPATH += $${VTKVTKPATH}/Common/Core
INCLUDEPATH += $${VTKVTKPATH}/Utilities
INCLUDEPATH += $${VTKVTKPATH}/Utilities/KWSys
INCLUDEPATH += $${VTKBINPATH}/Utilities/KWSys
#INCLUDEPATH += $${VTKBINPATH}/Utilities/KWSys/vtksys/ios
INCLUDEPATH += $${VTKVTKPATH}/Common/Core/Testing/Cxx
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/libproj4
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/libproj4/vtklibproj4
INCLUDEPATH += $${VTKBINPATH}/Utilities/DICOMParser
INCLUDEPATH += $${VTKVTKPATH}/Utilities/DICOMParser
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/freetype/vtkfreetype/include
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/freetype/vtkfreetype/include
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/netcdf/vtknetcdf
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/libsrc
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/libsrc4
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/include
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/libdispatch
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/exodusII/vtkexodusII/include
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/exodusII/vtkexodusII/include
INCLUDEPATH += $${VTKBINPATH}/Utilities/MaterialLibrary
INCLUDEPATH += $${VTKVTKPATH}/Utilities/MaterialLibrary
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/verdict/vtkverdict
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/verdict/vtkverdict
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/alglib
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/alglib
INCLUDEPATH += $${VTKVTKPATH}/Interaction/Style
INCLUDEPATH += $${VTKBINPATH}/Interaction/Style
INCLUDEPATH += $${VTKBINPATH}/ThirdParty/ftgl
INCLUDEPATH += $${VTKVTKPATH}/ThirdParty/ftgl/src
INCLUDEPATH += $${VTKVTKPATH}/Filters/Extraction
INCLUDEPATH += $${VTKBINPATH}/Filters/General

LIBS += $${CLAPACKPATH}/lapack_LINUX.a
LIBS += $${CLAPACKPATH}/blas_LINUX.a
LIBS += $${CLAPACKPATH}/F2CLIBS/libf2c.a
LIBS += -lm



LIBS += $${VTKBINPATH}/lib/libvtkFiltersGeometry-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsInfovis-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingLabel-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkInfovisLayout-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingHybrid-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsContext2D-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkInteractionWidgets-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingContext2D-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingOpenGL-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingAnnotation-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkInteractionStyle-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingVolumeOpenGL-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersGeneral-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingVolume-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingFreeType-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkftgl-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOXMLParser-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersExtraction-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonExecutionModel-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonDataModel-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOImage-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersSMP-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOParallelXML-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingContextOpenGL-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingLIC-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtktiff-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonTransforms-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkChartsCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonColor-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonComputationalGeometry-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonMath-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonMisc-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonSystem-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkDICOMParser-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkDomainsChemistry-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersAMR-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersFlowPaths-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersGeneric-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersHybrid-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersHyperTree-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersImaging-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersModeling-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersParallel-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersParallelImaging-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersProgrammable-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersSelection-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersSources-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersStatistics-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersTexture-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersVerdict-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkGeovisCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOAMR-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOEnSight-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOExodus-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOExport-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOGeometry-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOImport-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOInfovis-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOLSDyna-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOLegacy-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOMINC-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOMovie-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIONetCDF-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOPLY-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOParallel-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOSQL-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOVideo-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkIOXML-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingColor-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingFourier-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingGeneral-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingMath-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingMorphological-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingSources-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingStatistics-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingStencil-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkInfovisCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkInteractionImage-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkNetCDF-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkNetCDF_cxx-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkParallelCore-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingGL2PS-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingImage-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingLOD-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkTestingGenericBridge-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkTestingIOSQL-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkTestingRendering-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkalglib-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkexoIIc-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkexpat-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkfreetype-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkgl2ps-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkhdf5-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkhdf5_hl-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkjpeg-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkjsoncpp-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtklibxml2-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkmetaio-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkoggtheora-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkpng-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkproj4-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtksqlite-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtksys-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkverdict-6.3.a
LIBS += $${VTKBINPATH}/lib/libvtkzlib-6.3.a

## FOR LINUX system
INCLUDEPATH += /usr/include/GL
INCLUDEPATH += /usr/include/GL/internal
INCLUDEPATH += /usr/include/X11

LIBS += -ldl
LIBS += -lX11
LIBS += /usr/lib/x86_64-linux-gnu/libXt.so.6.0.0
LIBS += /usr/lib/x86_64-linux-gnu/libX11.a

## FOR MAC system
#INCLUDEPATH += /Developer/SDKs/MacOSX10.8.sdk/usr/X11/include/GL
#INCLUDEPATH += /Developer/SDKs/MacOSX10.8.sdk/usr/X11/include/GL/internal
#INCLUDEPATH += /Developer/SDKs/MacOSX10.8.sdk/usr/X11/include/X11
#LIBS += -framework Cocoa
#LIBS += -framework IOKit
#LIBS += -framework OpenGL
#LIBS += -framework AGL
#LIBS += -framework QtGui
#LIBS += -framework QtCore



#LIBS += -L/home/whuang/VTK/Bin-6.3.0/lib/ \
#-lvtkFiltersGeometry-6.3 \
#-lvtkViewsInfovis-6.3 \
#-lvtkRenderingLabel-6.3 \
#-lvtkInfovisLayout-6.3 \
#-lvtkImagingHybrid-6.3 \
#-lvtkViewsContext2D-6.3 \
#-lvtkInteractionWidgets-6.3 \
#-lvtkRenderingContext2D-6.3 \
#-lvtkRenderingOpenGL-6.3 \
#-lvtkRenderingAnnotation-6.3 \
#-lvtkInteractionStyle-6.3 \
#-lvtkRenderingVolumeOpenGL-6.3 \
#-lvtkViewsCore-6.3 \
#-lvtkFiltersGeneral-6.3 \
#-lvtkRenderingVolume-6.3 \
#-lvtkImagingCore-6.3 \
#-lvtkRenderingFreeType-6.3 \
#-lvtkftgl-6.3 \
#-lvtkRenderingCore-6.3 \
#-lvtkIOXMLParser-6.3 \
#-lvtkFiltersExtraction-6.3 \
#-lvtkFiltersCore-6.3 \
#-lvtkCommonExecutionModel-6.3 \
#-lvtkCommonDataModel-6.3 \
#-lvtkCommonCore-6.3 \
#-lvtkIOImage-6.3 \
#-lvtkFiltersSMP-6.3 \
#-lvtkIOParallelXML-6.3 \
#-lvtkRenderingContextOpenGL-6.3 \
#-lvtkRenderingLIC-6.3 \
#-lvtktiff-6.3 \
#-lvtkCommonTransforms-6.3 \
#-lvtkChartsCore-6.3 \
#-lvtkCommonColor-6.3 \
#-lvtkCommonComputationalGeometry-6.3 \
#-lvtkCommonMath-6.3 \
#-lvtkCommonMisc-6.3 \
#-lvtkCommonSystem-6.3 \
#-lvtkDICOMParser-6.3 \
#-lvtkDomainsChemistry-6.3 \
#-lvtkFiltersAMR-6.3 \
#-lvtkFiltersFlowPaths-6.3 \
#-lvtkFiltersGeneric-6.3 \
#-lvtkFiltersHybrid-6.3 \
#-lvtkFiltersHyperTree-6.3 \
#-lvtkFiltersImaging-6.3 \
#-lvtkFiltersModeling-6.3 \
#-lvtkFiltersParallel-6.3 \
#-lvtkFiltersParallelImaging-6.3 \
#-lvtkFiltersProgrammable-6.3 \
#-lvtkFiltersSelection-6.3 \
#-lvtkFiltersSources-6.3 \
#-lvtkFiltersStatistics-6.3 \
#-lvtkFiltersTexture-6.3 \
#-lvtkFiltersVerdict-6.3 \
#-lvtkGeovisCore-6.3 \
#-lvtkIOAMR-6.3 \
#-lvtkIOCore-6.3 \
#-lvtkIOEnSight-6.3 \
#-lvtkIOExodus-6.3 \
#-lvtkIOExport-6.3 \
#-lvtkIOGeometry-6.3 \
#-lvtkIOImport-6.3 \
#-lvtkIOInfovis-6.3 \
#-lvtkIOLSDyna-6.3 \
#-lvtkIOLegacy-6.3 \
#-lvtkIOMINC-6.3 \
#-lvtkIOMovie-6.3 \
#-lvtkIONetCDF-6.3 \
#-lvtkIOPLY-6.3 \
#-lvtkIOParallel-6.3 \
#-lvtkIOSQL-6.3 \
#-lvtkIOVideo-6.3 \
#-lvtkIOXML-6.3 \
#-lvtkImagingColor-6.3 \
#-lvtkImagingFourier-6.3 \
#-lvtkImagingGeneral-6.3 \
#-lvtkImagingMath-6.3 \
#-lvtkImagingMorphological-6.3 \
#-lvtkImagingSources-6.3 \
#-lvtkImagingStatistics-6.3 \
#-lvtkImagingStencil-6.3 \
#-lvtkInfovisCore-6.3 \
#-lvtkInteractionImage-6.3 \
#-lvtkNetCDF-6.3 \
#-lvtkNetCDF_cxx-6.3 \
#-lvtkParallelCore-6.3 \
#-lvtkRenderingGL2PS-6.3 \
#-lvtkRenderingImage-6.3 \
#-lvtkRenderingLOD-6.3 \
#-lvtkTestingGenericBridge-6.3 \
#-lvtkTestingIOSQL-6.3 \
#-lvtkTestingRendering-6.3 \
#-lvtkalglib-6.3 \
#-lvtkexoIIc-6.3 \
#-lvtkexpat-6.3 \
#-lvtkfreetype-6.3 \
#-lvtkgl2ps-6.3 \
#-lvtkhdf5-6.3 \
#-lvtkhdf5_hl-6.3 \
#-lvtkjpeg-6.3 \
#-lvtkjsoncpp-6.3 \
#-lvtklibxml2-6.3 \
#-lvtkmetaio-6.3 \
#-lvtkoggtheora-6.3 \
#-lvtkpng-6.3 \
#-lvtkproj4-6.3 \
#-lvtksqlite-6.3 \
#-lvtksys-6.3 \
#-lvtkverdict-6.3 \
#-lvtkzlib-6.3
