#-------------------------------------------------
#
# Project created by QtCreator 2010-07-11T13:11:32
#
#-------------------------------------------------

CONFIG += qt
QT += core gui widgets
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
    
#VTKBINPATH = /Users/GZhou/Desktop/Code-newversion/bin6.0.0
#VTKVTKPATH = /Users/GZhou/Desktop/Code-newversion/VTK6.0.0
#CLAPACKPATH = /Users/GZhou/Desktop/TreeScaper/clapack-3.2.1-CMAKE

VTKBINPATH = /Users/whuang/TreeScaper_dev/bin6.0.0
VTKVTKPATH = /Users/whuang/TreeScaper_dev/VTK6.0.0
CLAPACKPATH = /Users/whuang/TreeScaper_dev/clapack-3.2.1-CMAKE

INCPATH += $${CLAPACKPATH}/INCLUDE

INCPATH += $${VTKVTKPATH}/Common/Math
INCPATH += $${VTKBINPATH}/Common/Math
INCPATH += $${VTKVTKPATH}/Testing/Rendering
INCPATH += $${VTKBINPATH}/Testing/Rendering
INCPATH += $${VTKBINPATH}/Rendering/Volume
INCPATH += $${VTKBINPATH}/Filters/Hybrid
INCPATH += $${VTKBINPATH}/Interaction/Widgets
INCPATH += $${VTKVTKPATH}/Views/Context2D
INCPATH += $${VTKBINPATH}/Views/Context2D
INCPATH += $${VTKVTKPATH}/Views/Core
INCPATH += $${VTKBINPATH}/Views/Core
INCPATH += $${VTKVTKPATH}/Views/Infovis
INCPATH += $${VTKBINPATH}/Views/Infovis
INCPATH += $${VTKVTKPATH}/Infovis/Core
INCPATH += $${VTKVTKPATH}/Infovis/Layout
INCPATH += $${VTKVTKPATH}/Infovis/BoostGraphAlgorithms
INCPATH += $${VTKBINPATH}/Infovis/Core
INCPATH += $${VTKBINPATH}/Infovis/Layout
INCPATH += $${VTKBINPATH}/Rendering/OpenGL
INCPATH += $${VTKBINPATH}/Rendering/Context2D
INCPATH += $${VTKVTKPATH}/Rendering/Context2D
INCPATH += $${VTKVTKPATH}/Common/Transforms
INCPATH += $${VTKBINPATH}/Common/Transforms
INCPATH += $${VTKVTKPATH}/Filters
INCPATH += $${VTKBINPATH}
INCPATH += $${VTKBINPATH}/Common
INCPATH += $${VTKBINPATH}/Common/Core
INCPATH += $${VTKBINPATH}/Common/DataModel
INCPATH += $${VTKVTKPATH}/Common/DataModel
INCPATH += $${VTKVTKPATH}/Common/ExecutionModel
INCPATH += $${VTKBINPATH}/Common/ExecutionModel
INCPATH += $${VTKBINPATH}/Utilities
INCPATH += $${VTKBINPATH}/Rendering
INCPATH += $${VTKBINPATH}/Rendering/VolumeOpenGL
INCPATH += $${VTKBINPATH}/Charts/Core
INCPATH += $${VTKBINPATH}/ThirdParty/alglib
INCPATH += $${VTKVTKPATH}/IO
INCPATH += $${VTKVTKPATH}/Rendering/Annotation
INCPATH += $${VTKVTKPATH}/Geovis/Core
INCPATH += $${VTKVTKPATH}/Views
INCPATH += $${VTKVTKPATH}/IO/XML
INCPATH += $${VTKBINPATH}/IO/XML
INCPATH += $${VTKBINPATH}/IO/Geometry
INCPATH += $${VTKVTKPATH}/Filters/Geometry
INCPATH += $${VTKBINPATH}/Common/ExecutionModel
INCPATH += $${VTKBINPATH}/Filters/Geometry
INCPATH += $${VTKBINPATH}/Rendering/Core
INCPATH += $${VTKBINPATH}/Filters/Extraction
INCPATH += $${VTKBINPATH}/Filters/Statistics
INCPATH += $${VTKBINPATH}/IO/Image
INCPATH += $${VTKVTKPATH}/Filters/Sources
INCPATH += $${VTKBINPATH}/Filters/Sources
INCPATH += $${VTKBINPATH}/Filters/Imaging
INCPATH += $${VTKVTKPATH}/Imaging/Hybrid
INCPATH += $${VTKBINPATH}/Imaging/Hybrid
INCPATH += $${VTKVTKPATH}/Common/Misc
INCPATH += $${VTKBINPATH}/Common/Misc
INCPATH += $${VTKVTKPATH}/IO/Movie
INCPATH += $${VTKBINPATH}/IO/Movie
INCPATH += $${VTKVTKPATH}/IO/FFMPEG
INCPATH += $${VTKVTKPATH}/Imaging/Core
INCPATH += $${VTKBINPATH}/Imaging/Core
INCPATH += $${VTKVTKPATH}/Imaging/Sources
INCPATH += $${VTKBINPATH}/Imaging/Sources
INCPATH += $${VTKBINPATH}/Rendering/Annotation
INCPATH += $${VTKBINPATH}/IO/Export
INCPATH += $${VTKBINPATH}/Rendering/FreeType
INCPATH += $${VTKBINPATH}/Rendering/Label
INCPATH += $${VTKVTKPATH}/Interaction/Widgets
INCPATH += $${VTKVTKPATH}/Rendering/Core
INCPATH += $${VTKVTKPATH}/Charts/Core
INCPATH += $${VTKVTKPATH}/Rendering/OpenGL/Testing/Cxx
INCPATH += $${VTKVTKPATH}/IO/Core
INCPATH += $${VTKVTKPATH}/IO/Image
INCPATH += $${VTKVTKPATH}/Filters/General
INCPATH += $${VTKVTKPATH}/Filters/Modeling
INCPATH += $${VTKVTKPATH}/Filters/Core
INCPATH += $${VTKBINPATH}/Filters/Core
INCPATH += $${VTKVTKPATH}/Filters/Generic
INCPATH += $${VTKVTKPATH}/Common/Core
INCPATH += $${VTKVTKPATH}/Utilities
INCPATH += $${VTKVTKPATH}/Utilities/KWSys
INCPATH += $${VTKVTKPATH}/Common/Core/Testing/Cxx
INCPATH += $${VTKBINPATH}/ThirdParty/libproj4
INCPATH += $${VTKVTKPATH}/ThirdParty/libproj4/vtklibproj4
INCPATH += $${VTKBINPATH}/Utilities/DICOMParser
INCPATH += $${VTKVTKPATH}/Utilities/DICOMParser
INCPATH += $${VTKBINPATH}/ThirdParty/freetype/vtkfreetype/include
INCPATH += $${VTKVTKPATH}/ThirdParty/freetype/vtkfreetype/include
INCPATH += $${VTKBINPATH}/ThirdParty/netcdf/vtknetcdf
INCPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/libsrc
INCPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/libsrc4
INCPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/include
INCPATH += $${VTKVTKPATH}/ThirdParty/netcdf/vtknetcdf/libdispatch
INCPATH += $${VTKBINPATH}/ThirdParty/exodusII/vtkexodusII/include
INCPATH += $${VTKVTKPATH}/ThirdParty/exodusII/vtkexodusII/include
INCPATH += $${VTKBINPATH}/Utilities/MaterialLibrary
INCPATH += $${VTKVTKPATH}/Utilities/MaterialLibrary
INCPATH += $${VTKVTKPATH}/ThirdParty/verdict/vtkverdict
INCPATH += $${VTKVTKPATH}/ThirdParty/verdict/vtkverdict
INCPATH += $${VTKBINPATH}/ThirdParty/alglib
INCPATH += $${VTKVTKPATH}/ThirdParty/alglib
INCPATH += $${VTKVTKPATH}/Interaction/Style
INCPATH += $${VTKBINPATH}/Interaction/Style
INCPATH += $${VTKBINPATH}/ThirdParty/ftgl
INCPATH += $${VTKVTKPATH}/ThirdParty/ftgl/src
INCPATH += $${VTKVTKPATH}/Filters/Extraction
INCPATH += $${VTKBINPATH}/Filters/General
#INCPATH += /Users/whuang/TreeScaper_dev/vtkmpeg2encode

INCPATH += /Developer/SDKs/MacOSX10.8.sdk/usr/X11/include/GL
INCPATH += /Developer/SDKs/MacOSX10.8.sdk/usr/X11/include/GL/internal

#INCPATH += /usr/include/X11
INCPATH += /Developer/SDKs/MacOSX10.8.sdk/usr/X11/include/X11

#INCPATH += /usr/local/Cellar/boost/1.58.0/include

LIBS += $${CLAPACKPATH}/lapack_MAC.a
LIBS += $${CLAPACKPATH}/blas_MAC.a
LIBS += $${CLAPACKPATH}/F2CLIBS/libf2c.a
LIBS += -lm

LIBS += $${VTKBINPATH}/lib/libvtkChartsCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonColor-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonComputationalGeometry-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonDataModel-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonExecutionModel-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonMath-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonMisc-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonSystem-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkCommonTransforms-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkDICOMParser-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkDomainsChemistry-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersAMR-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersExtraction-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersFlowPaths-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersGeneral-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersGeneric-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersGeometry-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersHybrid-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersHyperTree-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersImaging-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersModeling-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersParallel-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersParallelImaging-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersProgrammable-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersSelection-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersSources-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersStatistics-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersTexture-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkFiltersVerdict-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkGeovisCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOAMR-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOEnSight-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOExodus-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOExport-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOGeometry-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOImage-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOImport-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOInfovis-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOLSDyna-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOLegacy-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOMINC-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOMovie-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIONetCDF-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOPLY-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOParallel-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOSQL-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOVideo-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOXML-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkIOXMLParser-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingColor-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingFourier-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingGeneral-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingHybrid-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingMath-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingMorphological-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingSources-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingStatistics-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkImagingStencil-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkInfovisCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkInfovisLayout-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkInteractionImage-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkInteractionStyle-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkInteractionWidgets-6.0.a
#LIBS += $${VTKBINPATH}/lib/libvtkLocalExample-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkNetCDF-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkNetCDF_cxx-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkParallelCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingAnnotation-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingContext2D-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingFreeType-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingFreeTypeOpenGL-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingGL2PS-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingHybridOpenGL-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingImage-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingLOD-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingLabel-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingOpenGL-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingVolume-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingVolumeAMR-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkRenderingVolumeOpenGL-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkTestingGenericBridge-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkTestingIOSQL-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkTestingRendering-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsContext2D-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsCore-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsGeovis-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkViewsInfovis-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkalglib-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkexoIIc-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkexpat-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkfreetype-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkftgl-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkgl2ps-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkhdf5-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkhdf5_hl-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkjpeg-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkjsoncpp-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtklibxml2-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkmetaio-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkoggtheora-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkpng-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkproj4-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtksqlite-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtksys-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtktiff-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkverdict-6.0.a
LIBS += $${VTKBINPATH}/lib/libvtkzlib-6.0.a


#LIBS += /usr/lib/libGL.so
#LIBS += /usr/lib/libglitz.a
#LIBS += /usr/lib/libglitz-glx.a
#LIBS += /usr/lib/caca/libgl_plugin.a
##LIBS += /Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib/libGL.dylib
##-fobjc-gc

#LIBS += /Library/Frameworks/QtGui.framework/Versions/4/QtGui_debug

#LIBS += /usr/lib/libXt.a
#LIBS += /usr/lib/libX11.a
#LIBS += /usr/lib/caca/libx11_plugin.a
#LIBS += /usr/lib/libXt.so
#LIBS += /Developer/SDKs/MacOSX10.4u.sdk/usr/X11R6/lib/libXt.a
#LIBS += /Developer/SDKs/MacOSX10.4u.sdk/usr/X11R6/lib/libX11.a

#LIBS += /Users/wenhuang/software/qt-everywhere-opensource-src-4.6.3/lib/QtGui.framework/Versions/4/QtGui
#LIBS += /Users/wenhuang/software/qt-everywhere-opensource-src-4.6.3/lib/QtCore.framework/Versions/4/QtCore

#LIBS += -framework Carbon
LIBS += -framework Cocoa
LIBS += -framework IOKit
#LIBS += -no-framework
LIBS += -framework OpenGL
LIBS += -framework AGL
LIBS += -framework QtGui
LIBS += -framework QtCore
