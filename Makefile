CC = g++

#WH HPC path
#ROOTPATH = /panfs/storage.local/genacc/home/wh08/CLAPACK-3.2.1
#INCDIRS = -I$(ROOTPATH)/INCLUDE
#CLAPLIB = $(ROOTPATH)/lapack_LINUX.a
#BLASLIB = $(ROOTPATH)/blas_LINUX.a
#F2CLIB = $(ROOTPATH)/F2CLIBS/libf2c.a
##BLWRLIB = $(ROOTPATH)/libcblaswr.a

#WH mac path
ROOTPATH = /home/whuang/CLAPACK-3.2.1
INCDIRS = -ID:\VTK\clapack-3.2.1-CMAKE\INCLUDE
CLAPLIB = D:\VTK\clapack\SRC\liblapack.a
BLASLIB = D:\VTK\clapack\BLAS\SRC\libblas.a
F2CLIB = D:\VTK\clapack\F2CLIBS\libf2c\libf2c.a

#common command
#LDLIBS  = $(CLAPLIB) $(BLASLIB) $(F2CLIB) -lm #$(BLWRLIB) -lm
LDLIBS  = $(CLAPLIB) $(BLASLIB) $(F2CLIB) -lm

CLVTreeScaper:
				$(CC) CLVmain.cpp -w randgen.cpp wstring.cpp warray.cpp wmapping.cpp wmix.cpp wfile.cpp wimport_form.cpp wDimEst.cpp wNLDR.cpp Trees.cpp TreeOPE.cpp Sparse_matrix.cpp greedy_louvain.cpp graph.cpp slicer.cpp label-map.cc community.cpp info.cpp hashfunc.cc hash.cc hungarian.c ClusterForest.cpp ClusterInstance.cpp Forest.cpp rspr.cpp SiblingPair.cpp SPRLCA.cpp SPRNode.cpp UndoMachine.cpp $(LDLIBS) $(INCDIRS) -DCOMMAND_LINE_VERSION -o CLVTreeScaper
				
CLVclean:
				rm -f CLVTreeScaper
