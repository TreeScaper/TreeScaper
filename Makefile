CC = g++

include makeCLVTreeScaper.inc

INCDIRS = -I$(CLAPPATH)/INCLUDE
CLAPLIB = $(CLAPPATH)/lapack$(PLAT).a
BLASLIB = $(CLAPPATH)/blas$(PLAT).a
F2CLIB = $(CLAPPATH)/F2CLIBS/libf2c.a

#common command
#LDLIBS  = $(CLAPLIB) $(BLASLIB) $(F2CLIB) -lm #$(BLWRLIB) -lm
LDLIBS  = $(CLAPLIB) $(BLASLIB) $(F2CLIB) -lm

all: CLVTreeScaper

.PHONY: generate_version
generate_version:
	@./generate_version.sh

CLVTreeScaper: generate_version
				$(CC) CLVmain.cpp -w randgen.cpp wstring.cpp warray.cpp wmapping.cpp wmix.cpp wfile.cpp wimport_form.cpp wDimEst.cpp wNLDR.cpp zdcommunity.cpp Trees.cpp TreeOPE.cpp Sparse_matrix.cpp greedy_louvain.cpp graph.cpp slicer.cpp label-map.cc community.cpp info.cpp hashfunc.cc hash.cc hungarian.c ClusterForest.cpp ClusterInstance.cpp Forest.cpp rspr.cpp SiblingPair.cpp SPRLCA.cpp SPRNode.cpp UndoMachine.cpp $(LDLIBS) $(INCDIRS) -DCOMMAND_LINE_VERSION -o CLVTreeScaper -DCOMMAND_LINE_VERSION

CLAPACK:
				make f2clib -C $(CLAPPATH)
				make blaslib -C $(CLAPPATH)
				make -C $(CLAPPATH)/INSTALL
				make -C $(CLAPPATH)/SRC
				
CLVclean:
				rm -f CLVTreeScaper
