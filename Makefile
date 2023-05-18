CC = g++

include ./makeCLVTreeScaper.inc

$(info CLAPACK path is "$(CLAPPATH)")

INCDIRS = -I $(CLAPPATH)/INCLUDE
LIBDIRS = -L $(CLAPPATH)/LIB

CLAPLIB = $(CLAPPATH)/lapack$(PLAT).a
BLASLIB = $(CLAPPATH)/blas$(PLAT).a
F2CLIB = $(CLAPPATH)/F2CLIBS/libf2c.a

LDLIBS = $(CLAPLIB) $(BLASLIB) $(F2CLIB) -lm


#CFLAG = -c $(INCDIRS) $(LIBDIRS) -lblas -lf2c -llapack
# CPPFLAG = -w -std=c++17 -fpermissive $(INCDIRS) -c
CPPFLAG = -w -g -std=c++17 -Wno-register -fpermissive $(INCDIRS) -c
#LFLAG = -std=c++17 -DNDEBUG $(INCDIRS) $(LIBDIRS) -lblas -lf2c -llapack
LFLAG = -w -std=c++17 -fpermissive -DNDEBUG $(INCDIRS) $(LDLIBS)


SRC = $(wildcard *.cpp)

OBJ := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

target = CLVTreeScaper2

all : $(target)

%.o:%.cpp
	$(CC) $(CPPFLAG) $< -o $@   

%.d:%.cpp
	@set -e; rm -f $@; $(CC) $(INCDIRS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

sinclude $(SRC:.cpp=.d)

$(target) : $(OBJ)
	$(CC) $^ -o $@  $(LFLAG)


CLAPACK:
	make f2clib -C $(CLAPPATH)
	make blaslib -C $(CLAPPATH)
	make -C $(CLAPPATH)/INSTALL
	make -C $(CLAPPATH)/SRC


.PHONY : clean
clean :
	rm -rf $(target) *.o *.d *.d.*