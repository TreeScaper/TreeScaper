CC = g++
CFLAG = -c -std=c++17 -Wall
LFLAG = -std=c++17 -Wall
LDFLAGS = -lcurl -lpugixml

OBJ = main.o array.o wstring.o sparse.o zdtree.o cra.o

all: CLVTreeScaper2

.PHONY: generate_version
generate_version:
	@./generate_version.sh

# There must be a "recipe" defined for this rule, or make may not detect
# that version.hpp was updated.
version.hpp: generate_version
	@echo Target generate_version run to generate version number.

CLVTreeScaper2 : $(OBJ)
	$(CC) $(LFLAG) $^ $(LDFLAGS) -o $@

main.o : array.hpp zdtree.hpp wstring.hpp zdtreeobj.hpp version.hpp
wstring.o : 
cra.cpp : cra.hpp
array.o : 
sparse.o : array.hpp
zdtree.o : array.hpp sparse.hpp
zdtreeobj.o : array.hpp sparse.hpp zdtree.hpp

.PHONY : clean
clean :
	rm CLVTreeScaper2 $(OBJ)
