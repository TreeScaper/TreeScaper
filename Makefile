CC = g++
CFLAG = -c -std=c++17 -Wall
LFLAG = -std=c++17 -Wall

OBJ = main.o wstring.o Sparse_matrix.o zdarray.o zdtree.o

all: CLVTreeScaper2

.PHONY: generate_version
generate_version:
	@./generate_version.sh

# There must be a "recipe" defined for this rule, or make may not detect
# that version.hpp was updated.
version.hpp: generate_version
	@echo Target generate_version run to generate version number.

CLVTreeScaper2 : $(OBJ)
	$(CC) $(LFLAG) $^ -o $@

main.o : zdarray.hpp zdtree.hpp wstring.hpp version.hpp
wstring.o : 
Sparse_matrix.o :
zdarray.o : 
zdtree.o : zdarray.hpp Sparse_matrix.hpp

.PHONY : clean
clean :
	rm CLVTreeScaper2 $(OBJ)
