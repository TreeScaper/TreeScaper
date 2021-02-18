CC = g++
CFLAG = -c -std=c++17 -Wall
LFLAG = -std=c++17 -Wall

OBJ = main.o wstring.o Sparse_matrix.o zdarray.o zdtree.o

all: generate_version CLVTreeScaper2

.PHONY: generate_version
generate_version:
	@./generate_version.sh

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
