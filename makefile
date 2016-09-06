# Declaration of variables

# File names
SRC_DIR =./tls
OBJ_DIR =./_obj
EXEC = ./_out/tls
EXEC_D = ./_out_d/tls
BOOST_DIR = /usr/local/include
#/hpc/software/boost/1.57.0/include
EIGEN_DIR = ../../Eigen
#/usr/include/Eigen
#/home/dza845/c++/include/Eigen
QI_DIR = ../../QIlib/QIlib
QI_LIB_DIR = ../../QIlib/lib
#COMPILER=gcc
COMPILER=g++-5

SOURCES = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/CMAES/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
LIBS = -L$(QI_LIB_DIR) -lrt -lQIlib


CC_FLAGS = -O3 -Ofast

setcc: 
	$(eval CC = $(COMPILER) -mavx -std=c++11 -fabi-version=10 -fopenmp -pthread -static-libgcc -Wno-deprecated-declarations -I $(QI_DIR) -I $(EIGEN_DIR) -I $(BOOST_DIR) -I /usr/include)

setcc_d: setcc
	$(eval CC += -g -Wall -Wno-unknown-pragmas -Wno-unused-local-typedefs) 

#OBJECTS = $(SOURCES:.cpp=.o)

# debug target -Werror
build_d: dirs setcc_d  $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC) $(LIBS)

# Main target
build: dirs setcc  $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC) $(LIBS)

# To obtain object files
$(OBJECTS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

dirs:
	mkdir -p $(OBJ_DIR)
	mkdir -p ./_out
	mkdir -p ./_out_d

remake_all:
	make clean
	make build
