# Apple toolchain (system libc++)
CXX=/usr/bin/clang++
CXXFLAGS=-std=c++17 -Wall -Wextra -pedantic -O3 -march=native \
         -Xpreprocessor -fopenmp \
         -I/opt/homebrew/opt/libomp/include \
         -Isrc
LDFLAGS=-L/opt/homebrew/opt/libomp/lib -lomp \
        -Wl,-rpath,/opt/homebrew/opt/libomp/lib

SRC := $(wildcard src/*.cpp src/*/*.cpp src/*/*/*.cpp)
OBJ := $(SRC:.cpp=.o)
TARGET := mc

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
