PROGRAM = run_hier
#PROGRAM = run_crosstalk



HEADERS = ./src
SRC_DIR := ./src
OBJ_DIR := ./obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

FLAGS =   -std=c++11 -O3 -MMD
CXX := g++

EXECUTABLE=./$(PROGRAM) 


all: $(EXECUTABLE) #run


$(PROGRAM): $(PROGRAM).o $(OBJ_FILES)
	$(CXX) $(FLAGS) -o  $(PROGRAM) $^ -I$(HEADERS) 

$(PROGRAM).o: $(PROGRAM).cpp 
	$(CXX) -c $(FLAGS) $(PROGRAM).cpp -I$(HEADERS)


# Compile source code
# see https://stackoverflow.com/questions/2908057/can-i-compile-all-cpp-files-in-src-to-os-in-obj-then-link-to-binary-in
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(FLAGS) -c -o $@ $<


run: $(EXECUTABLE)
	$(EXECUTABLE)

clean:
	rm -f *.o  $(OBJ_DIR)/*.o *.d $(OBJ_DIR)/*.d

# see https://stackoverflow.com/questions/2908057/can-i-compile-all-cpp-files-in-src-to-os-in-obj-then-link-to-binary-in
# see http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
#CXXFLAGS += -MMD
-include $(OBJ_FILES:.o=.d)

