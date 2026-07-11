# ==========================================
# Compiler
# ==========================================

CXX = g++

CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -fopenmp

INCLUDES = -Iinclude



# ==========================================
# Directories
# ==========================================

SRC_DIR = src
BUILD_DIR = build



# ==========================================
# Source files
# ==========================================

SOURCES = \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/Lattice.cpp \
	$(SRC_DIR)/Boundary.cpp \
	$(SRC_DIR)/Cavity.cpp



OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)



TARGET = cavity



# ==========================================
# Default target
# ==========================================

all: $(TARGET)



# ==========================================
# Linking
# ==========================================

$(TARGET): $(OBJECTS)

	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@



# ==========================================
# Compilation rule
# ==========================================

$(BUILD_DIR)/%.o: %.cpp

	mkdir -p $(dir $@)

	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@



# ==========================================
# Run
# ==========================================

run: $(TARGET)

	./$(TARGET)



# ==========================================
# Clean
# ==========================================

clean:

	rm -rf $(BUILD_DIR)
	rm -f $(TARGET)



# ==========================================
# Rebuild
# ==========================================

rebuild: clean all
