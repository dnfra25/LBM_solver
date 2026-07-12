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
# Solver sources
# ==========================================

SOLVER_SOURCES = \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/Lattice.cpp \
	$(SRC_DIR)/Boundary.cpp \
	$(SRC_DIR)/Cavity.cpp



SOLVER_OBJECTS = $(SOLVER_SOURCES:%.cpp=$(BUILD_DIR)/%.o)



TARGET = cavity



# ==========================================
# Postprocess sources
# ==========================================

POSTPROCESS_SOURCE = \
	$(SRC_DIR)/postprocess.cpp



POSTPROCESS_OBJECT = \
	$(BUILD_DIR)/src/postprocess.o



POSTPROCESS_TARGET = postprocess



# ==========================================
# Default
# ==========================================

all: $(TARGET)



# ==========================================
# Solver linking
# ==========================================

$(TARGET): $(SOLVER_OBJECTS)

	$(CXX) $(CXXFLAGS) $(SOLVER_OBJECTS) -o $@



# ==========================================
# Postprocess linking
# ==========================================

$(POSTPROCESS_TARGET): $(POSTPROCESS_OBJECT)

	$(CXX) $(CXXFLAGS) $(POSTPROCESS_OBJECT) -o $@



# ==========================================
# Compilation rule
# ==========================================

$(BUILD_DIR)/%.o: %.cpp

	mkdir -p $(dir $@)

	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@



# ==========================================
# Run solver
# ==========================================

run: $(TARGET)

	./$(TARGET)



# ==========================================
# Run postprocess
# ==========================================

post: $(POSTPROCESS_TARGET)

	./$(POSTPROCESS_TARGET)



# ==========================================
# Clean
# ==========================================

clean:

	rm -rf $(BUILD_DIR)

	rm -f $(TARGET)

	rm -f $(POSTPROCESS_TARGET)



# ==========================================
# Rebuild
# ==========================================

rebuild: clean all
