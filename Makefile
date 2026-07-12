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
# Cavity solver
# ==========================================

CAVITY_SOURCES = \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/Lattice.cpp \
	$(SRC_DIR)/Boundary.cpp \
	$(SRC_DIR)/Cavity.cpp


CAVITY_OBJECTS = \
	$(CAVITY_SOURCES:%.cpp=$(BUILD_DIR)/%.o)


CAVITY_TARGET = cavity



# ==========================================
# Postprocess
# ==========================================

POSTPROCESS_SOURCES = \
	$(SRC_DIR)/postprocess.cpp


POSTPROCESS_OBJECTS = \
	$(POSTPROCESS_SOURCES:%.cpp=$(BUILD_DIR)/%.o)


POSTPROCESS_TARGET = postprocess



# ==========================================
# Default: build both executables
# ==========================================

all: $(CAVITY_TARGET) $(POSTPROCESS_TARGET)



# ==========================================
# Link cavity
# ==========================================

$(CAVITY_TARGET): $(CAVITY_OBJECTS)

	$(CXX) $(CXXFLAGS) $(CAVITY_OBJECTS) -o $@



# ==========================================
# Link postprocess
# ==========================================

$(POSTPROCESS_TARGET): $(POSTPROCESS_OBJECTS)

	$(CXX) $(CXXFLAGS) $(POSTPROCESS_OBJECTS) -o $@



# ==========================================
# Compile rule
# ==========================================

$(BUILD_DIR)/%.o: %.cpp

	mkdir -p $(dir $@)

	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@



# ==========================================
# Run solver
# ==========================================

run: $(CAVITY_TARGET)

	./$(CAVITY_TARGET)



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

	rm -f $(CAVITY_TARGET)

	rm -f $(POSTPROCESS_TARGET)



# ==========================================
# Rebuild
# ==========================================

rebuild: clean all
