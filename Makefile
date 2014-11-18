# The compiler

FC = gfortran

# flags for debugging or for maximum performance, comment as necessary

#FCFLAGS = -g -fbounds-check 
FCFLAGS = -Wall

# flags forall (e.g. look for system .mod files, required in gfortran)

FCFLAGS += -Jmodules/mod -Imodules/mod

# List of executables to be built within the package

MODULES_SRC_DIR   = modules/src
MODULES_OBJ_DIR   = modules/obj

MODULES_F90  := $(wildcard $(MODULES_SRC_DIR)/*.f90)
MODULES_O    := $(patsubst $(MODULES_SRC_DIR)%.f90,$(MODULES_OBJ_DIR)%.o,$(MODULES_F90))

EXAMPLES_SRC_DIR   = examples/src
EXAMPLES_OBJ_DIR   = examples/obj
EXAMPLES_BIN_DIR   = examples/bin

EXAMPLES_F90  := $(wildcard $(EXAMPLES_SRC_DIR)/*.f90)
EXAMPLES_O    := $(patsubst $(EXAMPLES_SRC_DIR)%.f90,$(EXAMPLES_OBJ_DIR)%.o,$(EXAMPLES_F90))
EXAMPLES_EXE  := $(patsubst $(EXAMPLES_SRC_DIR)%.f90,$(EXAMPLES_BIN_DIR)%.exe,$(EXAMPLES_F90))

MY_PROJECT_SRC_DIR   = my_project/src
MY_PROJECT_OBJ_DIR   = my_project/obj
MY_PROJECT_BIN_DIR   = my_project/bin

MY_PROJECT_F90  := $(wildcard $(MY_PROJECT_SRC_DIR)/*.f90)
MY_PROJECT_O    := $(patsubst $(MY_PROJECT_SRC_DIR)%.f90,$(MY_PROJECT_OBJ_DIR)%.o,$(MY_PROJECT_F90))
MY_PROJECT_EXE  := $(patsubst $(MY_PROJECT_SRC_DIR)%.f90,$(MY_PROJECT_BIN_DIR)%.exe,$(MY_PROJECT_F90))
# "make" builds all
#TARGET = $(patsubst %.f90, %.exe, $(EXAMPLES))

OBJECTS  := $(EXAMPLES_F90:$(EXAMPLES_SRC_DIR)/%.f90=$(EXAMPLES_OBJ_DIR)/%.o)

# general rules

all: modules examples my_project
modules: $(MODULES_O)
examples:$(EXAMPLES_EXE) $(EXAMPLES_O)
my_project:$(MY_PROJECT_EXE) $(MY_PROJECT_O)

# rules for modules

$(MODULES_OBJ_DIR)/module_invasion_percolation.o: $(MODULES_OBJ_DIR)/module_disjoint_set.o $(MODULES_OBJ_DIR)/module_binary_tree.o $(MODULES_OBJ_DIR)/module_write_output_files.o $(MODULES_OBJ_DIR)/module_invasion_percolation_constants.o $(MODULES_OBJ_DIR)/module_trapping.o $(MODULES_OBJ_DIR)/module_random_media.o
$(MODULES_OBJ_DIR)/module_trapping.o: $(MODULES_OBJ_DIR)/module_invasion_percolation_constants.o $(MODULES_OBJ_DIR)/module_disjoint_set.o $(MODULES_OBJ_DIR)/module_label_clusters.o
$(MODULES_OBJ_DIR)/module_write_output_files.o: $(MODULES_OBJ_DIR)/module_trapping.o $(MODULES_OBJ_DIR)/module_cubic_indices.o $(MODULES_OBJ_DIR)/module_invasion_percolation_constants.o 
$(MODULES_O): $(MODULES_OBJ_DIR)/%.o:$(MODULES_SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# rules for examples

$(EXAMPLES_EXE): $(EXAMPLES_BIN_DIR)/%.exe : $(EXAMPLES_SRC_DIR)/%.f90 
	$(FC) $(FCFLAGS) -o $@ $^  $(MODULES_O)
	@echo "Linking complete!"

$(EXAMPLES_O): $(EXAMPLES_OBJ_DIR)/%.o : $(EXAMPLES_SRC_DIR)/%.f90 
	$(FC) $(FCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# rules for my_project

$(MY_PROJECT_EXE): $(MY_PROJECT_BIN_DIR)/%.exe : $(MY_PROJECT_SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -o $@ $^
	@echo "Linking complete!"

$(MY_PROJECT_O): $(MY_PROJECT_OBJ_DIR)/%.o : $(MY_PROJECT_SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean 

clean:
	rm -f modules/obj/*
	rm -f modules/mod/*
	rm -f examples/obj/*
	rm -f examples/bin/*
	rm -f my_project/obj/*
	rm -f my_project/bin/*