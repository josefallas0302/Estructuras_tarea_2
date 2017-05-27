####################################################
## Filename: Makefile
## Author: Esteban Zamora Alvarado
## PRIS-LAB
####################################################

#Compile/execute programs
#**************************
CC=gcc
CCX=g++

#Important directories
#**************************
TEST_DIR := tests
BIN_DIR := bin
OBJ_DIR := obj
DOC_DIR := doc
SRC_DIRS := src src/K_Means src/GMM_EM src/foo_dir
INC_DIRS := include src
REQ_DIRS := reqs    #Dependencies
LIB_DIRS := 

#Flags
#**************************
CPPFLAGS := -g -Wall -Wno-deprecated-declarations -std="c++11"
LDFLAGS :=

## include headers and libraries from project
INC_PARAMS := $(foreach d, $(INC_DIRS), -I$d)
LIB_PARAMS := $(foreach d, $(LIB_DIRS), -L$d)

#Rcpp/RInside Flags
#**************************
RCPP_RINSIDE = 0   #By default no R plotting test is generated

ifneq ($(RCPP_RINSIDE), 1)
INC_PARAMS += $(foreach d, $(REQ_DIRS), -I$d)
endif

ifeq ($(RCPP_RINSIDE), 1)

R_HOME := $(shell R RHOME)
CPPFLAGS += $(shell $(R_HOME)/bin/R CMD config CPPFLAGS)

## include headers and libraries for R 
RCPPFLAGS := $(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := $(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := $(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := $(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
#RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## include headers and libraries for RInside embedding classes
RINSIDEINCL := 	$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 	$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## RcppEigen headers
RCPPEIGENINCL := $(shell echo 'cat(paste("-I", system.file("include", package = "RcppEigen"), sep = ""))' | $(R_HOME)/bin/R --vanilla --slave)	

## add Rcpp/RInside headers and libraries to project
INC_PARAMS += $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(RCPPEIGENINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)
LIB_PARAMS += $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

endif

#Other Variables
#**************************
#.cpp files corresponding to executables (tests)
TESTS := $(wildcard $(TEST_DIR)/*.cpp)
PLOT_TESTS := $(wildcard $(TEST_DIR)/plot_tests/*.cpp)

#Executable files
BINS := $(addprefix $(BIN_DIR)/,$(notdir $(TESTS:.cpp=)))

ifeq ($(RCPP_RINSIDE), 1)
BINS += $(addprefix $(BIN_DIR)/,$(notdir $(PLOT_TESTS:.cpp=)))
endif

#.cpp source files in SRC_DIRS
SRCS := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cpp))

#Object files besides from TESTS are only .cpp files in SRC_DIRS
OBJS := $(addprefix $(OBJ_DIR)/,$(notdir $(SRCS:.cpp=.o))) 

#Append .hpp source files in SRC_DIRS to SRCS
SRCS := $(SRCS) $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.hpp))

#Directories that Make use to search for prerrequisites
VPATH := $(TEST_DIR) $(TEST_DIR)/plot_tests $(SRC_DIRS)

RM := rm -rf

#Rules
#**************************
all: $(BINS)

$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(OBJS) $(SRCS)
	$(CCX) $(LDFLAGS) $(LIB_PARAMS) -o $@ $< $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CCX) $(CPPFLAGS) $(INC_PARAMS) -c -o $@ $<

doc: $(DOC_DIR)/Doxyfile
	doxygen $(DOC_DIR)/Doxyfile

clean:
	$(RM) $(BIN_DIR)/* $(OBJ_DIR)/* $(DOC_DIR)/html $(DOC_DIR)/latex

.PHONY: doc clean

#Force Make to don't delete object files
.PRECIOUS: $(OBJ_DIR)/%.o
