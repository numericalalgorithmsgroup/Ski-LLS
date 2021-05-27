# This is a Makefile for Ski-LLS - a package for solving large-scale
# over-determined linear least square problems
# Call "make help" or refer to the documentation.

# Config files defines compilers, flags and the location of existing
# libraries. See config directory.
CONFIG ?= config/make_gcc.inc
include ${CONFIG}

# where to build
BUILDROOT := .
OBJDIR := $(BUILDROOT)/OBJ
LIBDIR := $(BUILDROOT)/LIB
EXEDIR := $(BUILDROOT)/EXE

# library source file
SRC := $(wildcard src/*.cpp) 
CSRC := hqrrp/lapack_compatible_sources/NoFLA_HQRRP_WY_blk_var4_skint.c

# header files
HPP = $(wildcard include/*.hpp)

# tests source files
TSTSRC := $(wildcard test/*.cpp) $(wildcard test/*.c)

# library objects
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o) $(CSRC:.c=.o)))

# library to be compiled
LIB := $(LIBDIR)/libski-lls.a 

# All test executables
EXE := $(addprefix $(EXEDIR)/, $(addsuffix .exe, $(notdir $(basename $(TSTSRC)))))
# Ski-LLS test executables - exclude the one for HQRRP
SKEXE := $(filter-out $(EXEDIR)/hqrrp%, $(EXE))

#$(info SRC = $(SRC))
#$(info CSRC = $(CSRC))
#$(info OBJ = $(OBJ))
#$(info LIB = $(LIB))
#$(info EXE = $(EXE))
#$(info SKEXE = $(SKEXE))

.SUFFIXES: 
.PHONY: all lib test hqrrp_test exe clean help first example $(notdir $(EXE))

VPATH = src : hqrrp/lapack_compatible_sources : test

INCLUDES = -I./include -I${SPARSE_INCLUDE} \
           -I$(BOOSTROOT) -I$(FFTW_INCLUDE)

LIBS = $(LIB) $(LIBS_SPARSE) $(LIBS_FFTW) $(LIBS_LAPACK) -lpthread -lm -ldl

# Show help if no target is specified
first: help

all: $(LIB) exe hqrrp_test test

lib: $(LIB)

$(OBJ) : $(HPP)

$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(INCLUDES) ${CXXFLAGS} -c $< -o $@

$(OBJDIR)/%.o : %.c
	@mkdir -p $(OBJDIR)
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $@

$(LIB):	$(OBJ)
	@mkdir -p $(LIBDIR)
	ar rc $(LIB) $(OBJ)

# all tests except HQRRP
$(SKEXE) : $(EXEDIR)/%.exe : $(OBJDIR)/%.o $(LIB)
	@mkdir -p $(EXEDIR)
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@

$(EXEDIR)/hqrrp_test.exe : $(OBJDIR)/hqrrp_test.o $(OBJDIR)/NoFLA_HQRRP_WY_blk_var4_skint.o 
	@mkdir -p $(EXEDIR)
	$(CC) $(CFLAGS) $^ $(LIBS_LAPACK) -lm -o $@

# Rule to allow executable files to be specified without directory
$(notdir $(EXE)) : %.exe : $(EXEDIR)/%.exe

# compile all executables
exe: $(EXE)

# run tests
hqrrp_test: $(EXEDIR)/hqrrp_test.exe
	$(EXEDIR)/hqrrp_test.exe

test: $(SKEXE)
	test/test_solvers.sh

example: EXE/ski-llsDenseExample.exe EXE/ski-llsSparseExample.exe
	@echo Running dense Ski-LLS example...
	EXE/ski-llsDenseExample.exe
	@echo =================================
	@echo Running sparse Ski-LLS example...
	EXE/ski-llsSparseExample.exe
	@echo =================================

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(EXEDIR)

help:
	@echo "==============================================================================="
	@echo "Ski-LLS - a package for solving large-scale over-determined"
	@echo "          linear least square problems"
	@echo "Please adapt your config makefile pointing to the location"
	@echo "of the dependencies. If you are not using the default one"
	@echo "call:  make CONFIG=config/my_config_file.inc lib [or other target]"
	@echo "==============================================================================="
	@echo "make help         : prints this message"
	@echo "make lib          : build the Ski-LLS library"
	@echo "                    $(LIB)"
	@echo "make exe          : build all the test executables"
	@echo "make hqrrp_test   : self check of HQRRP"
	@echo "make example      : build & run tiny dense and sparse Ski-LLS examples"
	@echo "make test         : run all tests of Ski-LLS on sample data"
	@echo "make all          : make the library & run the tests"
	@echo "make clean        : remove all build directories"
	@echo "Or call make on a single target"
	@echo "  $(notdir $(EXE))"
	@echo "==============================================================================="


