###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=PairwiseAlignment/WFA2-lib
FOLDER_LIB=PairwiseAlignment/WFA2-lib/lib

###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BUILD=PairwiseAlignment/WFA2-lib/build
FOLDER_BUILD_CPP=PairwiseAlignment/WFA2-lib/build/cpp

UNAME=$(shell uname)

CC:=$(CC)
CPP:=$(CXX)

CC_FLAGS=-w -Wall -g -fPIE

AR=ar
AR_FLAGS=-rsc

###############################################################################
# Configuration rules
###############################################################################
LIB_WFA=$(FOLDER_LIB)/libwfa.a
LIB_WFA_CPP=$(FOLDER_LIB)/libwfacpp.a
SUBDIRS=PairwiseAlignment/WFA2-lib/alignment \
        PairwiseAlignment/WFA2-lib/bindings/cpp \
        PairwiseAlignment/WFA2-lib/system \
        PairwiseAlignment/WFA2-lib/utils \
        PairwiseAlignment/WFA2-lib/wavefront

all: CC_FLAGS+=-O3 -march=native #-flto -ffat-lto-objects
all: setup build lib_wfa h4

# Debug target
debug: setup build

ASAN_OPT=-fsanitize=address -fsanitize=undefined -fsanitize=shift -fsanitize=alignment
ASAN_OPT+=-fsanitize=signed-integer-overflow -fsanitize=bool -fsanitize=enum
ASAN_OPT+=-fsanitize=pointer-compare -fsanitize=pointer-overflow -fsanitize=builtin

# AddressSanitizer target
asan: CC_FLAGS+=$(ASAN_OPT) -fno-omit-frame-pointer -fno-common
asan: setup build

###############################################################################
# Build rules
###############################################################################
# Ensure the setup step is performed first to create necessary directories
setup:
	@mkdir -p $(FOLDER_BUILD) $(FOLDER_BUILD_CPP) $(FOLDER_LIB)

# Build all subdirectories and libraries
build: $(SUBDIRS) lib_wfa

# Create the static libraries
lib_wfa: $(FOLDER_BUILD)/*.o $(FOLDER_BUILD_CPP)/*.o
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null
	$(AR) $(AR_FLAGS) $(LIB_WFA_CPP) $(FOLDER_BUILD)/*.o $(FOLDER_BUILD_CPP)/*.o 2> /dev/null

###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all

.PHONY: $(SUBDIRS)

###############################################################################
# Rules
###############################################################################
LIBS=-fopenmp -lm
ifeq ($(UNAME), Linux)
  LIBS+=-lrt 
endif
        
h4: *.cpp $(LIB_WFA) $(LIB_WFA_CPP)
	g++ $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) \
	./PairwiseAlignment/NeedlemanWunshReusable.cpp \
	./SuffixArray/parallel_import.cpp \
	./Utils/Arguments.cpp \
	./Utils/Fasta.cpp \
	./Utils/Graph.cpp \
	./Utils/Insertion.cpp \
	./Utils/NucleicAcidColumn.cpp \
	./Utils/Utils.cpp \
	./multi-thread/multi.cpp \
	./StarAlignment/StarAligner.cpp \
	stmsa.cpp -o halign4 -static-libstdc++ -std=c++17 -lpthread -lwfacpp $(LIBS)

# Clean target
clean: 
	rm -rf $(FOLDER_BUILD) $(FOLDER_BUILD_CPP) $(FOLDER_LIB) 2> /dev/null
	rm -rf $(FOLDER_TESTS)/*.alg $(FOLDER_TESTS)/*.log* 2> /dev/null
	rm -f halign4