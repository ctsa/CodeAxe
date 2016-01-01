#
#
# CodeAxe : phylogenetic analysis and simulation tools
#
#   http://www.phrap.org
#
#
# Copyright 2007 Christopher T Saunders (ctsa@u.washington.edu)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
#

# $Id: makefile 789 2007-09-06 01:06:59Z ctsa $

# requires lib paths for blas,lapack to be defined in the OS
#

BINTAG = CodeAxe

BIN_DIR = ../bin
OBJ_DIR = ../obj

DEF_INCLUDE_PATH = default_include
BOOST_PATH = $(DEF_INCLUDE_PATH)/boost-1.33.1-modified.subset

UTIL_DIR = util
UTIL_BIO_DIR = $(UTIL_DIR)/bio
UTIL_MATH_DIR = $(UTIL_DIR)/math
UTIL_GENERAL_DIR = $(UTIL_DIR)/general

ESSLLIB_PATH=$(HOME)/opt/code.lib/essl/lib





BINNAME = $(BIN_DIR)/$(BINTAG)

INCLUDE_DIRS = -I$(DEF_INCLUDE_PATH)
INCLUDE_DIRS += -I$(BOOST_PATH)

CXX_SRC = $(wildcard *.cc)
UTIL_BIO_SRC = $(wildcard $(UTIL_BIO_DIR)/*.cc)
UTIL_MATH_SRC = $(wildcard $(UTIL_MATH_DIR)/*.cc)
UTIL_GENERAL_SRC = $(wildcard $(UTIL_GENERAL_DIR)/*.cc)

OBJS  = $(patsubst %.cc,$(OBJ_DIR)/%.$(BUILD_TYPE).o,$(CXX_SRC))
OBJS += $(patsubst %.cc,$(OBJ_DIR)/%.$(BUILD_TYPE).o,$(UTIL_BIO_SRC))
OBJS += $(patsubst %.cc,$(OBJ_DIR)/%.$(BUILD_TYPE).o,$(UTIL_MATH_SRC))
OBJS += $(patsubst %.cc,$(OBJ_DIR)/%.$(BUILD_TYPE).o,$(UTIL_GENERAL_SRC))

DEPS  = $(patsubst %.cc,$(OBJ_DIR)/%.d,$(CXX_SRC))


GCCWARNFLAGS = -Wall -W -Wshadow -Wmissing-field-initializers -Wredundant-decls # -Wmissing-noreturn -Waggregate-return -Wunreachable-code -Wno-sign-compare
ICCWARNFLAGS = -Wall -w1 -wd411
GPP = g++
IXX = icpc

# shared defs:
CXX_DEP = g++
UNIFLAG =
CXXFLAGS = $(WARNFLAGS) $(OFLAGS) $(DFLAGS) $(UNIFLAG) $(INCLUDE_DIRS)
LDFLAGS = $(OFLAGS) $(DFLAGS) $(UNIFLAG)

LINUXLIBS=-llapack -lcblaswr -lcblas -latlas -lF77 # -lg2c

# gcc fast math is just not worth the trouble, specifically it:
# 1) causes undetectable nan's
# 2) is not consistent with the high res requirements of the probability model
# 3) seems to mess with the mkl
GCC_FAST_FLAGS=-O3 -DNDEBUG -finline-limit=6000 -fno-math-errno
GCC_FAST_FLAGS_32BIT=-malign-double
GCC_FAST_FLAGS_NON_PROF=-fomit-frame-pointer

ICC_FAST_FLAGS=-O3 -DNDEBUG -ip -Ob2

STRIP=touch
BLINK=echo
BLINK_NATIVE = (cd $(BIN_DIR); ln -fs $(BINTAG).$(BUILD_TYPE) $(BINTAG))

CXX = $(GPP)

ifeq ($(MAKECMDGOALS),gcc-linux-x86_64)
  BUILD_TYPE=gcc-linux-x86_64
  WARNFLAGS=$(GCCWARNFLAGS)
  OFLAGS = $(GCC_FAST_FLAGS) # -ffloat-store  #float-store prevents lapacks dlamch from getting stuck
  OFLAGS += -DUSE_BLAS -DUSE_LAPACK -DUSE_VML
  OFLAGS += $(GCC_FAST_FLAGS_NON_PROF)
#  OFLAGS += -pg   # gprof profiling
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = $(LINUXLIBS) -static
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lguide -lpthread -lvml # -static # static can occasionally cause crashes (threading issues?)
  STRIP = strip
  BLINK = $(BLINK_NATIVE)
endif

ifeq ($(MAKECMDGOALS),gcc-linux-i686)
  BUILD_TYPE=gcc-linux-i686
  WARNFLAGS=$(GCCWARNFLAGS)
  OFLAGS = $(GCC_FAST_FLAGS) # -ffloat-store  #float-store prevents lapacks dlamch from getting stuck
  OFLAGS += $(GCC_FAST_FLAGS_32BIT)
  OFLAGS += -DUSE_BLAS  -USE_LAPACK
  OFLAGS += $(GCC_FAST_FLAGS_NON_PROF)
#  OFLAGS += -pg       # gprof profiling
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = $(LINUXLIBS) -static
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lguide -lpthread # -static # static can occasionally cause crashes (threading issues?)
  STRIP = strip
endif

ifeq ($(MAKECMDGOALS),gcc-linux-p4)
  BUILD_TYPE=gcc-linux-p4
  WARNFLAGS=$(GCCWARNFLAGS)
  OFLAGS = $(GCC_FAST_FLAGS) -march=pentium4 -msse -msse2 -mfpmath=sse
  OFLAGS += $(GCC_FAST_FLAGS_32BIT)
  OFLAGS += -DUSE_BLAS  -DUSE_LAPACK -DUSE_VML
  OFLAGS += $(GCC_FAST_FLAGS_NON_PROF)
#  OFLAGS += -pg    # gprof profiling
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = $(LINUXLIBS)
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lguide -lpthread -lvml # -static # -lvml static can occasionally cause crashes (threading issues?) & so can dynamic :(
  STRIP=strip
  BLINK = $(BLINK_NATIVE)
endif

ifeq ($(MAKECMDGOALS),gcc-linux-debug)
  BUILD_TYPE=gcc-linux-debug
  WARNFLAGS=$(GCCWARNFLAGS)
  DFLAGS = -ggdb -fno-inline
  DFLAGS += -DDEBUG
  DFLAGS += -DUSE_BLAS -DUSE_LAPACK
#  DFLAGS += -pg
  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
  LDLIBS = $(LINUXLIBS)
#  LDFLAGS += -L$(MKL_LIB_PATH)
#  LDLIBS = -lmkl_lapack -lmkl -lguide -lpthread

  DFLAGS += -D_GLIBCXX_DEBUG
#  DFLAGS += -fmudflap
#  LDLIBS += -lmudflap
endif

ifeq ($(MAKECMDGOALS),gcc-linux-i686-debug)
  BUILD_TYPE=gcc-linux-i686-debug
  WARNFLAGS=$(GCCWARNFLAGS)
  DFLAGS = -ggdb -fno-inline
#  DFLAGS += -DDEBUG
  DFLAGS += -DUSE_BLAS -DUSE_LAPACK
#  DFLAGS += -pg
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = $(LINUXLIBS)
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lguide -lpthread

  DFLAGS += -D_GLIBCXX_DEBUG
#  DFLAGS += -fmudflap
#  LDLIBS += -lmudflap
endif

ifeq ($(MAKECMDGOALS),gcc-aix)
  BUILD_TYPE=gcc-aix
  WARNFLAGS=$(GCCWARNFLAGS)
  OFLAGS = -mcpu=powerpc64 -maix64 -O3 -DNDEBUG -ffast-math
  LDLIBS =
  STRIP = strip
endif


OSX_STD_DEF = -DDARWIN_HACK -DUSE_BLAS -DUSE_LAPACK
OSX_LDFLAGS = -framework veclib

ifeq ($(MAKECMDGOALS),gcc-osx)
  BUILD_TYPE=gcc-osx
  WARNFLAGS= $(GCCWARNFLAGS)
  OFLAGS =  $(GCC_FAST_FLAGS) # -fast -mcpu=7450     # -ffast-math -> this dies on tiger
  OFLAGS += $(GCC_FAST_FLAGS_32BIT)  # -- double alignment exposes malloc error, can't replicate in debug code so far
  OFLAGS += $(GCC_FAST_FLAGS_NON_PROF)
  OFLAGS +=  $(OSX_STD_DEF)
  LDFLAGS += $(OSX_LDFLAGS)
  STRIP = strip
  BLINK = $(BLINK_NATIVE)
endif

ifeq ($(MAKECMDGOALS),gcc-osx-debug)
  BUILD_TYPE=gcc-osx-debug
  WARNFLAGS= $(GCCWARNFLAGS)
  DFLAGS = -ggdb -fno-inline -D_GLIBCXX_DEBUG -DDEBUG # -I/System/Library/Frameworks/vecLib.framework/Headers
  DFLAGS +=  $(OSX_STD_DEF)
  LDFLAGS += $(OSX_LDFLAGS)
endif

ifeq ($(MAKECMDGOALS),gcc-osx-prescott)
  BUILD_TYPE=gcc-osx-prescott
  WARNFLAGS= $(GCCWARNFLAGS)
  OFLAGS =  $(GCC_FAST_FLAGS) $(GCC_FAST_FLAGS_32BIT) -march=prescott  -msse -msse2 -mfpmath=sse
  OFLAGS += $(GCC_FAST_FLAGS_NON_PROF)
  OFLAGS += $(OSX_STD_DEF)
  LDFLAGS += $(OSX_LDFLAGS)
  STRIP = strip
endif

ifeq ($(MAKECMDGOALS),icc-linux-it2)
  BUILD_TYPE=icc-linux-it2
  CXX = $(IXX)
  WARNFLAGS = $(ICCWARNFLAGS)
  OFLAGS = $(ICC_FAST_FLAGS) -mcpu=itanium2 -cxxlib-icc
#  OFLAGS += -p
  OFLAGS += -gcc-version=340
  OFLAGS += -DUSE_BLAS -DUSE_LAPACK -DUSE_VML
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = -llapack -lcblaswr -lcblas -latlas -lF77
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lvml -lguide -lpthread #-static
  STRIP = strip
endif

ifeq ($(MAKECMDGOALS),icc-linux-p4)
  BUILD_TYPE=icc-linux-p4
  CXX = $(IXX)
  WARNFLAGS = $(ICCWARNFLAGS)
  OFLAGS = $(ICC_FAST_FLAGS) -xN -cxxlib-icc -openmp -parallel
#  OFLAGS += -p
  OFLAGS += -gcc-version=340
  OFLAGS += -DUSE_BLAS -USE_LAPACK -DUSE_VML
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = -llapack -lcblaswr -lcblas -latlas -lF77
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lguide -lpthread #-static
  STRIP = strip
endif


ifeq ($(MAKECMDGOALS),icc-linux-x86_64)
  BUILD_TYPE=icc-linux-x86_64
  CXX = $(IXX)
  WARNFLAGS = $(ICCWARNFLAGS)
  OFLAGS = $(ICC_FAST_FLAGS) -xW -cxxlib-icc #-openmp -parallel
#  OFLAGS += -p
#  OFLAGS += -gcc-version=340
  OFLAGS += -DUSE_BLAS -DUSE_LAPACK -DUSE_VML
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = -llapack -lcblaswr -lcblas -latlas -lF77
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl_lapack -lmkl -lvml -lguide -lpthread #-static
  STRIP = strip
endif



ifeq ($(MAKECMDGOALS),icc-linux-debug)
  BUILD_TYPE=icc-linux-debug
  CXX = $(IXX)
  WARNFLAGS = $(ICCWARNFLAGS)
  DFLAGS = -DDEBUG -g -cxxlib-icc # -gcc-name=/usr/bin/gcc
  DFLAGS += -DUSE_BLAS -DUSE_LAPACK
#  CFLAGS += -ansi
#  LDFLAGS += -L$(CBLAS_LIB_PATH) -L$(CLAPACK_LIB_PATH)
#  LDLIBS = -llapack -lcblaswr -lcblas -latlas -lF77
  LDFLAGS += -L$(MKL_LIB_PATH)
  LDLIBS = -lmkl -lmkl_lapack -lguide -lpthread
endif

XLC_ANY_FLAGS=-qrtti=dynamiccast

ifeq ($(MAKECMDGOALS),xlc-aix)
  BUILD_TYPE=xlc-aix
  CXX = xlC
  OFLAGS = -DXLC_HACK -DNDEBUG
  OFLAGS += -O3 -qmaxmem=-1 -qstrict -qinline=300 -qarch=pwr4 -qtune=pwr4
#  OFLAGS += -qipa=level=2 -qipa=partition=large -qlibansi   # ipa opts
#  OFLAGS += -qthreaded -qsmp=auto                           # threading
  OFLAGS += -DUSE_BLAS -DUSE_LAPACK
  OFLAGS += -DBOOST_NO_INCLASS_MEMBER_INITIALIZATION -DBOOST_NO_INTEGRAL_INT64_T -DNO_BOOST_POOL_ALLOCATOR
#  OFLAGS += -g                                # profiling
  CXXFLAGS += $(XLC_ANY_FLAGS)
  LDFLAGS += -bmaxdata:0x70000000 -L$(CBLAS_LIB_PATH) #-L$(CLAPACK_LIB_PATH)
#  LDLIBS += -lptcblas -latlas -lpthread                    # threaded BLAS
  LDLIBS = -llapack -lcblaswr -lcblas -latlas -lF77 -lI77   # non-threaded BLAS
  LDLIBS += -L$(ESSLLIB_PATH) -lessl
  LDLIBS += -lpthread    # for boost::pool_allocator
  STRIP = strip
endif

ifeq ($(MAKECMDGOALS),xlc-aix-debug)
  BUILD_TYPE = xlc-aix-debug
  CXX = xlC
  DFLAGS = -DXLC_HACK -bmaxdata:0x70000000 -g -DDEBUG
  DFLAGS += -DUSE_BLAS
  CXXFLAGS += $(XLC_ANY_FLAGS)
  LDFLAGS += -L$(CBLAS_LIB_PATH) -DFLAGS
  DFLAGS += -DBOOST_NO_INCLASS_MEMBER_INITIALIZATION
#  LDLIBS = -llapack -lcblaswr -lF77
  LDLIBS = -lcblas -latlas
  LDLIBS += -L$(ESSLLIB_PATH) -lessl
  LDLIBS += -lpthread    # for boost::pool_allocator
endif


CREATE_OBJ_DIR=$(OBJ_DIR)/.dummy


#
# user targets
#
print_targets:
	$(user-targets-print)

define user-targets-print
@echo
@echo --make target options--
@echo
@echo gcc-linux-i686
@echo gcc-linux-x86_64
@echo gcc-linux-p4
@echo gcc-linux-debug
@echo gcc-linux-i686-debug
@echo gcc-aix
@echo gcc-osx
@echo gcc-osx-prescott
@echo gcc-osx-debug
@echo icc-linux-it2
@echo icc-linux-p4
@echo icc-linux-x86_64
@echo icc-linux-debug
@echo xlc-aix
@echo xlc-aix-debug
@echo
@echo clean
@echo
endef

gcc-linux-i686: build
gcc-linux-x86_64: build
gcc-linux-p4: build
gcc-linux-debug: build
gcc-linux-i686-debug: build
gcc-aix: build
gcc-osx: build
gcc-osx-prescott: build
gcc-osx-debug: build
icc-linux-it2: build
icc-linux-p4: build
icc-linux-x86_64: build
icc-linux-debug: build
xlc-aix: build
xlc-aix-debug: build

build: $(BINNAME).$(BUILD_TYPE)


$(CREATE_OBJ_DIR):
	mkdir -p $(OBJ_DIR)/$(UTIL_BIO_DIR)
	mkdir -p $(OBJ_DIR)/$(UTIL_MATH_DIR)
	mkdir -p $(OBJ_DIR)/$(UTIL_GENERAL_DIR)
	touch $(CREATE_OBJ_DIR)

$(BIN_DIR):
	mkdir -p $@


$(BINNAME).$(BUILD_TYPE): $(CREATE_OBJ_DIR) $(OBJS) $(BIN_DIR)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(BINNAME).$(BUILD_TYPE)
	$(STRIP) $(BINNAME).$(BUILD_TYPE)
	$(BLINK)


clean:
	$(RM) -r $(OBJ_DIR)
	$(RM) nfs*


$(OBJ_DIR)/%.$(BUILD_TYPE).o : %.cc
	@echo
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.d: $(CREATE_OBJ_DIR) %.cc
	@echo
	$(CXX_DEP) $(INCLUDE_DIRS) -MM -MT $(patsubst %.cc,$(OBJ_DIR)/%.'$$(BUILD_TYPE)'.o,$*.cc) $*.cc > $@
	echo "$@ : $*.cc" >> $@


ifneq ($(MAKECMDGOALS),clean)
include $(DEPS)
endif


