################################################################################
#
# Copyright 1993-2006 NVIDIA Corporation.  All rights reserved.
#
# NOTICE TO USER:   
#
# This source code is subject to NVIDIA ownership rights under U.S. and 
# international Copyright laws.  
#
# NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
# CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
# IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
# REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
# MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.   
# IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
# OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE 
# OR PERFORMANCE OF THIS SOURCE CODE.  
#
# U.S. Government End Users.  This source code is a "commercial item" as 
# that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of 
# "commercial computer software" and "commercial computer software 
# documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995) 
# and is provided to the U.S. Government only as a commercial end item.  
# Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
# 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
# source code with only those rights set forth herein.
#
################################################################################
#
# Build script for project
#
################################################################################

PROJECTDIR   := $(CURDIR)
ROOTDIR		 := /Developer/GPU\ Computing
CUDA_SDK_DIR := $(ROOTDIR)/C

SRCDIR     := $(PROJECTDIR)
ROOTBINDIR := $(PROJECTDIR)
BINDIR     := $(PROJECTDIR)
ROOTOBJDIR := $(ROOTBINDIR)
LIBDIR     := $(CUDA_SDK_DIR)/lib
COMMONDIR  := $(CUDA_SDK_DIR)/common
SHAREDDIR  := $(ROOTDIR)/shared
JRFXGLDIR  := /Users/josericardo/Projects/jrfxgl

LIB := -L$(JRFXGLDIR)/release -L/opt/local/lib -lJRFXGL -lSDL -lSDLmain -lGlew -framework cocoa
INCLUDES += -I$(JRFXGLDIR) -I/opt/local/include -I/opt/local/include/SDL -I/opt/local/include/FTGL \
            -I/opt/local/include/freetype2

verbose=1
x86_64=1
USEGLLIB=1

# add command line parameters so we can target multiple architectures
GENCODE_SM10 := -gencode=arch=compute_10,code=\"sm_10,compute_10\"
GENCODE_SM11 := -gencode=arch=compute_11,code=\"sm_11,compute_11\"
GENCODE_SM12 := -gencode=arch=compute_12,code=\"sm_12,compute_12\"
GENCODE_SM20 := -gencode=arch=compute_20,code=\"sm_20,compute_20\"

GENCODE_ARCH    := $(GENCODE_SM12)

#NVCCFLAGS := -D__RECURSION_GPU__
#NVCCFLAGS := --dryrun
NVCCFLAGS += -DCUDA -DNVCC -D_CUDA -D__CUDA__ -D__APPLE__
#CXXFLAGS := -D__APPLE__
CXXFLAGS += -D__CUDA__ -D_CUDA -DCUDA -D__APPLE__

# Add source files here
EXECUTABLE	:= FluidGrid

# CUDA source files (compiled with cudacc)
CUFILES		:= GPU_Kernels/K_Boundaries_Common.cu \
               GPU_Kernels/K_DeviceFunctions.cu GPU_Kernels/K_FluidSim.cu \
               GPU_Kernels/K_HashGrid.cu GPU_Kernels/K_MathUtils.cu \
            GPU_Kernels/K_MathUtils.cu

# CUDA dependency files
CU_DEPS		:= 

# C/C++ source files (compiled with gcc / g++)
CCFILES	 := ApplicationData.cpp AppMain.cpp Emitter.cpp FluidParticleRenderMaterial.cpp \
            FluidRenderer.cpp GPUSPHFluid.cpp HeterogeneousSimulator.cpp main.cpp
		    
	

################################################################################
# Rules and targets	

include common_custom.mk
    
