TARGET		= STATALTEXRamjet
# Makefile for SCRAMSPACE.mexa64

# HyPro
HyProIncl	= -I${HyPro}/include -I${HyPro}/test/${TARGET}
HyProLib    = -L${HyPro}/build -lHyPro

# OpenFoam
OpenFoamIncl = -I${FOAM_SRC}/OpenFOAM/lnInclude -I${FOAM_SRC}/finiteVolume/lnInclude -I${FOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude -I${FOAM_SRC}/OSspecific/POSIX/lnInclude -I${FOAM_SRC}/thermophysicalModels/specie/lnInclude
OpenFoamLib  = -L${FOAM_LIBBIN} -lspecie -lthermophysicalFunctions -lOpenFOAM -L${FOAM_LIBBIN}/${FOAM_MPI} -lPstream

# Overall
Include = ${HyProIncl} ${OpenFoamIncl}
Library = ${HyProLib} ${OpenFoamLib} ${GPCLib}

# Compiler flags.  Use -O to optimise
Definitions = -DHyProRescource=${HyPro}/resources -DSTXRescource=${HyPro}/test/STATALTEXRamjet/resources -DWM_DP -DNoRepository -D_GLIBCXX_DEBUG -Wno-unused-local-typedefs
CPPFLAGS	= ${Definitions} -O0 -g3 -Wall -c -fpermissive -std=c++0x -fPIC -MMD -MP

all: ${TARGET}.mexa64

clean:
	rm -f ${TARGET}.mexa64

${TARGET}.mexa64:
	mex ${Include} -g CXXFLAGS='${CXXFLAGS} ${CPPFLAGS}' ${Library} ${TARGET}.cpp
