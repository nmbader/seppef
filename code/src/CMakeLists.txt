#Add directories containing the header files for the code
include_directories(include ../include/base ../include/operators ../include/solvers)
include_directories (include ../../external/gnuplot_iostream/)

#Include genericIO directory (locally)
#include_directories(../../../external/include)

#Adding path to link genericIO libraries (locally)
#link_directories(../../../external/lib)

#set the directory to source code
set(SRC_DIR ./)

set (SOURCE_PEF2D ${SRC_DIR}/base/pef2D.cpp
		${SRC_DIR}/base/varPef2D.cpp
		${SRC_DIR}/base/functions.cpp
		${SRC_DIR}/base/floatHyperExt.cpp
		${SRC_DIR}/base/complexHyperExt.cpp
		${SRC_DIR}/base/doubleHyperExt.cpp
		${SRC_DIR}/operators/oper1D.cpp
		${SRC_DIR}/operators/oper2D.cpp
		${SRC_DIR}/operators/matrixMult.cpp
		${SRC_DIR}/operators/convolution2D.cpp
		${SRC_DIR}/operators/identity1D.cpp
		${SRC_DIR}/operators/identity2D.cpp
		${SRC_DIR}/operators/selection2D.cpp
		${SRC_DIR}/operators/weighting2D.cpp
		${SRC_DIR}/operators/zeroPhaseFiltering2D.cpp
		${SRC_DIR}/operators/fxTransform.cpp
		${SRC_DIR}/operators/fkTransform.cpp
		${SRC_DIR}/operators/nsPef2D.cpp
		${SRC_DIR}/operators/expNsPef2D.cpp
		${SRC_DIR}/operators/varDatConv2D.cpp
		${SRC_DIR}/operators/varConv2D.cpp
		${SRC_DIR}/operators/expVarConv2D.cpp
		${SRC_DIR}/solvers/lsolver.cpp
		${SRC_DIR}/solvers/lsolverReg.cpp
		${SRC_DIR}/solvers/cgls.cpp
		${SRC_DIR}/solvers/cglsReg.cpp
		${SRC_DIR}/solvers/cglsSuper.cpp
		${SRC_DIR}/solvers/sdls.cpp
	)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/../external/sep-ioLibs/)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../external/sep-ioLibs/external/genericIO/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external/sep-ioLibs/external/genericIO/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external/sep-ioLibs/external/genericIO/external/sepVector/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external/sep-ioLibs/external/genericIO/external/sepVector/external/hypercube/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external/sep-ioLibs/include)

#building the project library
add_library(pef2d STATIC ${SOURCE_PEF2D})
target_link_libraries(pef2d boost_iostreams boost_system genericCpp hypercube sepVector sep fftw3f ${TBB_LIBRARIES})

#add_executable(FILTER.x ${SRC_DIR}/main/filter.cpp)
#target_link_libraries(FILTER.x basic operators genericCpp sepVector hypercube)
#
#add_executable(TRANSPOSE.x ${SRC_DIR}/main/transpose.cpp)
#target_link_libraries(TRANSPOSE.x pef2d genericCpp sepVector hypercube)

add_executable(RESAMP.x ${SRC_DIR}/main/resamp.cpp)
target_link_libraries(RESAMP.x  pef2d genericCpp sepVector hypercube)

add_executable(FK_SPECTRUM.x ${SRC_DIR}/main/fk_spectrum.cpp)
target_link_libraries(FK_SPECTRUM.x pef2d genericCpp sepVector hypercube)

add_executable(FX_SPECTRUM.x ${SRC_DIR}/main/fx_spectrum.cpp)
target_link_libraries(FX_SPECTRUM.x pef2d genericCpp sepVector hypercube)

add_executable(FK_2_FV.x ${SRC_DIR}/main/fk_2_fv.cpp)
target_link_libraries(FK_2_FV.x pef2d genericCpp sepVector hypercube)

add_executable(CONVOLVE2D.x ${SRC_DIR}/main/convolve2D.cpp)
target_link_libraries(CONVOLVE2D.x pef2d  genericCpp sepVector hypercube)

add_executable(FK_FILTER.x ${SRC_DIR}/main/fk_filter.cpp)
target_link_libraries(FK_FILTER.x pef2d genericCpp sepVector hypercube)

add_executable(GENERATE_PLANE_WAVE.x ${SRC_DIR}/main/generate_plane_wave.cpp)
target_link_libraries(GENERATE_PLANE_WAVE.x pef2d genericCpp sepVector hypercube)

add_executable(GENERATE_RANDOM_NOISE.x ${SRC_DIR}/main/generate_random_noise.cpp)
target_link_libraries(GENERATE_RANDOM_NOISE.x pef2d genericCpp sepVector hypercube)

add_executable(GENERATE_HYPERBOLIC_WAVE.x ${SRC_DIR}/main/generate_hyperbolic_wave.cpp)
target_link_libraries(GENERATE_HYPERBOLIC_WAVE.x pef2d genericCpp sepVector hypercube)

add_executable(GENERATE_SPIKE.x ${SRC_DIR}/main/generate_spike.cpp)
target_link_libraries(GENERATE_SPIKE.x pef2d genericCpp sepVector hypercube)

add_executable(GENERATE_WAVELET.x ${SRC_DIR}/main/generate_wavelet.cpp)
target_link_libraries(GENERATE_WAVELET.x pef2d sep)

add_executable(GENERATE_MONOCHROMATIC_WAVE.x ${SRC_DIR}/main/generate_monochromatic_wave.cpp)
target_link_libraries(GENERATE_MONOCHROMATIC_WAVE.x pef2d sep)

add_executable(EXPAND.x ${SRC_DIR}/main/expand.cpp)
target_link_libraries(EXPAND.x pef2d sep)

#add_executable(CUT.x ${SRC_DIR}/main/cut.cpp)
#target_link_libraries(CUT.x pef2d sep)

#add_executable(PAD.x ${SRC_DIR}/main/pad.cpp)
#target_link_libraries(PAD.x pef2d sep)

add_executable(DIVERGE.x ${SRC_DIR}/main/diverge.cpp)
target_link_libraries(DIVERGE.x pef2d sep)

add_executable(DUPLICATE.x ${SRC_DIR}/main/duplicate.cpp)
target_link_libraries(DUPLICATE.x pef2d sep)

add_executable(STAT.x ${SRC_DIR}/main/stat.cpp)
target_link_libraries(STAT.x pef2d sep)

add_executable(ESTIMATE_PEF.x ${SRC_DIR}/main/estimate_pef.cpp)
target_link_libraries(ESTIMATE_PEF.x pef2d genericCpp)

add_executable(ESTIMATE_VAR_PEF.x ${SRC_DIR}/main/estimate_var_pef.cpp)
target_link_libraries(ESTIMATE_VAR_PEF.x  pef2d genericCpp)

add_executable(SOLVE_CONV.x ${SRC_DIR}/main/solve_conv.cpp)
target_link_libraries(SOLVE_CONV.x pef2d genericCpp)

add_executable(SOLVE_NS_PEF.x ${SRC_DIR}/main/solve_ns_pef.cpp)
target_link_libraries(SOLVE_NS_PEF.x pef2d genericCpp sepVector hypercube)

add_executable(SOLVE_EXP_NS_PEF.x ${SRC_DIR}/main/solve_exp_ns_pef.cpp)
target_link_libraries(SOLVE_EXP_NS_PEF.x pef2d genericCpp sepVector hypercube)

add_executable(SOLVE_VAR_CONV.x ${SRC_DIR}/main/solve_var_conv.cpp)
target_link_libraries(SOLVE_VAR_CONV.x pef2d genericCpp sepVector hypercube)

add_executable(SOLVE_EXP_VAR_CONV.x ${SRC_DIR}/main/solve_exp_var_conv.cpp)
target_link_libraries(SOLVE_EXP_VAR_CONV.x pef2d genericCpp sepVector hypercube)

add_executable(SOLVE_CONV_MASK.x ${SRC_DIR}/main/solve_conv_mask.cpp)
target_link_libraries(SOLVE_CONV_MASK.x pef2d genericCpp sepVector hypercube)

add_executable(SOLVE_PVAR_REG_PEXPVAR.x ${SRC_DIR}/main/solve_Pvar_reg_Pexpvar.cpp)
target_link_libraries(SOLVE_PVAR_REG_PEXPVAR.x pef2d genericCpp sepVector hypercube)

add_executable(SOLVE_PEXPVAR_REG_PVAR.x ${SRC_DIR}/main/solve_Pexpvar_reg_Pvar.cpp)
target_link_libraries(SOLVE_PEXPVAR_REG_PVAR.x pef2d  genericCpp sepVector hypercube)

add_executable(SOLVE_CONV_REG_IDENTITY.x ${SRC_DIR}/main/solve_conv_reg_identity.cpp)
target_link_libraries(SOLVE_CONV_REG_IDENTITY.x pef2d   genericCpp sepVector hypercube)

add_executable(SOLVE_CONV_REG_PEXPVAR.x ${SRC_DIR}/main/solve_conv_reg_Pexpvar.cpp)
target_link_libraries(SOLVE_CONV_REG_PEXPVAR.x pef2d  genericCpp sepVector hypercube)


#define the directories for the constructed libraries and executables
#set_target_properties(basic operators solvers pef2d
#						ESTIMATE_VAR_PEF.x SOLVE_VAR_CONV.x SOLVE_EXP_VAR_CONV.x MAIN.x
#						SOLVE_NS_PEF.x SOLVE_CONV.x SOLVE_EXP_NS_PEF.x
#						FILTER.x STAT.x
#						FX_SPECTRUM.x PAD.x CUT.x FK_FILTER.x GENERATE_SPIKE.x GENERATE_PLANE_WAVE.x GENERATE_HYPERBOLIC_WAVE.x
#						GENERATE_RANDOM_NOISE.x  PLOT1D.x TRANSPOSE.x
#						CONVOLVE2D.x FK_SPECTRUM.x PLOT2D.x EXPAND.x ESTIMATE_PEF.x PEF_LF_DEN_V1.x PEF_LF_DEN_V2.x PEF_LF_DEN_V3.x
#    PROPERTIES
#    ARCHIVE_OUTPUT_DIRECTORY "lib"
#    LIBRARY_OUTPUT_DIRECTORY "lib"
#    RUNTIME_OUTPUT_DIRECTORY "bin"
#)


#installing libraries and executables
install(TARGETS pef2d DESTINATION lib)
install(TARGETS
			GENERATE_RANDOM_NOISE.x GENERATE_PLANE_WAVE.x GENERATE_HYPERBOLIC_WAVE.x GENERATE_SPIKE.x GENERATE_MONOCHROMATIC_WAVE.x GENERATE_WAVELET.x
			ESTIMATE_PEF.x ESTIMATE_VAR_PEF.x
			SOLVE_CONV.x SOLVE_VAR_CONV.x SOLVE_EXP_NS_PEF.x SOLVE_NS_PEF.x SOLVE_EXP_VAR_CONV.x SOLVE_CONV_MASK.x SOLVE_PVAR_REG_PEXPVAR.x SOLVE_CONV_REG_IDENTITY.x SOLVE_CONV_REG_PEXPVAR.x SOLVE_PEXPVAR_REG_PVAR.x
			FK_SPECTRUM.x FX_SPECTRUM.x FK_2_FV.x
			#PAD.x CUT.x  
			FK_FILTER.x  CONVOLVE2D.x
			# MAIN.x
			# FILTER.x TRANSPOSE.x 
			RESAMP.x
			STAT.x
			EXPAND.x
			DIVERGE.x DUPLICATE.x
			DESTINATION bin)

#install Python executables
install( FILES ../python/seppyio.py PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION bin)
install( FILES ../python/pyPlot1D.py PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION bin)
install( FILES ../python/pyPlot2D.py PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION bin)
install( FILES ../python/pyPlot3D.py PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION bin)
