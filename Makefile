### Makefile to build the code and generate some examples: Low frequency de-noising using stationary and non-stationary prediction-error filters
### To run the examples, SEPlib must be installed a priori and included in the PATH

current_dir = $(shell pwd)

# Initialize local directories
clone:
	mkdir -p build
	mkdir -p local
	mkdir -p fig
	mkdir -p dat	

# Compile and build the library using CMake
install: clone
	cd ./build; cmake -DCMAKE_INSTALL_PREFIX=../local ../code/; make install -j8; cd ../
	export DATAPATH=${current_dir}/dat/
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${current_dir}/local/lib/

####################################################################

### Variables definition

# Executables directory
BINDIR = ./local/bin

# Figures directory
FIGDIR = ./fig

# Data directory
DATDIR = ./dat

B = ${BINDIR}
F = ${FIGDIR}
D = ${DATDIR}

################# Synthetic example, stationary data ###############

#filters
${D}/j1_lp_10.H:
	${B}/GENERATE_WAVELET.x nz=101 dz=0.01 high_cutoff=0.2 type=sinc phase=zero > $@

${D}/j1_lp_4.H:
	${B}/GENERATE_WAVELET.x nz=131 dz=0.01 high_cutoff=0.08 type=sinc phase=zero > $@

${D}/j1_lp_5.H:
	${B}/GENERATE_WAVELET.x nz=121 dz=0.01 high_cutoff=0.1 type=sinc phase=zero > $@

${D}/j1_hp_2p5.H:
	${B}/GENERATE_WAVELET.x nz=151 dz=0.01 low_cutoff=0.05 type=sinc phase=zero > $@

${D}/j1_butterworth.H:
	${B}/GENERATE_WAVELET.x nz=801 dz=0.01 type=butterworth low_cutoff=0.05 high_cutoff=0.2 half_order=6 shift=400 > $@

# data and noise
${D}/j1_signal_2p5_10_1sp.H: ${D}/j1_butterworth.H
	${B}/GENERATE_PLANE_WAVE.x n1=401 n2=401 o1=0 o2=0 d1=0.01 d2=10 z0=0.1 vel=2000 amp=10 interp=false | Filter \
	filter=$< pad=zero | ${B}/FK_FILTER.x kmin=-0.3 kmax=0.3 taper=0.1 > $@

${D}/j1_signal_2p5_4_1sp.H: ${D}/j1_signal_2p5_10_1sp.H ${D}/j1_lp_4.H
	Filter < $< filter=$(word 2,$^) pad=zero > $@

${D}/j1_signal_2p5_5_1sp.H: ${D}/j1_signal_2p5_10_1sp.H ${D}/j1_lp_5.H
	Filter < $< filter=$(word 2,$^) pad=zero > $@

${D}/j1_signal_5_10_1sp.H: ${D}/j1_signal_2p5_10_1sp.H ${D}/j1_signal_2p5_5_1sp.H
	Math file1=$< file2=$(word 2,$^) exp=file1-file2 > $@

${D}/j1_random_noise_2p5_4_1sp.H: ${D}/j1_hp_2p5.H ${D}/j1_lp_4.H
	${B}/GENERATE_RANDOM_NOISE.x n1=401 n2=401 o1=0 o2=0 d1=0.01 d2=10 min=-4 max=4 seed=1234 | Filter \
	filter=$< pad=zero | Filter filter=$(word 2,$^) pad=zero > $@

${D}/j1_monochromatic_2p5_1sp.H:
	${B}/GENERATE_MONOCHROMATIC_WAVE.x n1=401 n2=401 o1=0 o2=0 d1=0.01 d2=10 vel=-1500 amp=0.2 phase=0 w=0.05 > $@

${D}/j1_monochromatic_3p5_1sp.H:
	${B}/GENERATE_MONOCHROMATIC_WAVE.x n1=401 n2=401 o1=0 o2=0 d1=0.01 d2=10 vel=5000 amp=0.2 phase=0 w=0.07 > $@

${D}/j1_mono_noise_1sp.H: ${D}/j1_monochromatic_2p5_1sp.H ${D}/j1_monochromatic_3p5_1sp.H
	Add file1=$< file2=$(word 2,$^) > $@

${D}/j1_coherent_noise_2p5_4_1sp.H: ${D}/j1_butterworth.H ${D}/j1_lp_4.H
	${B}/GENERATE_PLANE_WAVE.x n1=601 n2=401 o1=-1 o2=0 d1=0.01 d2=10 z0=4 vel=-1000 amp=20 > temp1.H
	${B}/GENERATE_PLANE_WAVE.x n1=601 n2=401 o1=-1 o2=0 d1=0.01 d2=10 z0=-1 vel=1000 amp=20 > temp2.H
	Add file1=temp1.H file2=temp2.H | ${B}/FK_FILTER.x kmin=-0.3 kmax=0.3 taper=0.1 > temp.H
	Filter < temp.H filter=$< pad=zero | Filter filter=$(word 2,$^) pad=zero | Window \
	f1=101 n1=401 > $@
	Rm -f temp1.H temp2.H temp.H

${D}/j1_all_noise_2p5_4_1sp.H: ${D}/j1_random_noise_2p5_4_1sp.H ${D}/j1_mono_noise_1sp.H ${D}/j1_coherent_noise_2p5_4_1sp.H
	Add file1=$< file2=$(word 2,$^) file3=$(word 3,$^) > $@

${D}/j1_data_2p5_10_1sp.H: ${D}/j1_signal_2p5_10_1sp.H ${D}/j1_all_noise_2p5_4_1sp.H
	Add file1=$< file2=$(word 2,$^) > $@

${D}/j1_data_2p5_4_1sp.H: ${D}/j1_data_2p5_10_1sp.H ${D}/j1_lp_4.H
	Filter < $< filter=$(word 2,$^) pad=zero > $@

${D}/j1_data_2p5_5_1sp.H: ${D}/j1_data_2p5_10_1sp.H ${D}/j1_lp_5.H
	Filter < $< filter=$(word 2,$^) pad=zero > $@

${D}/j1_data_5_10_1sp.H: ${D}/j1_data_2p5_10_1sp.H ${D}/j1_data_2p5_5_1sp.H
	Math file1=$< file2=$(word 2,$^) exp=file1-file2 > $@

###################################################################
############################## PEF Prediction - stationary #########################
###################################################################

pef_zsize = 15
pef_xsize = 8
pef_lead = 7

# stationary PEF
${D}/j1_spef_5_10.H: ${D}/j1_data_5_10_1sp.H
	${B}/ESTIMATE_PEF.x < $< n1=${pef_zsize} n2=${pef_xsize} lead=${pef_lead} gap=0 niter=5 threshold=-0.001 cross_bounds=false > $@

${D}/j1_spef_5_10_exp.H: ${D}/j1_spef_5_10.H
	${B}/EXPAND.x < $< n1=true n2=true > $@

###################################################################
############################## PEF Prediction - non-stationary #########################
###################################################################

# Bank of filters
${D}/j1_nspef_5_10_S.H: ${D}/j1_data_5_10_1sp.H
	${B}/ESTIMATE_VAR_PEF.x < $< n1=${pef_zsize} n2=${pef_xsize} lead=${pef_lead} xinc=60 zinc=60 \
	niter=250 threshold=-0.001 obj_func=$@.func > $@

${D}/j1_nspef_5_10_D.H: ${D}/j1_data_2p5_4_1sp.H
	${B}/ESTIMATE_VAR_PEF.x < $< n1=${pef_zsize} n2=${pef_xsize} lead=${pef_lead} xinc=60 zinc=60 \
	niter=250 threshold=-0.001 obj_func=$@.func > $@

############### De-noising by Stationary PEF ###############
# solving S.m ~ 0 ; R.(dl - m) ~ 0 ; m_0 = d ; R = Id
# S: signal PEF (can be an expanded PEF from a clean HF data)
###################################################################

${D}/j1_data_2p5_4_1sp_cp.H: ${D}/j1_data_2p5_4_1sp.H
	Cp < $< >$@

${D}/j1_smod_2p5_4.H: ${D}/j1_data_2p5_4_1sp.H ${D}/j1_spef_5_10_exp.H ${D}/j1_data_2p5_4_1sp_cp.H
	${B}/SOLVE_CONV_REG_IDENTITY.x < $< fil=$(word 2,$^) prior=$(word 3,$^) cross_bounds=true \
	niter=20 threshold=-0.001 eps=0.05 solver=cgls \
	obj_func=$(subst mod,func,$@) mod_res=$(subst mod,mod_res,$@) dat_res=$(subst mod,dat_res,$@) output=$@

###################################################################
############### De-noising by 2D Non-Stationary PEF ###############
# solving  S.s=0 & eps.N.(s - d) = 0
# N: noise PEF
# S: signal PEF (can be an expanded PEF from a clean HF data)
# N could be replaced by D, the data PEF
# eps << : better noise removal
# eps >> : better signal preservation
###################################################################

${D}/j1_nsmod_2p5_4.H: ${D}/j1_data_2p5_4_1sp.H ${D}/j1_nspef_5_10_S.H ${D}/j1_nspef_5_10_D.H ${D}/j1_data_2p5_4_1sp_cp.H
	${B}/SOLVE_PEXPVAR_REG_PVAR.x < $< fil=$(word 2,$^) filb=$(word 3,$^) prior=$(word 4,$^) n1=${pef_zsize} n2=${pef_xsize} lead=${pef_lead} xinc=60 zinc=60 \
	niter=50 threshold=-0.001 eps=0.1 solver=cgls \
	obj_func=$@.func mod_res=$@.mres dat_res=$@.dres output=$@

######################### Figures #########################################

${F}/j1_signal_2p5_4_1sp.pdf: ${D}/j1_signal_2p5_4_1sp.H
	${B}/pyPlot2D.py input=$< xsize=6.66 zsize=6 xlabel="Location x (m)" zlabel="Time (s)" colormap=Greys min=-0.5 max=0.5 interpolation=sinc format=pdf colorbar=0 output=$@

${F}/j1_data_2p5_4_1sp.pdf: ${D}/j1_data_2p5_4_1sp.H
	${B}/pyPlot2D.py input=$< xsize=6.66 zsize=6 xlabel="Location x (m)" zlabel="Time (s)" colormap=Greys min=-0.5 max=0.5 interpolation=sinc format=pdf colorbar=0 output=$@

${F}/j1_smod_2p5_4.pdf: ${D}/j1_smod_2p5_4.H
	${B}/pyPlot2D.py input=$< xsize=6.66 zsize=6 xlabel="Location x (m)" zlabel="Time (s)" colormap=Greys min=-0.5 max=0.5 interpolation=sinc format=pdf colorbar=0 output=$@

${F}/j1_nsmod_2p5_4.pdf: ${D}/j1_nsmod_2p5_4.H
	${B}/pyPlot2D.py input=$< xsize=6.66 zsize=6 xlabel="Location x (m)" zlabel="Time (s)" colormap=Greys min=-0.5 max=0.5 interpolation=sinc format=pdf colorbar=0 output=$@

####################################################################
data: ${F}/j1_signal_2p5_4_1sp.pdf ${F}/j1_data_2p5_4_1sp.pdf ${F}/j1_smod_2p5_4.pdf ${F}/j1_nsmod_2p5_4.pdf

figures: ${F}/j1_signal_2p5_4_1sp.pdf ${F}/j1_data_2p5_4_1sp.pdf ${F}/j1_smod_2p5_4.pdf ${F}/j1_nsmod_2p5_4.pdf

burn:
	rm -f ${F}/*

clean:
	rm -f ${F}/*
	rm -f ${D}/*

uninstall:
	rm -rf ./build
	rm -rf ./local
	rm -rf ./fig
	rm -rf ./dat
