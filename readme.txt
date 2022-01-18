This folder contains Matlab (R2019a) and C codes for the manuscript: 
C. Huang, A. Pouget and B. Doiron (2020) Internally generated population activity in cortical networks hinders information transmission. bioRxiv. 

++++++++++++++++++++++++++++++++

To use, first compile all the C codes with the mex compiler in Matlab. Specifically, in the Matlab command line, run the following commands:
mex EIF1DRFfastslowSyn.c
mex spktime2count.c

Sim_Ori_gabor_L1.m simulates a two-layer network with different parameter sets. Simulation code for Figs 1-4, 6, 8.  Saves spike counts.   

Sim_Ori_gabor_multiL.m simulates multi-layer networks. Simulation code for Fig. 9. Saves spike counts.   

Sim_Ori_gabor_L1_Pts.m simulates a two-layer network with long-range tuning similarity dependent connections. Simulation code for Fig. 10C, D.  Saves spike counts.   

Sim_Ori_gabor_multiL_Pts.m simulates multi-layer networks with long-range tuning similarity dependent connections. Simulation code for Fig. 10E, F.  Saves spike counts.

Sim_Spont.m  Simulate spontaneous dynamics with uncorrelated and untuned inputs from L4. Related to Figs. 5C, 7C. 

CollectSpk.m collect spike counts from every 500 simulations into one file, which will be used to compute Fisher information.  

FIdecoder_cluster_L1.m computes fisher information by both bias-correction methods and training linear decoders. Calls function FIdecoder.m  

CollectData.m   Collects Fisher information vs. number of neuron data from different parameter sets.

SpatialFreq.m  computes spatial and temporal spectrum from simulated spontaneous spiking activity. Related to Figs. 5C, 7C. Run Sim_Spont.m first.  

neuralfield2D.m  computes eigenvalues per wavenumber of the neural field model (Figs. 5A,B, 7A,B) 

CollectTuning.m  computes tuning curves for each network.  Need to run Sim_Ori_gabor_L1.m or Sim_Ori_gabor_L1_Pts.m first and set the range of orientations to be 'testp.theta0=0.02:.02:1; '.    

RF2D3layer.m is the main simulation function. It contains default parameter values and uses the mex file EIF1DRFfastslowSyn.c for integration. 

demo.m is a demonstration code for two-layer network simulations. 

genXspk.m generate spike trains of input from Layer 1. 

gen_weights.m generates weight matrices without tuning-dependent connections 

gen_weights_tuning.m generates weight matrices with tuning-dependent connections 

ori_map.m generates a columnar orientation map. 

corr_d.m computes correlation as a function of distance. 

spktime2count.c converts spike time data to spike counts. 

raster2D_ani.m generates movie of spike rasters from 2D spatial networks. 


+++++++++++++++++++++++++++++++++++

One simulation of a two-layer network for 20 sec takes about 1.5 hours CPU time and under 3 gb memory. 
(Sim_Ori_gabor_L1.m)



