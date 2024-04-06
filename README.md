# SHGD
Symmetric projected gradient descent (SHGD) for spectral compressed sensing 



demo.m  
demo how to use SHGD , in 1-D signal case. 

Comparison_time_Montecarlo.m   
time comparisons to reach different recovery accuracies with Monte Carlo experiments,1-D signal case   

SHGD.m  
SHGD algorithm, 1D case

PGD.m  
PGD algorithm, 1D case

FIHT.m 
 FIHT algorithm, 1D case

generate_signal_1D.m
generate 1-D spectral sparse signal

--------------------------------------------------------------------------------
./2Dcase

SHGD_2D.m

SHGD algorithm for 2D case.

robust_recovery_2D_plot.m

recover  a 2D spectral sparse signal  robustly with noise and make super-resolution in frequencies via MUSIC

call_music.m

MUSIC for 2D  signal super-resolution

conv_fft.m

2D signal fast convolution via FFT

fhmvmultiply_2D.m

2-level block Hankel vector fast multiplication

generate_signal_2D.m

generate 2-D spectral sparse signal
