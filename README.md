# Fractional Brownian Motion (fBm) C Implementation Review
This project aims to validate the correctness of a C implementation of fractional Brownian motion (fBm) by checking its outputs against theoretical probability grounds. The implementation itself is included in the project, see `include/fbm.h`, and this review focuses solely on the validation process and results.

# Results summary
Major statistical properties were present in the generated fBm paths, but since a scaling factor was missing in the fft implementation (rendering the Fourier transform to be non-unitary) the produced trajectories were off scale in an assymptotical sense: this is shown by the iterated logarithm envelop on some of the plots (see the folder python/analysis_results_*).
