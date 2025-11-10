# Fractional Brownian Motion (fBM) C Implementation Review
This project aims to validate the correctness of a C implementation of fractional Brownian motion (fBM) by checking its outputs against theoretical probability grounds. The implementation itself is included in the project, see `include/fbm.h`, and this review focuses solely on the validation process and results.

# Results summary
Overall statistical properties were present, but since a scaling factor was missing in the fft implementation the produced trajectories were off scale: this is shown by the iterated logarithm envelop on some of the plots.
