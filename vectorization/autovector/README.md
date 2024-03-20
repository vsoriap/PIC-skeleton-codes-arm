# PIC Skeleton Codes:  Auto-Vectorization

These codes illustrate how to use OpenMP vectorization directives (e.g., #pragma omp simd).

For the 2D electrostatic:
  * no-vec = 35 nsec/particle/timestep
  * compiler vec = 18 nsec/particle/timestep
  * OpenMP vec = 

For the 2-1/2D electromagnetic:
  * no-vec = 29 nsec/particle/timestep
  * compiler vec = 26 nsec
  * OpenMP vec = 25 nsec

With SSE2 intrinsics one typically obtains about 3x speedup compared to no vectorization. Compiler vectorization achieves about 2x speedup. (All timings are on a 2.67GHz Intel Nehalem processor.)

1. 2D Parallel Electrostatic Spectral code:  vpic2
2. 2-1/2D Parallel Electromagnetic Spectral code:  vbpic2

### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


