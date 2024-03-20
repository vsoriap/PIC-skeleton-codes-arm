# PIC Skeleton Codes:  Vectorization

These codes illustrate how to use vectorization porting SSE2 to NEON.

For the 2D electrostatic:
  * no-vec = 35 nsec/particle/timestep
  * compiler vec = 18 nsec/particle/timestep
  * NEON = 12 nsec/particle/timestep

For the 2-1/2D electromagnetic:
  * no-vec = 100 nsec/particle/timestep
  * compiler vec = 60 nsec
  * NEON = 34 nsec

With SSE2 intrinsics one typically obtains about 3x speedup compared to no vectorization. Compiler vectorization achieves about 2x speedup. (All timings are on a 2.67GHz Intel Nehalem processor.)

### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


