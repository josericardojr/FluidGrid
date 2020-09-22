#ifndef __K_SimulationLib_CUH__

#define __K_SimulationLib_CUH__

#define USE_TEX 1

#if USE_TEX
#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCH(t, i) t[i]
#endif


#endif
