#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        USE_FFTW3\n"
"        ALLOW_IMBALANCED_GASPARTICLELOAD\n"
"        GALSF\n"
"        GRAVITY_ANALYTIC\n"
"        ADM\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        COOLING\n"
"        HYDRO_MESHLESS_FINITE_MASS\n"
"\n");
}
