% -g option of mex is necessary when you debug the code in C files
% PS: If you have any error during compilation please try again
mex computeKernelSizeFeature.c
mex -g computeKernelSizeFeature.c

mex computeKernelSizeShape.c
mex -g computeKernelSizeShape.c

mex Evolve_data.c
mex -g Evolve_data.c

mex Evolve_ver12L2.c
mex -g Evolve_ver12L2.c

mex FMM.c
mex -g FMM.c