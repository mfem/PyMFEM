## this is a place to put a test routine used for debugging.

## build command (nvcc)
> SITE=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
nvcc -ccbin mpic++  -I ${SITE}/mfem/external/par/include/ -I${SITE}/mfem/external/par/include/mfem -I ${SITE}/mfem/external/include/  -L ${SITE}/mfem/external/par/lib -lmfem -Xlinker "-rpath,${SITE}/mfem/external/par/lib" -expt-extended-lambda -O3 -DNDEBUG --generate-code=arch=compute_75,code=[compute_75,sm_75] -std=c++11 -MD  -x cu cpp/test_innerproduct.cpp

## build command (c++)
> SITE=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
> mpic++ cpp/test_innerproduct.cpp -I ${SITE}/mfem/external/par/include/ -I ${SITE}/mfem/external/include/  -L ${SITE}/mfem/external/par/lib -lmfem -Wl,-rpath,${SITE}/mfem/external/par/lib

## gdb with mpi
mpirun -np 2  xterm -e "gdb -ex "run" ./a.out"


## samples placed here

meshOptMWE.cpp  : minimum working example of mesh optimization (contribution by Ketan Mittel)





