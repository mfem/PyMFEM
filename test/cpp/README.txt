## this is a place to put a test routine used for debugging.

## build command (nvcc)
> SITE=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
> nvcc cpp/test_innerproduct.cpp -ccbin mpic++ -I ${SITE}/mfem/external/par/include/ -I ${SITE}/mfem/external/include/  -L ${SITE}/mfem/external/par/lib -lmfem -Xlinker "-rpath,${SITE}/mfem/external/par/lib"

## build command (c++)
> SITE=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
> mpic++ cpp/test_innerproduct.cpp -I ${SITE}/mfem/external/par/include/ -I ${SITE}/mfem/external/include/  -L ${SITE}/mfem/external/par/lib -lmfem -Wl,-rpath,${SITE}/mfem/external/par/lib
