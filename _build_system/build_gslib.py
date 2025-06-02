def make_gslib(serial=False):
    if bglb.verbose:
        print("Building gslib")

    path = os.path.join(extdir, 'gslib')
    if not os.path.exists(path):
        assert False, "gslib is not downloaded"

    pwd = chdir(path)
    make_call(['make', 'clean'])
    if serial:
        command = ['make', 'CC=' + cc_command, 'MPI=0', 'CFLAGS=-fPIC']
        make_call(command)
        command = ['make', 'MPI=0', 'DESTDIR=' + bglb.gslibs_prefix]
        make_call(command)
    else:
        command = ['make', 'CC=' + mpicc_command, 'CFLAGS=-O2 -fPIC']
        make_call(command)
        command = ['make', 'DESTDIR=' + bglb.gslibp_prefix]
        make_call(command)
    os.chdir(pwd)
