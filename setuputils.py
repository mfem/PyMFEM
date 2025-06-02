"""
Helper functions for setup.py
"""

import os
import sys
import configparser
from urllib import request
import itertools
import site
import subprocess
import multiprocessing
import ssl
import tarfile
from collections import namedtuple

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------

release = namedtuple('Release', ['version', 'hash', 'tarball'])
REPOS = dict(
    mfem = dict(
        url = "https://github.com/mfem/mfem.git",
        # version, hash, tarball
        releases = [
            release("4.7", "dc9128ef596e84daf1138aa3046b826bba9d259f", None),
            release("4.8", "a01719101027383954b69af1777dc828bf795d62", None),
        ]
    ),
    metis = dict(
        url = "https://github.com/KarypisLab/METIS",
        releases = [
            release("5.1.0", "94c03a6e2d1860128c2d0675cbbb86ad4f261256",
             "https://github.com/mfem/tpls/raw/gh-pages/metis-5.1.0.tar.gz"),
        ]
    ),
    gklib = dict(
        url = "https://github.com/KarypisLab/GKlib",
        releases = [
            release("5.1.1", "a7f8172703cf6e999dd0710eb279bba513da4fec",
             "https://github.com/KarypisLab/GKlib/archive/refs/tags/METIS-v5.1.1-DistDGL-0.5.tar.gz"),
        ]
    ),
    libceed = dict(
        url = "https://github.com/CEED/libCEED.git",
        releases = [
            release("0.12.0", None, "https://github.com/CEED/libCEED/archive/refs/tags/v0.12.0.tar.gz"),
        ]
    ),
    hypre = dict(
        url = None,
        releases = [
            release("2.28.0", None, "https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz"),
        ]
    ),
)

# ----------------------------------------------------------------------------------------
# Metadata
# ----------------------------------------------------------------------------------------

def read_mfem_tplflags(prefix):
    filename = os.path.join(prefix, 'share', 'mfem', 'config.mk')
    if not os.path.exists(filename):
        print("NOTE: " + filename + " does not exist.")
        print("returning empty string")
        return ""

    config = configparser.ConfigParser()
    with open(filename) as fp:
        config.read_file(itertools.chain(['[global]'], fp), source=filename)
    flags = dict(config.items('global'))['mfem_tplflags']
    return flags

# ----------------------------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------------------------
verbose = 1
dry_run = -1

def abspath(path):
    return os.path.abspath(os.path.expanduser(path))


def external_install_prefix(prefix, verbose=True):

    if hasattr(site, "getusersitepackages"):
        usersite = site.getusersitepackages()
    else:
        usersite = site.USER_SITE

    if verbose:
        print("running external_install_prefix with the following parameters")
        print("   sys.argv :", sys.argv)
        print("   sys.prefix :", sys.prefix)
        print("   usersite :", usersite)
        print("   prefix :", prefix)

    if '--user' in sys.argv:
        path = usersite
        if not os.path.exists(path):
            os.makedirs(path)
        path = os.path.join(path, 'mfem', 'external')
        return path

    else:
        # when prefix is given...let's borrow pip._internal to find the location ;D
        import pip._internal.locations
        path = pip._internal.locations.get_scheme(
            "mfem", prefix=prefix).purelib
        if not os.path.exists(path):
            os.makedirs(path)
        path = os.path.join(path, 'mfem', 'external')
        return path


def make_call(command, target='', force_verbose=False, env=None):
    '''
    call command
    '''
    print("calling ... " + " ".join(command))

    if dry_run:
        return
    kwargs = {'universal_newlines': True, 'env': env}
    if env is not None:
        env.update(os.environ)

    myverbose = verbose or force_verbose
    if not myverbose:
        kwargs['stdout'] = subprocess.DEVNULL
        kwargs['stderr'] = subprocess.DEVNULL

    p = subprocess.Popen(command, **kwargs)
    p.communicate()
    if p.returncode != 0:
        if target == '':
            target = " ".join(command)
        print("Failed when calling command: " + target)
        raise subprocess.CalledProcessError(p.returncode,
                                            " ".join(command))


def chdir(path):
    '''
    change directory to `path`; returns the previous directory
    '''
    pwd = os.getcwd()
    os.chdir(path)
    if verbose:
        print("Moving to a directory : " + path)
    return pwd


def remove_files(files):
    for f in files:
        if verbose:
            print("Removing : " + f)
        if dry_run:
            continue
        os.remove(f)


def make(target):
    '''
    make : add -j option automatically
    '''
    command = ['make', '-j', str(max((multiprocessing.cpu_count() - 1, 1)))]
    make_call(command, target=target, force_verbose=True)


def make_install(target, prefix=None):
    '''
    make install
    '''
    command = ['make', 'install']
    if prefix is not None:
        command.append('prefix='+prefix)
    make_call(command, target=target)


def download(xxx):
    '''
    download tar.gz from somewhere. xxx is name.
    url is given by repos above
    '''

    if os.path.exists(os.path.join(extdir, xxx)):
        print("Download " + xxx + " skipped. Use clean --all-exts if needed")
        return
    # Get the tarball for the latest release
    url = REPOS[xxx]["releases"][-1].tarball
    if url is None:
        raise RuntimeError(f"Could not find tarball URL for {xxx}")
    print("Downloading :", url)

    if use_unverifed_SSL:
        ssl._create_default_https_context = ssl._create_unverified_context

    ftpstream = request.urlopen(url)
    targz = tarfile.open(fileobj=ftpstream, mode="r|gz")
    targz.extractall(path=extdir)
    os.rename(os.path.join(extdir, targz.getnames()[0].split('/')[0]),
              os.path.join(extdir, xxx))


def gitclone(xxx, use_sha=False, branch='master'):
    cwd = os.getcwd()
    repo_xxx = os.path.join(extdir, xxx)
    if os.path.exists(repo_xxx):
        os.chdir(repo_xxx)
        command = ['git', 'checkout', branch]
        make_call(command)
        command = ['git', 'pull']
        make_call(command)
    else:
        repo = REPOS[xxx]["url"]
        if git_sshclone:
            repo = repo.replace("https://github.com/", "git@github.com:")

        os.chdir(extdir)
        command = ['git', 'clone', repo, xxx]
        make_call(command)

    if not dry_run:
        if not os.path.exists(repo_xxx):
            print(repo_xxx + " does not exist. Check if git clone worked")
        os.chdir(repo_xxx)

        if use_sha:
            sha = REPOS[xxx]["releases"][-1].hash
            command = ['git', 'checkout',  sha]
        else:
            command = ['git', 'checkout', branch]
        make_call(command)
    os.chdir(cwd)


def record_mfem_sha(mfem_source):
    pwd = chdir(mfem_source)
    command = ['git', 'rev-parse', 'HEAD']
    try:
        sha = subprocess.run(
            command, capture_output=True).stdout.decode().strip()
    except subprocess.CalledProcessError:
        print("subprocess failed to read sha...continuing w/o recording SHA")
        sha = None
    except BaseException:
        print("subprocess failed to read sha...continuing w/o recording SHA")
        sha = None

    chdir(pwd)

    sha_file = os.path.join('mfem', '__sha__.py')
    fid = open(sha_file, 'w')
    if sha is not None:
        fid.write('mfem = "' + sha + '"')
    fid.close()


def cmake(path, **kwargs):
    '''
    run cmake. must be called in the target directory
    '''
    command = ['cmake', path]
    for key, value in kwargs.items():
        command.append('-' + key + '=' + value)

    if osx_sysroot != '':
        command.append('-DCMAKE_OSX_SYSROOT=' + osx_sysroot)
    make_call(command)


def get_numpy_inc():
    command = ["python", "-c", "import numpy;print(numpy.get_include())"]
    try:
        numpyinc = subprocess.run(
            command, capture_output=True).stdout.decode().strip()
    except subprocess.CalledProcessError:
        assert False, "can not check numpy include directory"
    except BaseException:
        assert False, "can not check numpy include directory"
    return numpyinc


def get_mpi4py_inc():
    command = ["python", "-c", "import mpi4py;print(mpi4py.get_include())"]
    try:
        mpi4pyinc = subprocess.run(
            command, capture_output=True).stdout.decode().strip()
    except subprocess.CalledProcessError:
        assert False, "can not check numpy include directory"
    except BaseException:
        assert False, "can not check numpy include directory"
    return mpi4pyinc


def find_libpath_from_prefix(lib, prefix0):

    prefix0 = os.path.expanduser(prefix0)
    prefix0 = abspath(prefix0)

    soname = 'lib' + lib + dylibext
    aname = 'lib' + lib + '.a'

    path = os.path.join(prefix0, 'lib', soname)
    if os.path.exists(path):
        return path
    else:
        path = os.path.join(prefix0, 'lib64', soname)
        if os.path.exists(path):
            return path

    path = os.path.join(prefix0, 'lib', aname)
    if os.path.exists(path):
        return path
    else:
        path = os.path.join(prefix0, 'lib64', aname)
        if os.path.exists(path):
            return path
    print("Can not find library by find_libpath_from_prefix (continue)", lib, prefix0)

    return ''

# ----------------------------------------------------------------------------------------
# Build libraries
# ----------------------------------------------------------------------------------------

def cmake_make_hypre():
    '''
    build hypre
    '''
    if verbose:
        print("Building hypre")

    cmbuild = 'cmbuild'
    path = os.path.join(extdir, 'hypre', 'src', cmbuild)
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)

    pwd = chdir(path)

    cmake_opts = {'DBUILD_SHARED_LIBS': '1',
                  'DHYPRE_INSTALL_PREFIX': hypre_prefix,
                  'DHYPRE_ENABLE_SHARED': '1',
                  'DCMAKE_C_FLAGS': '-fPIC',
                  'DCMAKE_INSTALL_PREFIX': hypre_prefix,
                  'DCMAKE_INSTALL_NAME_DIR': os.path.join(hypre_prefix, 'lib'), }
    if verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    if enable_cuda and enable_cuda_hypre:
        # in this case, settitng CMAKE_C_COMPILER
        # causes "mpi.h" not found error. For now, letting CMAKE
        # to find MPI
        cmake_opts['DHYPRE_WITH_CUDA'] = '1'
        if cuda_arch != '':
            cmake_opts['DCMAKE_CUDA_ARCHITECTURES'] = cuda_arch
    else:
        cmake_opts['DCMAKE_C_COMPILER'] = mpicc_command

    cmake('..', **cmake_opts)

    make('hypre')
    make_install('hypre')

    os.chdir(pwd)


def make_metis_gklib(use_int64=False, use_real64=False):
    '''
    build GKlib/metis
    '''

    '''
    build/install GKlib
    '''
    if verbose:
        print("Building gklib")

    path = os.path.join(extdir, 'gklib')
    if not dry_run and not os.path.exists(path):
        assert False, "gklib is not downloaded"

    path = os.path.join(path, 'cmbuild')
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)
    pwd = chdir(path)

    cmake_opts = {'DBUILD_SHARED_LIBS': '1',
                  'DCMAKE_INSTALL_PREFIX': metis_prefix}
    if verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    cmake('..', **cmake_opts)
    make('gklib')
    make_install('gklib')
    os.chdir(pwd)

    '''
    build/install metis
    '''
    path = os.path.join(extdir, 'metis')
    if not dry_run and not os.path.exists(path):
        assert False, "metis is not downloaded"
    elif not os.path.exists(path):
        os.makedirs(path)
        os.makedirs(os.path.join(path, 'build'))

    pwd = chdir(path)

    gklibpath = os.path.dirname(find_libpath_from_prefix(
        'GKlib', metis_prefix))

    options = ['gklib_path='+metis_prefix]
    if use_int64:
        options.append('i64=1')

    if use_real64:
        options.append('r64=1')

    command = ['make', 'config', 'shared=1'] + options
    command = command + ['prefix=' + metis_prefix, 'cc=' + cc_command]
    make_call(command)

    chdir('build')
    cmake_opts = {'DGKLIB_PATH': metis_prefix,
                  'DSHARED': '1',
                  'DCMAKE_C_COMPILER': cc_command,
                  'DCMAKE_C_STANDARD_LIBRARIES': '-lGKlib',
                  'DCMAKE_INSTALL_RPATH': gklibpath,
                  'DCMAKE_BUILD_WITH_INSTALL_RPATH': '1',
                  'DCMAKE_INSTALL_PREFIX': metis_prefix}
    if verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    cmake('..', **cmake_opts)
    chdir(path)
    make('metis')
    make_install('metis')

    if platform == "darwin":
        command = ['install_name_tool',
                   '-id',
                   os.path.join(metis_prefix, 'lib', 'libGKlib.dylib'),
                   os.path.join(metis_prefix, 'lib', 'libGKlib.dylib'), ]
        make_call(command)
        command = ['install_name_tool',
                   '-id',
                   os.path.join(metis_prefix, 'lib', 'libmetis.dylib'),
                   os.path.join(metis_prefix, 'lib', 'libmetis.dylib'), ]
        make_call(command)
    os.chdir(pwd)


def make_metis(use_int64=False, use_real64=False):
    '''
    build metis
    '''
    if verbose:
        print("Building metis")

    path = os.path.join(extdir, 'metis')
    if not os.path.exists(path):
        assert False, "metis is not downloaded"

    pwd = chdir(path)

    if use_int64:
        pattern_int = "#define IDXTYPEWIDTH 32"
        replace_int = "#define IDXTYPEWIDTH 64"
    else:
        pattern_int = "#define IDXTYPEWIDTH 64"
        replace_int = "#define IDXTYPEWIDTH 32"
    with open("include/metis.h", "r") as metis_header_fid:
        metis_header_lines = metis_header_fid.readlines()
    with open("include/metis.h", "w") as metis_header_fid:
        for line in metis_header_lines:
            metis_header_fid.write(re.sub(pattern_int, replace_int, line))

    if use_real64:
        pattern_real = "#define REALTYPEWIDTH 32"
        replace_real = "#define REALTYPEWIDTH 64"
    else:
        pattern_real = "#define REALTYPEWIDTH 64"
        replace_real = "#define REALTYPEWIDTH 32"
    with open("include/metis.h", "r") as metis_header_fid:
        metis_header_lines = metis_header_fid.readlines()
    with open("include/metis.h", "w") as metis_header_fid:
        for line in metis_header_lines:
            metis_header_fid.write(re.sub(pattern_real, replace_real, line))

    command = ['make', 'config', 'shared=1',
               'prefix=' + metis_prefix,
               'cc=' + cc_command]
    make_call(command, env={'CMAKE_POLICY_VERSION_MINIMUM': '3.5'})
    make('metis')
    make_install('metis')

    if platform == "darwin":
        command = ['install_name_tool',
                   '-id',
                   os.path.join(metis_prefix, 'lib', 'libmetis.dylib'),
                   os.path.join(metis_prefix, 'lib', 'libmetis.dylib'), ]
        make_call(command)
    os.chdir(pwd)


def make_libceed(serial=False):
    if verbose:
        print("Building libceed")

    path = os.path.join(extdir, 'libceed')
    if not os.path.exists(path):
        assert False, "libceed is not downloaded"

    pwd = chdir(path)
    try:
        make_call(['make', 'clean'])
    except:
        pass

    if enable_cuda:
        command = ['make', 'configure', 'CUDA_DIR='+cuda_prefix]
        make_call(command)

    make('libceed')
    make_install('libceed', prefix=libceed_prefix)
    os.chdir(pwd)


def make_gslib(serial=False):
    if verbose:
        print("Building gslib")

    path = os.path.join(extdir, 'gslib')
    if not os.path.exists(path):
        assert False, "gslib is not downloaded"

    pwd = chdir(path)
    make_call(['make', 'clean'])
    if serial:
        command = ['make', 'CC=' + cc_command, 'MPI=0', 'CFLAGS=-fPIC']
        make_call(command)
        command = ['make', 'MPI=0', 'DESTDIR=' + gslibs_prefix]
        make_call(command)
    else:
        command = ['make', 'CC=' + mpicc_command, 'CFLAGS=-O2 -fPIC']
        make_call(command)
        command = ['make', 'DESTDIR=' + gslibp_prefix]
        make_call(command)
    os.chdir(pwd)


def cmake_make_mfem(serial=True):
    '''
    build MFEM
    '''
    cmbuild = 'cmbuild_ser' if serial else 'cmbuild_par'
    path = os.path.join(extdir, 'mfem', cmbuild)
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)

    ldflags = os.getenv('LDFLAGS') if os.getenv('LDFLAGS') is not None else ''
    metisflags = ''
    hypreflags = ''

    rpaths = []

    def add_rpath(p):
        if not p in rpaths:
            rpaths.append(p)

    cmake_opts = {'DBUILD_SHARED_LIBS': '1',
                  'DMFEM_ENABLE_EXAMPLES': '1',
                  'DMFEM_ENABLE_MINIAPPS': '0',
                  'DCMAKE_SHARED_LINKER_FLAGS': ldflags,
                  'DMFEM_USE_ZLIB': '1',
                  'DCMAKE_CXX_FLAGS': cxx11_flag,
                  'DCMAKE_BUILD_WITH_INSTALL_RPATH': '1'}

    if mfem_debug:
        cmake_opts['DMFEM_DEBUG'] = 'YES'

    if mfem_build_miniapps:
        cmake_opts['DMFEM_ENABLE_MINIAPPS'] = '1'

    if verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    if serial:
        cmake_opts['DCMAKE_CXX_COMPILER'] = cxx_command
        cmake_opts['DMFEM_USE_EXCEPTIONS'] = '1'
        cmake_opts['DCMAKE_INSTALL_PREFIX'] = mfems_prefix

        add_rpath(os.path.join(mfems_prefix, 'lib'))
        if enable_suitesparse:
            enable_metis = True
        else:
            enable_metis = False
    else:
        cmake_opts['DCMAKE_CXX_COMPILER'] = mpicxx_command
        cmake_opts['DMFEM_USE_EXCEPTIONS'] = '0'
        cmake_opts['DCMAKE_INSTALL_PREFIX'] = mfemp_prefix
        cmake_opts['DMFEM_USE_MPI'] = '1'
        cmake_opts['DHYPRE_DIR'] = hypre_prefix
        cmake_opts['DHYPRE_INCLUDE_DIRS'] = os.path.join(
            hypre_prefix, "include")

        add_rpath(os.path.join(mfemp_prefix, 'lib'))

        hyprelibpath = os.path.dirname(
            find_libpath_from_prefix(
                'HYPRE', hypre_prefix))

        add_rpath(hyprelibpath)

        hypreflags = "-L" + hyprelibpath + " -lHYPRE "

        if enable_strumpack:
            cmake_opts['DMFEM_USE_STRUMPACK'] = '1'
            cmake_opts['DSTRUMPACK_DIR'] = strumpack_prefix
            libpath = os.path.dirname(
                find_libpath_from_prefix("STRUMPACK", strumpack_prefix))
            add_rpath(libpath)
        if enable_pumi:
            cmake_opts['DMFEM_USE_PUMI'] = '1'
            cmake_opts['DPUMI_DIR'] = pumi_prefix
            libpath = os.path.dirname(
                find_libpath_from_prefix("pumi", strumpack_prefix))
            add_rpath(libpath)
        enable_metis = True

    if enable_metis:
        cmake_opts['DMFEM_USE_METIS_5'] = '1'
        cmake_opts['DMETIS_DIR'] = metis_prefix
        cmake_opts['DMETIS_INCLUDE_DIRS'] = os.path.join(
            metis_prefix, "include")
        metislibpath = os.path.dirname(
            find_libpath_from_prefix(
                'metis', metis_prefix))
        add_rpath(metislibpath)

        if use_metis_gklib:
            metisflags = "-L" + metislibpath + " -lmetis -lGKlib "
        else:
            metisflags = "-L" + metislibpath + " -lmetis "

    if ldflags != '':
        cmake_opts['DCMAKE_SHARED_LINKER_FLAGS'] = ldflags
        cmake_opts['DCMAKE_EXE_LINKER_FLAGS'] = ldflags

    if metisflags != '':
        cmake_opts['DMETIS_LIBRARIES'] = metisflags
    if hypreflags != '':
        cmake_opts['DHYPRE_LIBRARIES'] = hypreflags

    if enable_cuda:
        cmake_opts['DMFEM_USE_CUDA'] = '1'
        if cuda_arch != '':
            cmake_opts['DCMAKE_CUDA_ARCHITECTURES'] = cuda_arch

    if enable_libceed:
        cmake_opts['DMFEM_USE_CEED'] = '1'
        cmake_opts['DCEED_DIR'] = libceed_prefix
        libpath = os.path.dirname(
            find_libpath_from_prefix("ceed", libceed_prefix))
        add_rpath(libpath)

    if enable_gslib:
        if serial:
            cmake_opts['DMFEM_USE_GSLIB'] = '1'
            cmake_opts['DGSLIB_DIR'] = gslibs_prefix
        else:
            cmake_opts['DMFEM_USE_GSLIB'] = '1'
            cmake_opts['DGSLIB_DIR'] = gslibp_prefix

    if enable_suitesparse:
        cmake_opts['DMFEM_USE_SUITESPARSE'] = '1'
        if suitesparse_prefix != '':
            cmake_opts['DSuiteSparse_DIR'] = suitesparse_prefix

    if enable_lapack:
        cmake_opts['DMFEM_USE_LAPACK'] = '1'
    if blas_libraries != "":
        cmake_opts['DBLAS_LIBRARIES'] = blas_libraries
    if lapack_libraries != "":
        cmake_opts['DLAPACK_LIBRARIES'] = lapack_libraries

    cmake_opts['DCMAKE_INSTALL_RPATH'] = ":".join(rpaths)

    pwd = chdir(path)
    cmake('..', **cmake_opts)

    txt = 'serial' if serial else 'parallel'

    make('mfem_' + txt)
    make_install('mfem_' + txt)

    from shutil import copytree, rmtree

    print("copying mesh data for testing", "../data",
          cmake_opts['DCMAKE_INSTALL_PREFIX'])
    path = os.path.join(cmake_opts['DCMAKE_INSTALL_PREFIX'], "data")
    if os.path.exists(path):
        rmtree(path)
    copytree("../data", path)

    if do_bdist_wheel:
        ex_dir = os.path.join(cmake_opts['DCMAKE_INSTALL_PREFIX'], "examples")
        for x in os.listdir(ex_dir):
            path = os.path.join(ex_dir, x)
            command = ['chrpath', '-r', "$ORIGIN/../lib", path]
            make_call(command, force_verbose=True)

    os.chdir(pwd)
