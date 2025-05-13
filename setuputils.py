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


# meta data

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

# utilities

def abspath(path):
    return os.path.abspath(os.path.expanduser(path))


def external_install_prefix(verbose=True):

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
    url = repo_releases[xxx]
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

        # print("Deleting the existing " + xxx)
        # shutil.rmtree(repo_xxx)
    else:
        repo = repos[xxx]
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
            sha = repos_sha[xxx]
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
