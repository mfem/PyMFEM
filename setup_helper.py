import os
import subprocess
from subprocess import DEVNULL

def make_call_mp(command, target='', force_verbose=False, env=None, dry_run=False):
    '''
    call command
    '''
    print("calling ... " + " ".join(command))

    if dry_run:
        return
    kwargs = {'universal_newlines': True, 'env': env}
    if env is not None:
        env.update(os.environ)

    if not force_verbose:
        kwargs['stdout'] = DEVNULL
        kwargs['stderr'] = DEVNULL

    p = subprocess.Popen(command, **kwargs)
    p.communicate()
    if p.returncode != 0:
        if target == '':
            target = " ".join(command)
        print("Failed when calling command: " + target)
        raise subprocess.CalledProcessError(p.returncode,
                                            " ".join(command))

def get_numpy_inc():
    command = ["python", "-c", "import numpy;print(numpy.get_include())"]
    try:
        numpyinc = subprocess.run(command, capture_output=True).stdout.decode().strip()
    except subprocess.CalledProcessError:
        assert False, "can not check numpy include directory"
    except BaseException:
        assert False, "can not check numpy include directory"
    return numpyinc

def get_mpi4py_inc():
    command = ["python", "-c", "import mpi4py;print(mpi4py.get_include())"]
    try:
        mpi4pyinc = subprocess.run(command, capture_output=True).stdout.decode().strip()
    except subprocess.CalledProcessError:
        assert False, "can not check numpy include directory"
    except BaseException:
        assert False, "can not check numpy include directory"
    return mpi4pyinc
