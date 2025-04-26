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
