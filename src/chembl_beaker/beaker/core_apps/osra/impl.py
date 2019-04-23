__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
from subprocess import PIPE, Popen

# ----------------------------------------------------------------------------------------------------------------------


def _recogniseImage(img, osra, frmt=None, jaggy=False, adaptive=False, unpaper=0):
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    arguments = [osra]
    if jaggy:
        arguments.append('-j')
    if adaptive:
        arguments.append('-i')
    if unpaper and str(unpaper).isdigit():
        arguments.extend(['-u', str(unpaper)])
    if frmt and frmt in ('can', 'smi', 'sdf'):
        arguments.extend(['-f', frmt])
    arguments.append(fpath)
    p = Popen(arguments, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    result, err = p.communicate(input=img)
    os.remove(fpath)
    return result

# ----------------------------------------------------------------------------------------------------------------------


def _image2ctab(img, osra, **kwargs):
    kwargs['frmt'] = 'sdf'
    return _recogniseImage(img, osra, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


def _image2smiles(img, osra, **kwargs):
    return 'SMILES Name \n' + '\n'.join(filter(bool, _recogniseImage(img, osra, **kwargs).split('\n')))

# ----------------------------------------------------------------------------------------------------------------------
