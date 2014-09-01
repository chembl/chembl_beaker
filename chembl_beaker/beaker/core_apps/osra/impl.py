__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
from subprocess import PIPE, Popen

#-----------------------------------------------------------------------------------------------------------------------

def _recogniseImage(img, osra, frmt='smi'):
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    arguments = [osra, '-ij', '-f', frmt, fpath]
    print ' '.join(arguments)
    p = Popen(arguments, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    result, err = p.communicate(input=img)
    os.remove(fpath)
    return result

#-----------------------------------------------------------------------------------------------------------------------

def _image2ctab(img, osra):
    return _recogniseImage(img, osra, 'sdf')

#-----------------------------------------------------------------------------------------------------------------------

def _image2smiles(img, osra):
    return 'SMILES Name \n' + '\n'.join(filter(bool, _recogniseImage(img, osra).split('\n')))

#-----------------------------------------------------------------------------------------------------------------------
