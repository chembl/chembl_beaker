__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
from subprocess import PIPE, Popen
from rdkit import Chem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _getSDFString

#-----------------------------------------------------------------------------------------------------------------------

def _recogniseImage(img, osra):
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    p = Popen([osra, fpath], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a, err = p.communicate(input=img)
    os.remove(fpath)
    return filter(bool,a.split('\n'))

#-----------------------------------------------------------------------------------------------------------------------

def _image2ctab(img, osra):
    return _getSDFString(_apply(_recogniseImage(img, osra), Chem.MolFromSmiles))

#-----------------------------------------------------------------------------------------------------------------------
