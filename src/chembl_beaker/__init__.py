__author__ = 'mnowotka'

import rdkit

try:
    __version__ = '1.5.1'
except Exception as e:
    __version__ = 'development'


rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2019", "09", "2"]:
    raise ValueError("need an RDKit version >= 2019.09.2")