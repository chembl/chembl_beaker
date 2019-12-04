__author__ = 'mnowotka'

import rdkit

try:
    __version__ = __import__('pkg_resources').get_distribution('chembl_beaker').version
except Exception as e:
    __version__ = 'development'


rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2019", "09", "2"]:
    raise ValueError("need an RDKit version >= 2019.09.2")