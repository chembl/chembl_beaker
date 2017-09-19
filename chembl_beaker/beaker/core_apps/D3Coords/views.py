__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import AllChem
from chembl_beaker.beaker import app
from bottle import request
from chembl_beaker.beaker.core_apps.D3Coords.impl import _ctab23D, _smiles23D
from chembl_beaker.beaker.utils.io import _parseFlag
import base64

# ----------------------------------------------------------------------------------------------------------------------


def ctab23DView(data, params):
    kwargs = dict()
    kwargs['multi'] = int(params.get('multi', False))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['mmff'] = _parseFlag(params.get('mmff', False))

    return _ctab23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab23D/<ctab>', method=['OPTIONS', 'GET'], name="ctab23D")
def ctab23D(ctab):
    """
Generate 3D coordinates from molfile using Universal Force Field.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cuRL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab23D/$(cat no_coords.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab23D', method=['OPTIONS', 'POST'], name="ctab23D")
def ctab23D():
    """
Generate 3D coordinates from molfile using Universal Force Field.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @no_coords.mol ${BEAKER_ROOT_URL}ctab23D
    curl -X POST -F "file=@no_coords.mol" ${BEAKER_ROOT_URL}ctab23D
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles23DView(data, params):
    kwargs = dict()
    kwargs['multi'] = int(params.get('multi', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['mmff'] = _parseFlag(params.get('mmff', False))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles23D/<ctab>', method=['OPTIONS', 'GET'], name="smiles23D")
def smiles23D(ctab):
    """
Generate 3D coordinates from SMILES using Universal Force Field.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles23D/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles23D/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return smiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles23D', method=['OPTIONS', 'POST'], name="smiles23D")
def smiles23D():
    """
Generate 3D coordinates from SMILES using Universal Force Field.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles23D
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles23D
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles23D
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles23D
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

if hasattr(AllChem, 'MMFFOptimizeMolecule'):

    # TODO: this one is redundant (but you have to remember about conditional availability)
    def MMFFctab23DView(data, params):
        kwargs = dict()
        kwargs['multi'] = int(params.get('multi', False))
        kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
        kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
        kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
        kwargs['mmff'] = True

        return _ctab23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFctab23D/<ctab>', method=['OPTIONS', 'GET'], name="MMFFctab23D")
    def MMFFctab23D(ctab):
        """
Generate 3D coordinates from molfile using Merck Molecular Force Field.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}MMFFctab23D/$(cat no_coords.mol | base64 -w 0 | tr "+/" "-_")
        """

        data = base64.urlsafe_b64decode(ctab)
        return MMFFctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFctab23D', method=['OPTIONS', 'POST'], name="MMFFctab23D")
    def MMFFctab23D():
        """
Generate 3D coordinates from molfile using Merck Molecular Force Field.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @no_coords.mol ${BEAKER_ROOT_URL}MMFFctab23D
    curl -X POST -F "file=@no_coords.mol" ${BEAKER_ROOT_URL}MMFFctab23D    
        """

        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return MMFFctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    def MMFFsmiles23DView(data, params):
        kwargs = dict()
        kwargs['multi'] = int(params.get('multi', False))
        kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
        kwargs['delimiter'] = params.get('delimiter', ' ')
        kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
        kwargs['nameColumn'] = int(params.get('nameColumn', 1))
        kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
        kwargs['mmff'] = True

        if params.get('titleLine') is None and not data.startswith('SMILES Name'):
            kwargs['titleLine'] = False
        else:
            kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

        return _smiles23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFsmiles23D/<ctab>', method=['OPTIONS', 'GET'], name="MMFFsmiles23D")
    def MMFFsmiles23D(ctab):
        """
Generate 3D coordinates from SMILES using Merck Molecular Force Field.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}MMFFsmiles23D/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}MMFFsmiles23D/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")    
        """

        data = base64.urlsafe_b64decode(ctab)
        return MMFFsmiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFsmiles23D', method=['OPTIONS', 'POST'], name="MMFFsmiles23D")
    def MMFFsmiles23D():
        """
Generate 3D coordinates from SMILES using Merck Molecular Force Field.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}MMFFsmiles23D
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}MMFFsmiles23D
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}MMFFsmiles23D
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}MMFFsmiles23D
        """

        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return MMFFsmiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

if hasattr(AllChem, 'ETKDG'):

    # TODO: this one is redundant (but you have to remember about conditional availability)
    def ETKDGctab23DView(data, params):
        kwargs = dict()
        kwargs['multi'] = int(params.get('multi', False))
        kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
        kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
        kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
        kwargs['mmff'] = True

        return _ctab23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('ETKDGctab23D/<ctab>', method=['OPTIONS', 'GET'], name="ETKDGctab23D")
    def ETKDGctab23D(ctab):
        """
Generate 3D coordinates from molfile using Experimental-Torsion "basic Knowledge" Distance Geometry.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ETKDGctab23D/$(cat no_coords.mol | base64 -w 0 | tr "+/" "-_")
        """

        data = base64.urlsafe_b64decode(ctab)
        return ETKDGctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/ETKDGctab23D', method=['OPTIONS', 'POST'], name="ETKDGctab23D")
    def ETKDGctab23D():
        """
Generate 3D coordinates from molfile using Experimental-Torsion "basic Knowledge" Distance Geometry.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @no_coords.mol ${BEAKER_ROOT_URL}ETKDGctab23D
    curl -X POST -F "file=@no_coords.mol" ${BEAKER_ROOT_URL}ETKDGctab23D    
        """

        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return ETKDGctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    def ETKDGsmiles23DView(data, params):
        kwargs = dict()
        kwargs['multi'] = int(params.get('multi', False))
        kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
        kwargs['delimiter'] = params.get('delimiter', ' ')
        kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
        kwargs['nameColumn'] = int(params.get('nameColumn', 1))
        kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
        kwargs['mmff'] = True

        if params.get('titleLine') is None and not data.startswith('SMILES Name'):
            kwargs['titleLine'] = False
        else:
            kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

        return _smiles23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/ETKDGsmiles23D/<ctab>', method=['OPTIONS', 'GET'], name="ETKDGsmiles23D")
    def ETKDGsmiles23D(ctab):
        """
Generate 3D coordinates from SMILES using Experimental-Torsion "basic Knowledge" Distance Geometry.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ETKDGsmiles23D/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}ETKDGsmiles23D/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")    
        """

        data = base64.urlsafe_b64decode(ctab)
        return ETKDGsmiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/ETKDGsmiles23D', method=['OPTIONS', 'POST'], name="ETKDGsmiles23D")
    def ETKDGsmiles23D():
        """
Generate 3D coordinates from SMILES using Experimental-Torsion "basic Knowledge" Distance Geometry.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}ETKDGsmiles23D
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}ETKDGsmiles23D
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}ETKDGsmiles23D
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}ETKDGsmiles23D
        """

        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return ETKDGsmiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
if hasattr(AllChem, 'KDG'):

    # TODO: this one is redundant (but you have to remember about conditional availability)
    def KDGctab23DView(data, params):
        kwargs = dict()
        kwargs['multi'] = int(params.get('multi', False))
        kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
        kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
        kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
        kwargs['mmff'] = True

        return _ctab23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/KDGctab23D/<ctab>', method=['OPTIONS', 'GET'], name="KDGctab23D")
    def KDGctab23D(ctab):
        """
Generate 3D coordinates from molfile using "basic Knowledge" Distance Geometry.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}KDGctab23D/$(cat no_coords.mol | base64 -w 0 | tr "+/" "-_")
        """

        data = base64.urlsafe_b64decode(ctab)
        return KDGctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/KDGctab23D', method=['OPTIONS', 'POST'], name="KDGctab23D")
    def KDGctab23D():
        """
Generate 3D coordinates from molfile using "basic Knowledge" Distance Geometry.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @no_coords.mol ${BEAKER_ROOT_URL}KDGctab23D
    curl -X POST -F "file=@no_coords.mol" ${BEAKER_ROOT_URL}KDGctab23D    
        """

        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return KDGctab23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    def KDGsmiles23DView(data, params):
        kwargs = dict()
        kwargs['multi'] = int(params.get('multi', False))
        kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
        kwargs['delimiter'] = params.get('delimiter', ' ')
        kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
        kwargs['nameColumn'] = int(params.get('nameColumn', 1))
        kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
        kwargs['mmff'] = True

        if params.get('titleLine') is None and not data.startswith('SMILES Name'):
            kwargs['titleLine'] = False
        else:
            kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

        return _smiles23D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/KDGsmiles23D/<ctab>', method=['OPTIONS', 'GET'], name="KDGsmiles23D")
    def KDGsmiles23D(ctab):
        """
Generate 3D coordinates from SMILES using "basic Knowledge" Distance Geometry.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}KDGsmiles23D/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}KDGsmiles23D/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")    
        """

        data = base64.urlsafe_b64decode(ctab)
        return KDGsmiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

    @app.route('/KDGsmiles23D', method=['OPTIONS', 'POST'], name="KDGsmiles23D")
    def KDGsmiles23D():
        """
Generate 3D coordinates from SMILES using "basic Knowledge" Distance Geometry.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}KDGsmiles23D
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}KDGsmiles23D
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}KDGsmiles23D
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}KDGsmiles23D
        """

        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return KDGsmiles23DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
