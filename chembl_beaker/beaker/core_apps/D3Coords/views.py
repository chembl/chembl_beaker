__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.core_apps.D3Coords.impl import _ctab23D, _smiles23D
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab23D/<ctab>', method=['OPTIONS', 'GET'], name="ctab23D")
@app.route('/ctab23D/<ctab>/<multi>', method=['OPTIONS', 'GET'], name="ctab23D")
def ctab23D(ctab, multi=False):
    """
Generate 3D coordinates from molfile using Universal Force Field.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _ctab23D(data, int(multi), False)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab23D', method=['OPTIONS', 'POST'], name="ctab23D")
def ctab23D():
    """
Generate 3D coordinates from molfile using Universal Force Field.
CTAB is either single molfile or SDF file.
    """

    multi = int(request.forms.get('multi', False))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _ctab23D(data, multi, False)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles23D/<ctab>', method=['OPTIONS', 'GET'], name="smiles23D")
@app.route('/smiles23D/<ctab>/<multi>', method=['OPTIONS', 'GET'], name="smiles23D")
def smiles23D(ctab, multi=False):
    """
Generate 3D coordinates from SMILES using Universal Force Field.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _smiles23D(data, int(multi), False)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles23D', method=['OPTIONS', 'POST'], name="smiles23D")
def smiles23D():
    """
Generate 3D coordinates from SMILES using Universal Force Field.
CTAB is either single molfile or SDF file.
    """

    multi = int(request.forms.get('multi', False))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _smiles23D(data, multi, False)

#-----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import AllChem
if hasattr(AllChem, 'MMFFOptimizeMolecule'):

    @app.route('/MMFFctab23D/<ctab>', method=['OPTIONS', 'GET'], name="MMFFctab23D")
    @app.route('/MMFFctab23D/<ctab>/<multi>', method=['OPTIONS', 'GET'], name="MMFFctab23D")
    def MMFFctab23D(ctab, multi=False):
        """
    Generate 3D coordinates from molfile using Merck Molecular Force Field.
    CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
        """

        data = base64.urlsafe_b64decode(ctab)
        return _ctab23D(data, int(multi), True)

#-----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFctab23D', method=['OPTIONS', 'POST'], name="MMFFctab23D")
    def MMFFctab23D():
        """
    Generate 3D coordinates from molfile using Merck Molecular Force Field.
    CTAB is either single molfile or SDF file.
        """

        multi = int(request.forms.get('multi', False))
        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return _ctab23D(data, multi, True)

#-----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFsmiles23D/<ctab>', method=['OPTIONS', 'GET'], name="MMFFsmiles23D")
    @app.route('/MMFFsmiles23D/<ctab>/<multi>', method=['OPTIONS', 'GET'], name="MMFFsmiles23D")
    def MMFFsmiles23D(ctab, multi=False):
        """
    Generate 3D coordinates from SMILES using Merck Molecular Force Field.
    CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
        """

        data = base64.urlsafe_b64decode(ctab)
        return _smiles23D(data, int(multi), True)

#-----------------------------------------------------------------------------------------------------------------------

    @app.route('/MMFFsmiles23D', method=['OPTIONS', 'POST'], name="MMFFsmiles23D")
    def MMFFsmiles23D():
        """
    Generate 3D coordinates from SMILES using Merck Molecular Force Field.
    CTAB is either single molfile or SDF file.
        """

        multi = int(request.forms.get('multi', False))
        data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
        return _smiles23D(data, multi, True)

#-----------------------------------------------------------------------------------------------------------------------