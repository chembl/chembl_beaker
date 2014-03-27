__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request
from chembl_beaker.beaker.core_apps.conversions.impl import _ctab2smiles, _smiles2ctab, _inchi2ctab, _ctab2inchi, _inchi2inchiKey
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2smiles/<ctab>', method=['OPTIONS', 'GET'], name="ctab2smiles")
def ctab2smiles(ctab):
    """
Converts CTAB to SMILES format. CTAB is urlsafe_base64 encoded string containing single molfile or concatenation
of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _ctab2smiles(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2smiles', method=['OPTIONS', 'POST'], name="ctab2smiles")
def ctab2smiles():
    """
Converts CTAB to SMILES format. CTAB is either single molfile or SDF file.
    """

    data=request.body.read()
    return _ctab2smiles(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2ctab/<smiles>', method=['OPTIONS', 'GET'], name="smiles2ctab")
def smiles2ctab(smiles):
    """
Converts SMILES to CTAB. This method accepts urlsafe_base64 encoded string containing single or multiple SMILES
optionally containing header line, specific to *.smi format.
    """

    data = base64.urlsafe_b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    return _smiles2ctab(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2ctab', method=['OPTIONS', 'POST'], name="smiles2ctab")
def smiles2ctab():
    """
Converts SMILES to CTAB. This method accepts single or multiple SMILES or *.smi file.
    """

    data = request.body.read()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    return _smiles2ctab(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/inchi2ctab/<inchi>', method=['OPTIONS', 'GET'], name="inchi2ctab")
def inchi2ctab(inchi):
    """
Converts InChi to CTAB. This method accepts urlsafe_base64 encoded string containing one or multiple InChis.
    """

    inchis = base64.urlsafe_b64decode(inchi)
    return _inchi2ctab(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/inchi2ctab', method=['OPTIONS', 'POST'], name="inchi2ctab")
def inchi2ctab():
    """
Converts InChi to CTAB. This method accepts one or multiple InChis.
    """

    inchis = request.body.read()
    return _inchi2ctab(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2inchi/<ctab>', method=['OPTIONS', 'GET'], name="ctab2inchi")
def ctab2inchi(ctab):
    """
Converts CTAB to InChis. CTAB is urlsafe_base64 encoded string containing single molfile or concatenation
of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _ctab2inchi(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2inchi', method=['OPTIONS', 'POST'], name="ctab2inchi")
def ctab2inchi():
    """
Converts CTAB to InChis. CTAB is either single molfile or SDF file.
    """

    data=request.body.read()
    return _ctab2inchi(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/inchi2inchiKey/<inchi>', method=['OPTIONS', 'GET'], name="inchi2inchiKey")
def inchi2inchiKey(inchi):
    """
Converts InChis to InChiKeys. This method accepts urlsafe_base64 encoded string containing one or multiple InChis.
    """

    inchis = base64.urlsafe_b64decode(inchi)
    return _inchi2inchiKey(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/inchi2inchiKey', method=['OPTIONS', 'POST'], name="inchi2inchiKey")
def inchi2inchiKey():
    """
Converts InChis to InChiKeys. This method accepts one or multiple InChis.
    """

    inchis = request.body.read()
    return _inchi2inchiKey(inchis)

#-----------------------------------------------------------------------------------------------------------------------
