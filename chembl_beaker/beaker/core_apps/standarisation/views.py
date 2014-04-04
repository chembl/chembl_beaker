__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.core_apps.standarisation.impl import _break_bonds, _neutralise, _rules, _unsalt, _standardise
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/break_bonds/<ctab>', method=['OPTIONS', 'GET'], name="break_bonds")
def break_bonds(ctab):
    """
Break covalent bonds between oxygen or nitrogen atoms and Group I and II metal atoms.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01\_break\_bonds.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01_break_bonds.html)
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _break_bonds(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/break_bonds', method=['OPTIONS', 'POST'], name="break_bonds")
def break_bonds():
    """
Break covalent bonds between oxygen or nitrogen atoms and Group I and II metal atoms.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01\_break\_bonds.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01_break_bonds.html)
CTAB is either single molfile or SDF file.
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _break_bonds(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/neutralise/<ctab>', method=['OPTIONS', 'GET'], name="neutralise")
@app.route('/neutralise/<ctab>/<balance>', method=['OPTIONS', 'GET'], name="neutralise")
def neutralise(ctab, balance=False):
    """
Neutralise charges.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02\_neutralise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02_neutralise.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _neutralise(data, balance)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/neutralise', method=['OPTIONS', 'POST'], name="neutralise")
def neutralise():
    """
Neutralise charges.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02\_neutralise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02_neutralise.html)
CTAB is either single molfile or SDF file.
    """

    balance = int(request.forms.get('balance', False))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _neutralise(data, balance)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/unsalt/<ctab>', method=['OPTIONS', 'GET'], name="unsalt")
def unsalt(ctab):
    """
Remove counterions and solvate components.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04\_unsalt.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04_unsalt.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _unsalt(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/unsalt', method=['OPTIONS', 'POST'], name="unsalt")
def unsalt():
    """
Remove salt/solvates.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04\_unsalt.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04_unsalt.html)
CTAB is either single molfile or SDF file.
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _unsalt(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/rules/<ctab>', method=['OPTIONS', 'GET'], name="rules")
def rules(ctab):
    """
Apply structure-normalisation transformations.
List of rules and further documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03\_rules.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03_rules.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _rules(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/rules', method=['OPTIONS', 'POST'], name="rules")
def rules():
    """
Apply structure-normalisation transformations.
List of rules and further documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03\_rules.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03_rules.html)
CTAB is either single molfile or SDF file.
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _rules(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/standardise/<ctab>', method=['OPTIONS', 'GET'], name="standardise")
def standardise(ctab):
    """
Get standardised parent.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05\_standardise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05_standardise.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _standardise(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/standardise', method=['OPTIONS', 'POST'], name="standardise")
def standardise():
    """
Get standardised parent.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05\_standardise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05_standardise.html)
CTAB is either single molfile or SDF file.
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _standardise(data)

#-----------------------------------------------------------------------------------------------------------------------