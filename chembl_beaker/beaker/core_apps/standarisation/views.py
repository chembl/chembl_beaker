__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.utils.io import _parseFlag
from chembl_beaker.beaker.core_apps.standarisation.impl import _break_bonds, _neutralise, _rules, _unsalt, _standardise
import base64

#-----------------------------------------------------------------------------------------------------------------------

def breakBondsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    return _break_bonds(data, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/breakBonds/<ctab>', method=['OPTIONS', 'GET'], name="breakBonds")
def breakBonds(ctab):
    """
Break covalent bonds between oxygen or nitrogen atoms and Group I and II metal atoms.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01\_break\_bonds.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01_break_bonds.html)
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}breakBonds/$(cat breakBonds.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return breakBondsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/breakBonds', method=['OPTIONS', 'POST'], name="breakBonds")
def breakBonds():
    """
Break covalent bonds between oxygen or nitrogen atoms and Group I and II metal atoms.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01\_break\_bonds.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/01_break_bonds.html)
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @breakBonds.mol ${BEAKER_ROOT_URL}breakBonds
    curl -X POST -F "file=@breakBonds.mol" ${BEAKER_ROOT_URL}breakBonds
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return breakBondsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def neutraliseView(data, params):

    kwargs = dict()
    kwargs['balance'] = int(params.get('balance', False))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return _neutralise(data, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/neutralise/<ctab>', method=['OPTIONS', 'GET'], name="neutralise")
def neutralise(ctab):
    """
Neutralise charges.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02\_neutralise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02_neutralise.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}neutralise/$(cat neutralise.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return neutraliseView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/neutralise', method=['OPTIONS', 'POST'], name="neutralise")
def neutralise():
    """
Neutralise charges.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02\_neutralise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/02_neutralise.html)
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @neutralise.mol ${BEAKER_ROOT_URL}neutralise
    curl -X POST -F "file=@neutralise.mol" ${BEAKER_ROOT_URL}neutralise
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return neutraliseView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def unsaltView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    return _unsalt(data, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/unsalt/<ctab>', method=['OPTIONS', 'GET'], name="unsalt")
def unsalt(ctab):
    """
Remove counterions and solvate components.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04\_unsalt.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04_unsalt.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}unsalt/$(cat unsalt.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return unsaltView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/unsalt', method=['OPTIONS', 'POST'], name="unsalt")
def unsalt():
    """
Remove salt/solvates.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04\_unsalt.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/04_unsalt.html)
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @unsalt.mol ${BEAKER_ROOT_URL}unsalt
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}unsalt
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return unsaltView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def rulesView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    return _rules(data, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/rules/<ctab>', method=['OPTIONS', 'GET'], name="rules")
def rules(ctab):
    """
Apply structure-normalisation transformations.
List of rules and further documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03\_rules.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03_rules.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}rules/$(cat rules.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return rulesView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/rules', method=['OPTIONS', 'POST'], name="rules")
def rules():
    """
Apply structure-normalisation transformations.
List of rules and further documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03\_rules.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/03_rules.html)
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @rules.mol ${BEAKER_ROOT_URL}rules
    curl -X POST -F "file=@rules.mol" ${BEAKER_ROOT_URL}rules
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return rulesView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def standardiseView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    return _standardise(data, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/standardise/<ctab>', method=['OPTIONS', 'GET'], name="standardise")
def standardise(ctab):
    """
Get standardised parent.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05\_standardise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05_standardise.html)
CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}standardise/$(cat standardise.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return standardiseView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/standardise', method=['OPTIONS', 'POST'], name="standardise")
def standardise():
    """
Get standardised parent.
Examples and documentation: [https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05\_standardise.html](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/05_standardise.html)
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @standardise.mol ${BEAKER_ROOT_URL}standardise
    curl -X POST -F "file=@standardise.mol" ${BEAKER_ROOT_URL}standardise
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return standardiseView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------