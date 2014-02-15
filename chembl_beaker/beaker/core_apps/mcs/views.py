__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.mcs.impl import _mcs
from bottle import request
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/mcs/<ctab>', name="mcs")
def mcs(ctab):
    """
Returns Maximum Common Substructure of a compound. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _mcs(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/mcs', method=['OPTIONS', 'POST'], name="mcs")
def mcs():
    """
Returns Maximum Common Substructure of a compound. CTAB is either single molfile or SDF file.
    """

    data = request.body.getvalue()
    return _mcs(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------
