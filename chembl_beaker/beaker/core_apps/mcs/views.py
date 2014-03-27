__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.mcs.impl import _mcs
from bottle import request
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/mcs/<ctab>', method=['OPTIONS', 'GET'], name="mcs")
def mcs(ctab):
    """
Returns Maximum Common Substructure of a set of compounds. CTAB is urlsafe_base64 encoded string containing molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _mcs(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/mcs', method=['OPTIONS', 'POST'], name="mcs")
def mcs():
    """
Returns Maximum Common Substructure of a set of compounds. This method accepts compounds in SDF format.
    """

    data = request.body.read()
    return _mcs(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------
