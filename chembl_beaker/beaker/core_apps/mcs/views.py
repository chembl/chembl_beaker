__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.mcs.impl import _mcs
from chembl_beaker.beaker.utils.io import _parseFlag
from bottle import request
import base64

# ----------------------------------------------------------------------------------------------------------------------


def mcsView(data, params):

    kwargs = dict()
    kwargs['asSmiles'] = _parseFlag(params.get('asSmiles', '0'))
    kwargs['atomCompare'] = params.get('atomCompare', 'elements')
    kwargs['bondCompare'] = params.get('bondCompare', 'bondtypes')
    kwargs['ringMatchesRingOnly'] = _parseFlag(params.get('ringMatchesRingOnly', False))
    kwargs['completeRingsOnly'] = _parseFlag(params.get('completeRingsOnly', False))
    kwargs['threshold'] = params.get('threshold', None)
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['isomericSmiles'] = _parseFlag(params.get('isomericSmiles', False))
    kwargs['canonical'] = _parseFlag(params.get('canonical', True))
    kwargs['kekuleSmiles'] = _parseFlag(params.get('kekuleSmiles', False))

    return _mcs(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/mcs/<ctab>', method=['OPTIONS', 'GET'], name="mcs")
def mcs(ctab):
    """
Returns Maximum Common Substructure of a set of compounds. CTAB is urlsafe_base64 encoded string containing molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}mcs/$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return mcsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/mcs', method=['OPTIONS', 'POST'], name="mcs")
def mcs():
    """
Returns Maximum Common Substructure of a set of compounds. This method accepts compounds in SDF format.
cURL examples:

    curl -X POST --data-binary @mcs.sdf ${BEAKER_ROOT_URL}mcs
    curl -X POST -F "file=@mcs.sdf" ${BEAKER_ROOT_URL}mcs
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return mcsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
