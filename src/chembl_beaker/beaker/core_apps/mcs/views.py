__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from beaker import app
from beaker.core_apps.mcs.impl import _mcs
from beaker.utils.io import _parseFlag
from bottle import request

# ----------------------------------------------------------------------------------------------------------------------


def mcsView(data, params):

    kwargs = dict()
    kwargs['asSmiles'] = _parseFlag(params.get('asSmiles', '0'))
    kwargs['atomCompare'] = params.get('atomCompare', 'elements')
    kwargs['bondCompare'] = params.get('bondCompare', 'bondtypes')
    kwargs['ringMatchesRingOnly'] = _parseFlag(params.get('ringMatchesRingOnly', False))
    kwargs['completeRingsOnly'] = _parseFlag(params.get('completeRingsOnly', False))
    kwargs['threshold'] = params.get('threshold', None)
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    kwargs['isomericSmiles'] = _parseFlag(params.get('isomericSmiles', False))
    kwargs['canonical'] = _parseFlag(params.get('canonical', True))
    kwargs['kekuleSmiles'] = _parseFlag(params.get('kekuleSmiles', False))

    return _mcs(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/mcs', method=['OPTIONS', 'POST'], name="mcs")
def mcs():
    """
Returns Maximum Common Substructure of a set of compounds. This method accepts compounds in SDF format.
cURL examples:

    curl -X POST --data-binary @mcs.sdf ${BEAKER_ROOT_URL}mcs
    curl -X POST -F "file=@mcs.sdf" ${BEAKER_ROOT_URL}mcs
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return mcsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
