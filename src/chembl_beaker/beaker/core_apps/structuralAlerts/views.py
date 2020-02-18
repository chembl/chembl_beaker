__author__ = 'efelix'

#-----------------------------------------------------------------------------------------------------------------------

from beaker import app
from bottle import request, response
from beaker.utils.io import _parseFlag
from beaker.core_apps.structuralAlerts.impl import _get_alerts

#-----------------------------------------------------------------------------------------------------------------------


def getAlertsView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _get_alerts(data, **kwargs)


#-----------------------------------------------------------------------------------------------------------------------

@app.route('/structuralAlerts', method=['OPTIONS', 'POST'], name='structuralAlerts')
def get_alerts():
    """
Get structural alerts for a compound.
Examples and documentation: []()
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}structuralAlerts
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return getAlertsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------
