from beaker import app
from bottle import request
from beaker.core_apps.calculations.impl import _removeHs


# ----------------------------------------------------------------------------------------------------------------------


def removeHsView(data, params):
    return _removeHs(data)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/removeHs', method=['OPTIONS', 'POST'], name="removeHs")
def removeHs():
    """
Removes any hydrogens from the graph of a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @removeHs.mol ${BEAKER_ROOT_URL}removeHs
    curl -X POST -F "file=@removeHs.mol" ${BEAKER_ROOT_URL}removeHs
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return removeHsView(data, request.params)
