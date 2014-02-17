__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker import __version__ as version
from bottle import static_file
from chembl_beaker.beaker import STATIC_ROOT

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/status', name="status")
def status():
    """
Provides information about status of Beaker webservices. 200 means it's up and running.
Text message provides information about software version.
    """

    return "This is ChEMBL beaker, version %s" % version

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/static/<filename>')
def server_static(filename):
    return static_file(filename, root=STATIC_ROOT)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/favicon.ico')
def favicon():
    return static_file('favicon.ico', root=STATIC_ROOT)

#-----------------------------------------------------------------------------------------------------------------------
