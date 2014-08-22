#!/usr/bin/env python

__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from bottle import run
from optparse import OptionParser
from chembl_beaker.beaker import app, config
from chembl_beaker.beaker.plugins.enableCors import EnableCors

#-----------------------------------------------------------------------------------------------------------------------


parser = OptionParser()
parser.add_option("-c", "--config", dest="config_path",
              help="path to config file", default="beaker.conf")
(options, args) = parser.parse_args()
conf_path = options.config_path

config.load_config(conf_path)

#-----------------------------------------------------------------------------------------------------------------------

if config.get('enable_cors', 'true') != 'false':
    app.install(EnableCors())

#-----------------------------------------------------------------------------------------------------------------------

def main():
    run(app=app, host=config.get('bottle_host', 'localhost'), port=config.get('bottle_port', '8080'),
                                debug=config.get('debug', True), server=config.get('server_middleware', 'tornado'))

if __name__ == "__main__":
    main()

else:
    application = app

#-----------------------------------------------------------------------------------------------------------------------
