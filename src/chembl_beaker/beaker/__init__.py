__author__ = 'mnowotka'

import bottle
from bottle import Bottle
import beaker
from beaker.utils import import_class
import re
import os
import json
import rdkit


# ----------------------------------------------------------------------------------------------------------------------

try:
    __version__ = '1.5.1'
except Exception as e:
    __version__ = 'development'


rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2019", "09", "2"]:
    raise ValueError("need an RDKit version >= 2019.09.2")

# ----------------------------------------------------------------------------------------------------------------------

HTTP_CODES = bottle.HTTP_CODES.copy()
HTTP_CODES = dict((y, x) for x, y in list(HTTP_CODES.items()))
STATIC_ROOT = "chembl_beaker/static"

PARAM_REGEX = re.compile(r'<[^<>]+>')
DEFAULT_APPS = [
    "beaker",
    "beaker.core_apps.conversions",
    "beaker.core_apps.descriptors",
    "beaker.core_apps.marvin",
    "beaker.core_apps.mcs",
    "beaker.core_apps.osra",
    "beaker.core_apps.svgImages",
    "beaker.core_apps.standarisation",
    "beaker.core_apps.D2Coords",
    "beaker.core_apps.autoDocs",
    "beaker.core_apps.structuralAlerts",
    "beaker.core_apps.calculations"
    ]

DEFAULT_PLUGINS = [
    'beaker.plugins.enableCors.EnableCors',
    'beaker.plugins.restrictions.Restrictions',
    'beaker.plugins.throttling.Throttling',
    'beaker.plugins.caching.Caching',
]

# ----------------------------------------------------------------------------------------------------------------------


def loadPlugins(app, plugins):
    if not plugins:
        plugins = DEFAULT_PLUGINS
    for plugin in plugins:
        try:
            plugin_class = import_class(plugin)
            app.install(plugin_class())
        except Exception as e:
            print("Failed to load plugin %s because of error %s" % (plugin, str(e)))
            continue

# ----------------------------------------------------------------------------------------------------------------------


def loadApps(apps):
    if not apps:
        apps = DEFAULT_APPS
    for module in apps:
        try:
            __import__(module + ".views")
        except Exception as e:
            print("Loading module %s failed because of error: %s" % (module, str(e)))
            continue

# ----------------------------------------------------------------------------------------------------------------------

app = Bottle()
app.mount('/chembl/api/utils', app)
config = app.config

# ----------------------------------------------------------------------------------------------------------------------

if not getattr(config, 'load_config'):

        from configparser import ConfigParser
        from bottle import ConfigDict

        def load_config(self, filename):
            ''' Load values from an *.ini style config file.

                If the config file contains sections, their names are used as
                namespaces for the values within. The two special sections
                ``DEFAULT`` and ``bottle`` refer to the root namespace (no prefix).
            '''
            conf = ConfigParser()
            conf.read(filename)
            for section in conf.sections():
                for key, value in conf.items(section):
                    if section not in ('DEFAULT', 'bottle'):
                        key = section + '.' + key
                    self[key] = value
            return self

        ConfigDict.load_config = load_config

# ----------------------------------------------------------------------------------------------------------------------
