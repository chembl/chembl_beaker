__author__ = 'mnowotka'

import bottle
from bottle import Bottle
import chembl_beaker
from chembl_beaker.beaker.utils import import_class
import re
import os
import sys
import json

# ----------------------------------------------------------------------------------------------------------------------

HTTP_CODES = bottle.HTTP_CODES.copy()
HTTP_CODES = dict((y, x) for x, y in HTTP_CODES.iteritems())
STATIC_ROOT = os.path.join(os.path.split(chembl_beaker.__file__)[0], 'static')
PARAM_REGEX = re.compile(r'<[^<>]+>')
DEFAULT_APPS = [
    "chembl_beaker.beaker",
    "chembl_beaker.beaker.core_apps.calculations",
    "chembl_beaker.beaker.core_apps.conversions",
    "chembl_beaker.beaker.core_apps.descriptors",
    "chembl_beaker.beaker.core_apps.fingerprints",
    "chembl_beaker.beaker.core_apps.marvin",
    "chembl_beaker.beaker.core_apps.mcs",
    "chembl_beaker.beaker.core_apps.osra",
    "chembl_beaker.beaker.core_apps.rasterImages",
    "chembl_beaker.beaker.core_apps.svgImages",
    "chembl_beaker.beaker.core_apps.jsonImages",
    "chembl_beaker.beaker.core_apps.ringInfo",
    "chembl_beaker.beaker.core_apps.standarisation",
    "chembl_beaker.beaker.core_apps.D3Coords",
    "chembl_beaker.beaker.core_apps.similarityMaps",
    "chembl_beaker.beaker.core_apps.autoDocs",
    ]

DEFAULT_PLUGINS = [
    'chembl_beaker.beaker.plugins.enableCors.EnableCors',
    'chembl_beaker.beaker.plugins.restrictions.Restrictions',
    'chembl_beaker.beaker.plugins.throttling.Throttling',
    'chembl_beaker.beaker.plugins.caching.Caching',
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
            print "Failed to load plugin %s because of error %s" % (plugin, e.message)
            continue

# ----------------------------------------------------------------------------------------------------------------------


def loadApps(apps):
    if not apps:
        apps = DEFAULT_APPS
    for module in apps:
        try:
            __import__(module + ".views")
        except Exception as e:
            print "Loading module %s failed because of error: %s" % (module, e.message)
            continue

# ----------------------------------------------------------------------------------------------------------------------

app = Bottle()
config = app.config

# ----------------------------------------------------------------------------------------------------------------------

if not getattr(config, 'load_config'):

        py = sys.version_info
        py3k = py >= (3, 0, 0)

        if py3k:
            from configparser import ConfigParser
        else:
            from ConfigParser import SafeConfigParser as ConfigParser

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


