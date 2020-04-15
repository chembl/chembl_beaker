chembl_beaker
======

.. image:: https://img.shields.io/pypi/status/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: Development Status

.. image:: https://img.shields.io/pypi/l/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: License
    
What is Beaker?
--------

This is chembl_beaker package developed at `ChEMBL <https://www.ebi.ac.uk/chembl/>`_ group, `EMBL-EBI <https://www.ebi.ac.uk/>`_, Cambridge, UK.

This is wrapper for `RDKit <http://www.rdkit.org/>`_ and `OSRA <http://cactus.nci.nih.gov/osra/>`_, which exposes following methods:

 * `Format convertion <https://github.com/chembl/chembl_beaker/blob/master/src/chembl_beaker/beaker/core_apps/conversions/views.py>`_
 * `Compound recognition <https://github.com/chembl/chembl_beaker/blob/master/src/chembl_beaker/beaker/core_apps/osra/views.py>`_
 * `Vector image (SVG) generation, including Similarity Maps <https://github.com/chembl/chembl_beaker/blob/master/src/chembl_beaker/beaker/core_apps/svgImages/views.py>`_
 * `Descriptors <https://github.com/chembl/chembl_beaker/blob/master/src/chembl_beaker/beaker/core_apps/descriptors/views.py>`_
 * `Maximum Common Substructure <https://github.com/chembl/chembl_beaker/blob/master/src/chembl_beaker/beaker/core_apps/mcs/views.py>`_
 * `Structural alerts <https://github.com/chembl/chembl_beaker/blob/master/src/chembl_beaker/beaker/core_apps/structuralAlerts/views.py>`_
 * `ChEMBL Structure Pipeline <https://github.com/chembl/ChEMBL_Structure_Pipeline>`_
 * Marvin 4 JS compilant `webservices <https://marvin4js.chemaxon.com/marvin4js-latest/docs/dev/webservices.html>`_

As a portable, lightweight, `CORS <https://en.wikipedia.org/wiki/Cross-origin_resource_sharing>`_-ready, `REST <https://en.wikipedia.org/wiki/Representational_state_transfer>`_-speaking, `SPORE <https://github.com/SPORE/specifications>`_-documented webserver. This particular implementation wraps RDKit in `Bottle <http://bottlepy.org/docs/dev/>`_ on `Tornado <http://www.tornadoweb.org/en/stable/>`_.

Where is it used?
--------

It can be used as web service backend for `Marvin For Java Script <http://www.chemaxon.com/products/marvin/marvin-for-javascript/>`_ as it exposes methods compatible with it's webservice `specification <https://marvinjs-demo.chemaxon.com/latest/docs/dev/webservices.html>`_.
To do this you need to configure marvin sketcher instance:

::

    marvin.sketcherInstance = new marvin.Sketch("sketch");
    marvin.sketcherInstance.setServices(getDefaultServices({
        'clean2dws' : <url of  Beaker clean webservice>,
        'molconvertws' : <url of Beaker molExport webservice>,
        "stereoinfows" : <url of Beaker cipStereoInfo>
    }));

Usage
--------

The best way is to use Docker::

    docker run -p 5000:5000 chembl/beaker

open http://127.0.0.1:5000/docs

Configuration
--------
By default configuration is stored in ``beaker.conf`` file, located in current directory. You can specify location of
configuration file using ``--config (-c)`` parameter when running beaker. Configuration file format is standard ``*.ini``.
Beaker is distributed with example configuration file named ``beaker.conf.sample``.

 * **debug** - run bottle server in debug mode (True/False, default ``True``)
 * **bottle_port** - number of port on which Bottle server is listening for connections (integer, default ``8080``)
 * **bottle_host** - hostname of Bottle server (string, default ``localhost``)
 * **server_middleware** - networking middleware library used by Bottle (string, default ``tornado``)
 * **osra_binaries_location** - path to OSRA binary you want to use for compound recognition (string, default ``/usr/bin/osra``)
 * **enable_cors** - enable CORS plugin and respect all header settings below (True/False, default ``True``) 
 * **access_control_allow_origin** - content of 'Access-Control-Allow-Origin' header send with every response (string, default ``*``)
 * **access_control_allow_methods** - content of 'Access-Control-Allow-Methods' header send with every response (string, default ``GET, POST, PUT, OPTIONS``)
 * **installed_apps** - apps installed in beaker, default to [
    "beaker",
    "beaker.core_apps.autoDocs",
    "beaker.core_apps.conversions",
    "beaker.core_apps.descriptors",
    "beaker.core_apps.marvin",
    "beaker.core_apps.mcs",
    "beaker.core_apps.osra",
    "beaker.core_apps.rasterImages",
    "beaker.core_apps.svgImages",
    "beaker.core_apps.similarityMaps",
    "beaker.core_apps.standardisation",
    ]

Documentation
--------
Like every good software written in Python, beaker is self-documented. When you run beaker, open your browser and go to URL: ``BEAKER_ROOT/docs``. You will see live documentation genrated on the fly from all available webservices, exposed by beaker. You can immediately try them and see results return by the server. Every webservice should be documented - documentation is generated automatically as well, from docstring of every exposed webservice, interpreted as markdown.

Development - writing your own extentions
--------
Developing new app should be easy. The only required file is ``views.py`` where you should define your bottle ``routes``.
