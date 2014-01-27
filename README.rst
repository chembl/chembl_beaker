chembl_beaker
======

.. image:: https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Beaker.svg/200px-Beaker.svg.png
    :alt: logo

.. image:: https://pypip.in/v/chembl_beaker/badge.png
    :target: https://crate.io/packages/chembl_beaker/
    :alt: Latest PyPI version
    
What is Beaker?
--------

This is chembl_beaker package developed at `Chembl <https://www.ebi.ac.uk/chembl/>`_ group, `EMBL-EBI <https://www.ebi.ac.uk/>`_, Cambridge, UK.

This is wrapper for `RDKit <http://www.rdkit.org/>`_ and `OSRA <http://cactus.nci.nih.gov/osra/>`_, which exposes following RDKit functions:

 * Format convertion
 * Compound recognition
 * Image generation
 * Fingerprints
 * Descriptors
 * Marvin 4 JS compilant `webservices <https://marvin4js.chemaxon.com/marvin4js-latest/docs/dev/webservices.html>`_

As a portable, lightweight, `CORS <https://en.wikipedia.org/wiki/Cross-origin_resource_sharing>`_ ready webserver, speaking REST. This particular implementation wraps RDKit in `Bottle <http://bottlepy.org/docs/dev/>`_ on `Tornado <http://www.tornadoweb.org/en/stable/>`_.

Where is it used?
--------

Beaker is used in `Clippy <https://github.com/madgpap/chembl_clippy>`_ project but can be used as a standalone web server as well.

It can also be used as web service backend for `Marvin For Java Script <http://www.chemaxon.com/products/marvin/marvin-for-javascript/>`_ as it exposes methods compatible with it's webservice `specification <https://marvin4js.chemaxon.com/marvin4js-latest/docs/dev/webservices.html>`_. 
To do this you need to configure marvin sketcher instance:

::

    marvin.sketcherInstance = new marvin.Sketch("sketch");
    marvin.sketcherInstance.setServices(getDefaultServices({
        'clean2dws' : <url of  Beaker clean webservice>,
        'molconvertws' : <url of Beaker molExport webservice>,
        "stereoinfows" : <url of Beaker cipStereoInfo>
    }));

Software dependencies
--------

 * `RDKit <http://www.rdkit.org/>`_
 * `OSRA <http://cactus.nci.nih.gov/osra/>`_
 * `Bottle <http://bottlepy.org/docs/dev/>`_
 * `Tornado <http://www.tornadoweb.org/en/stable/>`_

Additional dependencies
--------

 * cairo/cairocffi (for SVG format support)
 * lxml (mrv file format)
 * matplotlib (generating similarity maps)

Installation
--------

The best way to install beaker is to use `PIP`:

    ``pip install chembl_beaker``
    
This command will install latest stable version with Bottle and Tornado. RDKit and OSRA must be installed separately.
You can of course clone development version from github but it's not guaranteed to be working.
If you want to install github version using `PIP`, invoke this command:

    ``sudo pip install git+https://github.com/mnowotka/chembl_beaker.git``

Configuration
--------
By default configuration is stored in ``beaker.conf`` file, located in current directory. You can specify location of
configuration file using ``--config (-c)`` parameter when running beaker. Configuration file format is standard ``*.ini``.
Beaker is distributed with example configuration file named ``beaker.conf.sample``.

 * **debug** - run bottle server in debug mode (True/False, default ``True``)
 * **bottle_port** - number of porn on which Bottle server is listening for connections (integer, default ``8080``)
 * **bottle_host** - hostname of Bottle server (string, default ``localhost``)
 * **server_middleware** - networking middleware library used by Bottle (string, default ``tornado``)
 * **osra_binaries_location** - path to OSRA binary you want to use for compound recognition (string, default ``/usr/bin/osra``)
 * **enable_cors** - enable CORS plugin and respect all header settings below (True/False, default ``True``) 
 * **access_control_allow_origin** - content of 'Access-Control-Allow-Origin' header send with every response (string, default ``*``)
 * **access_control_allow_methods** - content of 'Access-Control-Allow-Methods' header send with every response (string, default ``GET, POST, PUT, OPTIONS``)
 * **access_control_allow_headers** - content of 'Access-Control-Allow-Headers' header send with every response (string, default ``Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token``)

Running
--------
If you want to play with beaker run ``python run_beaker.py``
If you want to run beaker in production you should do this using virtualenv, uWSGI and NGINX as described `here <http://fclef.wordpress.com/2013/01/12/bottle-virtualenv-uwsgi-nginx-installation-on-ubuntu-12-04-1-lts/>`_. Other standard python deployment stacks will work as well.

More info and help
--------

More information can be found in `web based presentation <http://mnowotka.github.io/presentations/beaker>`_. You can always email the author: mmmnow@gmail.com
