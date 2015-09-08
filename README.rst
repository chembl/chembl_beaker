chembl_beaker
======

.. image:: https://img.shields.io/pypi/v/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: Latest Version

.. image:: https://img.shields.io/pypi/dm/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: Downloads

.. image:: https://img.shields.io/pypi/pyversions/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/status/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: Development Status

.. image:: https://img.shields.io/pypi/l/chembl_beaker.svg
    :target: https://pypi.python.org/pypi/chembl_beaker/
    :alt: License

.. image:: https://badge.waffle.io/chembl/chembl_beaker.png?label=ready&title=Ready 
 :target: https://waffle.io/chembl/chembl_beaker
 :alt: 'Stories in Ready'
    
What is Beaker?
--------

This is chembl_beaker package developed at `ChEMBL <https://www.ebi.ac.uk/chembl/>`_ group, `EMBL-EBI <https://www.ebi.ac.uk/>`_, Cambridge, UK.

This is wrapper for `RDKit <http://www.rdkit.org/>`_ and `OSRA <http://cactus.nci.nih.gov/osra/>`_, which exposes following methods:

 * `Format convertion <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/conversions/views.py>`_
 * `Compound recognition <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/osra/views.py>`_
 * `Raster image (PNG) generation <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/rasterImages/views.py>`_
 * `Vector image (SVG) generation <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/svgImages/views.py>`_
 * `HTML5 ready compound representation <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/jsonImages/views.py>`_
 * `Fingerprints <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/fingerprints/views.py>`_
 * `Descriptors <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/descriptors/views.py>`_
 * `Ring information <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/ringInfo/views.py>`_
 * `Maximum Common Substructure <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/mcs/views.py>`_
 * `Smiliarity maps <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/similarityMaps/views.py>`_
 * `ChEMBL standardisation process <https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/>`_, consisting of neutralisation, bond breaking, salt removal and applying various rules.
 * `3D coordinates generation, using Universal Force Field <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/D3Coords/views.py>`_
 * `Various other calculations (for example kekulisation) <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/calculations/views.py>`_
 * Marvin 4 JS compilant `webservices <https://marvin4js.chemaxon.com/marvin4js-latest/docs/dev/webservices.html>`_

As a portable, lightweight, `CORS <https://en.wikipedia.org/wiki/Cross-origin_resource_sharing>`_-ready, `REST <https://en.wikipedia.org/wiki/Representational_state_transfer>`_-speaking, `SPORE <https://github.com/SPORE/specifications>`_-documented webserver. This particular implementation wraps RDKit in `Bottle <http://bottlepy.org/docs/dev/>`_ on `Tornado <http://www.tornadoweb.org/en/stable/>`_.

Where is it used?
--------

Beaker is used in `Clippy <https://github.com/madgpap/chembl_clippy>`_ project but can be used as a standalone web server as well.

It can also be used as web service backend for `Marvin For Java Script <http://www.chemaxon.com/products/marvin/marvin-for-javascript/>`_ as it exposes methods compatible with it's webservice `specification <https://marvinjs-demo.chemaxon.com/latest/docs/dev/webservices.html>`_.
To do this you need to configure marvin sketcher instance:

::

    marvin.sketcherInstance = new marvin.Sketch("sketch");
    marvin.sketcherInstance.setServices(getDefaultServices({
        'clean2dws' : <url of  Beaker clean webservice>,
        'molconvertws' : <url of Beaker molExport webservice>,
        "stereoinfows" : <url of Beaker cipStereoInfo>
    }));

Beaker and myChEMBL
--------
Beaker is installed on myChEMBL Virtual Machine (currently version 19) so if you want to see how to deploy it for Apache or use it on your laptop or LAN without having to install RDKit and OSRA you can just grab a copy of myChEMBL.
The easiest way to do so is to install `vagrant <https://www.vagrantup.com/>`_ and type::

    vagrant init chembl/myChEMBL
    vagrant up

Software dependencies
--------

 * `RDKit <http://www.rdkit.org/>`_
 * `OSRA <http://cactus.nci.nih.gov/osra/>`_
 * `Bottle <http://bottlepy.org/docs/dev/>`_
 * `Tornado <http://www.tornadoweb.org/en/stable/>`_

Additional dependencies
--------

 * `pycairo <http://cairographics.org/pycairo/>`_/`cairocffi <https://github.com/SimonSapin/cairocffi>`_ (for `SVG <https://en.wikipedia.org/wiki/Scalable_Vector_Graphics>`_ format support)
 * `lxml <http://lxml.de/>`_ (`mrv <https://www.chemaxon.com/marvin/help/formats/mrv-doc.html>`_ file format)
 * `matplotlib <http://matplotlib.org/>`_ (generating similarity maps)
 * `standardiser <https://github.com/flatkinson/standardiser>`_ (Molecular standardisation tool used by Beaker standardisation app)

Installation
--------

The best way to install beaker is to use `PIP`:

    ``pip install chembl_beaker``
    
This command will install latest stable version with Bottle and Tornado. RDKit and OSRA must be installed separately.
You can of course clone development version from github but it's not guaranteed to be working.
If you want to install github version using `PIP`, invoke this command:

    ``sudo pip install git+https://github.com/mnowotka/chembl_beaker.git``

Full recipe for Mac users
--------

So I want to test it, I have a Mac and I don't know what rdkit, tornado and bottle is - how do I start?

First, install XQuartz from https://xquartz.macosforge.org/landing/, then::

      ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
      brew tap edc/homebrew-rdkit
      brew install rdkit
      export RDBASE=/usr/local/share/RDKit
      export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python2.7/site-packages
      export CFLAGS=-Qunused-arguments
      export CPPFLAGS=-Qunused-arguments
      sudo -E pip install cairocffi
      sudo -E pip install Pillow
      sudo -E pip install lxml
      sudo pip install standardiser
      sudo pip install chembl_beaker
      run_berker.py

Alternatively, you can use `this article <http://macinchem.org/reviews/cheminfo/cheminfoMac.php>`_ as an instllation guide.

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
    "chembl_beaker.beaker",
    "chembl_beaker.beaker.core_apps.calculations",
    "chembl_beaker.beaker.core_apps.conversions",
    "chembl_beaker.beaker.core_apps.descriptors",
    "chembl_beaker.beaker.core_apps.fingerprints",
    "chembl_beaker.beaker.core_apps.marvin",
    "chembl_beaker.beaker.core_apps.mcs",
    "chembl_beaker.beaker.core_apps.osra",
    "chembl_beaker.beaker.core_apps.rasterImages",
    "chembl_beaker.beaker.core_apps.ringInfo",
    "chembl_beaker.beaker.core_apps.svgImages",
    "chembl_beaker.beaker.core_apps.jsonImages",
    "chembl_beaker.beaker.core_apps.autoDocs",
    ]

Running
--------
If you want to play with beaker, type ``run_beaker``
If you want to run beaker in production, read section below .

Deploying on Apache/Nginx
--------
Beaker is a Bottle app so it's really easy to deploy it on Apache with mod_wsgi.
Only a few lines of code are required in your .wsgi file::

    from bottle import debug
    import json
    from chembl_beaker.beaker import app, config, loadPlugins, loadApps

    conf_path = "[path to config. file]"
    config.load_config(conf_path)

    apps = json.loads(config.get('installed_apps', '[]'))
    plugins = json.loads(config.get('plugins', '[]'))

    loadApps(apps)
    loadPlugins(app, plugins)

    debug(True)

    application = app

That's it! For details, refer to `this document <http://flask.pocoo.org/docs/deploying/mod_wsgi/>`_.
Everything that can be deployed on Apache with mod_wsgi, can be deployed on Nginx with uWSGI, details `here <http://fclef.wordpress.com/2013/01/12/bottle-virtualenv-uwsgi-nginx-installation-on-ubuntu-12-04-1-lts/>`_.

Documentation
--------
Like every good software written in Python, beaker is self-documented. When you run beaker, open your browser and go to URL: ``BEAKER_ROOT/docs``. You will see live documentation genrated on the fly from all available webservices, exposed by beaker. You can immediately try them and see results return by the server. Every webservice should be documented - documentation is generated automatically as well, from docstring of every exposed webservice, interpreted as markdown.

.. image:: https://dl.dropboxusercontent.com/u/10967207/static/docs.png
    :alt: docs screenshot

Development - writing your own extentions
--------
Developing new app should be easy. The only required file is ``views.py`` where you should define your botte ``routes``. Since your app is technically speaking a python module, ``__init__.py`` will be required as well.
You should wrap your module in ``PIP`` package and distribute via ``PyPi``. By doing so, a user who want to install your app has to install it via `PIP` and add it to ``installed_apps`` list.


More info and help
--------

More information can be found in `web based presentation <http://mnowotka.github.io/presentations/beaker>`_. You can always email the author: mmmnow@gmail.com
