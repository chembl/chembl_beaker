chembl_beaker
======

.. image:: https://dl.dropboxusercontent.com/u/10967207/static/logo.png
    :alt: logo

.. image:: https://pypip.in/v/chembl_beaker/badge.png
    :target: https://crate.io/packages/chembl_beaker/
    :alt: Latest PyPI version

.. image:: https://pypip.in/d/chembl_beaker/badge.png
    :target: https://crate.io/packages/chembl_beaker/
    :alt: Number of PyPI downloads
    
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
 * `Maximum Common Substructure <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/mcs/views.py>`_
 * `Smiliarity maps <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/similarityMaps/views.py>`_
 * `Various other calculation (for example kekulisation) <https://github.com/mnowotka/chembl_beaker/blob/master/chembl_beaker/beaker/core_apps/calculations/views.py>`_
 * Marvin 4 JS compilant `webservices <https://marvin4js.chemaxon.com/marvin4js-latest/docs/dev/webservices.html>`_

As a portable, lightweight, `CORS <https://en.wikipedia.org/wiki/Cross-origin_resource_sharing>`_-ready, `REST <https://en.wikipedia.org/wiki/Representational_state_transfer>`_-speaking, `SPORE <https://github.com/SPORE/specifications>`_-documented webserver. This particular implementation wraps RDKit in `Bottle <http://bottlepy.org/docs/dev/>`_ on `Tornado <http://www.tornadoweb.org/en/stable/>`_.

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

 * `pycairo <http://cairographics.org/pycairo/>`_/`cairocffi <https://github.com/SimonSapin/cairocffi>`_ (for `SVG <https://en.wikipedia.org/wiki/Scalable_Vector_Graphics>`_ format support)
 * `lxml <http://lxml.de/>`_ (`mrv <https://www.chemaxon.com/marvin/help/formats/mrv-doc.html>`_ file format)
 * `matplotlib <http://matplotlib.org/>`_ (generating similarity maps)

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
    "chembl_beaker.beaker.core_apps.svgImages",
    "chembl_beaker.beaker.core_apps.jsonImages",
    "chembl_beaker.beaker.core_apps.autoDocs",
    ]

Running
--------
If you want to play with beaker run ``python run_beaker.py``
If you want to run beaker in production you should do this using virtualenv, uWSGI and NGINX as described `here <http://fclef.wordpress.com/2013/01/12/bottle-virtualenv-uwsgi-nginx-installation-on-ubuntu-12-04-1-lts/>`_. Other standard python deployment stacks will work as well.

Documentation
--------
Like every good software written in Python, beaker is self-documented. When you run beaker, open your browser and go to URL: ``BEAKER_ROOT/docs``. You will see live documentation genrated on the fly from all available webservices, exposed by beaker. You can immediately try them and see results return by the server. Every webservice should be documented - documentation is generated automatically as well, from docstring of every exposed webservice, interpreted as markdown.

.. image:: https://dl.dropboxusercontent.com/u/10967207/static/docs.png
    :alt: docs screenshot

Development - writing your own extentions
--------
Developing new app should be easy. The only required file is ``views.py`` where you should define your botte ``routes``. Since your app is technically speaking a python module, ``__init__.py`` will be required as well.
You should wrap your module in ``PIP`` package and distribute via ``PyPi``. By doing so, a user who want to install your app has to install it via `PIP` and add it to ``installed_apps`` list.

FAQ
--------
1. How can I enable HTMl5 canvas support?

 - Current version of ``RDKit`` doesn't include ``JSONCanvas`` so there are two way to get this working. The first one requires a little time and almost no knowlegde, the second one requires little knowledge but it's fast.
     1. Checkout git branch with JSON Canvas and recompile RDKit:
     
::

         git clone -b JSONCanvas_Nov2013 https://github.com/rdkit/rdkit.git
         cd rdkit/External/INCHI-API/ & ./download-inchi.sh
         cd rdkit & mkdir build & cd build & cmake  -DRDK_BUILD_INCHI_SUPPORT=ON ..
         make install


2. Since JSON Canvas is implemented in Python there is no need to recompile, you can just apply a patch: https://github.com/rdkit/rdkit/compare/JSONCanvas_Nov2013 and this should do the trick.     
         

More info and help
--------

More information can be found in `web based presentation <http://mnowotka.github.io/presentations/beaker>`_. You can always email the author: mmmnow@gmail.com
