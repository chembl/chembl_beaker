chembl_beaker
======

What is Beaker?
--------

This is chembl_beaker package developed at `Chembl <https://www.ebi.ac.uk/chembl/>`_ group, `EMBL-EBI <https://www.ebi.ac.uk/>`_, Cambridge, UK.

This is wrapper for `RDKit <http://www.rdkit.org/>`_ and `OSRA <http://cactus.nci.nih.gov/osra/>`_, which exposes following RDKit functions:

 * Format convertion
 * Fingerprints
 * Descriptors

As a portable lightweight webserver, speaking REST. This particular implementation wraps RDKit in `Bottle <http://bottlepy.org/docs/dev/>`_ on `Tornado <http://www.tornadoweb.org/en/stable/>`_.

Where is it used?
--------

Beaker is used in `Clippy <https://github.com/madgpap/chembl_clippy>`_ project but can be used as a standalone web server as well.

Software dependencies
--------

 * `RDKit <http://www.rdkit.org/>`_
 * `OSRA <http://cactus.nci.nih.gov/osra/>`_
 * `Bottle <http://bottlepy.org/docs/dev/>`_
 * `Tornado <http://www.tornadoweb.org/en/stable/>`_

Configuration
--------
Right now all configuration is stored in chembl_beaker/settings.py. This file contains following configuration option, feel free to modify it to fit your environment:

 * **DEBUG** - run bottle server in debug mode (True/False)
 * **BOTTLE_PORT** - number of porn on which Bottle server is listening for connections (integer)
 * **BOTTLE_HOST** - hostname of Bottle server (string)
 * **SERVER_MIDDLEWARE** - networking middleware library used by Bottle
 * **OSRA_BINARIES_LOCATION** - dictionary of ``(version, path)`` pairs. Each pair describes location of certain OSRA version (for example ``{'2.0.0': '/usr/bin/osra'}``). Right now beaker will choose latest version.

Running
--------
If you want to play with beaker run ``python chembl_beaker/runserver.py``
If you want to run beaker in production you should do this using virtualenv, uWSGI and NGINX as described `here <http://fclef.wordpress.com/2013/01/12/bottle-virtualenv-uwsgi-nginx-installation-on-ubuntu-12-04-1-lts/>`_. Other standard python deployment stacks will work as well.

More info and help
--------

More information can be found in `web based presentation <https://github.com/mnowotka/beaker-presentation>`_. You can always mail the author: mmmnow@gmail.com
