__author__ = 'jfmosquera@ebi.ac.uk'

import os
import sys
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = os.path.abspath(os.path.join(SCRIPT_PATH, os.pardir))
BEAKER_PATH = os.path.join(SCRIPT_DIR, 'src')
sys.path.append(BEAKER_PATH)

import chembl_beaker
import chembl_beaker.run_beaker

chembl_beaker.__version__ += '-dev'

# WARNING: most methods including spore, and docs are cached in mongo
# Use -c <config_file> when running standalone
chembl_beaker.run_beaker.main()
