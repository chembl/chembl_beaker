__author__ = 'mnowotka'

import os
import unittest
import beaker
from run_beaker import app as beaker
import base64
from webtest import TestApp
import bottle
import re

# ----------------------------------------------------------------------------------------------------------------------


class TestServer(unittest.TestCase):

    def setUp(self):
        bottle.debug(True)
        self.app = TestApp(beaker)
        dr = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dr, "samples", "sample.sdf")) as f:
            self.sample_mol_data = f.read()

# ----------------------------------------------------------------------------------------------------------------------

    def test_status(self):
        status = self.app.get("/status")
        self.assertEqual(status.status_int, 200)

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2inchi(self):
        r = self.app.post("/ctab2inchi", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body.decode("utf-8"), 'InChI=1S/C16H17N5O2/c1-2-14(15(22)12-3-7-17-8-4-12)19-11-20-21-16(23)13-5-9-18-10-'
                                 '6-13/h2-10,14,19-20H,1,11H2,(H,21,23)\nInChI=1S/C15H19NO/c1-10-14-8-11-4-5-12(17)9'
                                 '-13(11)15(10,2)6-7-16(14)3/h4-5,9-10,14H,1,6-8H2,2-3H3/p+1\nInChI=1S/C17H16O3/c1-1'
                                 '2(17(19)20)11-16(18)15-9-7-14(8-10-15)13-5-3-2-4-6-13/h2-10,12H,11H2,1H3,(H,19,20)'
                                 '\n')

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2smiles(self):
        r = self.app.post("/ctab2smiles", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)

        allsmiles = r.body.splitlines()[1:]
        self.assertRegex(allsmiles[0], b'C=CC\(NCNNC\(=O\)c1ccncc1\)C\(=O\)c1ccncc1 [0]?')

        # If accepts either [CH2+]C1C2Cc3ccc(O)cc3C1(C)CCN2C or CN1CCC2(C)c3cc(O)ccc3CC1C2[CH4+]
        self.assertTrue(re.search(b'\[CH2\+\]C1C2Cc3ccc\(O\)cc3C1\(C\)CCN2C [1]?', allsmiles[1]) or
                        re.search(b'N1CCC2\(C\)c3cc\(O\)ccc3CC1C2\[CH4\+\] [1]?', allsmiles[1]))

        self.assertRegex(allsmiles[2], b'CC\(CC\(=O\)c1ccc\(-c2ccccc2\)cc1\)C\(=O\)O [2]?')
        self.assertEqual(allsmiles[3], b'C*CN(CC)CC#CC.Nc1cc(C2(C=O)CCCC2)ccc1N1CC1 3')

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2svg(self):
        r = self.app.post("/ctab2svg", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.content_type, 'image/svg+xml')
        r.mustcontain('svg')

# ----------------------------------------------------------------------------------------------------------------------


    def test_inchi2ctab(self):

        ctab = b"""
     RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""

        inchi = 'InChI=1S/CH4/h1H4'

        ctab_post = self.app.post("/inchi2ctab", inchi)
        self.assertEqual(ctab_post.status_int, 200)
        self.assertEqual(ctab_post.body, ctab)

# ----------------------------------------------------------------------------------------------------------------------

    def test_mcs(self):
        from rdkit import Chem
        from io import StringIO
        smis = ['c1ccccc1','c1ccccc1C','c1ccccc1O']
        sio = StringIO()
        w = Chem.SDWriter(sio)
        for smi in smis:
            w.write(Chem.MolFromSmiles(smi))
        w.flush()
        txt=sio.getvalue()
        ctab_post = self.app.post("/mcs", txt)
        self.assertEqual(ctab_post.status_int, 200)
        self.assertEqual(ctab_post.body.decode("utf-8"), '[#6]1=[#6]-[#6]=[#6]-[#6]=[#6]-1')

        ctab_post = self.app.post("/mcs?asSmiles=1", txt)
        self.assertEqual(ctab_post.status_int, 200)
        self.assertEqual(ctab_post.body, b'C1=CC=CC=C1')

# ----------------------------------------------------------------------------------------------------------------------

    def test_standardize(self):
        input_smiles = "[Na]OC(=O)Cc1ccc(C[NH3+])cc1.c1nnn[n-]1.O"
        r = self.app.post("/smiles2ctab", input_smiles)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/standardize", mol)
        self.assertEqual(r.status_int, 200)
        mol = r.json[0]['standard_molblock']
        r = self.app.post("/ctab2smiles", mol)
        self.assertEqual(r.status_int, 200)
        s_smiles = r.body
        self.assertEqual(s_smiles, b'SMILES Name \nNCc1ccc(CC(=O)[O-])cc1.O.[Na+].c1nnn[nH]1 0\n')

# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

# ----------------------------------------------------------------------------------------------------------------------
