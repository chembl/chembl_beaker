__author__ = 'mnowotka'

import os
import unittest
import chembl_beaker
from chembl_beaker.run_beaker import app as beaker
import base64
from webtest import TestApp
import bottle
import re

# ----------------------------------------------------------------------------------------------------------------------


class TestServer(unittest.TestCase):

    def setUp(self):
        bottle.debug(True)
        self.app = TestApp(beaker)
        dir = os.path.dirname(chembl_beaker.__file__)
        f = open(os.path.join(dir, 'samples', 'sample.sdf'))
        self.sample_mol_data = f.read()

# ----------------------------------------------------------------------------------------------------------------------

    def test_status(self):
        status = self.app.get("/status")
        self.assertEqual(status.status_int, 200)

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2inchi(self):
        r = self.app.post("/ctab2inchi", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, 'InChI=1S/C16H17N5O2/c1-2-14(15(22)12-3-7-17-8-4-12)19-11-20-21-16(23)13-5-9-18-10-'
                                 '6-13/h2-10,14,19-20H,1,11H2,(H,21,23)\nInChI=1S/C15H19NO/c1-10-14-8-11-4-5-12(17)9'
                                 '-13(11)15(10,2)6-7-16(14)3/h4-5,9-10,14H,1,6-8H2,2-3H3/p+1\nInChI=1S/C17H16O3/c1-1'
                                 '2(17(19)20)11-16(18)15-9-7-14(8-10-15)13-5-3-2-4-6-13/h2-10,12H,11H2,1H3,(H,19,20)'
                                 '\n')

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2smiles(self):
        r = self.app.post("/ctab2smiles", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)

        allsmiles = r.body.splitlines()[1:]
        self.assertRegexpMatches(allsmiles[0], 'C=CC\(NCNNC\(=O\)c1ccncc1\)C\(=O\)c1ccncc1 [0]?')

        # If accepts either [CH2+]C1C2Cc3ccc(O)cc3C1(C)CCN2C or CN1CCC2(C)c3cc(O)ccc3CC1C2[CH4+]
        self.assertTrue(re.search('\[CH2\+\]C1C2Cc3ccc\(O\)cc3C1\(C\)CCN2C [1]?', allsmiles[1]) or
                        re.search('N1CCC2\(C\)c3cc\(O\)ccc3CC1C2\[CH4\+\] [1]?', allsmiles[1]))

        self.assertRegexpMatches(allsmiles[2], 'CC\(CC\(=O\)c1ccc\(-c2ccccc2\)cc1\)C\(=O\)O [2]?')
        self.assertRegexpMatches(allsmiles[3], 'C\[\*\]CN\(CC\)CC#CC\.Nc1cc\(C2\(C=O\)CCCC2\)ccc1N1CC1 [3]?')

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2json(self):
        r = self.app.post("/ctab2json", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.content_type, 'application/json')
        self.assertTrue('path' in r, r)

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2xyz(self):
        r = self.app.post("/ctab2xyz", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        # it must contain the number of atoms of each of the molecules in the mol data
        # followed by an empty line
        r.mustcontain("40\n\n")
        r.mustcontain("37\n\n")
        r.mustcontain("36\n\n")
        r.mustcontain("60\n\n")

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2svg(self):
        r = self.app.post("/ctab2svg", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.content_type, 'image/svg+xml')
        r.mustcontain('<svg:svg')

# ----------------------------------------------------------------------------------------------------------------------

    def test_TPSA(self):
        r = self.app.post("/tpsa", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, '[96.01, 23.47, 54.37, 49.339999999999996]')

# ----------------------------------------------------------------------------------------------------------------------

    def test_MolWt(self):
        r = self.app.post("/molWt", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        weights = [float(w) for w in r.body.strip('[]').split(',')]
        self.assertAlmostEquals(weights[0], 311, delta=2.5)
        self.assertAlmostEquals(weights[1], 232, delta=2.5)
        self.assertAlmostEquals(weights[2], 268, delta=2.5)
        self.assertAlmostEquals(weights[3], 355, delta=2.5)

# ----------------------------------------------------------------------------------------------------------------------

    def test_logP(self):
        r = self.app.post("/logP", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        logps = [float(w) for w in r.body.strip('[]').split(',')]
        self.assertAlmostEquals(logps[0], 0.69559, places=4)
        self.assertAlmostEquals(logps[1], 2.6692, places=0)
        self.assertAlmostEquals(logps[2], 3.6471, places=4)
        self.assertAlmostEquals(logps[3], 3.5217, places=4)

# ----------------------------------------------------------------------------------------------------------------------

    def test_ctab2image(self):
        from PIL import Image
        import StringIO
        r = self.app.post("/ctab2image", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.content_type, 'image/png')
        buff = StringIO.StringIO(r.body)
        im = Image.open(buff)
        self.assertEqual(im.size, (800, 200))
        im.verify()

# ----------------------------------------------------------------------------------------------------------------------

    def test_inchi2ctab(self):
        inchi = 'InChI=1S/C19H20FNO3/c20-15-3-1-13(2-4-15)17-7-8-21-10-14(17)11-22-16-5-6-18-19(9-16)24-12-23-18/h1' \
                '-6,9,14,17,21H,7-8,10-12H2/t14?,17-/m0/s1/i1D,2D,3D,4D,8D2,11D2,12D2,14D,17D'

        ctab_post = self.app.post("/inchi2ctab", inchi)
        self.assertEqual(ctab_post.status_int, 200)

        ctab_get = self.app.get("/inchi2ctab/%s" % base64.urlsafe_b64encode(inchi))
        self.assertEqual(ctab_get.status_int, 200)

        self.assertEqual(ctab_post.body, ctab_get.body)

# ----------------------------------------------------------------------------------------------------------------------

    def test_mcs(self):
        from rdkit import Chem
        from StringIO import StringIO
        smis = ['c1ccccc1','c1ccccc1C','c1ccccc1O']
        sio = StringIO()
        w = Chem.SDWriter(sio)
        for smi in smis:
            w.write(Chem.MolFromSmiles(smi))
        w.flush()
        txt=sio.getvalue()        
        ctab_post = self.app.post("/mcs", txt)
        self.assertEqual(ctab_post.status_int, 200)
        self.assertTrue(ctab_post.body in ('[#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1', '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'))

        ctab_get = self.app.get("/mcs/%s" % base64.urlsafe_b64encode(txt))
        self.assertEqual(ctab_get.status_int, 200)
        self.assertEqual(ctab_get.body,ctab_post.body)

        ctab_post = self.app.post("/mcs?asSmiles=1", txt)
        self.assertEqual(ctab_post.status_int, 200)
        self.assertEqual(ctab_post.body,'c1ccccc1')

# ----------------------------------------------------------------------------------------------------------------------

    def test_break_bonds(self):
        input_smiles = "[Na]OC(=O)c1ccccc1"
        r = self.app.post("/smiles2ctab", input_smiles)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/breakBonds", mol)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/ctab2smiles", mol)
        self.assertEqual(r.status_int, 200)
        bb_smiles = r.body

        # accepts either [Na+].O=C([O-])c1ccccc1 or 'O=C([O-])c1ccccc1.[Na+]'
        self.assertTrue(re.match('SMILES Name \n\[Na\+\]\.O=C\(\[O\-\]\)c1ccccc1 [0]?\n', bb_smiles) or
                        re.match('SMILES Name \nO=C\(\[O\-\]\)c1ccccc1\.\[Na\+\] [0]?\n', bb_smiles),
                        msg="Response body doesn't match any of the valid regular expressions")

# ----------------------------------------------------------------------------------------------------------------------

    def test_neutralise(self):
        input_smiles = "C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])"
        r = self.app.post("/smiles2ctab", input_smiles)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/neutralise", mol)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/ctab2smiles", mol)
        self.assertEqual(r.status_int, 200)
        n_smiles = r.body
        self.assertRegexpMatches(n_smiles, 'SMILES Name \nNCC\(Cc1nn\[nH\]n1\)\(C\[N\+\]\(=O\)\[O\-\]\)C\(=O\)O [0]?\n')

# ----------------------------------------------------------------------------------------------------------------------

    def test_rules(self):
        input_smiles = "Oc1nccc2cc[nH]c(=N)c12"
        r = self.app.post("/smiles2ctab", input_smiles)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/rules", mol)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/ctab2smiles", mol)
        self.assertEqual(r.status_int, 200)
        r_smiles = r.body
        self.assertRegexpMatches(r_smiles, 'SMILES Name \nNc1nccc2cc\[nH\]c\(=O\)c12 [0]?\n')

# ----------------------------------------------------------------------------------------------------------------------

    def test_standardise(self):
        input_smiles = "[Na]OC(=O)Cc1ccc(C[NH3+])cc1.c1nnn[n-]1.O"
        r = self.app.post("/smiles2ctab", input_smiles)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/standardise", mol)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/ctab2smiles", mol)
        self.assertEqual(r.status_int, 200)
        s_smiles = r.body
        self.assertRegexpMatches(s_smiles, 'SMILES Name \nNCc1ccc\(CC\(=O\)O\)cc1 [0]?\n')

# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

# ----------------------------------------------------------------------------------------------------------------------
