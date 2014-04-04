__author__ = 'mnowotka'

import os
import unittest
import chembl_beaker
from chembl_beaker.run_beaker import app as beaker
import base64
from webtest import TestApp

#-----------------------------------------------------------------------------------------------------------------------

class TestServer(unittest.TestCase):

    def setUp(self):
        self.app = TestApp(beaker)
        dir = os.path.dirname(chembl_beaker.__file__)
        f = open(os.path.join(dir, 'samples', 'sample.sdf'))
        self.sample_mol_data = f.read()

#-----------------------------------------------------------------------------------------------------------------------

    def test_status(self):
        status = self.app.get("/status")
        self.assertEqual(status.status_int, 200)

#-----------------------------------------------------------------------------------------------------------------------

    def test_ctab2inchi(self):
        r = self.app.post("/ctab2inchi", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, 'InChI=1S/C16H17N5O2/c1-2-14(15(22)12-3-7-17-8-4-12)19-11-20-21-16(23)13-5-9-18-10-'
                                    '6-13/h2-10,14,19-20H,1,11H2,(H,21,23)\nInChI=1S/C15H19NO/c1-10-14-8-11-4-5-12(17)9'
                                    '-13(11)15(10,2)6-7-16(14)3/h4-5,9-10,14H,1,6-8H2,2-3H3/p+1\nInChI=1S/C17H16O3/c1-1'
                                    '2(17(19)20)11-16(18)15-9-7-14(8-10-15)13-5-3-2-4-6-13/h2-10,12H,11H2,1H3,(H,19,20)'
                                    '\n')

#-----------------------------------------------------------------------------------------------------------------------

    def test_ctab2smiles(self):
        r = self.app.post("/ctab2smiles", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, 'SMILES Name \nC=CC(NCNNC(=O)c1ccncc1)C(=O)c1ccncc1 \nCN1CCC2(C)c3cc(O)ccc3CC1C2[CH4+] '
                                 '\nCC(CC(=O)c1ccc(-c2ccccc2)cc1)C(=O)O \nC[*]CN(CC)CC#CC.Nc1cc(C2(C=O)CCCC2)ccc1N1CC1 '
                                 '\n')

#-----------------------------------------------------------------------------------------------------------------------

    def test_ctab2json(self):
        r = self.app.post("/ctab2json", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.content_type, 'application/json')
        r.mustcontain('M96.3780274809,109.985416295L83.8772809167,130.335345966')

#-----------------------------------------------------------------------------------------------------------------------

    def test_ctab2svg(self):
        r = self.app.post("/ctab2svg", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.content_type, 'image/svg+xml')
        r.mustcontain('<svg xmlns="http://www.w3.org/2000/svg"')

#-----------------------------------------------------------------------------------------------------------------------

    def test_TPSA(self):
        r = self.app.post("/tpsa", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, '[96.01, 23.47, 54.37, 49.339999999999996]')

#-----------------------------------------------------------------------------------------------------------------------

    def test_MolWt(self):
        r = self.app.post("/molWt", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, '[311.345, 232.347, 268.312, 355.5260000000001]')

#-----------------------------------------------------------------------------------------------------------------------

    def test_logP(self):
        r = self.app.post("/logP", self.sample_mol_data)
        self.assertEqual(r.status_int, 200)
        self.assertEqual(r.body, '[0.6955999999999999, 2.6692, 3.6471, 3.5217]')

#-----------------------------------------------------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------------------------------------------------

    def test_inchi2ctab(self):
        inchi = 'InChI=1S/C19H20FNO3/c20-15-3-1-13(2-4-15)17-7-8-21-10-14(17)11-22-16-5-6-18-19(9-16)24-12-23-18/h1' \
                '-6,9,14,17,21H,7-8,10-12H2/t14?,17-/m0/s1/i1D,2D,3D,4D,8D2,11D2,12D2,14D,17D'

        ctab_post = self.app.post("/inchi2ctab", inchi)
        self.assertEqual(ctab_post.status_int, 200)

        ctab_get = self.app.get("/inchi2ctab/%s" % base64.urlsafe_b64encode(inchi))
        self.assertEqual(ctab_get.status_int, 200)

        self.assertEqual(ctab_post.body, ctab_get.body)

#-----------------------------------------------------------------------------------------------------------------------

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
        self.assertEqual(ctab_post.body,'[#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1')

        ctab_get = self.app.get("/mcs/%s" % base64.urlsafe_b64encode(txt))
        self.assertEqual(ctab_get.status_int, 200)
        self.assertEqual(ctab_get.body,ctab_post.body)

        ctab_post = self.app.post("/mcs?asSmiles=1", txt)
        self.assertEqual(ctab_post.status_int, 200)
        self.assertEqual(ctab_post.body,'c1ccccc1')

#-----------------------------------------------------------------------------------------------------------------------

    def test_break_bonds(self):
        input_smiles = "[Na]OC(=O)c1ccccc1"
        r = self.app.post("/smiles2ctab", input_smiles)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/break_bonds", mol)
        self.assertEqual(r.status_int, 200)
        mol = r.body
        r = self.app.post("/ctab2smiles", mol)
        self.assertEqual(r.status_int, 200)
        bb_smiles = r.body
        self.assertEqual('SMILES Name \n[Na+].O=C([O-])c1ccccc1 \n', bb_smiles)

#-----------------------------------------------------------------------------------------------------------------------

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
        self.assertEqual('SMILES Name \nNCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O \n', n_smiles)

#-----------------------------------------------------------------------------------------------------------------------

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
        self.assertEqual('SMILES Name \nNc1nccc2cc[nH]c(=O)c12 \n', r_smiles)

#-----------------------------------------------------------------------------------------------------------------------

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
        self.assertEqual('SMILES Name \nNCc1ccc(CC(=O)O)cc1 \n', s_smiles)

#-----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

#-----------------------------------------------------------------------------------------------------------------------