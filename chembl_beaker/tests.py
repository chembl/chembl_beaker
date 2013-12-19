__author__ = 'mnowotka'

import os
import unittest
import requests
import chembl_beaker
import base64
from chembl_beaker import settings

#-----------------------------------------------------------------------------------------------------------------------

class TestServer(unittest.TestCase):

#-----------------------------------------------------------------------------------------------------------------------

    def test_status(self):
        status = requests.get("http://%s:%s/status" % (settings.BOTTLE_HOST, settings.BOTTLE_PORT))
        self.assertTrue(status.ok)

#-----------------------------------------------------------------------------------------------------------------------

    def test_ctab2inchi(self):
        dir = os.path.dirname(chembl_beaker.__file__)
        f = open(os.path.join(dir, 'samples', 'a.mol'))
        r = requests.post("http://%s:%s/ctab2inchi" % (settings.BOTTLE_HOST, settings.BOTTLE_PORT), data=f.read())
        self.assertTrue(r.ok)
        self.assertEqual(r.content, 'InChI=1S/C16H17N5O2/c1-2-14(15(22)12-3-7-17-8-4-12)19-11-20-21-16(23)13-5-9-18-10-'
                                    '6-13/h2-10,14,19-20H,1,11H2,(H,21,23)\nInChI=1S/C15H19NO/c1-10-14-8-11-4-5-12(17)9'
                                    '-13(11)15(10,2)6-7-16(14)3/h4-5,9-10,14H,1,6-8H2,2-3H3/p+1\nInChI=1S/C17H16O3/c1-1'
                                    '2(17(19)20)11-16(18)15-9-7-14(8-10-15)13-5-3-2-4-6-13/h2-10,12H,11H2,1H3,(H,19,20)'
                                    '\n')

#-----------------------------------------------------------------------------------------------------------------------

    def test_inchi2ctab(self):
        inchi = 'InChI=1S/C19H20FNO3/c20-15-3-1-13(2-4-15)17-7-8-21-10-14(17)11-22-16-5-6-18-19(9-16)24-12-23-18/h1' \
                '-6,9,14,17,21H,7-8,10-12H2/t14?,17-/m0/s1/i1D,2D,3D,4D,8D2,11D2,12D2,14D,17D'

        ctab_post = requests.post("http://%s:%s/inchi2ctab" % (settings.BOTTLE_HOST, settings.BOTTLE_PORT), data=inchi)
        self.assertTrue(ctab_post.ok)

        ctab_get = requests.get("http://%s:%s/inchi2ctab/%s" % (settings.BOTTLE_HOST, settings.BOTTLE_PORT,
                                                                                    base64.urlsafe_b64encode(inchi)))
        self.assertTrue(ctab_get.ok)

        self.assertEqual(ctab_post.content, ctab_get.content)

#-----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

#-----------------------------------------------------------------------------------------------------------------------
