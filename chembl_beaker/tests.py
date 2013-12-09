__author__ = 'mnowotka'

import os
import unittest
import requests
import chembl_beaker
from chembl_beaker import settings

class TestServer(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()