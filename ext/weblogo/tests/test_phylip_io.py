#!/usr/bin/env python

#  Copyright (c) 2006, The Regents of the University of California, through
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  (1) Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#
#  (2) Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and or other materials provided with the distribution.
#
#  (3) Neither the name of the University of California, Lawrence Berkeley
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors
#  may be used to endorse or promote products derived from this software
#  without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.

from io import StringIO

import unittest

from weblogo.seq import protein_alphabet
from weblogo.seq_io import phylip_io, plain_io, clustal_io
from . import data_stream


class test_phylip_io(unittest.TestCase):
    def test_read(self):
        f = data_stream("phylip_test_1.phy")
        seqs = phylip_io.read(f)
        # print seqs
        self.assertEqual(len(seqs), 10)
        self.assertEqual(seqs[0].name, "Cow")
        self.assertEqual(len(seqs[1]), 234)
        f.close()

    def test_iterseq(self):
        f = data_stream("phylip_test_1.phy")
        n = 0
        for seq in phylip_io.iterseq(f):
            n += 1
        assert n == 10
        f.close()

    def test_parse_plain_fail(self):
        # should fail with parse error
        f = StringIO(plain_io.example)
        self.assertRaises(ValueError,
                          phylip_io.read, f)

    def test_parse_phylip_test_2(self):
        f = data_stream('phylip_test_2.phy')
        seqs = phylip_io.read(f)
        self.assertEqual(len(seqs), 6)
        self.assertEqual(len(seqs[0]), 20)
        self.assertEqual(str(seqs[1]), "CGTTACTCGTTGTCGTTACT")
        self.assertEqual(seqs[1].name, "Hesperorni")
        f.close()

    def test_parse_clustal_fail(self):
        # should fail with parse error
        f = StringIO(clustal_io.example)
        self.assertRaises(ValueError,
                          phylip_io.read, f, protein_alphabet)

    def test_parse_phylip_test_3(self):
        f = data_stream('phylip_test_3.phy')
        seqs = phylip_io.read(f)
        self.assertEqual(len(seqs), 6)
        self.assertEqual(len(seqs[0]), 20)
        self.assertEqual(str(seqs[1]), "CGTTACTCGTTGTCGTTACT")
        self.assertEqual(seqs[1].name, "Hesperorni")
        f.close()

    def test_parse_phylip_test_4(self):
        f = data_stream('phylip_test_4.phy')
        seqs = phylip_io.read(f)
        self.assertEqual(len(seqs), 6)
        self.assertEqual(len(seqs[0]), 25)
        self.assertEqual(str(seqs[1]), "GTGGTGGTGGGCGCCGGCCGTGTGG")
        self.assertEqual(seqs[2].name, "ddrasa")
        f.close()

    def test_parse_phylip_test_5(self):
        f = data_stream('phylip_test_5.phy')
        seqs = phylip_io.read(f)
        self.assertEqual(len(seqs), 6)
        self.assertEqual(len(seqs[0]), 50)
        self.assertEqual(str(seqs[1]), "GTGGTGGTGGGCGCCGGCCGTGTGGGTGGTGGTGGGCGCCGGCCGTGTGG")
        self.assertEqual(seqs[2].name, "ddrasa")
        f.close()

    def test_parse_wrong_phylip_codes_1(self):
        f = data_stream('phylip_test_6.corrupt.phy')
        self.assertRaises(ValueError,
                          phylip_io.read, f, protein_alphabet)
        f.close()

    def test_parse_wrong_phylip_codes_2(self):
        f = data_stream('phylip_test_7.corrupt.phy')
        self.assertRaises(ValueError,
                          phylip_io.read, f, protein_alphabet)
        f.close()

    def test_parse_phylip_dna(self):
        f = data_stream('dna.phy')
        seqs = phylip_io.read(f)
        self.assertEqual(len(seqs), 10)
        self.assertEqual(len(seqs[0]), 705)
        self.assertEqual(str(seqs[1][0:10]), "ATGGCACACC")
        self.assertEqual(seqs[2].name, "Chicken")
        f.close()


if __name__ == '__main__':
    unittest.main()
