#!/usr/bin/env python

from io import StringIO
import unittest

import numpy as np

from weblogo import data
from weblogo.matrix import AlphabeticArray, Motif, SubMatrix
from weblogo.seq import protein_alphabet, Alphabet, unambiguous_protein_alphabet
from . import data_stream


class test_AlphabeticArray(unittest.TestCase):
    def test_create(self):
        matrix = AlphabeticArray((protein_alphabet, protein_alphabet))
        matrix['A', 'C'] = 10
        assert matrix[0, 1] == 10.0


class test_Motif(unittest.TestCase):
    def test_read_transfac_alphabet_superset(self):
        with data_stream("transfac_matrix.txt") as f:
            Motif.read_transfac(f, alphabet='TCGA')

        # Supplied alphabet can be superset of defacto alphabet.
        # Reverts to defacto alphabet
        with data_stream("transfac_matrix.txt") as f:
            Motif.read_transfac(f, alphabet='TCGAXYZ')

    def test_read_transfac(self):
        f = data_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)
        f.close()
        assert m[3, 'A'] == 0.0
        assert m[0, 'G'] == 2.0
        assert np.shape(m.array) == (12, 4)
        f.close()

        f = data_stream("transfac_matrix2.txt")
        m = Motif.read_transfac(f)
        f.close()
        assert m[3, 'A'] == 3.0
        assert m[0, 'G'] == 152.0
        assert np.shape(m.array) == (15, 4)

        # this one has extra Ps on start of each line
        f = data_stream("transfac_matrix3.txt")
        m = Motif.read_transfac(f)
        f.close()

    def test_reindex(self):
        f = data_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)
        f.close()
        m2 = m.reindex("TCGA")

        assert (str(m2.alphabet) == "TCGA")

        for k in range(0, 12):
            for i, a in enumerate("AGCT"):
                assert m[k, a] == m2[k, a]

    def test_reverse(self):
        f = data_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)
        f2 = data_stream("transfac_matrix.txt")
        m2 = Motif.read_transfac(f2)
        m2.reverse()

        (K, N) = np.shape(m2)
        for k in range(0, K):
            for n in range(0, N):
                assert (m[k, n] == m2[K - k - 1, n])

        f.close()
        f2.close()

    def test_complement(self):
        f = data_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)
        f2 = data_stream("transfac_matrix.txt")
        m2 = Motif.read_transfac(f2)
        m2.complement()

        (K, N) = np.shape(m2)
        for k in range(0, K):
            assert (m[k, 'A'] == m2[k, 'T'])
            assert (m[k, 'G'] == m2[k, 'C'])
            assert (m[k, 'C'] == m2[k, 'G'])
            assert (m[k, 'T'] == m2[k, 'A'])
        f.close()
        f2.close()

    def test_reverse_complement(self):
        f = data_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)

        f2 = data_stream("transfac_matrix.txt")
        m2 = Motif.read_transfac(f2)

        m.complement()
        m.reverse()

        m2.reverse_complement()

        assert (m.array == m2.array).all()
        f.close()
        f2.close()


class test_SubMatrix(unittest.TestCase):
    def test_create(self):
        ab = 'ABCD'
        ar = np.asarray([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        s = SubMatrix(ab, ar)

        assert s[0, 0] == 1
        assert s['A', 'A'] == 1
        assert s['B', 'C'] == 7
        s['B', 'C'] = -1
        assert s['B', 'C'] == -1

    def test_get(self):
        ab = Alphabet('ABCD')
        ar = np.asarray([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        s = SubMatrix(ab, ar)
        s1 = 'DCCBBBAAA'
        s2 = 'BA'
        v = s.index((s1, s2))
        # print v
        for m, i in enumerate(s1):
            for n, j in enumerate(s2):
                assert s[i, j] == v[m, n]

    def test_get_subMatrix(self):
        ab = Alphabet('ABCD')
        ar = np.asarray([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        mat = SubMatrix(ab, ar)

        mat2 = mat.reindex('ABC')
        assert np.all(mat2.array == np.asarray([[1, 2, 3], [5, 6, 7], [9, 10, 11]]))

        mat2 = mat.reindex('BA')
        assert np.all(mat2.array == np.asarray([[6, 5], [2, 1]]))

        mat2 = mat.reindex(Alphabet('BA'))
        assert np.all(mat2.array == np.asarray([[6, 5], [2, 1]]))

    def test_fail_get(self):
        ab = Alphabet('ABCD')
        ar = np.asarray([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        s = SubMatrix(ab, ar)

        self.assertRaises(IndexError, s.__getitem__, ('E', 'A'))
        self.assertRaises(IndexError, s.__getitem__, ('5', '6'))

        # FIXME
        self.assertRaises(IndexError, s.index, ('E', 'A'))

    def test_repr(self):
        ab = Alphabet('ABCD')
        ar = np.asarray([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        s = SubMatrix(ab, ar)

        repr(s)
        # print string

    def test_read(self):
        f = StringIO(test_matrix1)
        mat = SubMatrix.read(f)
        assert mat['a', 'a'] == 4

    def test_read_asymmetric_fail(self):
        f = StringIO(test_matrix4)
        self.assertRaises(ValueError, SubMatrix.read, f)

    def test_read_alphabets(self):

        # incompatable alphabets
        f = StringIO(test_matrix3)
        self.assertRaises(ValueError,
                          SubMatrix.read, f)

        f = StringIO(test_matrix3)
        SubMatrix.read(f, alphabet=Alphabet('ARNDCQEGHILKMFPSTWYV'))

        f2 = StringIO(test_matrix1)
        self.assertRaises(ValueError,
                          SubMatrix.read, f2, unambiguous_protein_alphabet)

    def test_read_corrupt(self):
        f = StringIO(test_matrix2)
        self.assertRaises(ValueError,
                          SubMatrix.read, f)

    def test_read_pam(self):
        f = data.data_stream("pam250")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[0, 0], 2.0)
        f.close()

        f = data.data_stream("pam120")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[4, 5], -7)
        f.close()

    def test_read_blosum(self):
        f = data.data_stream("blosum80")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[0, 10], -3)
        f.close()

        f = data.data_stream("blosum62")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[4, 5], -4)
        f.close()

    def test_read_blast(self):
        # New style blast matrices have letters at beginning of lines and a '*'
        f = data_stream("blosum35.blast.new")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[4, 5], -3)
        f.close()

        # Matrices formatted for old blast have a '*' (stop)
        # column and no letters at the beggining of lines
        f = data_stream("blosum35.blast")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[0, 10], -2)
        self.assertEqual(mat.array.shape, (23, 23))
        f.close()

        # For comparison, we'll also parse a matrix without '*'
        f = data_stream("pam250.mat")
        mat = SubMatrix.read(f)
        self.assertEqual(mat[4, 5], -5)
        f.close()


test_matrix1 = """# A Test Matrix
# More comments

# And blank line should be ignored

A    4  -2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  -2  -1  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3   3  -1  -1
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4   3  -1  -1
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5   0  -2  -2
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  -3   5  -2
Q   -1   1   0   0  -3   6   1  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3   0   1  -1
E   -1   0   0   1  -4   1   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4   0  -1  -1
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  -2  -2  -2
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3   0  -1  -1
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  -4  -2  -2
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  -4  -2  -2
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3   1  -1  -1
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  -3  -2  -1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  -3  -2  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  -2  -2  -2
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3   0  -1  -1
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  -1  -1  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  -3  -2  -1
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  -2  -2  -1
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  -4  -2  -2
B   -2   3   3   0  -3   0   0  -2   0  -4  -4   1  -3  -3  -2   0  -1  -3  -2  -4   3  -1  -1
Z   -1  -1  -1  -2   5   1  -1  -2  -1  -2  -2  -1  -2  -2  -2  -1  -1  -2  -2  -2  -1   3  -1
X   -1  -1  -1  -2  -2  -1  -1  -2  -1  -2  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1
"""

test_matrix2 = """# An invalid Test Matrix
# Its got a non-numerical value in it. Is the correct exception raised?

# And blank line should be ignored

A    4  -2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  -2  -1  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3   3  -1  -1
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4   3  -1  -1
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5   0  -2  -2
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  -3   5  -2
Q   -1   1   0   0  -3   6   1  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3   0   1  -1
E   -1   0   0   1  -4   x   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4   0  -1  -1
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  -2  -2  -2
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3   0  -1  -1
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  -4  -2  -2
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  -4  -2  -2
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3   1  -1  -1
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  -3  -2  -1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  -3  -2  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  -2  -2  -2
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3   0  -1  -1
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  -1  -1  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  -3  -2  -1
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  -2  -2  -1
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  -4  -2  -2
B   -2   3   3   0  -3   0   0  -2   0  -4  -4   1  -3  -3  -2   0  -1  -3  -2  -4   3  -1  -1
Z   -1  -1  -1  -2   5   1  -1  -2  -1  -2  -2  -1  -2  -2  -2  -1  -1  -2  -2  -2  -1   3  -1
X   -1  -1  -1  -2  -2  -1  -1  -2  -1  -2  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1
"""

test_matrix3 = """#
# This test matrix has a smaller alphabet
A    4  -2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0
Q   -1   1   0   0  -3   6   4  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3
E   -1   0   0   1  -4   4   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5
"""

test_matrix4 = """# This matrix is invalid because it is asymetric! (AR, RA)

A    4   2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  -2  -1  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3   3  -1  -1
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4   3  -1  -1
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5   0  -2  -2
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  -3   5  -2
Q   -1   1   0   0  -3   6   1  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3   0   1  -1
E   -1   0   0   1  -4   1   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4   0  -1  -1
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  -2  -2  -2
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3   0  -1  -1
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  -4  -2  -2
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  -4  -2  -2
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3   1  -1  -1
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  -3  -2  -1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  -3  -2  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  -2  -2  -2
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3   0  -1  -1
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  -1  -1  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  -3  -2  -1
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  -2  -2  -1
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  -4  -2  -2
B   -2   3   3   0  -3   0   0  -2   0  -4  -4   1  -3  -3  -2   0  -1  -3  -2  -4   3  -1  -1
Z   -1  -1  -1  -2   5   1  -1  -2  -1  -2  -2  -1  -2  -2  -2  -1  -1  -2  -2  -2  -1   3  -1
X   -1  -1  -1  -2  -2  -1  -1  -2  -1  -2  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1
"""

if __name__ == '__main__':
    unittest.main()
