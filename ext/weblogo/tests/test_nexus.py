#!/usr/bin/env python

import unittest

from weblogo.seq_io._nexus import Nexus
from . import data_stream


class test_nexus(unittest.TestCase):
    def test_create(self):
        n = Nexus()
        self.assertNotEqual(n, None)

    def test_parse_f0(self):
        f = data_stream("nexus/test_Nexus_input.nex")
        n = Nexus(f)
        # self.output_basics(n)

        expected = ['t1',
                    "t2 the name",
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t5", "t6", "t7", "t8", "t9"]
        taxa = n.taxlabels
        self.assertEqual(taxa, expected)
        f.close()

    def test_parse_protein(self):
        f = data_stream("nexus/protein.nex")
        Nexus(f)
        f.close()

    def test_parse_dna(self):
        f = data_stream("nexus/dna.nex")
        n = Nexus(f)

        taxa = n.taxlabels
        taxa.sort()
        self.assertEqual(len(taxa), 10)
        self.assertEqual(taxa[0], "Carp")
        self.assertEqual(taxa[-1], "Whale")
        f.close()

    def test_TreeTest1(self):
        """Test Tree module."""
        f = data_stream("nexus/test_Nexus_input.nex")
        n = Nexus(f)
        t3 = n.trees[2]
        n.trees[2]
        t3.root_with_outgroup(['t1', 't5'])

        # Return node_id of common ancestor if
        # taxon_list is monophyletic, -1 otherwise.
        self.assertEqual(t3.is_monophyletic(['t1', 't5']), 13)

        t3.split(parent_id=t3.search_taxon('t9'))
        f.close()


if __name__ == '__main__':
    unittest.main()
