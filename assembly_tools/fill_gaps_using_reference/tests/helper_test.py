#!/usr/bin/env python3

import sys
import os
import filecmp
import unittest
import copy
import pysam
import assembly_tools.fill_gaps_using_reference.helper as helper
from fastaq import sequences

modules_dir = os.path.dirname(os.path.abspath(helper.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestHelper(unittest.TestCase):
    def test_make_fasta_of_gap_flanks(self):
        '''Test make_fasta_of_gap_flanks()'''
        tmp_out = 'tmp.fa'
        gaps = {}
        helper.make_fasta_of_gap_flanks(os.path.join(data_dir, 'helper_test_to_be_filled.fa'), 3, tmp_out, gaps)
        self.assertTrue(filecmp.cmp(tmp_out, os.path.join(data_dir, 'helper_test_to_be_filled_gap_flanks.fa')))
        os.unlink(tmp_out)

    def test_paired_hit_samreader(self):
        '''Test paired_hit_samreader()'''
        samreader = helper.paired_hit_samreader(os.path.join(data_dir, 'helper_test_paired_hit_samreader.sam'))
        expected = [(['seq1:1-1.left'], ['seq1:1-1.right']),
                    (['seq1:6-6.left'] * 2, ['seq1:6-6.right'] * 2),
                    (['seq1:13-14.left'], ['seq1:13-14.right'] * 2),
                    (['seq1:9-9.left'] * 2, ['seq1:9-9.right'])]
        i = 0

        for left, right, samfile in samreader:
            left_names = [x.qname for x in left]
            right_names = [x.qname for x in right]
            self.assertListEqual(expected[i][0], left_names)
            self.assertListEqual(expected[i][1], right_names)
            i += 1

    def test_gap_flank_seqname_to_dict_key(self):
        '''Test gap_flank_seqname_to_dict_key()'''
        self.assertEqual(helper._gap_flank_seqname_to_dict_key('name:1-42.left'), ('name', (0, 41)))
        self.assertEqual(helper._gap_flank_seqname_to_dict_key('name:1-42.right'), ('name', (0, 41)))

        bad_names = [
            'name:a-42.left',
            'name:a-42.right',
            'name:1-a.leftt',
            'name:1-a.right',
            'name:1-42.x'
        ]

        for name in bad_names:
            with self.assertRaises(helper.Error):
                helper._gap_flank_seqname_to_dict_key(name)


    def test_parse_sam_file(self):
        '''Test test_parse_sam_file()'''
        # FIXME
        pass

if __name__ == '__main__':
    unittest.main()
