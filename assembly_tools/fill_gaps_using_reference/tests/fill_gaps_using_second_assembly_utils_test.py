#!/usr/bin/env python3

import sys
import os
import filecmp
import unittest
import fill_gaps_using_second_assembly_utils as utils

modules_dir = os.path.dirname(os.path.abspath(utils.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestUtils(unittest.TestCase):
    def test_make_fasta_of_gap_flanks(self):
        '''Check fasta of seqs flanking gaps made correctly'''
        tmp_out = 'tmp.fa'
        utils.make_fasta_of_gap_flanks('fill_gaps_using_second_assembly_utils_test_to_be_filled.fa', 3, tmp_out)
        self.assertTrue(filecmp.cmp(tmp_out, 'fill_gaps_using_second_assembly_utils_test_to_be_filled_gap_flanks.fa'))
        os.unlink(tmp_out)


    def test_gap_flank_seqname_to_gap(self):
        '''Test gap coords extracted OK from sequence name'''
        pass


    def test_paired_hit_samreader(self):
        '''Test paired_hit_samreader() works'''
        samreader = utils.paired_hit_samreader('fill_gaps_using_second_assembly_utils_test.sam')
        expected = [(['seq1:1-1.left'], ['seq1:1-1.right']),
                    (['seq1:6-6.left'] * 2, ['seq1:6-6.right'] * 2),
                    (['seq1:13-14.left'], ['seq1:13-14.right'] * 2),
                    (['seq1:9-9.left'] * 2, ['seq1:9-9.right'])]
        i = 0

        for left, right in samreader:
            left_names = [x.qname for x in left]
            right_names = [x.qname for x in right]
            self.assertListEqual(expected[i][0], left_names)
            self.assertListEqual(expected[i][1], right_names)
            i += 1

   
    def test_deal_with_hits(self):
        '''Test deal_with_hits works as expected'''
        pass

if __name__ == '__main__':
    unittest.main()
