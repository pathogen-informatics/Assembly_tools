#!/usr/bin/env python3

import sys
import os
import filecmp
import unittest
import copy
import pysam
import assembly_tools.fill_gaps_using_reference.gap as gap
from fastaq import sequences

modules_dir = os.path.dirname(os.path.abspath(gap.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestGap(unittest.TestCase):
    def test_construct_gap_and_str(self):
        '''Check gap constructor and string conversion'''
        expected_gap = gap.Gap()
        expected_gap.query_name = 'query_name'
        expected_gap.query_start = 9
        expected_gap.query_end = 41
        expected_gap.left_seq = 'LEFTSEQ'
        expected_gap.right_seq= 'RIGHTSEQ'
      
        line = '\t'.join(['query_name', '10', '42', 'LEFTSEQ', 'RIGHTSEQ' ])
        test_gap = gap.Gap(line)
        self.assertEqual(expected_gap, test_gap)
        expected_str = '\t'.join([ 'query_name', '10', '42', '0', '0', '*', '0', '0', '0', 'UNKNOWN' ])
       
        self.assertEqual(str(test_gap), expected_str)

        with self.assertRaises(gap.Error):
            gap.Gap('foo')

    def test_update_hits(self):
        '''Test hits updated OK'''
        def load_sam(filename):
            samfile = pysam.Samfile(filename, "r")
            left_hits = []
            right_hits = []
            for samrecord in samfile.fetch(until_eof=True):
                if samrecord.qname.endswith('.left'):
                    left_hits.append(samrecord)
                else:
                    right_hits.append(samrecord)
            return (left_hits, right_hits, samfile)

        samfiles = [
            'gap_test_update_hits.multiple.sam',
            'gap_test_update_hits.diff_seqs.sam',
            'gap_test_update_hits.same_seq_diff_strand.sam',
            'gap_test_update_hits.same_seq_diff_strand.2.sam',
            'gap_test_update_hits.same_seq_strand.sam',
            'gap_test_update_hits.same_seq_strand.2.sam',
            'gap_test_update_hits.same_seq_strand_overlap.sam',
            'gap_test_update_hits.same_seq_strand_wrong_order.sam',
            'gap_test_update_hits.same_seq_strand_wrong_order.2.sam',
            'gap_test_update_hits.unmapped.sam',
            'gap_test_update_hits.unmapped.2.sam',
        ]

        samfiles = [os.path.join(data_dir, x) for x in samfiles]
        expected_gap_types = [
            gap.MULTIPLE_HITS,
            gap.DIFF_SEQS,
            gap.SAME_SEQ_DIFF_STRAND,
            gap.SAME_SEQ_DIFF_STRAND,
            gap.SAME_SEQ_STRAND,
            gap.SAME_SEQ_STRAND,
            gap.SAME_SEQ_STRAND_OVERLAP,
            gap.SAME_SEQ_STRAND_WRONG_ORDER,
            gap.SAME_SEQ_STRAND_WRONG_ORDER,
            gap.UNMAPPED,
            gap.UNMAPPED,
        ]

        assert len(samfiles) == len(expected_gap_types)

        for i in range(len(samfiles)):
            test_gap = gap.Gap('\t'.join(['query_name', '10', '42', 'LEFTSEQ', 'RIGHTSEQ' ]))
            left, right, samfile = load_sam(samfiles[i])
            test_gap.update_hits(left, right, samfile)
            self.assertEqual(test_gap.gap_type, expected_gap_types[i])

        
    def test_left_or_right_fasta_name_prefix(self):
        '''Test _left_or_right_fasta_name_prefix()'''
        test_gap = gap.Gap('\t'.join(['query_name', '10', '42', 'LEFTSEQ', 'RIGHTSEQ' ]))
        self.assertEqual(test_gap._left_or_right_fasta_name_prefix(), 'query_name:10-42')


    def test_left_right_fasta(self):
        '''Test left_fasta() and right_fasta()'''
        test_gap = gap.Gap('\t'.join(['query_name', '10', '42', 'LEFTSEQ', 'RIGHTSEQ' ]))
        self.assertEqual(test_gap.left_fasta(), sequences.Fasta('query_name:10-42.left', 'LEFTSEQ'))
        self.assertEqual(test_gap.right_fasta(), sequences.Fasta('query_name:10-42.right', 'RIGHTSEQ'))


    def test_can_be_filled(self):
        '''Test can_be_filled()'''
        test_gap = gap.Gap('\t'.join(['query_name', '10', '42', 'LEFTSEQ', 'RIGHTSEQ' ]))
        gap_types = [x for x in gap.type_to_str.values() if x != 'SAME_SEQ_STRAND']
        for t in gap_types:
            test_gap.gap_type = t
            self.assertFalse(test_gap.can_be_filled())
        
        test_gap.gap_type = gap.SAME_SEQ_STRAND
        test_gap.ref_name = 'ref'
        test_gap.ref_start = 100
        test_gap.ref_end = 140
        self.assertTrue(test_gap.can_be_filled())
        self.assertFalse(test_gap.can_be_filled(abs_diff=0))
        self.assertFalse(test_gap.can_be_filled(relative_err=0))
        test_gap.gap_type = gap.SAME_SEQ_STRAND_OVERLAP
        self.assertTrue(test_gap.can_be_filled())


if __name__ == '__main__':
    unittest.main()
