#!/usr/bin/env python3

import sys
import filecmp
import os
import unittest
import copy
from assembly_tools.annotate_utrs_using_cufflinks import *
from assembly_tools.file_readers import gff
from fastaq import intervals

class Error (Exception): pass

class Test_transcript(unittest.TestCase):
    def setUp(self):
        self.gff_mRNA = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'mRNA', '42', '100', '.', '+', '.', 'ID=gene_id.1;Parent=gene_id']))
        self.gff_five_utr = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'five_prime_UTR', '42', '43', '.', '+', '.', 'ID=gene_id.1:5utr;Parent=gene_id.1']))
        self.gff_three_utr = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'three_prime_UTR', '95', '100', '.', '+', '.', 'ID=gene_id.1:3utr;Parent=gene_id.1']))
        self.gff_exon = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'CDS', '44', '53', '.', '+', '.', 'ID=gene_id.1:exon:1;Parent=gene_id.1']))
        self.gff_exon2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '60', '63', '.', '+', '.', 'ID=gene_id.1:exon:2;Parent=gene_id.1']))
        self.gff_exon3 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '90', '94', '.', '+', '.', 'ID=gene_id.1:exon:2;Parent=gene_id.1']))
        self.gff_pseudogenic_exon = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'pseudogenic_exon', '60', '63', '.', '+', '.', 'ID=gene_id.1:exon:3;Parent=gene_id.1']))
        self.gff_transcript = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'transcript', '42', '100', '.', '+', '.', 'ID=gene_id.1;Parent=gene_id']))
        self.gff_pseudogenic_transcript = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'pseudogenic_transcript', '42', '100', '.', '+', '.', 'ID=gene_id.1;Parent=gene_id']))
        self.gff_ncRNA = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'ncRNA', '42', '43', '.', '+', '.', 'ID=gene_id.1:ncRNA;Parent=gene_id.1']))
        self.gff_tRNA = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'tRNA', '42', '43', '.', '+', '.', 'ID=gene_id.1:tRNA;Parent=gene_id.1']))
        self.gff_snRNA = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'snRNA', '42', '43', '.', '+', '.', 'ID=gene_id.1:snRNA;Parent=gene_id.1']))
        self.gff_rRNA = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'rRNA', '42', '43', '.', '+', '.', 'ID=gene_id.1:rRNA;Parent=gene_id.1']))
        self.gff_other = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'eggs', '42', '4242', '.', '+', '.']))

        self.trans = transcript.Transcript(self.gff_mRNA)
 
    def test_init(self):
        '''Test __init__'''
        coords = intervals.Interval(42, 100)
        strand = '+'
        seqname = 'seqname'
        self.assertEqual(self.trans.coords, coords)
        self.assertEqual(self.trans.strand, strand)
        self.assertEqual(self.trans.seqname, seqname)

        for l in [self.trans.five_utr, self.trans.three_utr, self.trans.exons, self.trans.ncRNA, self.trans.rRNA, self.trans.tRNA, self.trans.snRNA]:
            self.assertEqual(len(l), 0)

        self.assertEqual(self.trans.mRNA, self.gff_mRNA)


    def test_add_gff_record(self):
        '''Test add_gff_record'''
        self.trans.add_gff_record(self.gff_five_utr)
        self.assertEqual(self.trans.five_utr, [self.gff_five_utr])

        self.trans.add_gff_record(self.gff_three_utr)
        self.assertEqual(self.trans.three_utr, [self.gff_three_utr])

        self.trans.add_gff_record(self.gff_exon)
        self.assertEqual(self.trans.exons, [self.gff_exon])

        self.trans.add_gff_record(self.gff_exon2)
        self.assertEqual(self.trans.exons, [self.gff_exon, self.gff_exon2])

        self.trans.add_gff_record(self.gff_pseudogenic_exon)
        self.assertEqual(self.trans.exons, [self.gff_exon, self.gff_exon2, self.gff_pseudogenic_exon])

        with self.assertRaises(transcript.Error):
            self.trans.add_gff_record(self.gff_transcript)
        with self.assertRaises(transcript.Error):
            self.trans.add_gff_record(self.gff_pseudogenic_transcript)

        self.trans.mRNA = None
        self.trans.add_gff_record(self.gff_transcript)
        self.assertEqual(self.trans.mRNA, self.gff_transcript)

        self.trans.mRNA = None
        self.trans.add_gff_record(self.gff_pseudogenic_transcript)
        self.assertEqual(self.trans.mRNA, self.gff_pseudogenic_transcript)

        self.trans.mRNA = None
        self.trans.add_gff_record(self.gff_mRNA)
        self.assertEqual(self.trans.mRNA, self.gff_mRNA)

        self.trans.add_gff_record(self.gff_ncRNA)
        self.assertEqual(self.trans.ncRNA, [self.gff_ncRNA])

        self.trans.add_gff_record(self.gff_tRNA)
        self.assertEqual(self.trans.tRNA, [self.gff_tRNA])

        self.trans.add_gff_record(self.gff_snRNA)
        self.assertEqual(self.trans.snRNA, [self.gff_snRNA])

        self.trans.add_gff_record(self.gff_rRNA)
        self.assertEqual(self.trans.rRNA, [self.gff_rRNA])

        self.trans.add_gff_record(self.gff_other)
        self.assertEqual(self.trans.other_gffs, [self.gff_other])

    def test_sort(self):
        '''Test sort'''
        unsorted_list = [
            gff.GFF_record('\t'.join(['x', 'x', 'x', '42', '43', '.', '.', '.'])),
            gff.GFF_record('\t'.join(['x', 'x', 'x', '12', '13', '.', '.', '.']))
        ]

        sorted_list = list(unsorted_list)
        sorted_list.sort()

        for l in [self.trans.five_utr, self.trans.three_utr, self.trans.exons, self.trans.ncRNA ,self.trans.rRNA, self.trans.tRNA, self.trans.snRNA]:
            l += unsorted_list

        self.trans._sort()
        for l in [self.trans.five_utr, self.trans.three_utr, self.trans.exons, self.trans.ncRNA ,self.trans.rRNA, self.trans.tRNA, self.trans.snRNA]:
            self.assertEqual(sorted_list, l)
        

    def test_set_coords(self):
        '''Test set_coords'''
        self.trans.coords = None
        self.trans._set_coords()
        self.assertEqual(self.trans.coords, intervals.Interval(42, 100))

        self.trans.add_gff_record(self.gff_exon)
        self.trans.coords = None
        self.trans._set_coords()
        self.assertEqual(self.trans.coords, intervals.Interval(44, 53))
        

    def test_total_exon_length(self):
        '''Test total_exon_length'''
        self.assertEqual(self.trans.total_exon_length(), 0)
        self.trans.add_gff_record(self.gff_exon)
        self.assertEqual(self.trans.total_exon_length(), 10)
        self.trans.add_gff_record(self.gff_exon2)
        self.assertEqual(self.trans.total_exon_length(), 14)


    def test_set_seqname(self):
        '''Test get_seqname'''
        self.trans.seqname = None
        self.trans._set_seqname()
        self.assertEqual(self.trans.seqname, 'seqname')

        self.trans.seqname = 'something else'
        with self.assertRaises(transcript.Error):
            self.trans._set_seqname()

        with self.assertRaises(transcript.Error):
            self.trans.add_gff_record(self.gff_exon)

        self.trans.seqname = 'seqname'
        self.trans.add_gff_record(self.gff_exon)
        self.trans.add_gff_record(self.gff_exon)
        self.trans.exons[0].seqname = 'something else'
        with self.assertRaises(transcript.Error):
            self.trans._set_seqname()
        
        
    def test_set_strand(self):
        '''Test set_strand'''
        self.trans.strand = None
        self.trans._set_strand()
        self.assertEqual(self.trans.strand, '+')

        self.trans.strand = '-'
        with self.assertRaises(transcript.Error):
            self.trans._set_strand()

        self.trans.strand = '+'
        self.trans.add_gff_record(self.gff_exon)
        self.trans.add_gff_record(self.gff_exon)
        self.trans.exons[0].strand = '-'
        with self.assertRaises(transcript.Error):
            self.trans._set_strand()
        
        transcript.lenient = True
        self.trans._set_strand()
        self.assertEqual(self.trans.strand, 'Inconsistent')
        transcript.lenient = False

    def test_intersects(self):
        '''Test intersects'''
        trans2 = copy.deepcopy(self.trans)
        not_intersects = [
            intervals.Interval(1,41),
            intervals.Interval(101,141),
        ]

        intersects = [
            intervals.Interval(1,42),
            intervals.Interval(42,50),
            intervals.Interval(50,60),
            intervals.Interval(50,100),
            intervals.Interval(100,142),
            intervals.Interval(20,424242),
        ]

        for i in not_intersects:
            trans2.coords = i
            self.assertFalse(self.trans.intersects(trans2))

        for i in intersects:
            trans2.coords = i
            self.assertTrue(self.trans.intersects(trans2))

    def test_might_extend(self):
        '''Test might_extend'''
        self.trans.add_gff_record(self.gff_exon)
        trans2 = copy.deepcopy(self.trans)
        self.assertFalse(self.trans.might_extend(trans2))
        trans2.coords = intervals.Interval(10, 150)
        self.assertTrue(self.trans.might_extend(trans2))

        trans2.seqname = 'other'
        self.assertFalse(self.trans.might_extend(trans2))
        trans2.seqname = 'seqname'

        trans2.strand = '-'
        self.assertFalse(self.trans.might_extend(trans2))
        trans2.strand = '+'

        trans2.coords = intervals.Interval(48, 50)
        self.assertFalse(self.trans.might_extend(trans2))
        self.assertTrue(trans2.might_extend(self.trans))

        might_extend = [
            intervals.Interval(1,41),
            intervals.Interval(41,42),
            intervals.Interval(1,42),
            intervals.Interval(1,50),
            intervals.Interval(1,150),
            intervals.Interval(90,142),
            intervals.Interval(100,142),
            intervals.Interval(101,142),
        ]

        self.trans.coords = intervals.Interval(42, 100)

        for i in might_extend:
            trans2.coords = i
            self.assertTrue(self.trans.might_extend(trans2))

        trans2.coords = intervals.Interval(40,42)
        self.assertTrue(self.trans.might_extend(trans2))
        self.assertTrue(self.trans.might_extend(trans2, min_extend=2))
        self.assertFalse(self.trans.might_extend(trans2, min_extend=3))
        self.assertFalse(self.trans.might_extend(trans2, min_extend=4))

        trans2.coords = intervals.Interval(50,102)
        self.assertTrue(self.trans.might_extend(trans2))
        self.assertTrue(self.trans.might_extend(trans2, min_extend=2))
        self.assertFalse(self.trans.might_extend(trans2, min_extend=3))
        self.assertFalse(self.trans.might_extend(trans2, min_extend=4))


    def test_can_extend_start(self):
        '''Test can_extend_start()'''
        self.trans.add_gff_record(self.gff_exon)
        trans2 = copy.deepcopy(self.trans)
        trans2.coords = intervals.Interval(42,44)
        self.assertTrue(self.trans.can_extend_start(trans2))
        self.assertTrue(self.trans.can_extend_start(trans2, min_extend=2))
        self.assertFalse(self.trans.can_extend_start(trans2, min_extend=3))
        self.assertFalse(self.trans.can_extend_start(trans2, min_extend=4))
        trans2.coords = intervals.Interval(44,50)
        self.assertFalse(self.trans.can_extend_start(trans2))

    def test_can_extend_end(self):
        '''Test can_extend_end()'''
        self.trans.add_gff_record(self.gff_exon)
        trans2 = copy.deepcopy(self.trans)
        trans2.coords = intervals.Interval(50,55)
        self.assertTrue(self.trans.can_extend_end(trans2))
        self.assertTrue(self.trans.can_extend_end(trans2, min_extend=2))
        self.assertFalse(self.trans.can_extend_end(trans2, min_extend=3))
        self.assertFalse(self.trans.can_extend_end(trans2, min_extend=4))

    def test_update_utrs(self):
        '''Test update_utrs()'''
        self.trans.add_gff_record(self.gff_exon)
        self.trans.add_gff_record(self.gff_exon2)
        self.trans.add_gff_record(self.gff_exon3)

        new_left_exon = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'CDS', '30', '50', '.', '+', '.']))
        new_left_exon2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'CDS', '10', '20', '.', '+', '.']))
        new_right_exon = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'CDS', '90', '100', '.', '+', '.']))
        new_right_exon2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'CDS', '110', '120', '.', '+', '.']))
        new_five_utr = gff.GFF_record('\t'.join(['seqname', 'UTR_updater', 'five_prime_UTR', '30', '43', '.', '+', '.', 'ID=gene_id.1:5utr;Parent=gene_id.1']))
        new_five_utr_2 = gff.GFF_record('\t'.join(['seqname', 'UTR_updater', 'five_prime_UTR', '30', '43', '.', '+', '.', 'ID=gene_id.1:5utr:2;Parent=gene_id.1']))
        new_five_utr_1 = gff.GFF_record('\t'.join(['seqname', 'UTR_updater', 'five_prime_UTR', '10', '20', '.', '+', '.', 'ID=gene_id.1:5utr:1;Parent=gene_id.1']))
        new_three_utr = gff.GFF_record('\t'.join(['seqname', 'UTR_updater', 'three_prime_UTR', '95', '100', '.', '+', '.', 'ID=gene_id.1:3utr;Parent=gene_id.1']))
        new_three_utr_1 = gff.GFF_record('\t'.join(['seqname', 'UTR_updater', 'three_prime_UTR', '95', '100', '.', '+', '.', 'ID=gene_id.1:3utr:1;Parent=gene_id.1']))
        new_three_utr_2 = gff.GFF_record('\t'.join(['seqname', 'UTR_updater', 'three_prime_UTR', '110', '120', '.', '+', '.', 'ID=gene_id.1:3utr:2;Parent=gene_id.1']))

        before_adding = [copy.deepcopy(self.trans) for x in range(5)]
        after_adding = [copy.deepcopy(self.trans) for x in range(5)]
        to_add = [copy.deepcopy(self.trans) for x in range(5)]

        # test extending the first exon
        to_add[0].exons.pop(0)
        to_add[0].add_gff_record(new_left_exon)
        after_adding[0].exons[0].coords.start = 44
        after_adding[0].add_gff_record(new_five_utr)

        # test adding a spliced UTR at the start
        to_add[1].exons.pop(0)
        to_add[1].add_gff_record(new_left_exon)
        to_add[1].add_gff_record(new_left_exon2)
        after_adding[1].exons[0].coords.start = 44
        after_adding[1].add_gff_record(new_five_utr_1)
        after_adding[1].add_gff_record(new_five_utr_2)

        # test extending the last exon
        to_add[2].exons.pop()
        to_add[2].add_gff_record(new_right_exon)
        after_adding[2].exons[-1].coords.end = 94
        after_adding[2].add_gff_record(new_three_utr)

        # test adding a spliced UTR at the end
        to_add[3].exons.pop()
        to_add[3].add_gff_record(new_right_exon)
        to_add[3].add_gff_record(new_right_exon2)
        after_adding[3].exons[-1].coords.end = 94
        after_adding[3].add_gff_record(new_three_utr_1)
        after_adding[3].add_gff_record(new_three_utr_2)
        
        # test trying to add too many new uts splices
        to_add[4].add_gff_record(new_left_exon)
        to_add[4].add_gff_record(new_left_exon2)
        to_add[4].add_gff_record(new_left_exon2)
        to_add[4].add_gff_record(new_left_exon2)

        # test adding a 5'UTR that overlaps an existing gene. Should do nothing
        t = copy.deepcopy(before_adding[0])
        t.update_utrs(to_add[0], exclude_coords = [intervals.Interval(10,31)])
        self.assertEqual(t, before_adding[0])
        
        # test adding a 3'UTR that overlaps an existing gene. Should do nothing
        t = copy.deepcopy(before_adding[2])
        t.update_utrs(to_add[2], exclude_coords = [intervals.Interval(100,110)])
        self.assertEqual(t, before_adding[2])

        for i in range(len(before_adding)):
            before_adding[i].update_utrs(to_add[i])
            self.assertEqual(before_adding[i], after_adding[i])

    def test_exon_splice_sites(self):
        self.assertEqual(self.trans.exon_splice_sites(), [])
        self.trans.add_gff_record(self.gff_exon)
        self.assertEqual(self.trans.exon_splice_sites(), [])
        self.trans.add_gff_record(self.gff_exon2)
        self.assertEqual(self.trans.exon_splice_sites(), [53, 60])
        self.trans.add_gff_record(self.gff_exon3)
        self.assertEqual(self.trans.exon_splice_sites(), [53, 60, 63, 90])
            
    def test_number_of_common_splice_site(self):
        trans2 = copy.deepcopy(self.trans)
        self.assertEqual(self.trans.number_of_common_splice_sites(trans2), 0)
        self.trans.add_gff_record(copy.deepcopy(self.gff_exon))
        trans2.add_gff_record(copy.deepcopy(self.gff_exon))
        self.assertEqual(self.trans.number_of_common_splice_sites(trans2), 0)
        self.trans.add_gff_record(copy.deepcopy(self.gff_exon2))
        trans2.add_gff_record(copy.deepcopy(self.gff_exon2))
        self.assertEqual(self.trans.number_of_common_splice_sites(trans2), 2)
        trans2.exons[-1].coords.start += 1
        self.assertEqual(self.trans.number_of_common_splice_sites(trans2), 1)
