#!/usr/bin/env python3

import sys
import filecmp
import os
import unittest
import copy
from assembly_tools.annotate_utrs_using_cufflinks import *
from assembly_tools.file_readers import gff
from pyfastaq import intervals

modules_dir = os.path.dirname(os.path.abspath(gene.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class Error (Exception): pass

class Test_gene(unittest.TestCase):
    def setUp(self):
        self.gff_mRNA = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'mRNA', '42', '100', '.', '+', '.', 'ID=gene_id.1;Parent=gene_id']))
        self.gff_gene = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'gene', '42', '100', '.', '+', '.', 'ID=gene_id']))
        self.gff_exon1 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '50', '60', '.', '+', '.', 'ID=gene_id.1:exon:1;Parent=gene_id.1']))
        self.gene = gene.Gene(self.gff_mRNA)
        self.gene.add_gff_record(self.gff_gene)
        self.gene.add_gff_record(self.gff_exon1)

    def test_set_seqname(self):
        self.gene.seqname = None
        self.gene._set_seqname()
        self.assertEqual(self.gene.seqname, 'seqname')
        self.gene.transcripts['gene_id.1'].seqname = 'something else'
        with self.assertRaises(gene.Error):
            self.gene._set_seqname()
     
    def test_set_strand(self):
        self.gene.strand = None
        self.gene._set_strand(self.gff_mRNA)
        self.assertEqual(self.gene.strand, '+')
        self.gff_mRNA.strand = '-'
        with self.assertRaises(gene.Error):
            self.gene._set_strand(self.gff_mRNA)
        
        gene.lenient = True
        self.gene._set_strand(self.gff_mRNA)
        self.assertTrue(self.gene.strand, 'Inconsistent')
        gene.lenient = False

    def test_set_coords(self):
        self.gene.coords = None
        self.gene._set_coords()
        self.assertEqual(self.gene.coords, self.gff_exon1.coords)
        
    def test_longest_transcript_by_exon_length(self):
        '''Test longest_transcript_by_exon_length'''
        self.assertEqual('gene_id.1', self.gene.longest_transcript_by_exon_length())
        gff_mRNA2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'mRNA', '42', '100', '.', '+', '.', 'ID=gene_id.2;Parent=gene_id']))
        gff_exon2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '50', '65', '.', '+', '.', 'ID=gene_id.2:exon:1;Parent=gene_id.2']))
        self.gene.add_gff_record(gff_mRNA2)
        self.gene.add_gff_record(gff_exon2)
        self.assertEqual('gene_id.2', self.gene.longest_transcript_by_exon_length())


    def test_remove_all_but_longest_transcript(self):
        '''Test remove_all_but_longest_transcript'''
        self.gene.remove_all_but_longest_transcript()
        self.assertEqual(['gene_id.1'], list(self.gene.transcripts.keys()))
        gff_mRNA2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'mRNA', '42', '100', '.', '+', '.', 'ID=gene_id.2;Parent=gene_id']))
        gff_exon2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '50', '65', '.', '+', '.', 'ID=gene_id.2:exon:1;Parent=gene_id.2']))
        self.gene.add_gff_record(gff_mRNA2)
        self.gene.add_gff_record(gff_exon2)
        self.gene.remove_all_but_longest_transcript()
        self.assertEqual(['gene_id.2'], list(self.gene.transcripts.keys()))


    def test_lt(self):
        gene2 = copy.deepcopy(self.gene)
        gene2.seqname = 'seqname'
        self.assertFalse(self.gene < gene2)
        gene2.coords.start -= 5
        self.assertTrue(gene2 < self.gene)
        gene2.seqname = 'foo'
        self.assertFalse(gene2 < self.gene)
        gene2.seqname = 'seqname'
        gene2.coords.end -= 2
        self.assertTrue(gene2 < self.gene)
        gene2.coords.start = 1
        gene2.coords.end = 10
        self.assertTrue(gene2 < self.gene)
       
        
    def test_intersects(self):
        gene2 = copy.deepcopy(self.gene)
        gene2.seqname = 'seqname'
        self.assertTrue(self.gene.intersects(gene2))
        gene2.seqname = 'foo'
        self.assertFalse(self.gene.intersects(gene2))
        gene2.seqname = 'seqname'
        gene2.coords.start += 1
        self.assertTrue(self.gene.intersects(gene2))
        gene2.coords.end -= 1
        self.assertTrue(self.gene.intersects(gene2))
        gene2.coords.start -= 3
        self.assertTrue(self.gene.intersects(gene2))
        gene2.coords.end += 3
        self.assertTrue(self.gene.intersects(gene2))

    def test_can_extend(self):
        new_left_exon = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '40', '55', '.', '+', '.', 'ID=gene_id.1:exon:1;Parent=gene_id.1']))
        new_right_exon = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '55', '70', '.', '+', '.', 'ID=gene_id.1:exon:1;Parent=gene_id.1']))
        gene2 = copy.deepcopy(self.gene)
        self.assertFalse(self.gene.can_extend(gene2))
        gene2.add_gff_record(new_left_exon)
        self.assertTrue(self.gene.can_extend(gene2))
        gene2.seqname = 'holyhandgrenade'
        self.assertFalse(self.gene.can_extend(gene2))
        gene2 = copy.deepcopy(self.gene)
        gene2.add_gff_record(new_right_exon)
        self.assertTrue(self.gene.can_extend(gene2))
      
        
      
