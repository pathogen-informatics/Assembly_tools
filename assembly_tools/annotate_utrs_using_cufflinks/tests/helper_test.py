#!/usr/bin/env python3

import sys
import filecmp
import os
import unittest
import copy
from assembly_tools.annotate_utrs_using_cufflinks import *
from assembly_tools.file_readers import gff
from pyfastaq import intervals

modules_dir = os.path.dirname(os.path.abspath(helper.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class Error (Exception): pass

class Test_helper(unittest.TestCase):
    def test_read_gff(self):
        level_records, other_records = helper.read_gff(os.path.join(data_dir, 'test_read_gff.gff'))
        expected_level_records = [
            [gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'gene', '42', '100', '.', '+', '.', 'ID=id']))],
            [gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'mRNA', '42', '100', '.', '+', '.', 'ID=id.1;Parent=id']))],
            [gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'exon', '52', '62', '.', '+', '.', 'ID=id.1:exon:1;Parent=id.1']))]
        ]
        expected_other_records = {
            'seqname': [gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'eggs', '150', '200', '.', '+', '.', 'ID=spam']))]
        }

        self.assertEqual(expected_level_records, level_records)
        self.assertEqual(expected_other_records, other_records)
 
    def test_update_other_records(self):
        other_records = {}
        gff1 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'eggs', '150', '200', '.', '+', '.', 'ID=spam']))
        gff2 = gff.GFF_record('\t'.join(['seqname', 'SOURCE', 'spam', '150', '200', '.', '+', '.', 'ID=eggs']))
        expected = {}
        
        helper.update_other_records(other_records, gff1)
        expected['seqname'] = [copy.deepcopy(gff1)]
        self.assertEqual(other_records, expected)

        helper.update_other_records(other_records, gff2)
        expected['seqname'].append(copy.deepcopy(gff2))
        self.assertEqual(other_records, expected)
        
    def test_load_ref_gff(self):
        #level_records, other_records = helper.read_gff(os.path.join(data_dir, 'test_get_genes_from_ref.gff'))
        genes, other_records = helper.load_ref_gff(os.path.join(data_dir, 'test_get_genes_from_ref.gff'))
        gene1_gene = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '42', '100', '.', '+', '.', 'ID=gene']))
        gene1_mRNA = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'mRNA', '42', '100', '.', '+', '.', 'ID=gene.1;Parent=gene']))
        gene1_exon = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'exon', '42', '62', '.', '+', '.', 'ID=gene.1:exon:1;Parent=gene.1']))
        gene1_exon2 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'exon', '92', '100', '.', '+', '.', 'ID=gene.1:exon:2;Parent=gene.1']))
        gene1 = gene.Gene(gene1_gene)
        gene1.add_gff_record(gene1_mRNA)
        gene1.add_gff_record(gene1_exon)
        gene1.add_gff_record(gene1_exon2)
       
        gene2_gene = gff.GFF_record('\t'.join(['seq2', 'SOURCE', 'gene', '1', '10', '.', '+', '.', 'ID=gene2']))
        gene2_mRNA = gff.GFF_record('\t'.join(['seq2', 'SOURCE', 'mRNA', '1', '10', '.', '+', '.', 'ID=gene2.1;Parent=gene2']))
        gene2_exon = gff.GFF_record('\t'.join(['seq2', 'SOURCE', 'exon', '1', '10', '.', '+', '.', 'ID=gene2.1:exon:1;Parent=gene2.1']))
        gene2 = gene.Gene(gene2_gene)
        gene2.add_gff_record(gene2_mRNA)
        gene2.add_gff_record(gene2_exon)

        expected_genes = {'seq': [gene1], 'seq2': [gene2]}
        expected_other = {'seq': [gff.GFF_record('\t'.join(['seq', 'SOURCE', 'eggs', '150', '200', '.', '+', '.', 'ID=spam']))]}
        self.assertEqual(expected_genes, genes)
        self.assertEqual(expected_other, other_records)
        

    def test_load_cufflinks_gtf(self):
        genes = helper.load_cufflinks_gtf(os.path.join(data_dir, 'test_load_cufflinks_gtf.gtf'))
        gene1_trans = gff.GFF_record('\t'.join(['seq', 'Cufflinks', 'transcript', '1', '100', '1000', '+', '.', 'gene_id "CUFF.1"; transcript_id "CUFF.1.1";']))
        gene1_exon = gff.GFF_record('\t'.join(['seq', 'Cufflinks', 'exon', '1', '100', '1000', '+', '.', 'gene_id "CUFF.1"; transcript_id "CUFF.1.1";']))
        gene2_trans = gff.GFF_record('\t'.join(['seq2', 'Cufflinks', 'transcript', '1', '200', '2000', '+', '.', 'gene_id "CUFF.2"; transcript_id "CUFF.2.1";']))
        gene2_exon1 = gff.GFF_record('\t'.join(['seq2', 'Cufflinks', 'exon', '1', '100', '2000', '+', '.', 'gene_id "CUFF.2"; transcript_id "CUFF.2.1";']))
        gene2_exon2 = gff.GFF_record('\t'.join(['seq2', 'Cufflinks', 'exon', '150', '200', '2000', '+', '.', 'gene_id "CUFF.2"; transcript_id "CUFF.2.1";']))
        gene1 = gene.Gene(gene1_trans)
        gene1.add_gff_record(gene1_exon)
        gene2 = gene.Gene(gene2_trans)
        gene2.add_gff_record(gene2_exon1)
        gene2.add_gff_record(gene2_exon2)
        expected_genes = {'seq': [gene1], 'seq2': [gene2]}
        self.assertEqual(expected_genes['seq2'], genes['seq2'])

        with self.assertRaises(helper.Error):
            helper.load_cufflinks_gtf(os.path.join(data_dir, 'test_load_cufflinks_gtf.no_parent.gtf'))

