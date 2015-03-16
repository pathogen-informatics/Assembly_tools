#!/usr/bin/env python3

import sys
import filecmp
import os
import unittest
import assembly_tools.file_readers.gff as gff
from pyfastaq import utils

modules_dir = os.path.dirname(os.path.abspath(gff.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class Error (Exception): pass

class Test_GFF_record(unittest.TestCase):
    def test_init_good_input(self):
        '''Test __init__ for good input'''
        test_input = [
            ['seq', 'SOURCE', 'gene', '42', '43', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '42', '43', '1.4', '-', '0'],
            ['seq', 'SOURCE', 'gene', '42', '43', '1.4', '+', '1'],
            ['seq', 'SOURCE', 'gene', '42', '43', '1.4', '+', '1', 'key=value'],
            ['seq', 'SOURCE', 'gene', '42', '43', '1.4', '+', '1', 'key=value;key2=value 2'],
            ['seq', 'Cufflinks', 'transcript', '42', '43', '.', '.', '.', 'key1 "val1";'],
            ['seq', 'Cufflinks', 'transcript', '42', '43', '.', '.', '.', 'key1 "val1";'],
            ['seq', 'Cufflinks', 'transcript', '42', '43', '.', '.', '.', 'key1 "val1"; key2 "val2";'],
        ]

        for l in test_input:
            record = gff.GFF_record('\t'.join(l))
            self.assertEqual('\t'.join(l), str(record))
   
    def test_init_bad_input(self):
        '''Test __init__ for bad input'''
        test_input = [
            ['seq', 'SOURCE', 'gene', '42', '43', '.', '.'],
            ['seq', 'SOURCE', 'gene', '42', '43', '.', '.', '.', '.', 'too many'],
            ['seq', 'SOURCE', 'gene', 'not_int', '43', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '42', 'not_int', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '42', '43', 'not_float', '.', '.'],
            ['seq', 'SOURCE', 'gene', '42', '43', '.', 'bad_strand', '.'],
            ['seq', 'SOURCE', 'gene', '42', '43', '.', '.', '.', 'attribute with no equals sign'],
        ]


        for l in test_input:
            with self.assertRaises(gff.Error):
                gff.GFF_record('\t'.join(l))

    def test_lenient_mode(self):
        '''Test __init__ when being lenient'''
        gff.lenient = True
        test_input = [
            ['seq', 'Cufflinks', 'transcript', '42', '43', '.', '.', '.', 'key1 "val1"; key2 "val2";'],
            ['seq', 'SOURCE', 'gene', '42', '43', '.', '.', '.', 'attribute with no equals sign'],
        ]

        for l in test_input:
            record = gff.GFF_record('\t'.join(l))
            self.assertEqual('\t'.join(l), str(record))

        gff.lenient = False


    def test_len(self):
        '''Test __len__'''
        g = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '42', '44', '.', '.', '.']))
        self.assertEqual(3, len(g))


    def test_less_than(self):
        '''Test less than operator'''
        gff_1 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '42', '44', '.', '.', '.']))
        gff_2 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '42', '43', '.', '.', '.']))
        gff_3 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '42', '45', '.', '.', '.']))
        gff_4 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '41', '42', '.', '.', '.']))
        gff_5 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '41', '45', '.', '.', '.']))
        gff_6 = gff.GFF_record('\t'.join(['different_seq', 'SOURCE', 'gene', '42', '45', '.', '.', '.']))

        self.assertFalse(gff_1 < gff_1)
        self.assertFalse(gff_1 < gff_2)
        self.assertTrue(gff_1 < gff_3)
        self.assertFalse(gff_1 < gff_4)
        self.assertFalse(gff_1 < gff_5)
        self.assertFalse(gff_1 < gff_6)

    def test_interscts(self):
        '''Test instersects'''
        gff_1 = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'gene', '42', '45', '.', '.', '.']))

        intersects = [
            ['seq', 'SOURCE', 'gene', '43', '44', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '42', '43', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '43', '45', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '41', '42', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '41', '43', '.', '.', '.']
        ]
        
        for l in intersects:
            record = gff.GFF_record('\t'.join(l))
            self.assertTrue(gff_1.intersects(record))

        not_intersects = [
            ['seq', 'SOURCE', 'gene', '40', '41', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '46', '50', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '1', '10', '.', '.', '.'],
            ['seq', 'SOURCE', 'gene', '100', '200', '.', '.', '.'],
            ['different_seq', 'SOURCE', 'gene', '43', '45', '.', '.', '.'],
        ]

        for l in not_intersects:
            record = gff.GFF_record('\t'.join(l))
            self.assertFalse(gff_1.intersects(record))

    def test_get_attribute(self):
        '''Test get_attribute'''
        gff.lenient = True
        gff_record = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'feature', '42', '43', '.', '.', '.', 'key1=val1;key2=val2;key3']))
        
        self.assertEqual(gff_record.get_attribute('key1'), 'val1')
        self.assertEqual(gff_record.get_attribute('key2'), 'val2')
        self.assertEqual(gff_record.get_attribute('key3'), None)

        with self.assertRaises(gff.Error):
            gff_record.get_attribute('killer rabbit')

        gff.lenient = False

    def test_set_attribute(self):
        '''Test get_attribute'''
        gff_record = gff.GFF_record('\t'.join(['seq', 'SOURCE', 'feature', '42', '43', '.', '.', '.']))
        gff_record.set_attribute('key1', '42')
        self.assertTrue(gff_record.get_attribute('key1'), '42')
        gff_record.set_attribute('key1', '43')
        self.assertTrue(gff_record.get_attribute('key1'), '43')

class Test_file_reader(unittest.TestCase):
    def test_read_write_gff_file(self):
        '''Test can read and write gff file OK'''
        outfile = 'tmp.out.gff'
        infiles = [
            os.path.join(data_dir, 'gff_io_test.gff'),
            os.path.join(data_dir, 'gff_io_test.cufflinks.gtf')
        ]

        for fname in infiles:
            f = utils.open_file_write(outfile)
            reader = gff.file_reader(fname)
            for record in reader:
                print(record, file=f)
            utils.close(f)
            self.assertTrue(filecmp.cmp(outfile, fname + '.out.gff'))

        os.unlink(outfile)
