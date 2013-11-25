#!/usr/bin/env python3

import argparse
import sys
import os
import pysam
import re
from fastaq import *
import assembly_tools.fill_gaps_using_reference


class Error (Exception): pass

def make_fasta_of_gap_flanks(fasta_in, flanking_bases, fasta_out, gaps):
    '''Makes a fasta file of the sequences flanking the gaps in a fasta/q file'''
    tmp_tabfile = fasta_out + '.tmp.tab'
    tasks.get_seqs_flanking_gaps(fasta_in, tmp_tabfile, flanking_bases, flanking_bases)
    fin = utils.open_file_read(tmp_tabfile)
    fout = utils.open_file_write(fasta_out)
    original_line_length = sequences.Fasta.line_length
    sequences.Fasta.line_length = 0

    for line in fin:
        if line.startswith('#'):
            continue

        gap = assembly_tools.fill_gaps_using_reference.gap.Gap(line)
        print(gap.left_fasta(), file=fout)
        print(gap.right_fasta(), file=fout)
        if gap.query_name not in gaps:
            gaps[gap.query_name] = {}
        gaps[gap.query_name][(gap.query_start, gap.query_end)] = gap


    utils.close(fin)
    utils.close(fout)
    os.unlink(tmp_tabfile)
    sequences.Fasta.line_length = original_line_length


def paired_hit_samreader(filename):
    '''Given a SAM file in read name order, yields a tuple of hits ([left_hits], [right_hits])'''
    samfile = pysam.Samfile(filename, "r")
    left_hits = []
    right_hits = []

    for samrecord in samfile.fetch(until_eof=True):
        if len(left_hits) == 0:
            if not samrecord.qname.endswith('.left'):
                raise Error('Expecting to get a "left" read in SAM but got this:' + samrecord.qname)
                sys.exit(1)
            left_hits.append(samrecord)
        elif len(right_hits) > 0 and samrecord.qname.endswith('.left'):
            yield left_hits, right_hits, samfile
            left_hits = [samrecord]
            right_hits = []
        elif samrecord.qname.endswith('.left'):
            left_hits.append(samrecord)
        elif samrecord.qname.endswith('.right'):
            right_hits.append(samrecord)
        else:
            raise Error('Unexpected error parsing SAM file. Cannot continue')
            sys.exit(1)

    yield left_hits, right_hits, samfile


def _gap_flank_seqname_to_dict_key(name):
    regex = re.compile('^(.*):(\d+)-(\d+)\.(?:left|right)')
    hits = regex.search(name)

    try:
        name = hits.group(1)
        gap_start = int(hits.group(2)) - 1
        gap_end = int(hits.group(3)) - 1
    except:
        raise Error('Error getting gap start/end coords from sequence with name ' + name)

    return name, (gap_start, gap_end)


def parse_sam_file(samfilename, gaps):
    samreader = paired_hit_samreader(samfilename)

    for left_hits, right_hits, samfile in samreader:
        # find gap corresponding to this pair of reads
        qry, coords = _gap_flank_seqname_to_dict_key(left_hits[0].qname)
        try:
            gaps[qry][coords].update_hits(left_hits, right_hits, samfile)
        except:
            raise Error('Error parsing line of SAM ' + left_hits[0])

