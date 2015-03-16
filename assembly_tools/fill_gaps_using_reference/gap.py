#!/usr/bin/env python3

import sys
import pysam
from pyfastaq import sequences

class Error (Exception): pass

UNKNOWN = 0
MULTIPLE_HITS = 1
UNMAPPED = 2
DIFF_SEQS = 3
SAME_SEQ_DIFF_STRAND = 4
SAME_SEQ_STRAND_WRONG_ORDER = 5
SAME_SEQ_STRAND = 6
SAME_SEQ_STRAND_OVERLAP = 7

type_to_str = {
    0: 'UNKNOWN',
    1: 'MULTIPLE_HITS',
    2: 'UNMAPPED',
    3: 'DIFF_SEQS',
    4: 'SAME_SEQ_DIFF_STRAND',
    5: 'SAME_SEQ_STRAND_WRONG_ORDER',
    6: 'SAME_SEQ_STRAND',
    7: 'SAME_SEQ_STRAND_OVERLAP'
}

class Gap:
    def __init__(self, gap_flanks_line=None):
        self.left_hits = []
        self.right_hits = []
        self.ref_name = None
        self.ref_start = -1  # start pos of sequence in ref matching gap in query
        self.ref_end = -1    # end pos of sequence in ref matching gap in query
        self.gap_type = UNKNOWN
        self.query_name = None
        self.query_start = None  # start coord of gap in query
        self.query_end = None    # end coord of gap in query
        self.left_seq = None
        self.right_seq = None
        self.query_replace_start = -1
        self.query_replace_end = -1
        self.reverse_hit = False


        if gap_flanks_line is None:
            return

        try:
            (self.query_name,
             self.query_start,  # start coord of gap in query
             self.query_end,    # end coord of gap in query
             self.left_seq,
             self.right_seq) = gap_flanks_line.rstrip().split('\t')

            # store everything zero-based
            self.query_start = int(self.query_start) - 1
            self.query_end = int(self.query_end) - 1
        except:
            raise Error('Error parsing this line of flanking seqs file:\n' + gap_flanks_line)

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        l = [
            self.query_name,
            self.query_start + 1,
            self.query_end + 1,
            self.query_replace_start + 1,
            self.query_replace_end + 1,
            self.ref_name,
            self.ref_start + 1,
            self.ref_end + 1,
            {True:1,False:0}[self.reverse_hit],
            type_to_str[self.gap_type],
        ]

        for i in range(len(l)):
            if l[i] is None:
                l[i] = '*'
            else:
                l[i] = str(l[i])

        return '\t'.join(l)


    def update_hits(self, left_hits, right_hits, samfile):
        # sanity check that we only have reads flanking one gap
        left_names = set([x.qname for x in left_hits])
        right_names = set([x.qname for x in right_hits])
        assert(len(left_names) == 1)
        assert(len(right_names) == 1)
        assert(left_names.pop().rsplit('.')[0] == right_names.pop().rsplit('.')[0])

        self.left_hits = left_hits
        self.right_hits = right_hits

        # For now, only do something where we have one hit either side of gap and the
        # hit positions are in same region of reference and in correct orientation
        if (len(left_hits) > 1 and len(right_hits) > 1):
            
            self.gap_type = MULTIPLE_HITS
            return

        self.left_hit = left_hits[0]
        self.right_hit = right_hits[0]
        left_hit = left_hits[0]
        right_hit = right_hits[0]

        if self.left_hit.is_unmapped or self.right_hit.is_unmapped:
            self.gap_type = UNMAPPED
        elif self.left_hit.tid != self.right_hit.tid:
            self.gap_type = DIFF_SEQS
        elif self.left_hit.is_reverse != self.right_hit.is_reverse:
            self.gap_type = SAME_SEQ_DIFF_STRAND
        elif (self.left_hit.is_reverse and self.left_hit.pos < self.right_hit.pos) \
          or ((not self.left_hit.is_reverse) and self.left_hit.pos > self.right_hit.pos):
            self.gap_type = SAME_SEQ_STRAND_WRONG_ORDER
            if self.left_hit.is_reverse:
                self.reverse_hit = True
        else:
            # if we're here then the hits are on the same sequence in the same
            # orientation and correct order, so now have to check the distance between them
            if not self.left_hit.is_reverse:
                self.ref_start = self.left_hit.aend + 1
                self.ref_end = self.right_hit.pos - 1
            else:
                self.ref_start = self.right_hit.aend + 1
                self.ref_end = self.left_hit.pos - 1

            self.ref_name = samfile.getrname(self.left_hit.tid)
            if self.ref_start <= self.ref_end:
                self.gap_type = SAME_SEQ_STRAND
                self.query_replace_start = self.query_start
                self.query_replace_end = self.query_end
            else:
                self.gap_type = SAME_SEQ_STRAND_OVERLAP
                
                if not self.left_hit.is_reverse:
                    self.ref_start = self.left_hit.pos
                    self.ref_end = self.right_hit.aend
                    self.query_replace_start = self.query_start - self.left_hit.qlen
                    self.query_replace_end = self.query_end + self.right_hit.qend
                else:
                    self.ref_start = self.right_hit.pos
                    self.ref_end = self.left_hit.aend
                    self.query_replace_start = self.query_start - self.right_hit.qlen
                    self.query_replace_end = self.query_end + self.right_hit.qend
                
            if self.left_hit.is_reverse:
                self.reverse_hit = True

                
    def can_be_filled(self, abs_diff=500, relative_err=3):
        if self.gap_type == SAME_SEQ_STRAND:
            gap_length_in_query = self.query_end - self.query_start + 1
            gap_length_in_ref = self.ref_end - self.ref_start + 1
            percent_diff = 1.0 * abs(gap_length_in_query - gap_length_in_ref) / gap_length_in_query
            return percent_diff < relative_err and abs(gap_length_in_query - gap_length_in_ref) < abs_diff
        else:
            return self.gap_type == SAME_SEQ_STRAND_OVERLAP
        

    def _left_or_right_fasta_name_prefix(self):
        return self.query_name + ":"  \
                 + str(self.query_start + 1) \
                 + '-' \
                 + str(self.query_end + 1)

    def left_fasta(self):
        return sequences.Fasta(self._left_or_right_fasta_name_prefix() + '.left', self.left_seq)

    def right_fasta(self):
        return sequences.Fasta(self._left_or_right_fasta_name_prefix() + '.right', self.right_seq)

    #def dict_key(self):
    #    return (self.query_name, self.query_start, self.query_end)

