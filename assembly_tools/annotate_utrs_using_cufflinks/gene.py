import sys
import copy
from assembly_tools.annotate_utrs_using_cufflinks import transcript
from fastaq import intervals

lenient = False
class Error (Exception): pass

feature_levels = [
    set(['gene', 'pseudogene']),
    set(['mRNA', 'ncRNA', 'rRNA', 'snRNA', 'tRNA', 'transcript', 'pseudogenic_transcript']),
    set(['five_prime_UTR', 'three_prime_UTR', 'CDS', 'exon', 'pseudogenic_exon'])
]


class Gene:
    def __init__(self, gff_record):
        self.gene_id = None
        self.transcripts = {}
        self.gene_record = None
        self.strand = None
        self.coords = None
        self.seqname = None
        self.add_gff_record(gff_record)

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def _set_seqname(self):
        names = set([t.seqname for t in self.transcripts.values()])

        if len(names) != 1:
            if self.gene_record is not None:
                names = set([self.gene_record.seqname])
            else:
                raise Error('Error getting seqname for gene - too many names. Names: ' + str(names))
             
        name = names.pop()
        if self.seqname is None:
            self.seqname = name
        elif self.seqname != name:
            raise Error('Error getting seqname for gene - too many names. Names:' + self.seqname + ', ' + str(names))

    def _set_strand(self, gff_record):
        if gff_record.strand != '.':
            if self.strand != None and self.strand != gff_record.strand:
                if lenient:
                    print('Warning: Strand inconsistency from this line of gff file:\n' + str(gff_record), file=sys.stderr)
                    self.strand = 'Inconsistent'
                else:
                    raise Error('Strand inconsistency from this line of gff file:\n' + str(gff_record))
            if self.strand == None:
                self.strand = gff_record.strand

    def _set_coords(self):
        try:
            start = min([t.coords.start for t in self.transcripts.values()])
            end = max([t.coords.end for t in self.transcripts.values()])
        except:
            return

        self.coords = intervals.Interval(start, end)

        if self.gene_record is not None:
            self.gene_record.coords = self.coords

    def add_gff_record(self, gff_record):
        gff_record = copy.deepcopy(gff_record)
        self._set_strand(gff_record)

        if gff_record.feature in feature_levels[0]:
            assert self.gene_record == None
            self.gene_id = gff_record.get_attribute('ID')
            self.gene_record = gff_record
        else:
            if gff_record.is_gtf:
                transcript_id = gff_record.get_attribute('transcript_id')
                gene_id = gff_record.get_attribute('gene_id')

            if gff_record.feature in feature_levels[1]:
                if not gff_record.is_gtf:
                    transcript_id = gff_record.get_attribute('ID')
                    gene_id = gff_record.get_attribute('Parent')
       
                self.gene_id = self.gene_id if self.gene_id is not None else gene_id
                if gene_id != self.gene_id:
                    raise Error('gene ID of the following line is not ' + str(self.gene_id) + '\n' + str(gff_record))
            elif gff_record.feature in feature_levels[2]:
                if not gff_record.is_gtf:
                    transcript_id = gff_record.get_attribute('Parent')
            else:
                raise Error('Error adding this line to gene information:\n' + str(gff_record))

            if transcript_id not in self.transcripts:
                self.transcripts[transcript_id] = transcript.Transcript(gff_record)
            else:
                self.transcripts[transcript_id].add_gff_record(gff_record)

        self._set_coords()
        self._set_seqname()

    def __lt__(self, other):
        return self.seqname == other.seqname and self.coords < other.coords

    def intersects(self, other):
       return self.seqname == other.seqname and self.coords.intersects(other.coords)

    def can_extend(self, other, min_extend=1):
        if self.seqname != other.seqname:
            return False

        for t in self.transcripts.values():
            for u in other.transcripts.values():
                if t.can_extend_start(u, min_extend=min_extend) or t.can_extend_end(u, min_extend=min_extend):
                    return True

        return False

    def extend(self, other, min_extend=1):
        if (not self.can_extend(other)) or other.coords.start + min_extend >= self.coords.start:
            return

        for t in self.transcripts.values():
           max_splices_in_common = -1
           max_splices_key = None
           for key, trans in other.transcripts.items():
               splices = t.number_of_common_splice_sites(trans)
               if splices > max_splices_in_common:
                   max_splices_in_common = splices
                   max_splices_key = key
           if max_splices_key is not None:
               t.update_utrs(other.transcripts[max_splices_key])
            
        self._set_coords()

    def __str__(self):
        return str(self.gene_id) + '\t\n' + str(self.gene_record) + '\n' + '\n'.join([str(x) for x in self.transcripts.values()]) + '\n----\n'


    def to_gff_list(self):
        l = []
        if self.gene_record is not None:
            l.append(self.gene_record)

        for t in self.transcripts.values():
            l += t.to_gff_list()
        return l
