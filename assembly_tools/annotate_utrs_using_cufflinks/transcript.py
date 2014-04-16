import copy
import sys
from fastaq import intervals
from assembly_tools import file_readers

lenient = False
class Error (Exception): pass

class Transcript:
    def __init__(self, gff_record):
        self.coords = None
        self.five_utr = []
        self.three_utr = []
        self.exons = []
        self.mRNA = None
        self.ncRNA = []
        self.rRNA = []
        self.tRNA = []
        self.snRNA = []
        self.polypeptide = None
        self.strand = None
        self.seqname = None
        self.other_gffs = []
        self.add_gff_record(gff_record)
        
    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def add_gff_record(self, gff_record):
        gff_record = copy.deepcopy(gff_record)
        if gff_record.feature == 'five_prime_UTR':
            self.five_utr.append(gff_record)
        elif gff_record.feature == 'three_prime_UTR':
            self.three_utr.append(gff_record)
        elif gff_record.feature in ['CDS', 'exon', 'pseudogenic_exon']:
            self.exons.append(gff_record)
        elif gff_record.feature in ['mRNA', 'transcript', 'pseudogenic_transcript']:
            if self.mRNA is not None:
                raise Error('\n'.join([
                    "Transcript already has this mRNA/transcript/pseudogenic_transcript:",
                    str(self.mRNA),
                    "Cannot add this one as well:",
                    str(gff_record)
                ]))
            self.mRNA = gff_record
        elif gff_record.feature == 'ncRNA':
            self.ncRNA.append(gff_record)
        elif gff_record.feature == 'rRNA':
            self.rRNA.append(gff_record)
        elif gff_record.feature == 'tRNA':
            self.tRNA.append(gff_record)
        elif gff_record.feature == 'snRNA':
            self.snRNA.append(gff_record)
        elif gff_record.feature == 'polypeptide':
            if self.polypeptide is not None:
                raise Error('\n'.join([
                    "Transcript already has this polypeptide:",
                    str(self.polypeptide),
                    "Cannot add this one as well:",
                    str(gff_record)
                ]))
            self.polypeptide = gff_record
        else:
            self.other_gffs.append(gff_record)

        self._set_coords()
        self._set_strand()
        self._set_seqname()
        self._sort()
        
    def _sort(self):
        for l in [self.five_utr, self.three_utr, self.exons, self.ncRNA ,self.rRNA, self.tRNA, self.snRNA, self.other_gffs]:
            l.sort()

    def _set_coords(self):
        try:
            start = min([t.coords.start for t in self.five_utr + self.three_utr + self.exons + self.ncRNA + self.rRNA + self.tRNA + self.snRNA])
            end = max([t.coords.end for t in self.five_utr + self.three_utr + self.exons + self.ncRNA + self.rRNA + self.tRNA + self.snRNA])
        except:
            if self.mRNA is not None:
                start = self.mRNA.coords.start
                end = self.mRNA.coords.end
            else:
                return

        self.coords = intervals.Interval(start, end)
        if self.mRNA is not None:
            self.mRNA.coords = self.coords
    

    def __str__(self):
        a = [
            'start-end\t' + str(self.coords),
            "5'UTR\n\t" + '\n\t'.join([str(x) for x in self.five_utr]),
            "3'UTR\n\t" + '\n\t'.join([str(x) for x in self.three_utr]),
            'mRNA\n\t' + str(self.mRNA),
            'ncRNA\n\t' + '\n\t'.join([str(x) for x in self.ncRNA]),
            'rRNA\n\t' + '\n\t'.join([str(x) for x in self.rRNA]),
            'tRNA\n\t' + '\n\t'.join([str(x) for x in self.tRNA]),
            'snRNA\n\t' + '\n\t'.join([str(x) for x in self.snRNA]),
            'exons\n\t' + '\n\t'.join([str(x) for x in self.exons]),
            'polypeptide\n\t' + str(self.polypeptide),
            'other\n\t' + '\n\t'.join([str(x) for x in self.other_gffs])
        ]

        return '\n'.join(a)
                   
    def total_exon_length(self):
        return sum([len(x) for x in self.exons])
      
            

    def _set_seqname(self):
        names = set([g.seqname for g in self.five_utr + self.three_utr + self.exons + self.ncRNA + self.rRNA + self.tRNA + self.snRNA + [self.mRNA] + self.other_gffs if g is not None])
        if len(names) != 1:
            raise Error('Error getting seqname for transcript. Cannot continue')
             
        name = names.pop()
        if self.seqname is None:
            self.seqname = name
        elif self.seqname != name:
            raise Error('Error getting seqname for transcript. Cannot continue')

    def _strand_error_message(self):
        return '*** Error getting strand info for transcript...\n' + str(self) + '\n***\n'

    def _set_strand(self):
        strands = set([g.strand for g in self.five_utr + self.three_utr + self.exons + self.ncRNA + self.rRNA + self.tRNA + self.snRNA + [self.mRNA] + self.other_gffs if g is not None])
        if len(strands) != 1:
            if lenient:
                self.strand = 'Inconsistent'
                print(self._strand_error_message(), file=sys.stderr)
            else:
                raise Error(self._strand_error_message())
             
        strand = strands.pop()
        if self.strand is None:
            self.strand = strand
        elif self.strand != strand:
            if lenient:
                self.strand = 'Inconsistent'
                print(self._strand_error_message(), file=sys.stderr)
            else:
                raise Error(self._strand_error_message())


    def intersects(self, other):
       return self.coords.intersects(other.coords)

    def might_extend(self, other, min_extend=1):
        coords = intervals.Interval(self.coords.start - min_extend + 1, self.coords.end + min_extend - 1)
        strands_ok = (self.strand == other.strand and self.strand not in ['.', 'Inconsistent']) \
                     or (self.strand in ['-', '+'] and other.strand == '.' and len(self.exons) == len(other.exons) == 1)
        return self.seqname == other.seqname \
             and strands_ok \
             and len(self.exons) * len(other.exons) != 0 \
             and (coords.intersects(other.coords) or other.coords.end + 1 == coords.start or coords.end + 1 == other.coords.start) \
             and (other.coords.start < coords.start or coords.end < other.coords.end)

    def can_extend_start(self, other, min_extend=1):
        return self.might_extend(other, min_extend=min_extend) \
                and other.coords.start + min_extend - 1 < self.coords.start

    def can_extend_end(self, other, min_extend=1):
        return self.might_extend(other, min_extend=min_extend) \
                and self.coords.end + min_extend - 1 < other.coords.end


    def exon_splice_sites(self):
       sites = []
       if len(self.exons) > 1:
           for e in self.exons:
               sites.extend([e.coords.start, e.coords.end])
           sites.pop(0)
           sites.pop()

       return sites
           
    def number_of_common_splice_sites(self, other):
        return len(set(self.exon_splice_sites()).intersection(set(other.exon_splice_sites())))

    def _add_utr_info(self, other, min_extend=1, max_new_utrs=3, extend_end=False, exclude_coords=[]):
        if (not extend_end and not self.can_extend_start(other, min_extend=min_extend)) \
          or (extend_end and not self.can_extend_end(other, min_extend=min_extend)):
            return

        if (self.strand == '+' and extend_end) or (self.strand == '-' and not extend_end):
            utr_list = self.three_utr
            gff_feature = 'three_prime_UTR'
            gff_id_suffix = ':3utr'
        else:
            utr_list = self.five_utr
            gff_feature = 'five_prime_UTR'
            gff_id_suffix = ':5utr'

        if len(utr_list):
            return

        if extend_end:
            i = len(other.exons) - 1
            def continue_while():
                return i > 0 and other.exons[i].coords.start > self.exons[-1].coords.end
            def update_i(i):
                return i - 1
        else:
            i = 0
            def continue_while():
                return i < len(other.exons) and other.exons[i].coords.end < self.exons[0].coords.start
            def update_i(i):
                return i + 1

        new_utrs = []

        while continue_while():
            new_gff = file_readers.gff.GFF_record('\t'.join([
                self.seqname,
                'UTR_updater',
                gff_feature,
                str(other.exons[i].coords.start),
                str(other.exons[i].coords.end),
                '.',
                self.strand,
                '.'
            ]))
            new_gff.set_attribute('ID', self.exons[0].attributes['Parent'] + gff_id_suffix)
            new_gff.set_attribute('Parent', self.exons[0].attributes['Parent'])
            intersects = False
            for c in exclude_coords:
                if new_gff.coords.intersects(c) or self.coords.intersects(c):
                    intersects = True
                    break
            if not intersects:
                new_utrs.append(new_gff)
            i = update_i(i)

        if len(new_utrs) > max_new_utrs:
            return

        if i < len(other.exons):
            if ((not extend_end) and self.exons[0].intersects(other.exons[i])) or \
               (extend_end and self.exons[-1].intersects(other.exons[i])):
                new_gff =  None

                try:
                    new_gff = file_readers.gff.GFF_record('\t'.join([
                        self.seqname,
                        'UTR_updater',
                        gff_feature,
                        str(self.exons[-1].coords.end + 1) if extend_end else str(other.exons[i].coords.start),
                        str(other.exons[i].coords.end) if extend_end else str(self.exons[0].coords.start - 1),
                        '.',
                        self.strand,
                        '.'
                    ]))
                except:
                    pass

                if new_gff is not None:
                    intersects = False
                    for c in exclude_coords:
                        if new_gff.coords.intersects(c) or self.coords.intersects(c):
                            intersects = True
                            break
                    if not intersects:
                        new_gff.set_attribute('ID', self.exons[0].attributes['Parent'] + gff_id_suffix)
                        new_gff.set_attribute('Parent', self.exons[0].attributes['Parent'])
                        new_utrs.append(new_gff)

        if len(new_utrs) > max_new_utrs:
            return

        counter = 1
        new_utrs.sort()

        for utr in new_utrs:
            if len(new_utrs) > 1:
                utr.set_attribute('ID', utr.get_attribute('ID') + ':' + str(counter))
                counter += 1
            self.add_gff_record(utr)


    def update_utrs(self, other, min_extend=1, max_new_utrs=3, exclude_coords=[]):
        self._add_utr_info(other, min_extend=min_extend, max_new_utrs=max_new_utrs, extend_end=False, exclude_coords=exclude_coords)
        self._add_utr_info(other, min_extend=min_extend, max_new_utrs=max_new_utrs, extend_end=True, exclude_coords=exclude_coords)
        self._set_coords()


    def to_gff_list(self):
        l = [x for x in self.five_utr + self.three_utr + self.exons + [self.mRNA] + self.ncRNA + self.rRNA + self.tRNA + self.snRNA + [self.polypeptide] + self.other_gffs if x is not None]
        l.sort()
        return l
