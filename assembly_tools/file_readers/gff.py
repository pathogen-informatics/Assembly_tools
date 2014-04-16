from fastaq import utils, intervals
import re
import sys

class Error (Exception): pass

is_gtf_re = re.compile('^(.*) "(.*)"$')
lenient = False
warnings = True

def file_reader(fname):
    f = utils.open_file_read(fname)
    for line in f:
        if line.startswith('##FASTA') or line.startswith('>'):
            break
        elif line.startswith('#'):
            continue
        else:
            yield GFF_record(line)

    utils.close(f)


class GFF_record:
    def __init__(self, line=''):
        data = line.rstrip().split('\t')
        if not (8 <= len(data) <= 9):
            raise Error('Error reading GFF file. The following line does not have 8 or 9 columns:\n' + line)

        self.seqname = data[0]
        self.source= data[1]
        self.feature = data[2]

        try:
            start = int(data[3])
            end = int(data[4])
        except:
            raise Error('Error reading GFF file. The following line\'s start or end coord is not an integer:\n' + line)
        self.coords = intervals.Interval(start, end)

        if (data[5] == '.'):
            self.score = None
        else:
            try:
                self.score = float(data[5])
            except:
                raise Error('Error reading GFF file. The following line\'s score does not appear to be a number:\n' + line)

        self.strand = data[6]
        if self.strand not in ['-', '+', '.']:
            raise Error('Error reading GFF file. The following line\'s frame is not +,- or .:\n' + line)

        self.frame = data[7]

        self.attributes = {}
        self.attribute_keys = []
        self.is_gtf = False

        if len(data) == 9:
            for att in data[8].rstrip(';').split(';'):
                hits = is_gtf_re.search(att.strip())
                if hits is not None:
                    key, val = hits.group(1), hits.group(2)
                    self.is_gtf = True
                else:
                    try:
                        (key, val) = att.split('=', 1)
                    except ValueError:
                        error_message = 'Error splitting into key/value pair:\n' +  att +  '\nfrom GFF file line:\n' + line
                        if lenient:
                            key = att
                            val = None
                            if warnings:
                                print(error_message, file=sys.stderr)
                        else:
                            raise Error(error_message)

                self.attributes[key] = val
                self.attribute_keys.append(key)


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        if self.score is None:
            score = '.'
        else:
            score = str(self.score)

        s = '\t'.join([self.seqname, self.source, self.feature, str(self.coords.start), str(self.coords.end), score, self.strand, self.frame])

        if (len(self.attributes) == 0):
            return s
        else:
            atts = []
            for k in self.attribute_keys:
                if self.attributes[k] is None:
                    atts.append(k)
                elif self.is_gtf:
                    atts.append(k + ' "' + self.attributes[k] + '"')
                else:
                    atts.append(k + '=' + self.attributes[k])

            if self.is_gtf:
                return s + '\t' + '; '.join(atts) + ';'
            else:
                return s + '\t' + ';'.join(atts)

    def __len__(self):
        return len(self.coords)

    def __lt__(self, other):
        return self.seqname == other.seqname and self.coords < other.coords

    def intersects(self, other):
        return self.seqname == other.seqname and self.coords.intersects(other.coords)

    def get_attribute(self, key):
        try:
            return self.attributes[key]
        except:
            raise Error('Attribute "' + key + '" not found from:\n' + str(self))

    def set_attribute(self, key, value):
        if key not in self.attribute_keys:
            self.attribute_keys.append(key)

        self.attributes[key] = value

