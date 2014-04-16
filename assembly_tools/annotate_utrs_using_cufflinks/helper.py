class Error (Exception): pass

import sys
from fastaq import *
from assembly_tools import file_readers
from assembly_tools.annotate_utrs_using_cufflinks import gene


def read_gff(filename):
    level_records = [[], [], []]
    other_records = {}
    gff_reader = file_readers.gff.file_reader(filename)
    for g in gff_reader:
        for i in range(len(gene.feature_levels)):
            if g.feature in gene.feature_levels[i]:
                level_records[i].append(g)
                break
        else:
            update_other_records(other_records, g)

    return level_records, other_records


def update_other_records(other_records, gff_record):
    if gff_record.seqname not in other_records:
        other_records[gff_record.seqname] = []
    other_records[gff_record.seqname].append(gff_record)


def get_genes_from_ref(level_records, other_records):
    level2_to_level1 = {}
    genes = {}

    while len(level_records[0]):
        g = level_records[0].pop()
        g_id = g.get_attribute('ID')
        if g.seqname not in genes:
            genes[g.seqname] = {}
        genes[g.seqname][g_id] = gene.Gene(g)

    # add the transcript info to each gene
    while len(level_records[1]):
        g = level_records[1].pop()
        g_id = g.get_attribute('ID')

        if g_id is None:
            raise Error('Error getting ID/gene_id from GFF line\n' + str(g))

        parent_id = g.get_attribute('Parent')

        if g.seqname not in genes:
            raise Error('No parent in sequence "' + g.seqname + '" for this:\n' + str(g))

        if parent_id not in genes[g.seqname]:
            update_other_records(other_records, g)
            print('Warning: Parent with id "' + parent_id + '" not found for this feature:', str(g), sep='\n\t', file=sys.stderr)
            continue

        genes[g.seqname][parent_id].add_gff_record(g)
        level2_to_level1[g_id] = parent_id

    # add the exon/CDSs to each gene
    while len(level_records[2]):
        g = level_records[2].pop()
        if g.feature ==  'polypeptide':
            parent_id = g.get_attribute('Derives_from')
        else:
            parent_id = g.get_attribute('Parent')

        try:
            gene_id = level2_to_level1[parent_id]
        except:
            update_other_records(other_records, g)
            print('Warning: Parent of "' + parent_id + '" not found, originating from this line:', str(g), sep='\n\t', file=sys.stderr)
            continue

        genes[g.seqname][gene_id].add_gff_record(g)


    return genes, other_records


def load_ref_gff(filename):
    level_records, other_records = read_gff(filename)
    genes, other_records = get_genes_from_ref(level_records, other_records)
    sort_gene_dict_values(genes)
    return genes, other_records


def initialize_genes_dict_cufflinks(records):
    genes = {}
    for gff_record in records:
        gene_id = gff_record.get_attribute('gene_id')
        if gene_id is None:
            raise Error('Error! No gene_id from this cufflinks gff line:\n' + str(gff_record))

        if gff_record.seqname not in genes:
            genes[gff_record.seqname] = {}

    return genes


def get_genes_from_cufflinks(level_records):
    genes = initialize_genes_dict_cufflinks(level_records[1])

    # add the transcript info to each gene
    while len(level_records[1]):
        g = level_records[1].pop()
        transcript_id = g.get_attribute('transcript_id')
        gene_id = g.get_attribute('gene_id')

        if g.seqname not in genes:
            raise Error('No parent in sequence "' + g.seqname + '" for this:\n' + str(g))

        if gene_id not in genes[g.seqname]:
            genes[g.seqname][gene_id] = gene.Gene(g)
        else:
            genes[g.seqname][gene_id].add_gff_record(g)

    # add the exon/CDSs to each gene
    while len(level_records[2]):
        g = level_records[2].pop()
        gene_id = g.get_attribute('gene_id')
        if g.seqname not in genes or gene_id not in genes[g.seqname]:
            raise Error('No parent for this:\n' + str(g))
        genes[g.seqname][gene_id].add_gff_record(g)

    return genes


def load_cufflinks_gtf(filename):
    level_records, x = read_gff(filename)
    genes = get_genes_from_cufflinks(level_records)
    sort_gene_dict_values(genes)
    return genes

def gene_dict_to_sorted_list(d):
    l = []
    while(len(d)):
        g_id, g = d.popitem()
        l.append(g)
    l.sort()
    return l

def sort_gene_dict_values(d):
    for k in d:
        d[k] = gene_dict_to_sorted_list(d[k])

