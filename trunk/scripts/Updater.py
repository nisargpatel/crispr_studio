#!/usr/bin/env python

"""

"""

__author__ = "Ethan Ward, Ruben Acuna"
__copyright_ = "Copyright(c) 2011, ASU iGEM Team"

from utilities import UtilWeb
import datetime
import pickle

DATA_PATH = "data/"

CRISPRDB_PATH = DATA_PATH + "crisprDB/"
CRISPRDB_SPACER_RAW_PATH = CRISPRDB_PATH + "spacer_info.txt"
CRISPRDB_REPEAT_RAW_PATH = CRISPRDB_PATH + "repeat_info.txt"
CRISPRDB_SPACER_URL = "http://crispr.u-psud.fr/crispr/BLAST/Spacer/Spacerdatabase"
CRISPRDB_REPEAT_URL = "http://crispr.u-psud.fr/crispr/BLAST/DR/DRdatabase"
REPEAT_DB_PATH = DATA_PATH + "repeatDB.txt"
SPACER_DB_PATH = DATA_PATH + "spacerDB.txt"
ARRAY_DB_PATH = DATA_PATH + "arrayDB.txt"

NCBI_PATH = DATA_PATH + "ncbi/"
LPROKS_0_PATH = NCBI_PATH + "lproks_0.txt"
LPROKS_1_PATH = NCBI_PATH + "lproks_1.txt"
LPROKS_2_PATH = NCBI_PATH + "lproks_2.txt"
TAX_INFO_PATH = DATA_PATH + "tax_info.txt"
TAX_INFO_FILTERED_PATH = DATA_PATH + "tax_info_filtered.txt"
GENOMEPRJ_URL = ["ftp.ncbi.nih.gov", "genomes/genomeprj/"]

UPDATE_PATH = DATA_PATH + "last_update.txt"

def update_tax_info():
    def update_source():
        
        UtilWeb.fetchFileFTP(GENOMEPRJ_URL[0], GENOMEPRJ_URL[1], 'lproks_0.txt', LPROKS_0_PATH)            
        UtilWeb.fetchFileFTP(GENOMEPRJ_URL[0], GENOMEPRJ_URL[1], 'lproks_1.txt', LPROKS_1_PATH)            
        UtilWeb.fetchFileFTP(GENOMEPRJ_URL[0], GENOMEPRJ_URL[1], 'lproks_2.txt', LPROKS_2_PATH)            

    def process_source():
        lproks_0 = [a.split('\t') for a in open(LPROKS_0_PATH).readlines()][2:]
        lproks_1 = [a.split('\t') for a in open(LPROKS_1_PATH).readlines()][2:]

        ref_seq = []
        tax_id = []
        species_name = []
        accession = []

        for line in lproks_0:
            ref_seq.append(line[0])
            tax_id.append(line[2])
            species_name.append(line[3])
            accession.append('')
        for line in lproks_1:
            if line[2] in tax_id:
                accession[tax_id.index(line[2])] = line[13]

        result = []
        for a in range(len(tax_id)):
            result.append([ref_seq[a], tax_id[a], species_name[a], accession[a]])
            
        wfile = open(TAX_INFO_PATH, "w")
        pickle.dump(result, wfile)
        wfile.close()
    update_source()
    process_source()
        
def update_repeats():
    def update_source():
        UtilWeb.fetchFileHTTP(CRISPRDB_REPEAT_URL, CRISPRDB_REPEAT_RAW_PATH)
        
    def process_source():
        DR_database = [a.replace('\n', '') for a in open(CRISPRDB_REPEAT_RAW_PATH).readlines()]

        result = []
        for a in range(len(DR_database) - 1):
            if a % 2 == 0:
                DR_cons = DR_database[a + 1]
                for a_number in DR_database[a][1:].split('|'):
                    result.append([DR_cons, a_number])

        wfile = open(REPEAT_DB_PATH, "w")
        pickle.dump(result, wfile)
        wfile.close()
    update_source()
    process_source()

def update_spacers():
    def update_source():
        UtilWeb.fetchFileHTTP(CRISPRDB_SPACER_URL, CRISPRDB_SPACER_RAW_PATH)

    def process_source():
        spacer_db = [a.replace('\n', '') for a in open(CRISPRDB_SPACER_RAW_PATH).readlines()]

        result = []
        for a in range(len(spacer_db) - 1):
            if a % 2 == 0:
                spacer = spacer_db[a + 1]
                for a_number in spacer_db[a][1:].split('|'):
                    result.append([spacer, a_number])
        wfile = open(SPACER_DB_PATH, "w")
        pickle.dump(result, wfile)
        wfile.close()
    update_source()
    process_source()

def create_array_db():
    repeat_db = pickle.load(open(REPEAT_DB_PATH))
    spacer_db = pickle.load(open(SPACER_DB_PATH))
    tax_info = pickle.load(open(TAX_INFO_PATH))

    locus_tax_id_dict = {}

    for a in tax_info:
        if a[3] != '':
            for b in a[3].split(','):
                locus_tax_id_dict[b] = a[1]

    # {base locus: {i_locus: [repeat, [s1, s2...]], i_locus2: ...
    final_dict = {}
    
    for a in repeat_db:
        locus_base = a[1][::-1][a[1][::-1].index('_') + 1:][::-1]
        if locus_tax_id_dict.has_key(locus_base):
            if final_dict.has_key(locus_tax_id_dict[locus_base]):
                final_dict[locus_tax_id_dict[locus_base]][a[1]] = [a[0], []]
            else:
                final_dict[locus_tax_id_dict[locus_base]] = {a[1]: [a[0], []]}

    for a in spacer_db:
        repeat_locus = a[1][::-1][a[1][::-1].index('_') + 1:][::-1]
        locus_base = repeat_locus[::-1][repeat_locus[::-1].index('_') + 1:][::-1]
        if locus_tax_id_dict.has_key(locus_base):
            if final_dict.has_key(locus_tax_id_dict[locus_base]):
                final_dict[locus_tax_id_dict[locus_base]][repeat_locus][1].append(a[0])

    wfile = open(ARRAY_DB_PATH, "w")
    pickle.dump(final_dict, wfile)
    wfile.close()

def filter_tax_info():
    # filter tax_info to use only organisms in arraydb
    tax_info = pickle.load(open(TAX_INFO_PATH))
    arrayDB = pickle.load(open(ARRAY_DB_PATH))

    tax_info_filtered = [a for a in tax_info if a[1] in arrayDB.keys()]
    wfile = open(TAX_INFO_FILTERED_PATH, "w")
    pickle.dump(tax_info_filtered, wfile)
    wfile.close()
    
def update_date():
    f = open(UPDATE_PATH, "w")
    f.write(str(datetime.datetime.now()))
    f.close()

def main():
    print "Building taxon data."
    update_tax_info()
    print "Building repeat data."
    update_repeats()
    print "Building spacer data."
    update_spacers()
    print "Building sequence dictionary."
    create_array_db()
    filter_tax_info()
    print "Successfully updated all files"
    update_date()
