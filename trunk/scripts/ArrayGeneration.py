#!/usr/bin/env python

"""
Implements basic functionality for generating CRISPR spacers.
"""

__author__ = "Ruben Acuna"
__copyright_ = "Copyright(c) 2011, ASU iGEM Team"

################################### IMPORTS ####################################
from Bio import Restriction
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from ftplib import FTP
import os.path
import sys
import pickle
import Config
import re

################################## CONSTANTS ###################################
DISALLOWED_RE = pickle.load(open(Config.SETTINGS))[0]
BLAST_FILE_DB = "blast_fasta.txt"
BLAST_FILE_INPUT = "blast_input.txt"
BLAST_FILE_OUTPUT = "blast_output.xml"
BLAST_EVALUE = 10 #web interface default
HAIRPIN_MAX = 8

################################## FUNCTIONS ###################################
def getPossibleSpacers(sequence, spacerLength):
    """Given a sequence and a spacerLength, returns all continuous subsets of spacerLength length."""
    # PAM: [required, [name, seq, pos]]
    
    settings_pickle = pickle.load(open(Config.SETTINGS))
    PAM = settings_pickle[1][1][settings_pickle[1][0]]
    PAM_seq_re = PAM[1][1].lower().replace('n', '.')
    PAM_seq_re_comp = re.compile(PAM_seq_re)
    possibleSpacers = []

    #Y Guido no implement real recursion????
    for i in xrange(0, len(sequence) - spacerLength + 1):
        
        if PAM[0] == 1:
            if PAM[1][2] == 1:
                if i + spacerLength + len(PAM_seq_re) <= len(sequence):
                    test_seq = str(sequence[i + spacerLength:i + spacerLength + len(PAM_seq_re)]).lower()
                    re_match = re.match(PAM_seq_re_comp, test_seq)
                    if re_match != None:
                        if re_match.group() == test_seq:
                            possibleSpacers.append(sequence[i:i+spacerLength])
            else:
                if i >= len(PAM_seq_re):
                    test_seq = str(sequence[i - len(PAM_seq_re):i]).lower()
                    re_match = re.match(PAM_seq_re_comp, test_seq)
                    if re_match != None:
                        if re_match.group() == test_seq:
                            possibleSpacers.append(sequence[i:i+spacerLength])
        else:
            possibleSpacers.append(sequence[i:i+spacerLength])

    return possibleSpacers

def filterRE(sequences):
    """Given a list of sequences, returns a subset containing only those sequences lacking disallowed RE sites."""

    #make sure this up to date
    settings_pickle = pickle.load(open(Config.SETTINGS))
    DISALLOWED_RE = settings_pickle[0][1]
    if settings_pickle[0][0] == False:
        return [sequence for sequence in sequences if not [re for re in DISALLOWED_RE if re.search(sequence)]]
    else:
        return sequences

def filterHomology(sequences, other):
    """Given a list of sequences and another sequence, returns a subset containing only those sequences which not do
       have homology with the single sequence."""

    SeqIO.write([SeqRecord(sequence, id="seq%s" % i) for i, sequence in enumerate(sequences)], open(BLAST_FILE_INPUT, "wb"), "fasta") 
    SeqIO.write(SeqRecord(other, id="other"), open(BLAST_FILE_DB, "wb"), "fasta")
    
    dbGenerate(BLAST_FILE_DB)
    
    commandBLASTN = NcbiblastnCommandline(query=BLAST_FILE_INPUT, db=BLAST_FILE_DB, evalue=BLAST_EVALUE, outfmt=5, out=BLAST_FILE_OUTPUT)
    commandBLASTN()

    blastRecords = list(NCBIXML.parse(open(BLAST_FILE_OUTPUT)))  
    nonhomologousSequences = []
    
    for i, record in enumerate(blastRecords):
        if len(record.alignments) == 0:
            nonhomologousSequences += [sequences[i]]
    
    return nonhomologousSequences

def filterHomologyDatabase(sequences, database):
    """Given a list of sequences and NCBI BLAST+ database name, returns a subset containing only those sequences which
    not do have homology with the database."""

    SeqIO.write([SeqRecord(sequence, id="seq%s" % i) for i, sequence in enumerate(sequences)], open(BLAST_FILE_INPUT, "wb"), "fasta")
    
    commandBLASTN = NcbiblastnCommandline(query=BLAST_FILE_INPUT, db=database, evalue=BLAST_EVALUE, outfmt=5, out=BLAST_FILE_OUTPUT)
    commandBLASTN()    

    blastRecords = list(NCBIXML.parse(open(BLAST_FILE_OUTPUT)))
    nonhomologousSequences = []
    
    for i, record in enumerate(blastRecords):
        if len(record.alignments) == 0:
            nonhomologousSequences += [sequences[i]]           
            
    return nonhomologousSequences
    
def filterHairpinning(sequences):
    """Given a list of sequences, returns a subset containing only those sequences which have less than the required
       number of paired bases in secondary structure."""

    return [sequence for sequence in sequences if pairedBases(sequence) <= HAIRPIN_MAX]

def pairedBases(seq):
    """Implements the algorithm from Tardos and Kleinberg (2005) for detecting RNA secondary structure. Returns maximum
       paired bases in structure."""

    M = [[0]*len(seq) for i in xrange(len(seq))]

    for k in xrange(4, len(seq) - 1):
        for i in xrange(0, len(seq) - k - 1):

            j = i + k + 1

            #cannot hairpin if less than four residues away
            if i >= j - 4:
                #print "skip"
                M[i][j] = 0
            #can hairpin
            elif comp(seq[i]) == seq[j]:
                #print "pairing ", seq[i], " with ", seq[j]
                M[i][j] = 1 + max([M[i][t - 1] + M[t + 1][j - 1] for t in xrange(i, j - 4)])
            #cannot hairpin
            else:
                #print "not pairing ", seq[i], " with ", seq[j]
                M[i][j] = M[i][j - 1]

    return M[0][len(seq) - 1]

#borrowed from Ethan
def comp(n):
    return ['A', 'T', 'C', 'G'][['T', 'A', 'G', 'C'].index(n)]

def dbGenerate(filename):
    """Given the file name of a FASTA file, creates a BLAST DB with the same name. Returns database name."""

    if not (os.path.exists(filename + ".nin") and os.path.exists("data/ncbi/" + filename + ".nsq") and os.path.exists("data/ncbi/" + filename + ".nhr")):
        os.system("makeblastdb -in " + filename + " -dbtype nucl")
        
    return filename

def dbFetchNCBI(uid):
    """ Given a NCBI genome uid, downloads the associated FASTA and stores it locally using the default name. Returns file name. """

    genomeServer = FTP("ftp.ncbi.nih.gov")
    genomeServer.login()    
    genomeServer.cwd("genomes/Bacteria")
    genomeServer.cwd([file for file in genomeServer.nlst() if str(uid) in file][0])

    filename = [file for file in genomeServer.nlst() if ".fna" in file][0]
    
    if not os.path.exists(filename):
        fileDBFASTA = open(filename, "wb")
        genomeServer.retrbinary("RETR " + filename, fileDBFASTA.write)
        fileDBFASTA.close()
    
    genomeServer.close()
    
    return filename

def getValidSpacers(sequence, repeat, organismUID, spacerLength, filterHairPinnning):
    """Given a target sequence, a repeat sequence, a ncbi genonme refseq, and a spacer length, returns a set of spacers
       meeting the required criteria:
         1) No BioBrick restriction sites.
         2) No homology with repeat sequence
         3) Minimal hairpinning
         4) No homology with host genome"""

    sys.stdout.write("    Finding possible spacers...")
    spacers = getPossibleSpacers(sequence, spacerLength)
    sys.stdout.write(" found %s.\n" % len(spacers))

    sys.stdout.write("    Removing spacers with REs...")
    spacers = filterRE(spacers)
    sys.stdout.write(" %s remaining spacers.\n" % len(spacers))

    sys.stdout.write("    Removing repeat homologous spacers...")
    spacers = filterHomology(spacers, repeat)
    sys.stdout.write(" %s remaining spacers.\n" % len(spacers))

    if filterHairPinnning:
        sys.stdout.write("    Removing hairpinning spacers...")
        spacers = filterHairpinning(spacers)
        sys.stdout.write(" %s remaining spacers.\n" % len(spacers))

    sys.stdout.write("    Removing genome homologous spacers...")
    spacers = filterHomologyDatabase(spacers, dbGenerate(dbFetchNCBI(organismUID)))
    sys.stdout.write(" %s remaining spacers.\n" % len(spacers))

    return spacers
