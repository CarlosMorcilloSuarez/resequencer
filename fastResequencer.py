#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
fastResequencer.py
    Generates a fastq file simulating a resequencing of an individual with
    a genome identical to the reference genome

    It's a faster and less memory consuming version of resequencer.py but
    doesn't allow for SVs (structural variants) introduction.
'''

__author__ = "Carlos Morcillo Suarez"
__license__ = "GPL"
__version__ = "1.4"

import sys
import re
import random

import getopt
import numpy as np
import gzip
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser


def usage():
    print '''
        fastResequencer.py

        Generates a simulated fastq file with reads from a NGS experiment
        from a reference genome in fasta format

        Use

            fastResequencer.py [-r,--reference] genome.fa [options]

        Options

            -o, --output-file
                name of the fastaq generated file (default: "sample.fq")
                with *.gz the output will be gzipped

            -c, --coverage
                read coverage of the genome (default: 15)

            -e, --error-rate
                error rate of reads respect to the reference genome
                (default: 0)

            -l, --length
                length of reads to be generated
                (default: 100)

            -s, --seed
                a seed for the random generator
                
            -f, --full
                fragments the genome in reads with the given overlap instead 
                of stochastically generating reads at random places
                --full 10
                    Next read begins at position 10 of previous one.
                --full  0
                    Contiguous fragments without overlap 
                    
                If overlap is bigger than length, gaps between reads will
                result.

    '''


# Process command line options
def proccessCommandLine(argv):
    try:
        opts, remainder = getopt.getopt(
                                argv,
                                "r:o:c:e:s:l:f:",
                                ["reference=", "output-file=", "coverage=",
                                 "error-rate=", "seed=",
                                 "length=","full="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Checks that mandatory options are present
    if '-r' not in [options for options,values in opts]:
        if '--reference' not in [options for options,values in opts]:
            usage()
            sys.exit(2)

    # Process read options
    for opt, arg in opts:
        if opt in ("-r", "--reference"):
            global referenceGenomeFileName
            referenceGenomeFileName = arg
        elif opt in  ("-o", "--output-file"):
            global outputFileName
            outputFileName = arg
        elif opt in ("-c", "--coverage"):
            global coverage
            coverage = float(arg)
        elif opt in ("-e", "--error-rate"):
            global errorRate
            errorRate = float(arg)
        elif opt in ("-l", "--length"):
            global readLength
            readLength = int(arg)
        elif opt in ("-s", "--seed"):
            global seed
            seed = arg
        elif opt in ("-f", "--full"):
            global overlap
            global full
            full = True
            overlap = int(arg)        


def kmerGenerator(referenceGenomeFileName,length,overlap):
    '''
        Generator that returns all kmers from the reference sequence
        of the given length

        It navigates the chromosome structure of the fasta file
        It returns kmers shorter than the given length when it arrives at
        the end of a chromosome
        
        Returns tuple:
            chromosome,
            init position of fragment (BED coordinates, i.e begins wihth 0),
            end position of fragment  (BED coordinates, i.e end position
                                       not included),
            sequence of fragment
    '''
    with open(referenceGenomeFileName,"r") as referenceGenomeFile:
        for chr,seq in SimpleFastaParser(referenceGenomeFile):
            currentPosition = 0
            while currentPosition < len(seq):
                yield ( chr,
                        currentPosition,
                        currentPosition+length,
                        seq[currentPosition:currentPosition+length]
                        )
                currentPosition += overlap


if __name__ == '__main__':

    referenceGenomeFileName = ''
    outputFileName = 'sample.fq'
    coverage = 15
    readLength = 100
    errorRate = 0
    seed = None
    full = False
    overlap = 0

    # Process command line options
    proccessCommandLine(sys.argv[1:])

    # If there is seed defined, initializes the random machine
    if seed:
        random.seed(seed)
        np.random.seed(int(seed))

    # Selects file open function (normal / gzipped) to create output file
    fileOpen = open
    if re.search('.gz$',outputFileName):
        fileOpen = gzip.open

    # Creates reads and writes them to fastaq file
    # For each kmer in the fasta file generates a number of reads
    # Following a poisson distribution
    readNumber = 0
    with fileOpen(outputFileName, "w") as fastqFile:
        
        # Consecutive fragmenting of all genome
        if full:
            for chr, initPosition, endPosition, kmer in kmerGenerator(
                            referenceGenomeFileName,
                            readLength,
                            overlap):
                readNumber += 1
                fastqFile.write("@%s.%s:%d-%d\n" % (outputFileName,
                                                    chr,
                                                    initPosition,
                                                    endPosition)
                                )
                fastqFile.write("%s\n" % (kmer))
                fastqFile.write("+%s.%s:%d-%d\n" % (outputFileName,
                                                    chr,
                                                    initPosition,
                                                    endPosition)
                                )
                fastqFile.write("%s\n" % (re.sub(".","A",kmer)))                
        
        # Stochastic generation of reads
        else:
            for chr, initPosition, endPosition, kmer in kmerGenerator(
                            referenceGenomeFileName,
                            readLength,
                            1):
                # Determines how many times this kmer will appear in the
                # output file. Uses a random Poisson distribution
                repeats = np.random.poisson(float(coverage)/readLength)
                for repeat in range(repeats):
                    readNumber += 1
                    read = kmer

                    # half of the reads are from + strand, half from - strand
                    if (random.choice("+-") == "-"):
                        read = str(Seq(read).reverse_complement())

                    # Genotyping errors
                    if errorRate != 0:
                        tmpList = list(read)
                        for i,nucleotide in enumerate(tmpList):
                            if random.random() < errorRate:
                                choices = re.sub(tmpList[i],'','ACGT')
                                tmpList[i] = random.choice(choices)
                        read = ''.join(tmpList)


                    fastqFile.write("@%s.%d\n" % (outputFileName,readNumber))
                    fastqFile.write("%s\n" % (read))
                    fastqFile.write("+%s.%d\n" % (outputFileName,readNumber))
                    fastqFile.write("%s\n" % (re.sub(".","A",read)))
