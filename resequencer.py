#!/usr/bin/env python2
# -*- coding: latin-1 -*-
'''
resequencer.py

    Generates a fastq file simulating a resequencing of an individual with
    the reference genome
'''

__author__ = "Carlos Morcillo Suarez"
__license__ = "GPL"
__version__ = "1.0"

import sys
import re
import random
import getopt

# Globals -----------------------------------------
referenceGenomeFileName = ''
cnvFileName = ''
sampleName = 'sample'
coverage = 15
averageLength = 100
errorRate = 0.000
seed = None

# Functions -----------------------------------------

def usage():
    print '''
        resequencer.py

        Generates a simulated fastq file with reads from a NGS experiment
        from a reference genome in fasta format

        Use

            resequencer.py [-r,--reference] genome.fa [options]

        Options

            -n, --name
                name of the fastaq generated file (default: "sample")

            -c, --coverage
                read coverage of the genome (default: 15)

            -e, --error-rate
                error rate of reads respect to the reference genome
                (default: 0)

            -l, --length
                length of reads to be generated

            -s, --seed
                a seed for the random generator,

            -d, --duplications
                a file describing the CNVs that will be randomly introduced

    '''

# Process command line options
def proccessCommandLine(argv):
    try:
        opts, remainder = getopt.getopt(
                                argv,
                                "r:n:c:e:s:d:l:",
                                ["reference=", "name=", "coverage=",
                                    "error-rate=", "seed=","duplications=",
                                    "length="])

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
        elif opt in  ("-n", "--name"):
            global sampleName
            sampleName = arg
        elif opt in ("-c", "--coverage"):
            global coverage
            coverage = float(arg)
        elif opt in ("-e", "--error-rate"):
            global errorRate
            errorRate = float(arg)
        elif opt in ("-l", "--length"):
            global averageLength
            averageLength = int(arg)
        elif opt in ("-s", "--seed"):
            global seed
            seed = arg
        elif opt in ("-d", "--duplications"):
            global cnvFileName
            cnvFileName = arg

# Returns complementary strand 'AATCC' -> 'GGATT'
def changeReadStrand(sequence):
    equiv = { "A":"T",
              "C":"G",
              "G":"C",
              "T":"A",
              "N":"N",
              "a":"t",
              "c":"g",
              "g":"c",
              "t":"a",
              "n":"n",
            }
    return "".join([equiv[char] for char in sequence][::-1])

def getRandomPosition(sequence):
    return(random.randint(0,len(sequence)-1))

if __name__ == '__main__':

    # Process command line options
    proccessCommandLine(sys.argv[1:])

    # If there is seed defined, initializes the random machine
    if seed:
        random.seed(seed)

    # Reads reference genome
    genome = []
    chromosome = ""
    chromosomeNumber = -1
    with open(referenceGenomeFileName,"r") as referenceGenomeFile:
        for line in referenceGenomeFile:
            if re.search("^>",line):
                chromosomeNumber += 1
                genome.append(chromosome)
            else:
                genome[chromosomeNumber] += line.strip()

    # All the reference sequence of the organism as a string
    # with '*'s separating chromosomes
    genomeSequence = "*".join(genome)


    # If there is a CNV configuration file, CNV variation is introduced
    if cnvFileName != '':
        cnvConfig = []
        with open(cnvFileName,'r') as cnvFile:
            # Uploads CNV configuration
            for line in cnvFile:
                line = line.strip()
                if re.search('^#',line) or len(line) == 0:
                    continue
                cnvConfig.append([int(value) for value in line.split()])

            # Inserts CNVs into genomeSequence
            for length, repeats in cnvConfig:
                randomPosition = getRandomPosition(genomeSequence)
                cnvSequence = genomeSequence[randomPosition:randomPosition+length]
                for repeat in range(repeats):
                    insertionPoint = getRandomPosition(genomeSequence)
                    genomeSequence = genomeSequence[:insertionPoint]+\
                                     cnvSequence+\
                                     genomeSequence[insertionPoint:]


    # Creates random reads
    numberOfReads = coverage * len(genomeSequence) / averageLength
    fastqFileName = sampleName+".fastq"
    with open(fastqFileName, "w") as fastqFile:
        for readNumber in range(int(numberOfReads)):
            randomPosition = getRandomPosition(genomeSequence)
            read = genomeSequence[ randomPosition : randomPosition+averageLength ]

            # if sequence contains '*' (chromosome separator) removes
            # the tail of the read
            read = re.sub("\*.*$","",read)

            # if the read is empty, discard it
            # (this would happen if randomPosition fell into a positon
            # containing '*')
            if len(read) == 0:
                continue

            # half of the reads are from + strand, half from - strand
            if (random.choice("+-") == "-"):
                read = changeReadStrand(read)

            # Genotyping errors
            tmpList = list(read)
            for i,nucleotide in enumerate(tmpList):
                if random.random() < errorRate:
                    choices = re.sub(tmpList[i],'','ACGT')
                    tmpList[i] = random.choice(choices)
            read = ''.join(tmpList)


            fastqFile.write("@%s.%d\n" % (sampleName,readNumber))
            fastqFile.write("%s\n" % (read))
            fastqFile.write("+%s.%d\n" % (sampleName,readNumber))
            fastqFile.write("%s\n" % (re.sub(".","A",read)))
