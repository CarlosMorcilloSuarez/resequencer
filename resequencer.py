#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
resequencer.py

    Generates a fastq file simulating a resequencing of an individual with
    the reference genome
'''

__author__ = "Carlos Morcillo Suarez"
__license__ = "GPL"
__version__ = "1.1"

import sys
import re
import random
import getopt
import numpy as np
from Bio.Seq import Seq

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
                length of reads to be generated (default: 100)

            -s, --seed
                a seed for the random generator,

            -d, --duplications <file>
                a file describing the CNVs that will be randomly introduced

            -p, --pair-ended <meanFragmentLength>
                generates two *fastq files corresponding to pair ended readings
                from DNA fragments of meand <meanFragmentLength> length.

    '''


# Process command line options
def proccessCommandLine(argv):
    try:
        opts, remainder = getopt.getopt(
                                argv,
                                "r:n:c:e:s:d:l:p:",
                                ["reference=", "name=", "coverage=",
                                    "error-rate=", "seed=","duplications=",
                                    "length=", "pair-ended="])

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
            global averageReadLength
            averageReadLength = int(arg)
        elif opt in ("-p", "--pair-ended"):
            global averageFragmentLength
            averageFragmentLength = int(arg)
        elif opt in ("-s", "--seed"):
            global seed
            seed = arg
        elif opt in ("-d", "--duplications"):
            global cnvFileName
            cnvFileName = arg


def getRandomPosition(sequence):
    return(random.randint(0,len(sequence)-1))


def introduceGenotypingErrors(sequence,errorRate):
    tmpList = list(sequence)
    for i,nucleotide in enumerate(tmpList):
        if random.random() < errorRate:
            choices = re.sub(tmpList[i],'','ACGT')
            tmpList[i] = random.choice(choices)
    return ''.join(tmpList)


if __name__ == '__main__':

    referenceGenomeFileName = ''
    cnvFileName = ''
    sampleName = 'sample'
    coverage = 15
    averageReadLength = 100
    errorRate = 0.000
    seed = None
    averageFragmentLength = 0

    # Process command line options
    proccessCommandLine(sys.argv[1:])

    # If there is seed defined, initializes the random machine
    if seed:
        random.seed(seed)
        np.random.seed(int(seed))

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
    numberOfReads = coverage * len(genomeSequence) / averageReadLength

    # For single file fastq file
    if averageFragmentLength == 0:
        fastqFileName = sampleName+".fastq"
        with open(fastqFileName, "w") as fastqFile:
            for readNumber in range(int(numberOfReads)):
                randomPosition = getRandomPosition(genomeSequence)
                read = genomeSequence[ randomPosition : randomPosition+averageReadLength ]

                # if sequence contains '*' (chromosome separator) removes
                # the tail of the read
                read = re.sub("\*.*$","",read)

                # if the read is empty, discard it
                # (this would happen if randomPosition fell into a position
                # containing '*')
                if len(read) == 0:
                    continue

                # half of the reads are from + strand, half from - strand
                if (random.choice("+-") == "-"):
                    read = str(Seq(read).reverse_complement())

                # Genotyping errors
                read = introduceGenotypingErrors(read,errorRate)


                fastqFile.write("@%s.%d\n" % (sampleName,readNumber))
                fastqFile.write("%s\n" % (read))
                fastqFile.write("+%s.%d\n" % (sampleName,readNumber))
                fastqFile.write("%s\n" % (re.sub(
                                            "^AAA",
                                            "77A",
                                            re.sub(".","A",read)
                                            )
                                          )
                                )

    # for pair ended fastq files
    else:
        fastqFileName1 = sampleName+"_1.fastq"
        fastqFileName2 = sampleName+"_2.fastq"
        with open(fastqFileName1, "w") as fastqFile1:
            with open(fastqFileName2, "w") as fastqFile2:
                for readNumber in range(int(numberOfReads)):
                    randomPosition = getRandomPosition(genomeSequence)
                    fragmentLength = int(np.random.normal(
                                            averageFragmentLength
                                            ))
                    fragment = genomeSequence[
                        randomPosition : randomPosition+fragmentLength
                        ]

                    # If fragment is too short (pair-end reads overlap) it
                    # is discarded. This may happen when the fragment comes
                    # from the end of the genome
                    if len(fragment) < averageReadLength*2.8:
                        continue

                    # if fragment contains '*' (chromosome separator) discards
                    # the fragment
                    if re.search('\*',fragment):
                        continue

                    readPlus = fragment[:averageReadLength]
                    readMinus = str(
                              Seq(fragment[-averageReadLength:]
                              ).reverse_complement()
                            )

                    # Half of the reads are from + strand, half from - strand
                    if (random.choice("+-") == "+"):
                        read1 = readPlus
                        read2 = readMinus
                    else:
                        read2 = readPlus
                        read1 = readMinus

                    # Genotyping errors
                    read1 = introduceGenotypingErrors(read1,errorRate)
                    read2 = introduceGenotypingErrors(read2,errorRate)


                    fastqFile1.write("@%s_1.%d\n" % (sampleName,readNumber))
                    fastqFile1.write("%s\n" % (read1))
                    fastqFile1.write("+%s_1.%d\n" % (sampleName,readNumber))
                    fastqFile1.write("%s\n" % (re.sub(
                                                "^AAA",
                                                "77A",
                                                re.sub(".","A",read1)
                                                )
                                              )
                                    )

                    fastqFile2.write("@%s_2.%d\n" % (sampleName,readNumber))
                    fastqFile2.write("%s\n" % (read2))
                    fastqFile2.write("+%s_2.%d\n" % (sampleName,readNumber))
                    fastqFile2.write("%s\n" % (re.sub(
                                                "^AAA",
                                                "77A",
                                                re.sub(".","A",read2)
                                                )
                                              )
                                    )
