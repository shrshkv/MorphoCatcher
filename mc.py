#Dear users of the MorphoCatcher script tool,
#this program is free only for academic users.
#If your company will commercializes the primer sets obtained by the tool,
#please, contact with the developer, Fedor V. Shirshikov.
#Contact e-mail: shrshkv@ya.ru

#!/usr/bin/env python3

import sys
import os
import argparse
import collections
import re

import zipfile
import numpy as np
import matplotlib
matplotlib.use('agg')
import pylab
import seaborn


def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def get_max_length(input_filename):
    with open(input_filename) as input_file:
        lengths = list()
        for line in input_file.readlines():
            if line.strip():
                length = line.rsplit(maxsplit=1)[-1]
                if length.isdigit():
                    lengths.append(int(length))
        return max(lengths)


def get_label(filename):
    label = filename.rstrip('.divis.txt')
    if label.count('.') == 1:
        label = label.replace('.', ' ', 1)
    else:
        label = re.sub(r'\.(\S+)\.', r' [\1] ', label, 1)
    return label

# class for each species
class Species():
    def __init__(self, name='new species'):
        self.name = name # name of species
        self.strainList = [] # list of strains for this species
        self.consensus = [] # output consensus for this species
        self.checkCons = [] # secondary consensus

# class for each species
class Strain():
    def __init__(self, name, start):
        self.name = name # str, name of strain
        self.start = start # int, start position of the strain sequence
        self.seq = [] # list, strain sequence

# return length of the input sequence
def getSeqLen(seq):
    return len(seq) - seq.count('-')

# return data from each alignment string of input file
def parseString(line):
    splitted = line.split()
    species = splitted[0].split('|', 1)[0]
    strain = splitted[0].split('|', 1)[1]
    seq = list(splitted[1])
    position = splitted[2]
    return species, strain, seq, position


# return name of the output consensus file
def getSmallOutFileName(species, args, suffix=''):
    filename = '%s.%s.%stxt' % (species.name.replace('|', '.'), species.strainList[0].name.replace('|', '.'), suffix)
    return os.path.join(args.outputDir, filename)

# return name of the output consensus file without dir for debug
def getSmallOutFileName_without_dir(species, args, suffix=''):
    filename = '%s.%s.%stxt' % (species.name.replace('|', '.'), species.strainList[0].name.replace('|', '.'), suffix)
    return filename

# check that there is no gap in position of consensus
def checkNoGapInSmallOutput(species, k):
    noGap = True
    if species.strainList[0].seq[k] == '-':
        noGap = False
    return noGap

# print alignment block for consensus output file
def printBlockInSmallOutput(species, outFile, start, end):
    for i in range(start, end):
        if checkNoGapInSmallOutput(species, i):
            outFile.write(species.strainList[0].seq[i])
    outFile.write('\n')
    for i in range(start, end):
        if checkNoGapInSmallOutput(species, i):
            outFile.write(species.consensus[i])
    outFile.write('\n')

# print data for consensus output file
def getSmallOutput(species, outFile, args):
    seqLen = len(species.strainList[0].seq)
    i = 0
    for i in range(args.blockLen, seqLen, args.blockLen):
        printBlockInSmallOutput(species, outFile, i - args.blockLen, i)
        outFile.write('\n')
    printBlockInSmallOutput(species, outFile, i, seqLen)

# print each number of gaps in divis output file
def printBlockInGapCountOutput(species, outFile, args, start, end, tab=0):
    tabBetween = 0
    if tab > 0:
        outFile.write(''.ljust(tab))
        tabBetween = args.gapCountLen
    i = start
    for i in range(start + args.gapCountLen, end, args.gapCountLen):
        outFile.write(str(species.consensus[i - args.gapCountLen:i].count('-')).ljust(tabBetween))
        if tab == 0:
            outFile.write('\n')
    outFile.write(str(species.consensus[i:end].count('-')).ljust(tabBetween))
    outFile.write('\n')

# print data for divis output file
def getGapCountOutput(species, outFile, args):
    seqLen = len(species.strainList[0].seq)
    i = 0
    for i in range(args.blockLen, seqLen, args.blockLen):
        printBlockInGapCountOutput(species, outFile, args, i - args.blockLen, i)
    printBlockInGapCountOutput(species, outFile, args, i, seqLen)

# print data for alignment output file
def getBigOutput(speciesList, outFile, args):
    seqLen = len(speciesList[0].strainList[0].seq)

    # get maximal length of description (name of species + name of strain) in order to proper alignment
    maxNameLen = 0
    for species in speciesList:
        for strain in species.strainList:
            nameLen = len('%s|%s' % (species.name, strain.name))
            if nameLen > maxNameLen:
                maxNameLen = nameLen

    # print blocks of output alignment
    i = 0
    for i in range(args.blockLen, seqLen, args.blockLen):
        printBlockInBigOutput(speciesList, outFile, maxNameLen, args, i - args.blockLen, i)
    printBlockInBigOutput(speciesList, outFile, maxNameLen, args, i, seqLen)

# print block of alignment for alignment output file
def printBlockInBigOutput(speciesList, outFile, maxNameLen, args, start, end):
    # set spaces between columns
    tabFirst = 6
    tabSecond = 3

    for species in speciesList:
        for strain in species.strainList:
            combinedName = '%s|%s' % (species.name, strain.name)
            position = strain.start + getSeqLen(strain.seq[0:end])
            outFile.write('%s%s%s%s%d\n' % (combinedName.ljust(maxNameLen), ''.ljust(tabFirst), ''.join(strain.seq[start:end]), ''.ljust(tabSecond), position))
        outFile.write('%s%s%s\n' % (''.ljust(maxNameLen), ''.ljust(tabFirst), ''.join(species.consensus[start:end])))
        printBlockInGapCountOutput(species, outFile, args, start, end, tab=tabFirst + maxNameLen)
    outFile.write('\n\n')

# create consensus output file
def generateSmallOutputFile(speciesList, args):
    for species in speciesList:
        outFile = open(getSmallOutFileName(species, args), 'w', newline = args.lineBreakFormat)
        getSmallOutput(species, outFile, args)
        outFile.close()

# create divis output file
def generateGapCountFile(speciesList, args):
    for species in speciesList:
        outFile = open(getSmallOutFileName(species, args,  suffix='divis.'), 'w', newline = args.lineBreakFormat)
        getGapCountOutput(species, outFile, args)
        outFile.close()

# create alignment output file (Alignment.txt)
def generateBigOutputFile(speciesList, args):
    output_filename = os.path.join(args.outputDir, args.outFileName)
    outFile = open(output_filename, 'w', newline = args.lineBreakFormat)
    getBigOutput(speciesList, outFile, args)
    outFile.close()

# create debugging output file (All.data.txt)
def generateDebugFile(speciesList, args):
    debug_filename = os.path.join(args.outputDir, args.debugFileName)
    outFile = open(debug_filename, 'w', newline=args.lineBreakFormat)

    # write small data
    for species in speciesList:
        outFile.write('%s\n' % (getSmallOutFileName_without_dir(species, args)))
        getSmallOutput(species, outFile, args)
        outFile.write('\n\n\n\n\n')

    # write big data
    getBigOutput(speciesList, outFile, args)

    outFile.close()

# return list of species from alignment of input file
def parseInputFile(args):
    speciesList = []
    inFile = open(args.inFileName)
    for line in inFile:
        line = line.strip()
        if not (line.startswith('CLUSTAL') or line.startswith('*') or line == ''):
            curSpecies, curStrain, seq, position = parseString(line)

            # check if species already in the list, then define index of found species or create new species
            speciesIndex = -1
            for i, species in enumerate(speciesList):
                if curSpecies == species.name:
                    speciesIndex = i
            if speciesIndex == -1:
                newSpecies = Species(curSpecies)
                speciesList.append(newSpecies)
                speciesIndex = len(speciesList) - 1

            # check if strain already in the list, then define index of found strain or create new strain
            strainIndex = -1
            for i, strain in enumerate(speciesList[speciesIndex].strainList):
                if curStrain == strain.name:
                    strainIndex = i
            if strainIndex == -1:
                start = int(position) - getSeqLen(seq)
                newStrain = Strain(curStrain, start)
                speciesList[speciesIndex].strainList.append(newStrain)
                strainIndex = len(speciesList[speciesIndex].strainList) - 1

            # add seq to strain
            speciesList[speciesIndex].strainList[strainIndex].seq.extend(seq)

    return speciesList

# return secondary consensus for each position of species
def strainBpConsensus(bpList):
    if '-' in bpList:
        consensusBp = '>'
    else:
        firstBp = bpList[0]
        consensusBp = firstBp
        for bp in bpList:
            if bp != firstBp:
                consensusBp = 'N'
    return consensusBp

# check is this strain in this position could be distinguished from other strains
def getDiff(diffList):
    diff = ''
    if len(diffList) == 0:
        diff = '?'
    elif False in diffList:
        diff = '*'
    else:
        diff = '-'
    return diff

# check if there is at least two pair of the identical letters in position of secondary alignment
def checkTwoPairs(bpList):
    commonList = collections.Counter(bpList).most_common(2)
    if len(commonList) == 2 and commonList[0][1] >= 2 and commonList[1][1] >= 2:
        atLeastTwoPairs = True
    else:
        atLeastTwoPairs = False
    return atLeastTwoPairs

# return consensus for each position of species
def speciesBpConsensus(bpList):
    bpListLen = len(bpList)
    consensusBpList = ['0' for x in range(bpListLen)]
    if 'N' in bpList:
        consensusBpList = ['*' for x in range(bpListLen)]
    elif checkTwoPairs(bpList):
        consensusBpList = ['?' for x in range(bpListLen)]
    else:
        for i in range(bpListLen):
            diffList = []
            for k in range(bpListLen):
                if k != i:
                    if bpList[i] != bpList[k]:
                        diffList.append(True)
                    else:
                        diffList.append(False)
            consensusBpList[i] = getDiff(diffList)
    return consensusBpList

# add consensus strings for list of species
def getSpeciesConsensus(speciesList):
    # get subsidiary string with secondary strain consensus for each species
    seqLen = len(speciesList[0].strainList[0].seq)
    for species in speciesList:
        for i in range(seqLen):
            bpList = [strain.seq[i] for strain in species.strainList]
            species.checkCons += strainBpConsensus(bpList)

    # get species consensus for each species
    for i in range(seqLen):
        bpList = [species.checkCons[i] for species in speciesList]
        consensusBpList = speciesBpConsensus(bpList)
        for k, species in enumerate(speciesList):
            species.consensus += consensusBpList[k]

    # replace '*' to '>' in case of gaps in strain consensus
    for i in range(seqLen):
        for species in speciesList:
            if species.checkCons[i] == '>':
                species.consensus[i] = '>'

    return speciesList

def generate_diagrams(args):
    output_files_lists = [[], [], []]
    for root, dirs, files in os.walk(args.outputDir):
        for filename in files:
            if filename in ('Alignment.txt', 'Data.txt'):
                output_files_lists[0].append(filename)
            elif filename[-9:] == 'divis.txt':
                output_files_lists[2].append(filename)
            else:
                output_files_lists[1].append(filename)
    output_files_lists[2].sort()
    seaborn.set(style='white', font_scale=2.4)
    matplotlib.rcParams['lines.linewidth'] = 2
    matplotlib.rc('axes', edgecolor='black', linewidth=1.5)
    figure = pylab.figure(figsize=(12, 10), dpi=200)
    ax = figure.add_subplot(111)
    for filename in output_files_lists[2]:
        with open(os.path.join(args.outputDir, filename), 'r') as f:
            l = [float(line) for line in f.readlines()]
            label = get_label(filename)
            ax.plot(range(160, len(l)*20 - 120, 20), moving_average(l, 15), label=label)
    start, end = ax.get_xlim()
    max_length = get_max_length(args.inFileName)

    # set some initial variables to auto generate tick scale
    nmajor = 100
    minticks = 1
    addticks = 2
    majticks = 5

    # auto scale to have less than 15 major ticks on the map
    while nmajor > 8 :
        nmajor = minticks
        minticks = addticks
        addticks = majticks
        majticks = nmajor * 10
        nmajor = int(max_length/majticks)
    ax.xaxis.set_ticks(np.arange(0, end, majticks))
    ax.xaxis.set_ticks(np.arange(0, end, minticks), minor=True)
    ax.xaxis.grid(which='minor', linestyle='dotted', linewidth=1.5)
    ax.xaxis.grid(which='major', linestyle='dashed', linewidth=1.5)

    ax.yaxis.grid(which='major', linestyle='dashed', linewidth=1.5)
    ax.tick_params(axis='both', direction='out', length=5, width=1.5, color='black')
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(0, end, 1))
    ax.set_xlabel('Nucleotide position', labelpad=20)
    ax.set_ylabel('Average mutation index', labelpad=20)
    ax.legend()
    figure.tight_layout()
    figure.savefig(os.path.join(args.outputDir, 'Diagram.png'))
    figure.savefig(os.path.join(args.outputDir, 'Diagram.tiff'))
    figure.savefig(os.path.join(args.outputDir, 'Diagram.svg'))
    figure.savefig(os.path.join(args.outputDir, 'Diagram.eps'))
    figure.savefig(os.path.join(args.outputDir, 'Diagram.pdf'))


    return output_files_lists

def main():
    #[options]
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                      '--inFileName',
                      type = str,
                      help = ' [str], no default value')
    parser.add_argument('-o',
                      '--outFileName',
                      type = str,
                      default = 'Alignment.txt',
                      help = ' [str], default value = "Alignment.txt"')
    parser.add_argument('-d',
                      '--debugFileName',
                      type = str,
                      default = 'All.data.txt',
                      help = ' [str], default value = "All.data.txt"')
    parser.add_argument('-b',
                      '--blockLen',
                      type = int,
                      default = 60,
                      help = 'Length of output alignment block [int], default value = 60')
    parser.add_argument('-g',
                      '--gapCountLen',
                      type = int,
                      default = 20,
                      help = 'Length of alignment block in which we count gaps [int], default value = 20')
    parser.add_argument('-f',
                      '--lineBreakFormat',
                      type = str,
                      default = 'u',
                      help = 'Type of line breaks: "u" for unix "\\n", "w" for windows "\\r\\n" [str], default value = "u"')
    parser.add_argument('-D',
                      '--outputDir',
                      type = str,
                      default = 'result',
                      help = 'output dir name')
    args = parser.parse_args()

    if args.outputDir:
        try:
            if not os.path.exists(args.outputDir):
                os.makedirs(args.outputDir)
        except:
            print('Please, enter correct output directory name')
            sys.exit(1)

    if args.inFileName == None:
        print('Please, define --inFileName parameter. To read help use -h. Program is broken.')
        sys.exit(1)

    if args.lineBreakFormat == 'u':
        args.lineBreakFormat = '\n'
    elif args.lineBreakFormat == 'w':
        args.lineBreakFormat = '\r\n'
    else:
        print('Type of line breaks: "u" for unix "\\n", "w" for windows "\\r\\n". Please, set another --lineBreakFormat parameter. To read help use -h. Program is broken.')
        sys.exit(1)

    if args.gapCountLen > args.blockLen:
        print('gapCountLen should be less or equal to blockLen. Please, define another --gapCountLen parameter. To read help use -h. Program is broken.')
        sys.exit(1)

    # END_OF [options]

    speciesList = parseInputFile(args)
    speciesList = getSpeciesConsensus(speciesList)
    generateSmallOutputFile(speciesList, args)
    generateGapCountFile(speciesList, args)
    generateBigOutputFile(speciesList, args)
    generateDebugFile(speciesList, args)
    generate_diagrams(args)

    return 0
# def main

if __name__ == '__main__':
    sys.exit(main())
