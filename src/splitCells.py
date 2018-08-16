#!/usr/bin/env python

#SBATCH --time=72:00:00
#SBATCH --mem=32g

import os
import gzip
from argparse import ArgumentParser
import sys

parser = ArgumentParser(description='split fastq to individual cells')
parser.add_argument('-i', '--inputFile', help='input file', required=True)
parser.add_argument('-o', '--outputDir', help='output directory', required=True)
parser.add_argument('-s', '--sampleName', help='sample name', required=True)
parser.add_argument('-b', '--cellIDFile', help='output bed/bed.gz file', required=True)

options = parser.parse_args()

# input parsing
fastqInput = options.inputFile
fastqDir = options.outputDir + "/fastq"
sampleName = options.sampleName
cellIDFile = options.cellIDFile

print('fastqCombinedFile: ' + fastqInput)
print('fastqDir: ' + fastqDir)
print('sampleName: ' + sampleName)
print('barcodeFile: ' + cellIDFile)

# Open fastq for reading
# os.path.isfile(fastqInput)
# os.path.isfile(barcodeFile)
# sys.exit("Inputs are good")

def parse_fastq(f):
    """parse every 4 lines of fastq file"""
    name = f.readline().strip()
    read = f.readline().strip()
    plus = f.readline().strip()
    qual = f.readline().strip()
    return [name,read,plus,qual]


with open(cellIDFile) as f:
    barcodesKeep = f.read().splitlines()
fastqFile = gzip.open(fastqInput, 'rb')
cellFiles = { k : gzip.open(fastqDir + "/" + k + "." + sampleName + ".fastq.gz" , "wb") \
    for k in  barcodesKeep}

while True:
    lines = parse_fastq(fastqFile)
    if "" in lines: break
    barcode = lines[0][1:33]
    if barcode in barcodesKeep:
        _ = cellFiles[barcode].write(lines[0] + "\n")
        _ = cellFiles[barcode].write(lines[1] + "\n")
        _ = cellFiles[barcode].write(lines[2] + "\n")
        _ = cellFiles[barcode].write(lines[3] + "\n")


fastqFile.close()

for f in cellFiles.values():
    f.close()

