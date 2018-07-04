#!usr/bin/env python3

import sys
    
    
def read_meta(fp_fna):
     """
     input: filepath of fna
     output: accessions and fellowing sequences
     """
     with open(fp_fna) as filehandle:
        accession = None
        sequence = ""
        for line in filehandle:
            # removes newline character from the end of the line
            line = line.strip()
            if line.startswith(">"):
                # will be True if accession!=None
                if accession:
                    items = accession.split()
                    yield (items[0], sequence, items[4].split(':')[1].split('S')[0])
                accession = line
                sequence = ""
            else:
                sequence += line.strip()
        if accession:
            items = accession.split()
            yield (items[0], sequence, items[4].split(':')[1].split('S')[0])

def read_fna(fp_fna):
     """
     input: filepath of fna
     output: accessions and fellowing sequences
     """
     with open(fp_fna) as filehandle:
        accession = None
        sequence = ""
        for line in filehandle:
            # removes newline character from the end of the line
            line = line.strip()
            if line.startswith(">"):
                # will be True if accession!=None
                if accession:
                    yield (accession[1:], sequence)
                accession = line
                sequence = ""
            else:
                sequence += line.strip()
        if accession:
            yield (accession, sequence)

def read_labels(fp_labels):
    labels = []
    with open(fp_labels,'r') as fl:
        for line in fl:
            yield line.split('\t')
    fl.close()

def read_gff3(fp):
   """
   input: file path of gff3
   output: a generator of rRNA locations which
                index 0, species id
                index 3, start position
                index 4, end position
                index 6, strand(0 is positive, 16 is negitive)
                index 8, description of rRNA
   """
   with open(fp, 'r') as fl:
      fl.readline()
      for line in fl:
         items = line.split('\t')
         if items[6] == '+':
            items[6] = '0'
         else:
            items[6] = '16'
         items[8]=items[8][0:-1]
         yield items
      fl.close
 
def read_sam(fp):
    """
    input : path of sam file
    output : a generator of reads which
                   index 0, seq id
                   index 1, strand
                   index 2, specie id
                   index 3, start position
                   index 4, flag
                   ....(only use these much information)
    """
    with open(fp) as fl:
        li = []
        for line in fl:     
            if (line[0]!='@'):
                read = line.split('\t')
                for associate in li:
                    if associate[0] == read[2]:
                        read.append(associate[1])
                        break
                yield read
            elif (line[0]=='@'):
                sn = line.split('\t')[1].split(',')[0].split(' ',1)
                seqid = sn[0].split(':')[1]
                species = sn[-1]
                li.append((seqid,species))
    fl.closed

if __name__ == "__main__":
    for i in read_sam2(sys.argv[1]):
        print(i)
