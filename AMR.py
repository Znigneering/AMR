import json
import pandas as pd
import readFile
import numpy as np
import pysam as ps
import os
from datetime import datetime

class CARD(object):
    def __init__(self):
        data = json.load(open('card/card.json'))
        
        resistance_variants = {}
        rRNAs = {}
        table = {}
        for index in data:
            if index[0] != '_':
                record = data.get(str(index))
                if record.get('model_type_id') == '40295':
                    #snps
                    snps = record.get('model_param').get('snp').get('param_value')
                    #species_name with rRNA_type
                    name = record.get('ARO_name').split()
                    for i in range(0,len(name)):
                        if name[i] == 'rRNA':
                            name = ' '.join(name[0:i])
                            break
                    #sequence
                    accession,fmin,fmax,strand,sequence = next(iter(
                            record.get('model_sequences').get('sequence').values())).get('dna_sequence').values()
                    #ARO
                    ARO_accession = record.get('ARO_accession')
                    #check whether snps are valid
                    item = resistance_variants.get(accession)
                    """
                    isVaild = True
                    for i in snps:
                        item = snps.get(i)
                        nt = item[0]
                        pos = int(item[1:-1]) - 1
                        if sequence[pos] != nt and nt != 'U':
                            isVaild == False
                            break
                    #add snps into the dictionary if they're vaild
                    if isVaild:
                    """
                    #add rRNA into list
                    if item is None:
                        resistance_variants[accession] = CARD.parser_snps(snps,0,ARO_accession,fmin)
                        rRNAs[accession] = [sequence]
                        #fetch fasta from entrez
                        fn = '_'.join(name.split()) + '.fna'
                        os.system('esearch -query '+accession+' -db nuccore | efetch -format fasta > genome/'+fn)
                    else:
                        seqid = 0
                        notFound = True
                        for rRNA in rRNAs[accession]:
                            if rRNA == sequence:
                                notFound = False
                                break
                            seqid += 1
                        if notFound:
                            rRNAs[accession].append(sequence)
                        resistance_variants[accession] = resistance_variants[accession].append(
                            CARD.parser_snps(snps,seqid,ARO_accession,fmin), ignore_index=True)
                    #add matched accession and species name and rRNA type into dictionary
                    if table.get(accession) == None:
                        table[accession] = name
            self.resistance_variants = resistance_variants
            self.rRNAs = rRNAs
            self.table = table
    def parser_snps(snps,seqid,ARO_accession,offset):
        #parser snps into prev nt, position of nt, curr nt, seqid, ARO, count
        df = []
        for i in snps:
            i = snps.get(i)
            df.append([i[0],int(i[1:-1])-1,i[-1],seqid,ARO_accession,0])
        return pd.DataFrame(data = df,columns=['prev','pos','curr','seqid','ARO_accession','total'])

def load_dataset(fp_mt):
    """
    input: file path of mt.seq
    output: a np array have accession, dna sequence, and rRNA type
    """
    return np.array([ i for i in readFile.read_meta('test_files/mt.seq')],
                    [('accession','<U50'),('sequence','<U400'),('type','i')])

def bwa_sequences(prefix,fp_mt):
    """
    input: 
        prefix of index
        file path of meta.seq
    output:
        bwa.sam which has the result of alignment
    """
    os.system('bwa mem bwa_index/'+prefix+' '+fp_mt+' > alignments/bwa.sam')

def build_bwa_index(method,rRNAs,prefix):
    """
    input: 
        method selection
        a list of rRNAs sequences
        prefix of index
    output:
        bwa index which use to do bwa alignment
    """
    if (method == 0):
        reference_file = open('bwa_index/'+prefix+'.fna','w')
        for i in rRNAs.items():
            for index in range(0,len(i[1])):
                reference_file.write('>'+i[0]+':'+str(index)+'\n')
                reference_file.write(i[1][index]+'\n')
        reference_file.close()
        os.system('bwa index -p bwa_index/'+prefix+' bwa_index/'+prefix+'.fna')
    elif (method == 1):
        os.system('cat genome/* > bwa_index/'+prefix+'.fna')
        os.system('bwa index -p bwa_index/'+prefix+' bwa_index/'+prefix+'.fna')

def count_ARO(fp_sam,card):
    """
    input: 
        file path of bwa.sam
        card database
    output:
        count number of ARO in all the sequence
    """
    sam = ps.AlignmentFile('alignments/bwa.sam','r')
    references = [j for i in card.rRNAs.items() for j in i[1]]
    for i in sam:
        if i.reference_id!=-1 and i.flag == 0:
            qseq = i.query_alignment_sequence
            rseq = references[i.reference_id]
            rname = i.reference_name.split(':')[0]
            # 1-based?
            start = i.reference_start
            end = i.reference_end
            if card.table.get(rname) == None:
                variant = card.resistance_variants[rname.split('.')[0]]
            else:
                variant = card.resistance_variants[rname]
            for j in range(0,len(variant)):
                pos = variant.loc[j,'pos']
                if pos < end and pos >= start:
                    if variant.loc[j,'curr'] == 'U' or variant.loc[j,'curr'] == qseq[pos-start]:
                        variant.loc[j,'total'] +=1

if __name__ == '__main__':
    #arguments
    fp_dataset = 'test_files/sim.fna'
    prefix = 'card'
    
    timer = datetime.now()
    #loading card database
    print('Loading card database......')
    start_time = datetime.now()
    card = CARD()
    print('Loading completed ({}) '.format(datetime.now() - start_time ))
    #filter rRNA from dataset
    print('Filtering rRNA......')
    start_time = datetime.now()
    #os.system('bash metaThing.sh '+fp_dataset)
    print('Filtering completed ({}) '.format(datetime.now() - start_time ))
    #build BWA index
    print('Building bwa index......')
    start_time = datetime.now()
    build_bwa_index(1,card.rRNAs,prefix)
    print('Building completed ({}) '.format(datetime.now() - start_time ))
    #align rRNA reads
    print('Aligning sequences(direct alignment)......')
    start_time = datetime.now()
    bwa_sequences(prefix,'test_files/mt.seq')
    print('Alignment completed ({}) '.format(datetime.now() - start_time ))
    #count the number of snps in the reads
    print('Counting snps......')
    start_time = datetime.now()
    count_ARO('alignments/bwa.sam',card)
    print('Counting completed ({}) '.format(datetime.now() - start_time ))
    #save into snps.count
    with open('snps.count','w') as fw:
        for i in card.resistance_variants:
            for j in card.resistance_variants.get(i).iterrows():
                if j[1].total != 0:
                    fw.write(j[1].prev+str(j[1].pos)+j[1].curr+'\t'+j[1].ARO_accession+'\t'+str(j[1].total)+"\n")
    fw.close()
    print('Result is saved to snps.count file(total cost: {})'.format(datetime.now() - timer ))
