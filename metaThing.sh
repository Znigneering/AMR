#!/bin/bash

# stage 1, use metaRNA to predict ribosome RNAs 
split -l 1000000 $1 MetaRNA/rRNA_hmm_fs_wst_v0/input/;
for i in `ls MetaRNA/rRNA_hmm_fs_wst_v0/input/`;
	 
  do MetaRNA/rRNA_hmm_fs_wst_v0/rna_hmm3.py -k bac -i MetaRNA/rRNA_hmm_fs_wst_v0/input/$i -L MetaRNA/rRNA_hmm_fs_wst_v0/HMM3/ -o MetaRNA/rRNA_hmm_fs_wst_v0/output/$i;
   
done;
cat MetaRNA/rRNA_hmm_fs_wst_v0/output/*.seq > test_files/mt.seq;
rm MetaRNA/rRNA_hmm_fs_wst_v0/output/*;
rm MetaRNA/rRNA_hmm_fs_wst_v0/input/*;





