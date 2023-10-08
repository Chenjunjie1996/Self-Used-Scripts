### Split human and mouse barcode
- usage: python split_hs_mu.py {analysis directory}
- example: /SGRNJ06/randd/USER/cjj/celedev/rna/20230516hs_mu

### Output invalid reads
- usage: python invalid_reads.py rawdata_path
- example: /SGRNJ06/randd/USER/cjj/celedev/rna/20230627/3/Mus_0614PZ_LGOT_1_3lib/01.barcode

### Filter bam file
- usage: python filter_bam.py --bam3 --bam5
- example: /SGRNJ06/randd/USER/cjj/celedev/rna/filter_bam/20230801

### human mouse 排污
- usage: python run.py -c 
- example: /SGRNJ06/randd/USER/cjj/celedev/rna/human_mouse/20230815

### lncRNA umi 占比 boxplot 
- usage python lncRNA_boxplot.py --matrix_files
- example: /SGRNJ06/randd/USER/cjj/celedev/rna/lncRNA/script_test

### read count 统计
- usage: python read_count.py --bam3 --bam5 --cb
- example: /SGRNJ06/randd/USER/cjj/celedev/rna/20231007read_count/