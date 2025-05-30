### Downsample rawdata
- usage: python downsample.py --fq_path /SGRNJ06/DATA04/23_02/2023_02_22/RD22110301/A0214_4_1_CCRF_adapter1_less3bp_101010_3W_Drop_noT7_1.0chunhua/2023-02-22-83
- example: /SGRNJ06/randd/USER/cjj/celedev/atac/20230224downsample

### atac mkref (customized)
- usage: /SGRNJ06/randd/USER/cjj/ATAC/cellranger-atac mkref --config=config
- example: /SGRNJ06/randd/USER/cjj/celedev/atac/20230328guoying

### Downsample 10X rawdata
- usage: python downsample10X.py --fq_path /SGRNJ06/DATA04/23_04/2023_04_04/P22062806/B21-A/2023-04-04-174 --reads_num 8000000
- example: /SGRNJ06/randd/PROJECT/scATAC/20230427_downsamle_sc

### FragDisPlot
- usage: python fragdisplot.py
- example: /SGRNJ06/randd/PROJECT/scATAC/self_pipe/20240407_K562_3T3_824_sc_rr

### tssPlot
- usage: python tssplot.py mouse/human
- example: /SGRNJ06/randd/PROJECT/scATAC/self_pipe/20240407_K562_3T3_824_sc_rr

### convert(multi_omics)
- usage: python convert.py --atac_path --rna_path
- example: /SGRNJ06/randd/USER/cjj/celedev/atac/20250206multi_omics