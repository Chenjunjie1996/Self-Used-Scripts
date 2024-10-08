### make pie-plot of CDR3aa
---
### cell calling for flv_cr result
- usage: python cell_calling.py --cr_dir --match_dir
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230321beida_cellcall

对is_cell=False的细胞 calling回组装出多条高可信度的productive链的barcode不被判定为细胞，条件:
1. 满足productive=True
2. high_confidence=True
3. groupby barcode chain，取最高umi的链，仅call双链细胞
4. calling回的细胞合并至filtered annotation文件
5. outdir: 07.cell_calling 
---
### cell calling more cells for flv_cr result
- usage: python cell_calling_more.py --cr_dir --match_dir --method --coeff --seqtype
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/tmp/20230331/T1-wyz-T, /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230309ningxiaBCR_downsample/7g/SS_02
- Added to celescope https://github.com/singleron-RD/CeleScope/pull/234

对 is_cell=False的细胞 calling部分主要由于组装出多条productive链而不被判定为细胞的barcode。
条件:
1. 满足productive=True
2. 若满足heavy-light唯一配对，保留
3. 对heavy or light 链>=2条的barcode，进行背景过滤(SNR, AUTO, NOT_FILTER)
4. outdir: 07.cell_calling 
---
### Human V region + mouse C region
- example: /SGRNJ06/randd/USER/cjj/celedev/jicuiyaokang/
### Human V region + Rat C region
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230423hsV_ratC
### Rat reference
- path: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230423hsV_ratC/ref/vdj_IMGT_rattus
---
### Downsample rawdata
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230427/downsample20g/downsample_script
---
### 计算单链，双链，多链数量
- usage: python cal_chain_cell.py T/BCR
---
### Add 33bp to converted r1
- usage: python addr1.py --mapfile
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230612last33bp
### 统计因多条链被过滤的barcode count
- usage: python count_multi_chain.py
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20230614
### 检测flvCR注释文件区域是否有空值
- usage: python region.py
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20231108by_script/Script
### 过滤clonotype文件，仅保留双链组合(IGH+IGK/L)
- usage: python pair_cl.py
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20231116pair_cl
### clonotype_id	cdr3s_aa	cdr3s_nt	barcode	contig_id	chain	fl_seq
- usage: python combine.py combine_multi.py(single or multi chain clonotypes)
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj10x/20231120/