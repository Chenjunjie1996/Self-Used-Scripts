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

对 is_cell=False的细胞 calling部分主要由于组装出多条productive链而不被判定为细胞的barcode。
条件:
1. 满足productive=True
2. 若满足heavy-light唯一配对，保留
3. 对heavy or light 链>=2条的barcode，进行背景过滤(SNR, AUTO, NOT_FILTER)
4. outdir: 07.cell_calling 
---
### Human V region + mouse C region
- example: /SGRNJ06/randd/USER/cjj/celedev/jicuiyaokang/
