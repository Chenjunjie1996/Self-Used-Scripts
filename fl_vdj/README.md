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
