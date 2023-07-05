### Split index for bulk_tcr 
- usage: python split_index --mapfile mapfile --whitelist whitelist
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230116multi/hs
---
### Calculate clonotypes intersection of two samples
- usage: python /SGRNJ03/randd/cjj/Script/bulk_tcr/clonotypes_intersection R230201001E_clonetypes.csv R230201001F_clonetypes.csv
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230207multi_index/intersection_clonotypes
---
### vnplot online
- http://jvenn.toulouse.inra.fr/app/example.html
### barplot and heatmap for 96 index
- usage: python /SGRNJ03/randd/cjj/Script/bulk_tcr/plot96index.py --path ./
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230403_96/
### calculate clonotypes(total number and types)
- usage: python /SGRNJ03/randd/cjj/Script/bulk_tcr/calcu_clone.py
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/2023041496tcr
### split mapfile
- usage: python /SGRNJ03/randd/cjj/Script/bulk_tcr/split_mapfile --mapfile mapfile
- example: /SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230424bcr_96/R230414020 
### use plotly
- 20230704plot
