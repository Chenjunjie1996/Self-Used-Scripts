def convert_gff_to_gtf(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            if line.startswith("#") or line.strip() == "":
                continue  # 跳过注释和空行

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue  # 确保数据完整

            # feature_type = parts[2]

            # # 如果特征类型是 CDS, five_prime_UTR, 或 three_prime_UTR，则将其改为 exon
            # if feature_type in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
            #     parts[2] = 'exon'

            # 提取 GFF 属性并将其转换为 GTF 格式
            attributes = parts[8]
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=')
                    attr_dict[key] = value

            # 构建 GTF 格式的属性列
            gtf_attributes = []
            if 'Parent' in attr_dict:
                gtf_attributes.append(f'gene_id "{attr_dict["Parent"]}"; transcript_id "{attr_dict["Parent"]}";')
            if 'ID' in attr_dict:
                gtf_attributes.append(f'gene_id "{attr_dict["ID"]}"; transcript_id "{attr_dict["ID"]}";')

            parts[8] = " ".join(gtf_attributes)
            
            if parts[2] == "mRNA":
                result1, result2, result3 = parts.copy(), parts.copy(), parts.copy()
                result1[2] = "gene"
                result2[2] = "transcript"
                result3[2] = "exon"
                outfile.write("\t".join(result1) + "\n")
                outfile.write("\t".join(result2) + "\n")
                outfile.write("\t".join(result3) + "\n")
            elif parts[2] == "CDS":
                result1, result2 = parts.copy(), parts.copy()
                result1[2] = "exon"
                outfile.write("\t".join(result1) + "\n")
                outfile.write("\t".join(result2) + "\n")
            else:
                outfile.write("\t".join(parts) + "\n")

# 使用函数转换文件
convert_gff_to_gtf("Squaliobarbus_Curriculus.coding.gene.V1.20200416.gff", "Squaliobarbus_Curriculus.gtf")
