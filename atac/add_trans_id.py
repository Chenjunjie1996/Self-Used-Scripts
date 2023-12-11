import sys


def main():
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            gene_name = parts[-1].split('"')[1]
            del parts[-1]
            parts.append(f'''gene_id "{gene_name}"; transcript_id "{gene_name}"; exon_number "1"; gene_name "{gene_name}"; gene_biotype "protein_coding"; transcript_name "{gene_name}"; transcript_biotype "protein_coding"; exon_id "{gene_name}";''')
            outfile.write("\t".join(parts) + "\n")

if __name__ == '__main__':
    input_filename, output_filename = sys.argv[1], sys.argv[2]
    main()