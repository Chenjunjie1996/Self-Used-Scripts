import itertools

BARCODES = """TTACTCGC
TCGTTAGC
TACCGAGC
TGTTCTCC
TTCGCACC
TTGCGTAC
TCTACGAC
TGACAGAC
TAGAACAC
TCATCCTA
TGCTGATA
TAGACGGA
TGTGAAGA
TCTCTTCA
TTGTTCCA
TGAAGCCA
TACCACCA
TGCGTGAA
GGTGAGTT
GATCTCTT
GTGTCCTT
GACGGATT
GCAACATT
GGTCGTGT
GAATCTGT
GTACATCT
GAGGTGCT
GCATGGCT
GTTAGCCT
GTCGCTAT
GGAATGAT
GAGCCAAT
GCTCCTTG
GTAAGGTG
GAGGATGG
GTTGTCGG
GGATTAGG
GATAGAGG
GTGTGTCG
GCAATCCG
GACCTTAG
GCCTGTTC
GCACTGTC
GCTAACTC
GATTCATC
GTCTTGGC"""

def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def main():
    out = open('/SGRNJ06/randd/USER/cjj/celedev/atac/20240328hammingindex/res.txt', 'w')
    
    a = BARCODES.split('\n')
    c = set(itertools.combinations(a, 8))
    res = []
    for i in c:
        tmp = list(itertools.combinations(i, 2))
        for j in tmp:
            if hamming_distance(j[0], j[1])<5:
                break
        res.append(i)
    
    for i in res:
        out.write(f'{i}\n')


if __name__ == '__main__':
    main()