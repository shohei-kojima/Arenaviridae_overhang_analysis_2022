#!/usr/bin/env python


import os,sys,glob,gzip
import string


start='CGCA'
start_pos_threshold=8
length_threshold=1000
len_start=len(start)


def parse_fasta(path_to_file):
    tmp={}
    seq=[]
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=''.join(seq)
                header=line.strip().replace('>', '')
                seq=[]
            elif '>' in line and not seq:
                header=line.strip().replace('>', '')
            else:
                seq.append(line.strip())
        tmp[header]=''.join(seq)
    return tmp

def complement(seq):
    return seq.translate(str.maketrans('ATGCatgc', 'TACGtacg'))[::-1]

def format_overhang(seq):
    for pos in range(len(seq) - len_start + 1):
        if seq[pos:pos + len_start] == start:
            overhang=seq[:pos]
            break
    return ['%s|%s' % (overhang, seq[pos:]), '%s|%s' % (overhang, seq[pos:pos + len_start]), overhang]


f='Arenaviridae_210107_nuclotide.fa'
fa=parse_fasta(f)
out=['refseq\tfasta_header\tseq(5...3)\t5end\t3end(complement)\t5overhang\t3overhang(complement)\t5overhang_length\t3overhang_length\n']
for h in fa:
    seq=fa[h].upper().replace('U', 'T')
    if len(seq) < length_threshold:
        continue
    cseq=complement(seq)
    if start in seq[:start_pos_threshold] and start in cseq[:start_pos_threshold]:
        l_format_full,l_format_trimmed,l_overhang=format_overhang( seq[:start_pos_threshold])
        r_format_full,r_format_trimmed,r_overhang=format_overhang(cseq[:start_pos_threshold])
        r_format_full=complement(r_format_full)
        out.append('%s\t%s\t%s...%s\t%s\t%s\t%s\t%s\t%d\t%d\n' % (h.split(' ')[0], h, l_format_full, r_format_full, l_format_trimmed, r_format_trimmed, l_overhang, r_overhang, len(l_overhang), len(r_overhang)))

print(len(out))
with open('arenaviridae_overhangs.tsv', 'w') as outfile:
    outfile.write(''.join(out))
