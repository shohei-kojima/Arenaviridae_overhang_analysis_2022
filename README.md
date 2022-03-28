# Analyze overhangs of Arenaviruses  
  
### Analysis details  
The sequences of arenaviruses in NCBI refseq database was downloaded by the following command on 2022 Jan 07. Then, the overhangs in the downloaded sequences were analyzed by a custom Python script.  
```
# download sequences
esearch -db nucleotide -query "Viruses[Organism] AND srcdb_refseq[PROP] NOT cellular organisms[Organism] AND Arenaviridae[Organism] OR Arenaviridae[All Fields]" | efetch -format fasta > Arenaviridae_210107_nuclotide.fa

# analyze presence of overhangs
python batch_220107_233000.py
```
  
  
### Software versions  
- Python 3.7.4  


# Analyze genome length of Nairoviridae
### Analysis details  
The sequences of nairoviruses in NCBI refseq database was downloaded by the following command on 2022 Jan 16. Then, the overhangs in the downloaded sequences were analyzed by a custom Python script.  
```
# download sequences
esearch -db nucleotide -query "Viruses[Organism] AND srcdb_refseq[PROP] NOT cellular organisms[Organism] AND Nairoviridae[Organism] OR Nairoviridae[All Fields]" | efetch -format fasta > Nairoviridae_210116_nuclotide.fa

# analyze presence of overhangs
python batch_220116_175252.py
```
  
