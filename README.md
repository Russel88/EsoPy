# EsoPy
Esoteric python scripts

## mask_fasta
Replace specific parts of a fasta file with N's to mask it in a BLAST (or BLAST-like) search

```sh
./mask_fasta.py input.fasta mask.csv output.fasta
```
The mask.csv file should be a comma-delimeted file with information on where to mask formatted as: fastaID,start,end. Start and end positions are inclusive in the masking

## cdhit_cluster_parse
Parse the clustering output from CDHIT to a tab-delimted table
```sh
./cdhit_cluster_parse.py cdhit.clstr
```

## crispr_repeat_cluster
Cluster BLAST matches of repeat sequences to putative CRISPR arrays
```sh
./crispr_repeat_cluster.py -h
```

## minced_parser
Parse output from minced to get all spacers in a fasta, all repeats in a fasta, or just a tab-delimeted files with positions of crisprs
```sh
./minced_parser.py minced.out tab
./minced_parser.py minced.out repeats
./minced_parser.py minced.out spacers
```
