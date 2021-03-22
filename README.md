# NG-MMAT
# Next-generation Mitochondrial Metagenomics Analysis Toolkit
A pipeline for assembly, annotation, and phylogenomic analysis of mitochondrial sequences derived from next-generation sequencing data.  
The program is designed for metagenomic data but is equally useful for data from individual animals.  

![alt text](https://github.com/joseph7e/NG-MMAT/blob/main/img/diagram-NGMMAT.png?raw=true)

## Installation and Dependencies

```
git clone https://github.com/Joseph7e/NG-MMAT.git
```


* python3 (tested with v3.6.9)  
* SPAdes assembler (tested with v3.13.1)
* BLAST (v2.9.0+)
* BWA-MEM
* Mitos CL (https://gitlab.com/Bernt/MITOS)
* ART-illumina (only required for read simulation)
* PALADIN (for read filtering or assembly-free analysis)


## Usage:

```
usage: ng-mmat.py [-h] -i INPUT [INPUT ...] -o OUTPUT [-c REFERENCE] [-s SAM]
                  [-ass ASSEMBLY] [-r RANK] [-q QUALITY] [-a ALIGN_OPTS] [-f]
                  [-sim]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...]  input sequence(s), or accessions
  -o OUTPUT             output directory
  -c REFERENCE          custom reference, genbank format
  -s SAM                existing SAM
  -ass ASSEMBLY         existing (meta)genome assembly
  -r RANK               taxonomic rank (default: species)
  -q QUALITY            mapping quality filter (default: none)
  -a ALIGN_OPTS         additional aligner options (default: '-t 24')
  -f                    force rebuild of reference
  -sim                  input is assumed as list of headers for read
                        simulation

```

## Mitochondrial-Sequence-Database
By default, the program will construct a reference database of complete mitochondrial genomes and inddividual gene sequences automatically, assuming no current datbase exists (with a last update within a month). The database is stored in ~/.mitobin/.  
The user can override this database construction and use custom databases if chosen.


### Database file descriptions
```bash
#database location
ls ~/.mitobin/
```
reference.faa -> amino acid translations of the protein-coding genes.  
reference.fasta -> full length source sequences (typically complete mitochondrial genomes).  
reference.fna -> nucleotide sequences for all genes with annotations in source files.  
reference.tsv -> metadata and descriptions of database sequences (genetic codes, gene lengths, etc.)  
simulated.fasta -> mitochondrial genomes sequences used to construct read datasets for pipeline testing.  


## Sequence format
NG-MMATT expects the FASTA files in the refreence database to adhere to a specific format.


## Mito-genome annotations

### Supported genetic codes https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
(transl_table=2) - The Vertebrate Mitochondrial Code  
(transl_table=4) - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (Cnidaria/Ctenophora)  
(transl_table=5) - The Invertebrate Mitochondrial Code  
(transl_table=9) - The Echinoderm and Flatworm Mitochondrial Code (Echinodermata, Platyhelminthes)  
(transl_table=13) - The Ascidian Mitochondrial Code (Tunicates)  
(transl_table=14) - The Alternative Flatworm Mitochondrial Code (Platyhelminthes)  
(transl_table=21) - Trematode Mitochondrial Code (Platyhelminthes)    
(transl_table=24) - Rhabdopleuridae Mitochondrial Code (Hemichordata)  
(transl_table=33) - Cephalodiscidae Mitochondrial UAA-Tyr Code (Hemichordata)  
  
