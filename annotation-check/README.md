# Tutorial on correction of gene annotations

This tutorial is designed as a python juptyer notebook. Instruction on how to install and run a juypter notebook can be found here.  
https://test-jupyter.readthedocs.io/en/latest/install.html

Running the program will require little-to-no python experience but you are required to edit certain parts of the code to adjust gene coordinates..



For batch submissions I recommend using tbl2asn to create the submission files.

https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/

## Step 1: Create the submission template file.
 - This includes info about the BioProject, sequence authors, and biosamples (if you do one at a time)

## Step 2: Gather the fasta and tbl annotation files.
The FASTA and tbl need to have the same sequence IDs in the header.
The FASTA should have source identifiers added on.
example; JLN-Belize-A1_1 [organism=Ototyphlonemertes erneba] [location=mitochondrion] [moltype=genomic DNA] [gcode=5]

## Step 3: Validate and create sqn files. (It expects to find a tbl with the same name as the fasta).

tbl2asn -i <FASTA> -r ./ -t template.sbt -V vb -a s -j "[mgcode=5]"

## Step 4: Check the val file for errors and correct. Also see the genbank file.

See the python notebook for a method on how to do this. Start by writing down the genes that need to be addressed.
Running prokka (prokka 4A.fasta -o prokka_4A --gcode 5) can help determine appropriate start/stop

## Step 5: Repeat Step 3 and Step 4 until the TBL is fixed and tbl2asn runs without errors in the validation file.

Move onto submission with bankit.


 
 # Running 






To start you will need to edit the first few lines of the file, these specify the name of the genome and associated gff file.


