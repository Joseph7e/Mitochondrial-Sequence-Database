#! /usr/bin/env python3
# Authors: Joseph Sevigny and Anthony Westbrook
# Copyright 2021, University of New Hampshire

import argparse
import gzip
import multiprocessing
import os
import resource
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from urllib import request
from Bio import SeqIO

LOG_ERROR, LOG_WARN, LOG_INFO = range(3)
PALADIN_EXEC = "paladin"
ART_EXEC = "art_illumina"
ARTFIX_EXEC = "artFullHeaders.py"
SPADES_EXEC = "spades.py"
RANK_LIST = ["superkingdom", "kingdom", "phylum", "class", "order", "suborder", "infraorder", "superfamily", "family", "genus", "species"]
RANK_ORDER = {rank: order for (order, rank) in enumerate(RANK_LIST)}
# All ranks used by NCBI: ["superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "superorder", "order", "parvorder", "suborder", "infraorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma"]
MITO_GENES = ['12S', '16S', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'NAD1', 'NAD2', 'NAD3', 'NAD4', 'NAD4L', 'NAD5', 'NAD6']

class FileStore:
    """ Manage temporary and cache files """
    FTYPE_TEMP, FTYPE_CACHE, FTYPE_OUTPUT, FTYPE_USER = range(4)
    FOPT_NORMAL, FOPT_GZIP, FOPT_GZIP_DECOMPRESS, FOPT_TAR = range(4)

    # Entries and system generated temporary path
    _entries = dict()
    _output_base = ""
    _cache_path = ""
    _temp_path = ""
    _output_path = ""
    _expire_age = 0
    _close_target = 0
    _open_count = 0

    @classmethod
    def init(cls, output_base, cache_path, output_path, expire_age, close_target):
        """ Initialize the file store  """
        cls._output_base = os.path.expanduser(output_base)
        cls._cache_path = os.path.expanduser(cache_path)
        cls._temp_path = tempfile.mkdtemp(prefix=output_base)
        cls._output_path = os.path.expanduser(output_path)
        cls._expire_age = expire_age
        cls._close_target = close_target

        # Create normal directories
        if not os.path.exists(cls._output_path):
            os.makedirs(cls._output_path)
        if not os.path.exists(cls._cache_path):
            os.makedirs(cls._cache_path)

        # Populate entries
        cls._populate()

    @classmethod
    def destroy(cls):
        """ Destroy the file store """
        # Close all open files
        for group in cls._entries.values():
            for entry in group.values():
                if entry.handle:
                    entry.handle.close()
                    cls._open_count -= 1

        # Delete temporary directorty
        shutil.rmtree(cls._temp_path)

    @classmethod
    def get_entry(cls, name, group=None):
        """ Get a file entry from the store """
        if group:
            return cls._entries[group][name]
        else:
            return cls._entries[name][name]

    @classmethod
    def get_group(cls, group):
        """ Get all entries for a group """
        return cls._entries[group].values()

    @classmethod
    def check_max(cls):
        """ Check if maximum number of files are open, and if so, close percentage of files """
        if cls._open_count > resource.getrlimit(resource.RLIMIT_NOFILE)[0] - 5:
            target_count = int(cls._open_count * cls._close_target)

            # Close open write files until target count has been reached
            for group in cls._entries:
                for name in cls._entries[group]:
                    entry = cls.get_entry(name, group)

                    if entry.mode and ("r" in entry.mode or "+" in entry.mode):
                        continue

                    # Close file and reduce count
                    if entry.handle:
                        entry.handle.close()
                        entry.handle = None
                        cls._open_count -= 1

                    if cls._open_count < target_count:
                        break

    @classmethod
    def _populate(cls):
        """ Add all entries to the file store """
        # group, name, filename, url, info
        FileStore("input", "input", "input.fq", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("forward", "forward", "forward.fq", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("reverse", "reverse", "reverse.fq", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("assembly", "assembly", "ngmmat-assembly.fasta", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("sim_for", "sim_for", "simulated-1.altered.fq", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("sim_rev", "sim_rev", "simulated-2.altered.fq", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("sim_ref", "sim_ref", "simulated.fasta", None, cls.FTYPE_CACHE, cls.FOPT_NORMAL)
        FileStore("reference_faa", "reference_faa", "reference.faa", None, cls.FTYPE_CACHE, cls.FOPT_NORMAL)
        FileStore("reference_fna", "reference_fna", "reference.fna", None, cls.FTYPE_CACHE, cls.FOPT_NORMAL)
        FileStore("reference_fasta", "reference_fasta", "reference.fasta", None, cls.FTYPE_CACHE, cls.FOPT_NORMAL)
        FileStore("reference_genes", "reference_genes", "reference.tsv", None, cls.FTYPE_CACHE, cls.FOPT_NORMAL)
        FileStore("sam", "sam", "{0}.sam".format(cls._output_base), None, cls.FTYPE_OUTPUT, cls.FOPT_NORMAL)
        FileStore("synonyms", "synonyms", "synonyms.tsv", None, cls.FTYPE_USER, cls.FOPT_NORMAL)
        FileStore("lineage", "lineage", "taxdump.tar.gz", "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", cls.FTYPE_TEMP, cls.FOPT_TAR)
        FileStore("lineage-names", "lineage-names", "names.dmp", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("lineage-nodes", "lineage-nodes", "nodes.dmp", None, cls.FTYPE_TEMP, cls.FOPT_NORMAL)
        FileStore("mito-gbk", "mito-gbk1", "mitochondrion.1.genomic.gbff.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz", cls.FTYPE_TEMP, cls.FOPT_GZIP)
        FileStore("mito-gbk", "mito-gbk2", "mitochondrion.2.genomic.gbff.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.genomic.gbff.gz", cls.FTYPE_TEMP, cls.FOPT_GZIP)
        # FileStore("mito-gff", "mito-gff1", "mito1.gpff.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.gpff.gz", cls.FTYPE_TEMP, cls.FOPT_GZIP)
        # FileStore("mito-gff", "mito-gff2", "mito2.gpff.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.protein.gpff.gz", cls.FTYPE_TEMP, cls.FOPT_GZIP)
        # FileStore("mito-seq", "mito-seq1", "mito1.faa.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.faa.gz", cls.FTYPE_TEMP, cls.FOPT_GZIP)
        # FileStore("mito-seq", "mito-seq2", "mito2.faa.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.protein.faa.gz", cls.FTYPE_TEMP, cls.FOPT_GZIP)


    def __init__(self, group, fid, name, url, ftype, options):
        self.fid = fid
        self.url = url
        self.ftype = ftype
        self.options = options
        self.handle = None
        self.mode = None

        # Set full path
        self.set_path(name, ftype)

        # Add current entry to file store
        FileStore._entries.setdefault(group, dict())
        FileStore._entries[group][fid] = self

    def get_path(self):
        if self.options != FileStore.FOPT_GZIP_DECOMPRESS:
            return self.path

        # Check if file already decompressed (partial/interrupted possible)
        decompressed = self.path[:-3]
        if os.path.exists(decompressed):
            return decompressed
        else:
            return self.path

    def set_path(self, name, ftype):
        """ Calculate and set the full path associated with the file """
        if ftype == FileStore.FTYPE_TEMP:
            self.path = os.path.join(FileStore._temp_path, name)
        if ftype == FileStore.FTYPE_CACHE:
            self.path = os.path.join(FileStore._cache_path, name)
        if ftype == FileStore.FTYPE_OUTPUT:
            self.path = os.path.join(FileStore._output_path, name)
        if ftype == FileStore.FTYPE_USER:
            self.path = name

    def get_handle(self, mode=None):
        """ Return file handle, and (re)open if mode specified """
        # Check if we've reached max number of open files
        FileStore.check_max()

        # Check if handle never opened and no mode specified
        if not mode and not self.mode:
            return None

        # Check for existing handle
        if not mode and self.handle:
            return self.handle

        # (Re)open file, change previous write to append
        cur_mode = mode if mode else self.mode.replace("w", "a")

        # Close if previously open
        if self.handle:
            self.handle.close()
            FileStore._open_count -= 1

        if self.options == FileStore.FOPT_GZIP:
            self.handle = gzip.open(self.get_path(), cur_mode)
        else:
            self.handle = open(self.get_path(), cur_mode)

        # Record mode and increase open count
        self.mode = cur_mode
        FileStore._open_count += 1
        
        return self.handle

    def prepare(self):
        """ Download and/or decompress the file """
        # Download file if URL present
        if self.url:
            request.urlretrieve(self.url, self.path)

        # Decompress if gzip marked for decompression
        if self.options == FileStore.FOPT_GZIP_DECOMPRESS:
            with gzip.open(self.path, "rb") as handle_in:
                with open(self.path[:-3], "wb") as handle_out:
                    handle_out.write(handle_in.read())

            # Remove old file, update path to decompressed file
            os.remove(self.path)
            self.path = self.path[:-3]

        # Extract if an archive
        if self.options == FileStore.FOPT_TAR:
            handle = tarfile.open(self.path, "r")
            handle.extractall(os.path.dirname(self.path))

    def exists(self):
        """ Check if file exists """
        return os.path.exists(self.get_path())

    def expired(self):
        """ Check if file is expired """
        mtime = os.path.getmtime(self.get_path())
        age = time.time() - mtime
        return (age / 86400) > FileStore._expire_age


def log(message, level):
    """ Render log messages """
    pre = ["ERROR", "WARN", "INFO"][level]
    stamp = time.strftime("%Y%m%d %H:%M:%S", time.localtime())

    print("{0} [{1}]: {2}".format(pre, stamp, message), flush=True)

    # Quit if error
    if level == 0:
        sys.exit(1)


def parse_args():
    """ Parse command line arguments """
    # Calculate default alignment options
    cores = multiprocessing.cpu_count()
    align_default = "-t {0}".format(cores)

    # Parse arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-i", dest="input", type=str, nargs="+", required=True, help="input sequence(s), or accessions")
    arg_parser.add_argument("-o", dest="output", type=str, required=True, help="output directory")
    arg_parser.add_argument("-c", dest="reference", type=str, help="custom reference, genbank format")
    arg_parser.add_argument("-s", dest="sam", type=str, help="existing SAM")
    arg_parser.add_argument("-ass", dest="assembly", type=str, help="existing (meta)genome assembly")
    arg_parser.add_argument("-r", dest="rank", type=int, default=0, help="taxonomic rank (default: species)")
    arg_parser.add_argument("-q", dest="quality", type=int, default=0, help="mapping quality filter (default: none)")
    arg_parser.add_argument("-a", dest="align_opts", type=str, default=align_default, help="additional aligner options (default: '{0}')".format(align_default))
    arg_parser.add_argument("-f", dest="force", action="store_true", help="force rebuild of reference")
    arg_parser.add_argument("-sim", dest="simulate", action="store_true", help="input is assumed as list of headers for read simulation")

    return arg_parser.parse_args()


def prepare_workflow(args):
    """ Prepare values dependent on workflow path """
    # Existing reference skips reference generation, existing SAM skips all both
    work_ref, work_align, work_assemble  = True, False, True

    if args.sam:
        work_ref, work_align = False, False
    if args.reference:
        work_ref = False
    if args.assembly:
        work_assemble = False

    # Directly use input file if a single input
    if len(args.input) == 1:
        FileStore.get_entry("input").set_path(args.input[0], FileStore.FTYPE_USER)

    # Redirect reference to custom path if specified
    if args.reference:
        FileStore.get_entry("reference_fasta").set_path(args.reference, FileStore.FTYPE_USER)

    # Redirect SAM to existing path if specified
    if args.sam:
        FileStore.get_entry("sam").set_path(args.sam, FileStore.FTYPE_USER)

    # Redirect reference to custom path if specified
    if args.assembly:
        FileStore.get_entry("assembly").set_path(args.assembly, FileStore.FTYPE_USER)



    return work_ref, work_align, work_assemble


def download_dbs(force):
    """ Download sequence, annotation, and taxonomy data from NCBI """
    ref_entry = FileStore.get_entry("reference_fasta")

    if not force:
        # Check if reference exists and is still fresh
        if ref_entry.exists() and not ref_entry.expired():
            log("Fresh reference still exists, skipping preparation", LOG_INFO)
            return False

    log("Downloading databases", LOG_INFO)

    # Download every file with a URL
    for group in FileStore._entries:
        for entry in FileStore.get_group(group):
            entry.prepare()

    return True


def populate_synonyms():
    """ Create {synonym: official gene name} lookup, eg {CO1: COX1} """
    synonym_lookup = dict()

    for line in FileStore.get_entry("synonyms").get_handle("rt"):
        synonym, standard = line.rstrip().split("\t")
        synonym_lookup[synonym] = standard

    return synonym_lookup


def populate_accession():
    """ Create {accession: (gene name, tax ID)} lookup, eg {YP_003058231: (ND1, 123685) """
    accession_lookup = dict()

    log("Populating accession lookup", LOG_INFO)

    for entry in FileStore.get_group("mito-gff"):
        handle = entry.get_handle("rt")

        # Parse the annotations for features of interest
        for record in SeqIO.parse(handle, "genbank"):
            acc_id = record.annotations['accessions'][0]
            gene, tax_id = "unknown", None

            # Search for CDS and source features within record
            for feature in record.features:
                if feature.type == "CDS":
                    # Check for valid gene qualifier
                    if "gene" not in feature.qualifiers: continue
                    gene = feature.qualifiers["gene"][0]

                if feature.type == "source":
                    # Check for valid cross-reference qualifier
                    if "db_xref" not in feature.qualifiers: continue
                    for db in feature.qualifiers["db_xref"]:
                        if db.startswith("taxon:"):
                            tax_id = db.split(":")[1]

            # Ensure tax_id was found
            if tax_id is None: continue

            accession_lookup[acc_id] = (gene, tax_id)

    return accession_lookup


def populate_taxonomy():
    """ Create {taxid: (parent_id, rank, name)} lookup, eg {6: (335928, "genus", "Azorhizobium")} """
    taxonomy_lookup = dict()
    node_handle = FileStore.get_entry("lineage-nodes").get_handle("rt")
    name_handle = FileStore.get_entry("lineage-names").get_handle("rt")

    log("Populating taxonomy lookup", LOG_INFO)

    # Populate lineage lookup
    for node_line in node_handle:
        node_fields = [x.strip() for x in node_line.split("|")]

        # Find corresponding entry in names table
        while True:
            name_line = name_handle.readline()
            name_fields = [x.strip() for x in name_line.split("|")]

            if int(name_fields[0]) > int(node_fields[0]):
                log("Corresponding name not found for taxonomy node {0}".format(node_fields[0]), LOG_ERROR)

            if (name_fields[0] == node_fields[0]) and (name_fields[3] == "scientific name"):
                break

        taxonomy_lookup[node_fields[0]] = tuple(node_fields[1:3] + [name_fields[1].replace(" ", "_")])

    return taxonomy_lookup


def generate_lineage(tax_id, taxonomy_lookup):
    """ Generate a lineage starting at the requested taxonomy ID """
    lineage = ['unknown_' + x for x in RANK_LIST]
    current_id = tax_id

    # Traverse up the lineage
    while True:
        # If current ID is unknown, return lineage as is
        if current_id not in taxonomy_lookup:
            return lineage

        # Check for top of tree
        entry = taxonomy_lookup[current_id]
        if entry[0] == "1":
            break

        # Record name at supported rank
        if entry[1] in RANK_ORDER:
            lineage[RANK_ORDER[entry[1]]] = entry[2]

        # Iterate to parent
        current_id = entry[0]

    return lineage


def write_reference(taxonomy_lookup, synonym_lookup):
    """ Create modified reference from original protein fasta and new headers """
    output_handle = FileStore.get_entry("reference_fasta").get_handle("wt")
    output_handle_faa = FileStore.get_entry("reference_faa").get_handle("wt")
    output_handle_fna = FileStore.get_entry("reference_fna").get_handle("wt")
    output_handle_genes = FileStore.get_entry("reference_genes").get_handle("wt")

    log("Parsing GenBank File and Creating References (see ~/.mitobin/)", LOG_INFO)
    genome_lookup = dict()
    for entry in FileStore.get_group("mito-gbk"):
        handle = entry.get_handle("rt")
        for record in SeqIO.parse(handle, "genbank"):
            accession = record.id
            info = record.annotations['source']
            full_sequence = record.seq
            tax_id = None

            # Identify taxonomic lineage
            for feature in record.features:
                if feature.type == "source":
                    # Check for valid cross-reference qualifier
                    if "db_xref" not in feature.qualifiers: continue
                    for db in feature.qualifiers["db_xref"]:
                        if db.startswith("taxon:"):
                            tax_id = db.split(":")[1]

            #Ensure tax_id was found
            if tax_id is None: continue
            lineage = "|".join(generate_lineage(tax_id, taxonomy_lookup))
            if not "|Metazoa|" in lineage: continue

            # Write full genomes. Header format "AccID|GeneID|superkingdom|kingdom|...|species supplementary"
            output_handle.write(">{0}|{1}|{2} {3}\n".format(accession, "source", lineage, info))
            output_handle.write(str(full_sequence) + '\n')

            # Search for CDS and source features within record
            for feature in record.features:
                if feature.type == "CDS":
                    # Check for valid gene qualifier
                    if "gene" not in feature.qualifiers: continue
                    gene = feature.qualifiers["gene"][0]
                    fixed_gene = synonym_lookup.get(gene, "unknown")
                    acc_id = record.annotations['accessions'][0]
                    strand = feature.location.strand
                    start_loc = feature.location.start.position
                    end_loc = feature.location.start.position
                    nuc_sequence = feature.extract(record.seq)
                    translation = feature.qualifiers['translation'][0]
                    #translation_custom = nuc_sequence.translate(table=5, cds=True)

                    # write FAA reference file
                    output_handle_faa.write(">{0}|{1}|{2} {3}\n".format(acc_id, fixed_gene, lineage, info))
                    output_handle_faa.write(translation + '\n')

                    # write FNA nucleotide reference file
                    output_handle_fna.write(">{0}|{1}|{2} {3}\n".format(acc_id, fixed_gene, lineage, info))
                    output_handle_fna.write(str(nuc_sequence) + '\n')

                    # write gene lookup for whole genome
                    # header = accession, gene, location, direction, info
                    header_match = "{}|{}|{}".format(accession, "source", lineage)
                    output_handle_genes.write("{}\t{}\t{}\t{}\t{}\n".format(header_match, acc_id, fixed_gene, start_loc, end_loc, strand))

    return genome_lookup

def simulateReads(coverage='10'):
    '''Simulate reads using an input of accessions, pulls from reference fasta'''
    log("Simulating reads based on list of input headers", LOG_INFO)
    sim_output_handle = FileStore.get_entry("sim_ref").get_handle("wt")

    # open reference and store in dictionary
    reference_lookup = SeqIO.to_dict(SeqIO.parse(FileStore.get_entry("reference_fasta").path, "fasta"))

    # parse wanted accessions and write new fasta file
    for line in FileStore.get_entry("input").get_handle("rt"):
        cur_head = line.rstrip().split(' ')[0]
        if cur_head in reference_lookup.keys():
            sim_output_handle.write(">{}\n{}\n".format(cur_head, reference_lookup[cur_head].seq))
    sim_output_handle.close()

    # run art with new reference
    art_command = "{} -i {} -o simulated- -f {} -s 150 -p -l 250 -sam -m 650".format(ART_EXEC, FileStore.get_entry("sim_ref").path, coverage)
    subprocess.run(art_command, shell=True)

    log("Adding gene location info to ARTs fastq files", LOG_INFO)
    # fix fastq header
    for file in ["simulated-1.fq", "simulated-2.fq"]:
        fix_command = "python3 {} {} {} {}".format(ARTFIX_EXEC, file, "simulated-.sam", FileStore.get_entry("reference_genes").path)
        subprocess.run(fix_command, shell=True)

    # clean up unwanted files
    os.remove("simulated-1.aln"); os.remove("simulated-2.aln")
    os.remove("simulated-2.fq"); os.remove("simulated-1.fq")
    os.remove("simulated-.sam")

    # change the input file to the newly constructed reads
    FileStore.get_entry("forward").set_path("simulated-1.altered.fq", FileStore.FTYPE_USER)
    FileStore.get_entry("reverse").set_path("simulated-2.altered.fq", FileStore.FTYPE_USER)

def assemble_reads(forward, reverse, assemblers=["spades"]):
    """ run assemblies  """
    log("Starting assembly workflow", LOG_INFO)
    assemblies = []
    if "spades" in assemblers:
        log("Starting SPADES assembly", LOG_INFO)
        spades_command = "{} -1 {} -2 {} -o ngmmatt-assembly-spades/".format(SPADES_EXEC, forward, reverse)
        subprocess.run(spades_command, shell=True)

        # check if assembly exists,
        if os.path.exists("ngmmatt-assembly-spades/contigs.fasta"):
            assemblies.append("ngmmatt-assembly-spades/contigs.fasta")

    if "metaspades" in assemblers:
        log("Starting SPADES assembly", LOG_INFO)
        spades_command = "{} -1 {} -2 {} -o ngmmatt-assembly-metaspades/".format(SPADES_EXEC, forward, reverse)
        subprocess.run(spades_command, shell=True)

        # check if assembly exists,
        if os.path.exists("ngmmatt-assembly-metaspades/contigs.fasta"):
            assemblies.append("ngmmatt-assembly-metaspades/contigs.fasta")


    if len(assemblies) > 1:
        print ("MERGING THIS SHIT")
    else:
        assembly = assemblies[0]


    FileStore.get_entry("assembly").set_path(assembly, FileStore.FTYPE_USER)


def identify_mito_contigs(assembly, reference):
    subprocess.run("echo " + assembly, shell=True)


def index_reference():
    """ Create aligner index for the reference """
    command = "{0} index -r 3 {1}".format(PALADIN_EXEC, FileStore.get_entry("reference").path)
    subprocess.run(command, shell=True)


def combine_inputs(input_paths):
    """ Concatenate input sequences into single file """
    # Exit if only a single file
    if len(input_paths) == 1:
        return

    # Get a handle to the temp file to be used as the concatenated input
    output_handle = FileStore.get_entry("input").get_handle("wt")

    # Iterate through all lines of each input
    for input_file in input_paths:
        file_name = os.path.split(input_file)[1]
        FileStore("inputs", file_name, input_file, None, FileStore.FTYPE_USER, FileStore.FOPT_NORMAL)
        input_entry = FileStore.get_entry(file_name, "inputs")

        if not input_entry.exists():
            log("Input file does not exist: " + input_file, LOG_ERROR)

        # Write concatenated output to temp handle
        input_handle = input_entry.get_handle("rt")
        for line in input_handle:
            output_handle.write(line)


def exec_alignment(options):
    """ Align reads to reference """
    # 1. Execute aligner for the given reads, reference, options and output
    ref_path = FileStore.get_entry("reference").path
    input_path = FileStore.get_entry("input").path
    sam_path = FileStore.get_entry("sam").path

    # Build and execute command
    log("Executing sequence alignment:", LOG_INFO)
    command = "{0} align {1} {2} {3} > {4}".format(PALADIN_EXEC, ref_path, input_path, options, sam_path)
    subprocess.run(command, shell=True)


def populate_bins(quality, rank):
    """ Filter SAM and group aligned reads into bins """
    bin_lookup = dict()
    handle = FileStore.get_entry("sam").get_handle("rt")

    for line in handle:
        sam_fields = line.split("\t")

        # Skip headers, unmapped, and low quality
        if line.startswith("@"):
            continue

        if sam_fields[1] == "4":
            continue

        if int(sam_fields[4]) < quality:
            continue

        # Confine rank to maximum if necessary
        ref_fields = sam_fields[2].split("|")
        rank_count = len(ref_fields) - 3
        if rank > rank_count:
            rank = rank_count

        # Save to lookup
        read_id = ":".join(sam_fields[0].split(":")[3:])
        bin_lookup[read_id] = (ref_fields[1], ref_fields[len(ref_fields) - rank - 1])

    return bin_lookup


def write_bins(bin_lookup):
    """ Separate and write reads into appropriate fastq and SAM file """
    input_handles = dict()
    input_handles["fq"] = FileStore.get_entry("input").get_handle("rt")
    input_handles["sam"] = FileStore.get_entry("sam").get_handle("rt")

    for input_type in input_handles:
        for line in input_handles[input_type]:
            write_handle = None

            # Skip SAM headers
            if input_type == "sam" and line.startswith("@"):
                continue

            # Parse read ID
            if input_type == "fq":
                read_id = line[1:].split(" ")[0].rstrip()
            else:
                read_id = ":".join(line.split()[0].split(":")[3:])

            # Strip paired end strand ID
            if len(read_id) > 1 and read_id[-2] == "/":
                read_id = read_id[:-2]

            # Obtain appropriate handle and write
            if read_id in bin_lookup:
                bin_id = "{0}-{1}.".format(*bin_lookup[read_id]) + input_type

                # Add to FileStore and open if new
                if bin_id not in FileStore._entries:
                    write_entry = FileStore(bin_id, bin_id, bin_id, None, FileStore.FTYPE_OUTPUT, FileStore.FOPT_NORMAL)
                    write_entry.get_handle("wt")

                # Write header
                write_handle = FileStore.get_entry(bin_id).get_handle()
                write_handle.write(line)

            # Write/skip subsequent 3 lines for FASTQ files
            if input_type == "fq":
                for _ in range(3):
                    remaining = input_handles["fq"].readline()
                    if write_handle:
                        write_handle.write(remaining)


# Parse arguments
args = parse_args()

# Initialize the file store
FileStore.init("mitobin", "~/.mitobin", args.output, 30, 0.5)

# Prepare workflow steps
work_ref, work_align, work_assembly = prepare_workflow(args)

try:
    # Generate reference (optional)
    if work_ref:
        if download_dbs(args.force):
            #accession_lookup = populate_accession()
            taxonomy_lookup = populate_taxonomy()
            synonym_lookup = populate_synonyms()
            write_reference(taxonomy_lookup, synonym_lookup)
    if args.simulate:
        simulateReads('10')

    if work_assembly:
        assemble_reads(FileStore.get_entry("forward").path, FileStore.get_entry("reverse").path, ['spades', 'metaspades'])

    identify_mito_contigs(FileStore.get_entry("assembly").path, FileStore.get_entry("reference_fasta").path)

    #         index_reference()
    #
    # # Combine inputs (mandatory)
    # combine_inputs(args.input)
    #
    # # Align reads (optional)
    # if work_align:
    #     exec_alignment(args.align_opts)
    #
    # # Populate and write bins (mandatory)
    # bin_lookup = populate_bins(args.quality, args.rank)
    # write_bins(bin_lookup)

finally:
    # Destroy the file store on exit
    FileStore.destroy()
