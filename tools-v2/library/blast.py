import os
import re
import json
import glob
import random
import shutil
import hashlib
import datetime
from functools import lru_cache
from subprocess import call, PIPE, STDOUT

import polars as pl

from Bio import SeqIO
from Bio.Seq import Seq

"""
# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti
# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt

# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastp_application_options
# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options
# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.tblastn_application_options
# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options

# apt-get update && apt-get install -y ncbi-blast+

import polars as pl
import blast

fasta_file = "data/ensembl/Homo_sapiens.GRCh38.pep.all.fa"
db_prefix = "data/blastdb.prot/output"
blast.get_blast_db_from_fasta(fasta_file, db_prefix, dbtype="prot", overwrite=True)

fasta_file = "data/ensembl/Homo_sapiens.GRCh38.cds.all.fa"
db_prefix = "data/blastdb.nucl/output"
blast.get_blast_db_from_fasta(fasta_file, db_prefix, dbtype="nucl", overwrite=True)

# blastp
db_prefix = "data/blastdb.prot/output"
output_dir = "data/test-blast"
sequences = [
    "WYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVC",
    "KRGRIYLKAEVADEKLHVTVRDAKNLIPMDPNGLSDPYVKLKLIPDPKNESKQKTKTIRSTLNPQWNESFTFKLKPSDKDRRLSVEIWDWDRTTRNDFMGSLSFGVSELMKMPASGWY",
    "PYWTRPERMDKKLLAVPAANTVRFRCPAAGNPTPSISWLKNGREFRGEHRIGGIKLRHQQWSLVMESVVPSDRGNYTCVVENKFGSIRQTYTLD",
    ]
hsps_parquet, hits_parquet, querys_parquet = blast.get_hits(sequences, output_dir, db_prefix, program="blastp")
merge = blast.get_mapped_hits_with_context(hsps_parquet, hits_parquet, querys_parquet)
display(merge)

# blastn
db_prefix = "data/blastdb.nucl/output"
sequences = [
    "GAAGGATCTACAAAGTTTATCTTCAGAAATGAACTTTTCTGGGGTGGCACTGGAAATGCCTAAGCTGCATAAGATTTTGGGGCTTGGAACCTTCTGTACGGTAACTTGGGCCAAGGCATTCGAGTGATCGGTGAGAATTTCATAAGGAAAAGTGTGCTTGGGAGGGGAGATC",
    "CAGTGTGGGATCTCCCCTCCCAAGCACACTTTTCCTTATGAAATTCTCACCGATCACTCGAATGCCTTGGCCCAAGTTACCGTACAGAAGGTTCCAAGCCCCAAAATCTTATGCAGCTTAGGCATTTCCAGTGCCACCCCAGAAAAGTTCATTTCTGAAGATAAACTTTGTAGATCCTTC",
    "AATCTGAGCAGAAACAGGTTCAGGAGATCTGAGACCAACAGGTACAAATACTGGAGCTTTGTTTACTTCTTCGGGAGCAAGATGATAGCTTTCTTCCTCATAGTCAGACTCAAAATCGTCAGCATCATCATCTTCTTCCTGGGAAGCCCGAAGAGTCACCAGCATGGCCTTGGAAGCTCTAGGTTCCTTCTTAGAGGTCTTACCCCGAACTCGTTTT"
    "GATCAAAACGAGTTCGGGGTAAGACCTCTAAGAAGGAACCTAGAGCTTCCAAGGCCATGCTGGTGACTCTTCGGGCTTCCCAGGAAGAAGATGATGATGCTGACGATTTTGAGTCTGACTATGAGGAAGAAAGCTATCATCTTGCTCCCGAAGAAGTAAACAAAGCTCCAGTATTTGTACCTGTTGGTCTCAGATCTCCTGAACCTGTTTCTGCTCAGATTGAGGAA",
    "TCTGCTCACCCATGTGGACGTCCTGTTCAGCGACACCTTCACCTCCGCCGGCCTCGACCCTGCAGGCCGCTGCCTGCTCCCCAGGCCCAAGTCCCTTGCGGGCAGCTGCCCCTCCACCCGCCTGCTGACGCTGGAGGAAGCCCAGGCACGCACCCAGGGCCGGCTGGGGACGCCCACGGAGCCCACAACTCCCAAGGCCCCGGCCTCACCTGCGGAAAGGAGGA",
    ]
hsps_parquet, hits_parquet, querys_parquet = blast.get_hits(sequences, output_dir, db_prefix, program="blastn")
merge = blast.get_mapped_hits_with_context(hsps_parquet, hits_parquet, querys_parquet)
display(merge)

# tblastn
db_prefix = "data/blastdb.nucl/output"
sequences = [
    "WYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVC",
    "KRGRIYLKAEVADEKLHVTVRDAKNLIPMDPNGLSDPYVKLKLIPDPKNESKQKTKTIRSTLNPQWNESFTFKLKPSDKDRRLSVEIWDWDRTTRNDFMGSLSFGVSELMKMPASGWY",
    "PYWTRPERMDKKLLAVPAANTVRFRCPAAGNPTPSISWLKNGREFRGEHRIGGIKLRHQQWSLVMESVVPSDRGNYTCVVENKFGSIRQTYTLD",
    ]
hsps_parquet, hits_parquet, querys_parquet = blast.get_hits(sequences, output_dir, db_prefix, program="tblastn")
merge = blast.get_mapped_hits_with_context(hsps_parquet, hits_parquet, querys_parquet)
display(merge)

# blastx
db_prefix = "data/blastdb.prot/output"
sequences = [
    "GAAGGATCTACAAAGTTTATCTTCAGAAATGAACTTTTCTGGGGTGGCACTGGAAATGCCTAAGCTGCATAAGATTTTGGGGCTTGGAACCTTCTGTACGGTAACTTGGGCCAAGGCATTCGAGTGATCGGTGAGAATTTCATAAGGAAAAGTGTGCTTGGGAGGGGAGATC",
    "CAGTGTGGGATCTCCCCTCCCAAGCACACTTTTCCTTATGAAATTCTCACCGATCACTCGAATGCCTTGGCCCAAGTTACCGTACAGAAGGTTCCAAGCCCCAAAATCTTATGCAGCTTAGGCATTTCCAGTGCCACCCCAGAAAAGTTCATTTCTGAAGATAAACTTTGTAGATCCTTC",
    "AATCTGAGCAGAAACAGGTTCAGGAGATCTGAGACCAACAGGTACAAATACTGGAGCTTTGTTTACTTCTTCGGGAGCAAGATGATAGCTTTCTTCCTCATAGTCAGACTCAAAATCGTCAGCATCATCATCTTCTTCCTGGGAAGCCCGAAGAGTCACCAGCATGGCCTTGGAAGCTCTAGGTTCCTTCTTAGAGGTCTTACCCCGAACTCGTTTT"
    "GATCAAAACGAGTTCGGGGTAAGACCTCTAAGAAGGAACCTAGAGCTTCCAAGGCCATGCTGGTGACTCTTCGGGCTTCCCAGGAAGAAGATGATGATGCTGACGATTTTGAGTCTGACTATGAGGAAGAAAGCTATCATCTTGCTCCCGAAGAAGTAAACAAAGCTCCAGTATTTGTACCTGTTGGTCTCAGATCTCCTGAACCTGTTTCTGCTCAGATTGAGGAA",
    "TCTGCTCACCCATGTGGACGTCCTGTTCAGCGACACCTTCACCTCCGCCGGCCTCGACCCTGCAGGCCGCTGCCTGCTCCCCAGGCCCAAGTCCCTTGCGGGCAGCTGCCCCTCCACCCGCCTGCTGACGCTGGAGGAAGCCCAGGCACGCACCCAGGGCCGGCTGGGGACGCCCACGGAGCCCACAACTCCCAAGGCCCCGGCCTCACCTGCGGAAAGGAGGA",
    ]
hsps_parquet, hits_parquet, querys_parquet = blast.get_hits(sequences, output_dir, db_prefix, program="blastx")
merge = blast.get_mapped_hits_with_context(hsps_parquet, hits_parquet, querys_parquet)
display(merge)

"""

######### build

def get_blast_db_from_sequences(sequences, db_prefix, dbtype="prot", makeblastdb="makeblastdb", overwrite=True):

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    fasta_file = f"{temp_name}/input.fasta"

    output_dir = os.path.dirname(db_prefix)
    if os.path.exists(output_dir):
        if overwrite:
            shutil.rmtree(output_dir)
        else:
            raise Exception("")
    
    is_sequence  = lambda x: bool(re.search("^([A-Za-z*]+|[A-Za-z-.]+)$", x))
    sequence_to_unaligned = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "-."))
    assert all(map(is_sequence, sequences))
    with open(fasta_file, 'w') as w:
        for n, sequence in enumerate(sequences):
            w.write(f">{n}\n{sequence_to_unaligned(sequence)}\n")

    ### subprocess
    cmd = f"{makeblastdb} -parse_seqids -dbtype {dbtype} -in {fasta_file} -out {db_prefix}"
    print(cmd)
    call(cmd.split(), stdout=PIPE, stderr=STDOUT)

    ### output : db
    shutil.rmtree(temp_name)
    return db_prefix

def get_blast_db_from_fasta(fasta_file, db_prefix, dbtype="prot", makeblastdb="makeblastdb", overwrite=True):

    ### staging
    output_dir = os.path.dirname(db_prefix)
    if os.path.exists(output_dir):
        if overwrite:
            shutil.rmtree(output_dir)
        else:
            raise Exception("")
    
    ### subprocess
    cmd = f"{makeblastdb} -parse_seqids -dbtype {dbtype} -in {fasta_file} -out {db_prefix}"
    print(cmd)
    call(cmd.split(), stdout=PIPE, stderr=STDOUT)

    ### output : db
    return db_prefix

######### retrieve

def get_sequences_from_db(db_prefix, identifiers, blastdbcmd="blastdbcmd"):
    
    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    input_file = f"{temp_name}/input.txt"
    output_fasta = f"{temp_name}/output.fasta"
    
    with open(input_file, 'w') as w:
        w.write("\n".join(identifiers))
    
    ### subprocess
    cmd = f"{blastdbcmd} -db {db_prefix} -entry_batch {input_file} -out {output_fasta}"
    print(cmd)
    call(cmd.split(), stdout=PIPE, stderr=STDOUT)

    ### parse output
    parser = SeqIO.FastaIO.SimpleFastaParser(open(output_fasta))
    output = list(zip(*[(*header.split(None, 1), sequence) for header, sequence in parser]))
    table = pl.DataFrame(output, schema=pl.Schema({
        "accession" : pl.String,
        "title"     : pl.String,
        "sequence"  : pl.String,
        }))
    shutil.rmtree(temp_name)
    return table

######### search

def get_hits(sequences, output_dir, db_prefix, program="blastn", options=""):

    ### parameters
    output_dir = os.path.normpath(output_dir)
    db_extension = set(i.split('.')[-1][0] for i in glob.glob(f"{db_prefix}*"))
    assert 1 == len(db_extension)

    match program:
        case "blastn":  # nucl query - nucl db
            assert "n" == db_extension.pop()
        case "blastp":  # prot query - prot db
            assert "p" == db_extension.pop()
        case "tblastn": # prot query - nucl db
            assert "n" == db_extension.pop()
        case "blastx":  # nucl query - prot db
            assert "p" == db_extension.pop()

    ### staging
    temp_name = f'{output_dir}/temp'
    os.makedirs(temp_name, exist_ok=True)

    input_fasta = f"{temp_name}/query.fasta"
    output_json = f"{temp_name}/output.json"
    hsps_parquet = f"{output_dir}/hsps.parquet"
    hits_parquet = f"{output_dir}/hits.parquet"
    querys_parquet = f"{output_dir}/querys.parquet"

    get_hash = lambda x: hashlib.sha256((x).encode()).hexdigest()
    is_sequence  = lambda x: bool(re.search("^([A-Za-z*]+|[A-Za-z-.]+)$", x))
    sequence_to_unaligned = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "-."))
    assert all(map(is_sequence, sequences))
    with open(input_fasta, 'w') as w:
        for sequence in sequences:
            w.write(f">{get_hash(sequence)}\n{sequence_to_unaligned(sequence)}\n")

    ### subprocess
    cmd = f"{program} -query {input_fasta} -db {db_prefix} -out {output_json} -outfmt 15 {options.strip()}".rstrip()
    print(cmd)
    call(cmd.split(), stdout=PIPE, stderr=STDOUT)

    ### output : hsps
    with open(output_json, 'r') as r:
        data = json.load(r)
    
    table_hsps = []
    for results in data["BlastOutput2"]:
        search = results["report"]["results"]["search"]
        for hit in search["hits"]:
            table_hsps += [{
                "program"     : results["report"]["program"],
                "version"     : results["report"]["version"],
                "query_id"    : search["query_id"],
                "query_title" : search["query_title"],
                "description" : hit["description"],
                "hsps"        : hit["hsps"],           
                }]

    table_hsps = (pl.LazyFrame(table_hsps, 
        schema=pl.Schema({
            "program"          : pl.String,
            "version"          : pl.String,
            "query_id"         : pl.String,
            "query_title"      : pl.String,
            "description"      : pl.List(pl.Struct({
                'id'           : pl.String,
                'accession'    : pl.String,
                'title'        : pl.String,
                })),
            "hsps"             : pl.List(pl.Struct({
                "num"          : pl.Int64,
                "bit_score"    : pl.Float64,
                "score"        : pl.Int64,
                "evalue"       : pl.Float64,
                "identity"     : pl.Int64,
                "positive"     : pl.Int64,
                "query_from"   : pl.Int64,
                "query_to"     : pl.Int64,
                "query_strand" : pl.String,
                "hit_from"     : pl.Int64,
                "hit_to"       : pl.Int64,
                "hit_strand"   : pl.String,
                "align_len"    : pl.Int64,
                "gaps"         : pl.Int64,
                "qseq"         : pl.String,
                "hseq"         : pl.String,
                "midline"      : pl.String,
                })),
            }))
        .select(
            pl.col("program"),
            pl.col("version"),
            pl.col("query_id"),
            pl.col("query_title"),
            pl.col("description").list.eval(
                pl.element().name.prefix_fields("hit_")),
            pl.col("hsps").list.eval(
                pl.element().name.prefix_fields("hsp_")),
            )
        .explode("description")
        .unnest("description")
        .explode("hsps")
        .unnest("hsps")
        .collect()
        .write_parquet(hsps_parquet)
        )
    
    ### output : hits
    hit_accessions = pl.read_parquet(hsps_parquet, columns=["hit_accession"])["hit_accession"].unique().to_list()
    (get_sequences_from_db(db_prefix, hit_accessions)
        .write_parquet(hits_parquet)
        )
    
    ### output : query
    parser = SeqIO.FastaIO.SimpleFastaParser(open(input_fasta))
    (pl.DataFrame(list(zip(*parser)), schema=pl.Schema({
            "title": pl.String, 
            "sequence": pl.String
            }))
        .write_parquet(querys_parquet)
        )

    ### output
    shutil.rmtree(temp_name)
    return hsps_parquet, hits_parquet, querys_parquet

######### post-processing

def get_mapped_hits_with_context(hsps_parquet, hits_parquet, querys_parquet):

    @lru_cache
    def _get_reverse_complement(sequence):
        return str(Seq(sequence).reverse_complement())

    merge = (pl.scan_parquet(hsps_parquet)
        .join(pl.scan_parquet(hits_parquet) # join context : hits
            .select(
                pl.col("accession").alias("hit_accession"),
                pl.col("sequence").alias("$hit_sequence"),
                ),
            on="hit_accession", how="left")
        .join(pl.scan_parquet(querys_parquet) # join context : querys
            .select(
                pl.col("title").alias(("query_title")),
                pl.col("sequence").alias(("$query_sequence")),
                ),
            on="query_title", how="left")
        .with_columns(
            pl.col("hsp_identity") # percent identity
                .truediv(pl.col("hsp_align_len"))
                .cast(pl.Float64)
                .alias("$hsp_identity"),
            pl.col("hsp_positive") # percent similarity
                .truediv(pl.col("hsp_align_len"))
                .cast(pl.Float64)
                .alias("$hsp_similarity"),
            pl.when(pl.col("program").is_in({"blastn", "blastx"})) # query type
                .then(pl.lit("nucl"))
                .when(pl.col("program").is_in({"blastp", "tblastn"}))
                .then(pl.lit("prot"))
                .otherwise(None)
                .alias("$query_type"),
            pl.when(pl.col("program").is_in({"blastn", "tblastn"})) # hit type
                .then(pl.lit("nucl"))
                .when(pl.col("program").is_in({"blastp", "blastx"}))
                .then(pl.lit("prot"))
                .otherwise(None)
                .alias("$hit_type"),
            )
        .with_columns(
            pl.min_horizontal( # query start
                pl.col("hsp_query_from"),
                pl.col("hsp_query_to"))
                .alias("$hsp_query_start"),
            pl.max_horizontal( # query end
                pl.col("hsp_query_from"),
                pl.col("hsp_query_to"))
                .alias("$hsp_query_end"),
            pl.when(pl.col("$query_type").eq("prot")) # query strand
                .then(pl.lit("."))
                .when(pl.col("$query_type").eq("nucl"))
                .then(pl
                    .when(pl.col("hsp_query_from").lt(pl.col("hsp_query_to")))
                    .then(pl.lit("+"))
                    .when(pl.col("hsp_query_from").gt(pl.col("hsp_query_to")))
                    .then(pl.lit("-"))
                    .otherwise(None))
                .otherwise(None)
                .alias("$hsp_query_strand"),
            pl.min_horizontal( # hit start
                pl.col("hsp_hit_from"),
                pl.col("hsp_hit_to"))
                .alias("$hsp_hit_start"),
            pl.max_horizontal( # hit end
                pl.col("hsp_hit_from"),
                pl.col("hsp_hit_to"))
                .alias("$hsp_hit_end"),
            pl.when(pl.col("$hit_type").eq("prot")) # hit strand
                .then(pl.lit("."))
                .when(pl.col("$hit_type").eq("nucl"))
                .then(pl
                    .when(pl.col("hsp_hit_from").lt(pl.col("hsp_hit_to")))
                    .then(pl.lit("+"))
                    .when(pl.col("hsp_hit_from").gt(pl.col("hsp_hit_to")))
                    .then(pl.lit("-"))
                    .otherwise(None))
                .otherwise(None)
                .alias("$hsp_hit_strand"),
            )
        .with_columns(
            pl.col("$query_sequence") # query segment (intermediate)
                .str.slice(
                    pl.col("$hsp_query_start").sub(1),
                    pl.col("$hsp_query_end").add(1).sub(pl.col("$hsp_query_start")),
                ).alias("$hsp_query_segment"),
            pl.col("$hit_sequence") # hit segment (intermediate)
                .str.slice(
                    pl.col("$hsp_hit_start").sub(1),
                    pl.col("$hsp_hit_end").add(1).sub(pl.col("$hsp_hit_start")),
                ).alias("$hsp_hit_segment"),
            )
        .with_columns(
            pl.when(pl.col("$hsp_query_strand").is_in({"+", "."})) # query segment
                .then(pl.col("$hsp_query_segment"))
                .when(pl.col("$hsp_query_strand").is_in({"-"}))
                .then(pl.col("$hsp_query_segment").map_elements(_get_reverse_complement, return_dtype=pl.String))
                .otherwise(None)
                .alias("$hsp_query_segment"),
            pl.when(pl.col("$hsp_hit_strand").is_in({"+", "."})) # hit segment
                .then(pl.col("$hsp_hit_segment"))
                .when(pl.col("$hsp_hit_strand").is_in({"-"}))
                .then(pl.col("$hsp_hit_segment").map_elements(_get_reverse_complement, return_dtype=pl.String))
                .otherwise(None)
                .alias("$hsp_hit_segment"),
            )
        .with_columns(
            pl.col("hsp_qseq") # validate mapping : query
                .str.to_uppercase()
                .str.replace_all("-", "", literal=True)
                .eq(pl.col("$hsp_query_segment")
                    .str.to_uppercase())
                .alias("$check1"),
            pl.col("hsp_hseq") # validate mapping : hit
                .str.to_uppercase()
                .str.replace_all("-", "", literal=True)
                .eq(pl.col("$hsp_hit_segment")
                    .str.to_uppercase())
                .alias("$check2"),
            
            )
        .with_columns(
            pl.col("hsp_qseq") # validate mapping : query
                .str.to_uppercase()
                .str.replace_all("-", "", literal=True)
                .eq(pl.col("$hsp_query_segment")
                    .str.to_uppercase())
                .alias("$check1"),
            pl.col("hsp_hseq") # validate mapping : hit
                .str.to_uppercase()
                .str.replace_all("-", "", literal=True)
                .eq(pl.col("$hsp_hit_segment")
                    .str.to_uppercase())
                .alias("$check2"),
            pl.col("hsp_midline") # validate length
                .str.len_chars()
                .eq(pl.col("hsp_align_len"))
                .alias("$check3"),
            )
        .with_columns(
            pl.concat_list( # merge validations
                    pl.col("$check1"),
                    pl.col("$check2"),
                    pl.col("$check3"))
                .list.all()
                .alias("$valid"),
            )
    .select(
            pl.col("program"),
            pl.col("version"),
            pl.col("$query_type").alias("query_type"),
            pl.col("$query_sequence").alias("query_sequence"),
            pl.col("$hit_type").alias("hit_type"),
            pl.col("hit_accession"),
            pl.col("hit_title"),
            pl.col("$hit_sequence").alias("hit_sequence"),
            pl.col("hsp_bit_score"),
            pl.col("hsp_score"),
            pl.col("hsp_evalue"),
            pl.col("$hsp_identity").alias("hsp_identity"),
            pl.col("$hsp_similarity").alias("hsp_similarity"),
            pl.col("$hsp_query_start").alias("hsp_query_start"),
            pl.col("$hsp_query_end").alias("hsp_query_end"),
            pl.col("$hsp_query_strand").alias("hsp_query_strand"),
            pl.col("$hsp_hit_start").alias("hsp_hit_start"),
            pl.col("$hsp_hit_end").alias("hsp_hit_end"),
            pl.col("$hsp_hit_strand").alias("hsp_hit_strand"),
            pl.col("hsp_qseq"),
            pl.col("hsp_hseq"),
            pl.col("hsp_midline"),
            pl.col("$valid").alias("valid"),
            )
        )

    merge = merge.collect()
    # assert merge["valid"].to_numpy().all()
    return merge
