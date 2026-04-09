import os
import sys
import gzip
import shutil
import hashlib

import polars as pl

from Bio import SeqIO

import pyhmmer

#########

{
    ### Options controlling reporting thresholds:
    "E"           : 10.0  , # The per-target E-value threshold for reporting a hit.
    "T"           : None  , # The per-target bit score threshold for reporting a hit. If given, takes precedence over E.
    "domE"        : 10.0  , # The per-domain E-value threshold for reporting a domain hit.
    "domT"        : None  , # The per-domain bit score threshold for reporting a domain hit. If given, takes precedence over domE.

    ### Options controlling inclusion (significance) thresholds:
    "incE"        : 0.01  , # The per-target E-value threshold for including a hit in the resulting TopHits.
    "incT"        : None  , # The per-target bit score threshold for including a hit in the resulting TopHits. If given, takes precedence over incE.
    "incdomE"     : 0.01  , # The per-domain E-value threshold for including a domain in the resulting TopHits.
    "incdomT"     : None  , # The per-domain bit score thresholds for including a domain in the resulting TopHits. If given, takes precedence over incdomE.

    ### Options for model-specific thresholding:
    "bit_cutoffs" : None  , # The model-specific thresholding option to use for reporting hits. With None (the default), use global pipeline options; otherwise pass one of "noise", "gathering" or "trusted" to use the appropriate cutoffs.

    ### Options controlling acceleration heuristics:
    "F1"          : 0.02  , # The MSV filter threshold.
    "F2"          : 1e-3  , # The Viterbi filter threshold.
    "F3"          : 1e-5  , # The uncorrected Forward filter threshold.
    "bias_filter" : True  , # Whether or not to enable composition bias filter. Defaults to True.

    ### Other expert options:
    "null2"       : True  , # Whether or not to compute biased composition score corrections. Defaults to True.
    "Z"           : None  , # The effective number of comparisons done, for E-value calculation. Leave as None to auto-detect by counting the number of sequences queried.
    "domZ"        : None  , # The number of significant sequences found, for domain E-value calculation. Leave as None to auto-detect by counting the number of sequences reported.
    "seed"        : 42    , # The seed to use with the random number generator. Pass 0 to use a one-time arbitrary seed, or None to keep the default seed from HMMER.
    "cpus"        : 2     ,
    }

######### utils

class BufferedParquetWriter:

    def __init__(self, dir_name, n_rows=3, mode="w"):

        self.dir_name = dir_name
        self.n_rows = n_rows
        self.dataframe = None

        match mode:
            case "w":
                if os.path.exists(self.dir_name):
                    shutil.rmtree(self.dir_name)
                os.makedirs(self.dir_name, exist_ok=True)
                self.partition_num = 0
                self.schema = None
            case "a":
                if os.path.exists(self.dir_name):
                    self.partition_num = max(int(i.split(".")[0]) for i in os.listdir(dir_name))
                    self.schema = pl.scan_parquet(self.dir_name).collect_schema()
                else:
                    os.makedirs(self.dir_name, exist_ok=True)
                    self.partition_num = 0
                    self.schema = None
            case _:
                raise Exception()

    def get_partition_name(self):

        self.partition_num += 1
        return f"{self.dir_name}/{self.partition_num:0>8}.parquet"
    
    def add_dataframe(self, dataframe):

        if self.schema is None:
            self.schema = dataframe.schema
        if self.dataframe is None:
            self.dataframe = dataframe
        else:
            assert dataframe.schema == self.schema
            self.dataframe = pl.concat([self.dataframe, dataframe])
        
        *partitions, self.dataframe = self.dataframe.iter_slices(self.n_rows)
        for partition in partitions:
            partition.write_parquet(self.get_partition_name())

    def dump_buffer(self):

        if self.dataframe.shape[0] != 0:
            self.dataframe.write_parquet(self.get_partition_name())
            self.dataframe = self.schema.to_frame()

######### hmmsearch

def get_size_from_fasta(fasta_file):

    is_gzipped = lambda x: open(x, "rb").read(2) == b'\x1f\x8b'
    handle = gzip.open(fasta_file, "rt") if is_gzipped(fasta_file) else open(fasta_file, "r")
    for n, _ in enumerate(SeqIO.FastaIO.SimpleFastaParser(handle), 1):
        if n % 1000 == 0:
            sys.stdout.write(f"{n}\r")
    sys.stdout.write(f"{n}\r")
    return n

def get_chunks_from_fasta(fasta_file, max_size=100_000):

    max_size = max(max_size, 1)
    buffer, size = [], 0
    alphabet = pyhmmer.easel.Alphabet.amino()
    is_gzipped = lambda x: open(x, "rb").read(2) == b'\x1f\x8b'
    handle = gzip.open(fasta_file, "rt") if is_gzipped(fasta_file) else open(fasta_file, "r")

    for header, sequence in SeqIO.FastaIO.SimpleFastaParser(handle):
        sequence_hash = hashlib.sha256((header + sequence).encode()).hexdigest()
        buffer += [(sequence_hash, header, sequence,
            pyhmmer.easel.TextSequence(name=sequence_hash.encode(), sequence=sequence).digitize(alphabet))]
        size += len(sequence)
        if size > max_size:
            yield buffer
            buffer, size = [], 0
    if size != 0:
        yield buffer

def get_hits_from_chunk(chunk, hmms, e_value = 1e-2, db_size = None, cpus = 2):

    *sequences, digital_sequences = zip(*chunk)
    sequences = pl.DataFrame(sequences, schema=pl.Schema({
        "$hash"    : pl.String,
        "header"   : pl.String,
        "sequence" : pl.String,
        }))
    
    buffer = []
    for hits in pyhmmer.hmmer.hmmsearch(hmms, digital_sequences, cpus=cpus, E=e_value, Z=db_size):
        for hit in hits:
            for domain in hit.domains:
                buffer += [{
                    "$hash"                   : domain.alignment.target_name,

                    "hmm_name"                : domain.alignment.hmm_name,                # The name of the query HMM.
                    "hmm_accession"           : domain.alignment.hmm_accession,           # The accession of the query, or its name if it has none.

                    "target_sequence"         : domain.alignment.target_sequence,         # The sequence of the target sequence in the alignment.
                    "hmm_sequence"            : domain.alignment.hmm_sequence,            # The sequence of the query HMM in the alignment.
                    "identity_sequence"       : domain.alignment.identity_sequence,       # The identity sequence between the query and the target.
                    "posterior_probabilities" : domain.alignment.posterior_probabilities, # Posterior probability annotation of the alignment.

                    "strand"                  : domain.strand,                            # The strand where the domain is located.
                    "target_from"             : domain.alignment.target_from,             # The start coordinate of the alignment in the target sequence.
                    "target_to"               : domain.alignment.target_to,               # The end coordinate of the alignment in the target sequence.
                    "target_length"           : domain.alignment.target_length,           # The length of the target sequence in the alignment.

                    "env_from"                : domain.env_from,                          # The start coordinate of the domain envelope.
                    "env_to"                  : domain.env_to,                            # The end coordinate of the domain envelope.

                    "hmm_from"                : domain.alignment.hmm_from,                # The start coordinate of the alignment in the query HMM.
                    "hmm_to"                  : domain.alignment.hmm_to,                  # The end coordinate of the alignment in the query HMM.
                    "hmm_length"              : domain.alignment.hmm_length,              # The length of the query HMM in the alignment.

                    "i_evalue"                : domain.i_evalue,                          # The independent e-value for the domain.
                    "c_evalue"                : domain.c_evalue,                          # The conditional e-value for the domain.
                    "score"                   : domain.score,                             # The overall score in bits, null2-corrected.
                    "pvalue"                  : domain.pvalue,                            # The p-value of the domain bitscore.
                    "bias"                    : domain.bias,                              # The null2 score contribution to the domain score.
                    "correction"              : domain.correction,                        # The null2 score when calculating a per-domain score.
                    "envelope_score"          : domain.envelope_score,                    # The forward score in the envelope, without null2 correction.

                    "included"                : domain.included,                          # Whether this domain is marked as included.
                    "reported"                : domain.reported,                          # Whether this domain is marked as reported.
                    }]
    
    hits = pl.DataFrame(buffer, schema=pl.Schema({
        "$hash"                   : pl.String,
        "hmm_name"                : pl.String,
        "hmm_accession"           : pl.String,
        "target_sequence"         : pl.String,
        "hmm_sequence"            : pl.String,
        "identity_sequence"       : pl.String,
        "posterior_probabilities" : pl.String,
        "strand"                  : pl.String,
        "target_from"             : pl.Int64,
        "target_to"               : pl.Int64,
        "target_length"           : pl.Int64,
        "env_from"                : pl.Int64,
        "env_to"                  : pl.Int64,
        "hmm_from"                : pl.Int64,
        "hmm_to"                  : pl.Int64,
        "hmm_length"              : pl.Int64,
        "i_evalue"                : pl.Float64,
        "c_evalue"                : pl.Float64,
        "score"                   : pl.Float64,
        "pvalue"                  : pl.Float64,
        "bias"                    : pl.Float64,
        "correction"              : pl.Float64,
        "envelope_score"          : pl.Float64,
        "included"                : pl.Boolean,
        "reported"                : pl.Boolean,
        }))

    hits = (sequences
        .join(
            hits
                .group_by("$hash", maintain_order=True)
                .agg(pl.struct(pl.all().exclude("$hash")).alias("hits"))
            , on="$hash", how="right")
        .select(
            pl.col("header"),
            pl.col("sequence"),
            pl.col("hits"),
            )
        )
    
    return hits

def get_hmm_hits_parquets(
    input_fasta,
    profile_hmm,
    output_parquet,
    fasta_buffer_size = 300000,
    parquet_partition_size = 1000,
    e_value = 1e-2,
    db_size = None,
    cpus = 2,
    ):

    ### parameters
    db_size = get_size_from_fasta(input_fasta) if db_size == None else db_size
    hmms = [hmm for hmm in pyhmmer.plan7.HMMFile(profile_hmm)]
    assert len(hmms) == len({hmm.name for hmm in hmms})

    ### staging
    output_prefix = f"{output_parquet}.TEMP"
    os.makedirs(output_prefix, exist_ok=True)
    
    ### output parts
    n_sequences, n_hits, n_domains = 0, 0, 0
    writer = BufferedParquetWriter(output_prefix, n_rows=parquet_partition_size, mode="w")
    for chunk in get_chunks_from_fasta(input_fasta, max_size=fasta_buffer_size):
        table = get_hits_from_chunk(chunk, hmms, e_value = e_value, db_size = db_size, cpus = cpus)
        writer.add_dataframe(table)

        n_sequences += len(chunk)
        n_hits += table.shape[0]
        n_domains += table.explode("hits").shape[0]
        sys.stdout.write(f"{db_size} total; {n_sequences} searched; {n_hits} sequences ({n_domains} hits)\r")
    
    sys.stdout.write(f"\n")
    writer.dump_buffer()

    ### merge parts
    pl.read_parquet(output_prefix).write_parquet(output_parquet)
    shutil.rmtree(output_prefix)

def get_mapped_hits_with_context(hits_parquet, extend_by_gaps=False, extend_terms=(5, 5)):

    merge = (pl.scan_parquet(hits_parquet)
        .explode("hits")
        .unnest("hits")
        .with_columns(
            pl.concat_list( # a2m sequence
                pl.lit("-").repeat_by(pl.col("hmm_from").sub(1)),
                pl.col("target_sequence"),
                pl.lit("-").repeat_by(pl.col("hmm_length").sub(pl.col("hmm_to"))),
                ).list.join("").alias("$a2m"),
            pl.col("sequence") # nterm sequence
                .str.slice(
                    0, 
                    pl.col("target_from").sub(1))
                .alias("$nterm"),
            pl.col("sequence") # unaligned hit
                .str.slice(
                    pl.col("target_from").sub(1),
                    pl.col("target_to").sub(pl.col("target_from")).add(1))
                .alias("$hit"),
            pl.col("sequence") # cterm sequence
                .str.slice(
                    pl.col("target_to"), 
                    pl.col("sequence").str.len_chars())
                .alias("$cterm"),
            )
        .with_columns(
            pl.concat_list( # validate full mapping
                    pl.col("$nterm"),
                    pl.col("$hit"),
                    pl.col("$cterm"))
                .list.join("")
                .eq(pl.col("sequence"))
                .alias("$check1"),
            pl.col("$a2m") # validate a2m mapping
                .str.to_uppercase()
                .str.replace_all("-", "", literal=True)
                .eq(pl.col("$hit"))
                .alias("$check2"),
            pl.when(extend_by_gaps) # nterm extension length
                .then(pl.col("$a2m").str.extract("^(-*)").str.len_chars())
                .otherwise(0)
                .add(extend_terms[0])
                .clip(0, pl.col("$nterm").str.len_chars())
                .cast(pl.Int64)
                .alias("$nlen"),
            pl.when(extend_by_gaps) # cterm extension length
                .then(pl.col("$a2m").str.extract("(-*)$").str.len_chars())
                .otherwise(0)
                .add(extend_terms[1])
                .clip(0, pl.col("$cterm").str.len_chars())
                .cast(pl.Int64)
                .alias("$clen"),
            )
        .with_columns(
            pl.concat_str( # extended a2m
                pl.col("$nterm")
                    .str.slice(
                        pl.col("$nlen").sub(pl.col("$nterm").str.len_chars()).abs(),
                        pl.col("$nlen"))
                    .str.to_lowercase(),
                pl.col("$a2m"),
                pl.col("$cterm")
                    .str.slice(
                        0,
                        pl.col("$clen"))
                    .str.to_lowercase(),
                ).alias("$a2m_extended")
            )
        .with_columns(
            pl.col("$a2m_extended") # occupied positions count
                .str.count_matches("[A-Z]", literal=False)
                .cast(pl.Int64)
                .alias("$n_occ"),
            pl.col("$a2m_extended") # aligned positions count
                .str.count_matches("[A-Z-]", literal=False)
                .cast(pl.Int64)
                .alias("$n_pos"),
            )
        .with_columns(
            pl.col("$n_occ").cast(pl.Float64) # occupancy
                .truediv(pl.col("$n_pos").cast(pl.Float64))
                .alias("$occ"),
            pl.col("target_from") # start position
                .sub(pl.col("$nlen"))
                .alias("$start"),
            pl.col("target_to") # end position
                .add(pl.col("$clen"))
                .alias("$end")            
            )
        .with_columns(
            pl.col("sequence") # validate extended hit
                .str.slice(
                    pl.col("$start").sub(1),
                    pl.col("$end").sub(pl.col("$start")).add(1))
                .eq(pl.col("$a2m_extended")
                    .str.to_uppercase()
                    .str.replace_all("-", "", literal=True))
                .alias("$check3"),
            pl.col("$n_pos") # validate hmm size
                .eq(pl.col("hmm_length"))
                .alias("$check4")
            )
        .with_columns(
            pl.concat_list( # merge validations
                    pl.col("$check1"),
                    pl.col("$check2"),
                    pl.col("$check3"),
                    pl.col("$check4"),
                    )
                .list.all()
                .alias("$valid"),
            )
        .select(
            pl.col("header"),
            pl.col("sequence"),
            pl.col("hmm_name"),
            pl.col("hmm_accession"),
            pl.col("hmm_length"),
            pl.col("$a2m_extended").alias("a2m"),
            pl.col("$start").alias("start"),
            pl.col("$end").alias("end"),
            pl.col("$occ").alias("occupancy"),
            pl.col("i_evalue"),
            pl.col("c_evalue"),
            pl.col("$valid").alias("valid"),
            )
        )
    
    merge = merge.collect()
    assert merge["valid"].to_numpy().all()
    assert merge.group_by("hmm_name").agg(pl.col("hmm_length")).select(pl.col("hmm_length").list.n_unique().eq(1)).to_numpy().all()
    return merge.drop("valid")
