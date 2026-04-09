import os
import sys
from subprocess import Popen, PIPE
import polars as pl


INPUT_FASTA    = os.environ.get("INPUT_FASTA", None)
OUTPUT_DIR     = os.environ.get("OUTPUT_DIR", "output_cdd")
EVALUE_CUTOFF  = os.environ.get("EVALUE_CUTOFF", 0.01)
DATA_MODE      = os.environ.get("DATA_MODE", "rep") # {'rep', 'std', 'full'}
TARGET_DATA    = os.environ.get("TARGET_DATA", "both") # {'doms', 'feats', 'both'}
NUM_THREADS    = os.environ.get("NUM_THREADS", 1)

RPSBPROC_DIR   = "/tmp/rpsbproc"
OUTPUT_RAW     = f"{OUTPUT_DIR}/rpsbproc.txt"
OUTPUT_DOMAINS = f"{OUTPUT_DIR}/cdd_domains.parquet"
OUTPUT_SITES   = f"{OUTPUT_DIR}/cdd_sites.parquet"

sys.stdout.write(f"\n")
sys.stdout.write(f"INPUT_FASTA    = {INPUT_FASTA}\n")
sys.stdout.write(f"OUTPUT_DIR     = {OUTPUT_DIR}\n")
sys.stdout.write(f"EVALUE_CUTOFF  = {EVALUE_CUTOFF}\n")
sys.stdout.write(f"DATA_MODE      = {DATA_MODE}\n")
sys.stdout.write(f"TARGET_DATA    = {TARGET_DATA}\n")
sys.stdout.write(f"NUM_THREADS    = {NUM_THREADS}\n")
sys.stdout.write(f"\n")
sys.stdout.write(f"RPSBPROC_DIR   = {RPSBPROC_DIR}\n")
sys.stdout.write(f"OUTPUT_RAW     = {OUTPUT_RAW}\n")
sys.stdout.write(f"OUTPUT_DOMAINS = {OUTPUT_DOMAINS}\n")
sys.stdout.write(f"OUTPUT_SITES   = {OUTPUT_SITES}\n")
sys.stdout.write(f"\n")

assert INPUT_FASTA is not None
os.makedirs(OUTPUT_DIR, exist_ok=True)


cmd1 = f"rpsblast+ -query {INPUT_FASTA} -db {RPSBPROC_DIR}/db/Cdd -evalue {EVALUE_CUTOFF} -num_threads {NUM_THREADS} -outfmt 11"
cmd2 = f"rpsbproc -o {OUTPUT_RAW} -e {EVALUE_CUTOFF} -m {DATA_MODE} -t {TARGET_DATA} -d {RPSBPROC_DIR}"
proc1 = Popen(cmd1.split(), stdout=PIPE, text=True)
proc2 = Popen(cmd2.split(), stdin=proc1.stdout, stdout=PIPE, text=True)
proc2.communicate()

cdd_info = f"{RPSBPROC_DIR}/data/cdd.info"
if os.path.isfile(cdd_info):
    line = next(open(cdd_info))
    cdd_version = line.strip().split()[-1]
    assert line.startswith("cdd version")
else:
    cdd_version = None

cddid_tbl = f"{RPSBPROC_DIR}/data/cddid.tbl"
cddid = pl.read_csv(cddid_tbl, comment_prefix="#", separator="\t", quote_char=None, schema=pl.Schema([
    ("PSSM-Id"     , pl.String),
    ("accession"   , pl.String),
    ("short-name"  , pl.String),
    ("description" , pl.String),
    ("PSSM-Length" , pl.String),
    ]))


def parse_rpsbproc(output_file):

    session = []
    query = []
    domains = []
    sites = []

    is_data = False
    is_session = False
    is_query = False
    is_domains = False
    is_sites = False

    with open(output_file, "r") as r:
        for line in r:
            section = None if line.isspace() else next(iter(line.split()))
            match section:
                case "DATA":
                    is_data = True
                case "ENDDATA":
                    is_data = False
                case "SESSION":
                    is_session = True
                    session += [line.strip().split("\t")[1:]]
                case "ENDSESSION":
                    is_session = False
                case "QUERY":
                    is_query = True
                    query += [line.strip().split("\t")[1:]]
                case "ENDQUERY":
                    is_query = False
                case "DOMAINS":
                    is_domains = True
                case "ENDDOMAINS":
                    is_domains = False
                case "SITES":
                    is_sites = True
                case "ENDSITES":
                    is_sites = False
                case _:
                    if is_domains:
                        assert all([is_data, is_session, is_query])
                        domains += [line.strip().split("\t")]
                    if is_sites:
                        assert all([is_data, is_session, is_query])
                        sites += [line.strip().split("\t")]
    
    query = pl.DataFrame(query, orient="row", schema=pl.Schema([
        ("query-id"        , pl.String),
        ("seq-type"        , pl.String),
        ("seq-length"      , pl.String),
        ("definition-line" , pl.String),
        ]))
    
    domains = pl.DataFrame(domains, orient="row", schema=pl.Schema([
        ("session-ordinal"        , pl.String),
        ("query-id[readingframe]" , pl.String),
        ("hit-type"               , pl.String),
        ("PSSM-ID"                , pl.String),
        ("from"                   , pl.String),
        ("to"                     , pl.String),
        ("E-Value"                , pl.String),
        ("bitscore"               , pl.String),
        ("accession"              , pl.String),
        ("short-name"             , pl.String),
        ("incomplete"             , pl.String),
        ("superfamily PSSM-ID"    , pl.String),
        ]))
    
    sites = pl.DataFrame(sites, orient="row", schema=pl.Schema([
        ("session-ordinal"        , pl.String),
        ("query-id[readingframe]" , pl.String),
        ("annot-type"             , pl.String),
        ("title"                  , pl.String),
        ("residue(coordinates)"   , pl.String),
        ("complete-size"          , pl.String),
        ("mapped-size"            , pl.String),
        ("source-domain"          , pl.String),
        ]))
    
    return query, domains, sites

query, domains, sites = parse_rpsbproc(OUTPUT_RAW)

(query
    .join(domains, how="inner", left_on="query-id", right_on="query-id[readingframe]")
    .join(cddid
        .select(
            pl.col("PSSM-Id").alias("superfamily PSSM-ID"),
            pl.col("accession").alias("superfamily accession"),
            pl.col("short-name").alias("superfamily short-name"),
            ),
        how="left", on="superfamily PSSM-ID")
    .select(
        # pl.col("query-id"),
        pl.col("seq-type"),
        pl.col("seq-length").cast(pl.Int64),
        pl.col("definition-line"),
        # pl.col("session-ordinal").cast(pl.Int64),
        pl.col("hit-type"),
        pl.col("PSSM-ID"),
        pl.col("accession"),
        pl.col("short-name"),
        pl.col("superfamily PSSM-ID").replace("-", None).fill_null(pl.col("PSSM-ID")),
        pl.col("superfamily accession").fill_null(pl.col("accession")),
        pl.col("superfamily short-name").fill_null(pl.col("short-name")),
        pl.col("from").cast(pl.Int64),
        pl.col("to").cast(pl.Int64),
        pl.col("E-Value").cast(pl.Float64),
        pl.col("bitscore").cast(pl.Float64),
        pl.col("incomplete"),
        pl.lit(cdd_version).alias("cdd-version"),
        )
    .write_parquet(OUTPUT_DOMAINS)
    )

(query
    .join(sites, how="inner", left_on="query-id", right_on="query-id[readingframe]")
    .select(
        # pl.col("query-id"),
        pl.col("seq-type"),
        pl.col("seq-length"),
        pl.col("definition-line"),
        # pl.col("session-ordinal"),
        pl.col("annot-type"),
        pl.col("title"),
        pl.col("residue(coordinates)"),
        pl.col("complete-size").cast(pl.Int64),
        pl.col("mapped-size").cast(pl.Int64),
        pl.col("source-domain"),
        pl.lit(cdd_version).alias("cdd-version"),
        )
    .write_parquet(OUTPUT_SITES)
    )
