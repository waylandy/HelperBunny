import sys
import pandas as pd
import pyhmmer

"""

HMMforage.to_text(
    'database.fa',
    'profiles.hmm',
    'hits.csv',            # hits, includes a2m alignment
    'fullseq.fa',          # full sequence of hits
    e_value = 0.01,        # cutoff e-value
    threads = 6,           # number of threads
    chunk_size = 10000000, # sequence chars per chunk
    db_size = None)        # number of sequences, for calculating e-value

"""

class utils:
    
    @staticmethod
    def get_fasta_size(fa_file):
        db_size = 0
        with open(fa_file, 'r') as r:
            for line in r:
                if line.startswith('>'):
                    db_size += 1
                    if db_size % 200 == 0:
                        sys.stderr.write(f'{db_size}\r')
        sys.stderr.write(f'{db_size}\n')
        return db_size
    
    @staticmethod
    def fasta_chunker(seq_file, chunk_size=0):
        with pyhmmer.easel.SequenceFile(seq_file, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seqs:
            if chunk_size <= 1:
                _chunk_seqs = list(seqs)
                yield len(_chunk_seqs), _chunk_seqs
                return

            _chunk_size, _chunk_seqs = 0, []

            for n, seq in enumerate(seqs):
                _chunk_size += seq.sequence.shape[0]
                _chunk_seqs += [seq]

                if _chunk_size > chunk_size:                
                    yield n + 1 , _chunk_seqs
                    _chunk_size, _chunk_seqs = 0, []

            yield n + 1 , _chunk_seqs

class HMMforage:
    
    @staticmethod
    def hmmsearch_dispatch(hmms_obj, seqs_obj, **params):
        dispatcher = pyhmmer.hmmer.hmmsearch(hmms_obj, seqs_obj, **params)
        output = []
        for tophits in dispatcher:
            n_pos = tophits.searched_nodes
            for hit in tophits:
                seq_name = hit.name.decode()
                seq_description = '' if hit.description is None else hit.description.decode()
                evalue = hit.evalue
                domain = hit.best_domain
                hmm_name = domain.alignment.hmm_name.decode()
                hmm_accession = domain.alignment.hmm_accession.decode()
                hmm_from, hmm_to = domain.alignment.hmm_from, domain.alignment.hmm_to
                # ali_from, ali_to = domain.alignment.target_from, domain.alignment.target_to
                env_from, env_to = domain.env_from, domain.env_to
                a2m = ((hmm_from - 1)*'-') + domain.alignment.target_sequence + ('-'*(n_pos - hmm_to))
                output += [(
                    hmm_name,
                    hmm_accession,
                    seq_name,
                    seq_description,
                    evalue,
                    env_from,
                    env_to,
                    a2m,
                )]
        columns = ['hmm_name', 'hmm_accession', 'seq_name', 'seq_description', 'evalue', 'env_from', 'env_to', 'a2m']
        return pd.DataFrame(output, columns=columns)
    
    @staticmethod
    def retrieve_sequences(seqs_obj, names_descriptions):
        output = []
        for seq in seqs_obj:
            name = seq.name.decode()
            description = seq.description.decode()
            if tuple([name, description]) in names_descriptions:
                output += [[name, description, seq.textize().sequence]]
        columns = ['name', 'description', 'sequence']
        return pd.DataFrame(output, columns=columns)
    
    @staticmethod
    def to_text(
        input_fa,
        input_hmm,
        output_csv,            # hits, includes a2m alignment
        output_fa,             # full sequence of hits
        e_value = 0.01,        # cutoff e-value
        threads = 6,           # number of threads
        chunk_size = 10000000, # sequence chars per chunk
        db_size = None,        # number of sequences, for calculating e-value
        ):

        if db_size == None:
            db_size = utils.get_fasta_size(input_fa)

        w_hits = open(output_csv, 'w')
        w_fullseq = open(output_fa, 'w')

        columns = ['hmm_name', 'hmm_accession', 'seq_name', 'seq_description', 'evalue', 'env_from', 'env_to', 'a2m']
        w_hits.write(','.join(columns)+'\n')

        for n, seqs in utils.fasta_chunker(input_fa, chunk_size=chunk_size):
            with pyhmmer.plan7.HMMFile(input_hmm) as hmms:
                hits = HMMforage.hmmsearch_dispatch(hmms, seqs, cpus=threads, E=e_value, Z=db_size)
                hits['evalue'] = hits['evalue'].apply(lambda x: f'{x:.1e}')

                names_descriptions = set(map(tuple,hits[['seq_name', 'seq_description']].values))
                fullseqs = HMMforage.retrieve_sequences(seqs, names_descriptions)

                w_hits.write(hits.to_csv(None, index=False, header=False))
                w_fullseq.write(''.join(f'>{h} {d}\n{s}\n\n' for h, d, s in fullseqs.values))

                sys.stderr.write(f'{n}\r')
        sys.stderr.write(f'\n')

        w_hits.close()
        w_fullseq.close()

    @staticmethod
    def to_db(
        input_fa,
        input_hmm,
        output_db,
        e_value = 0.01,
        threads = 6,
        chunk_size = 10000000,
        db_size = None,
        ):
        return # do smarter way to store it or somethin


