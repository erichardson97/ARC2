import tempfile
from Bio import SearchIO, Seq
from Bio.SeqRecord import SeqRecord

record = SeqRecord(Seq.Seq(sequence), id='Nevermind')


clf = SeqClassifier()


with tempfile.NamedTemporaryFile(mode="w") as hmm_out:
    receptor, chain_type = None, None
    clf.run_hmmscan(record, hmm_out)
    hmmer_query = SearchIO.read(hmm_out.name, "hmmer3-text")
hit_table, top_descriptions = clf.parse_hmmer_query(hmmer_query, bit_score_threshold=bit_score_threshold)
clf.prototype_chain_type(top_descriptions)