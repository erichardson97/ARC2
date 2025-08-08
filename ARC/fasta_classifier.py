"""
Classifies input sequences into BCR, TCR, or MHC.
Specifies chain type including constant regions.
Contains code from ANARCI.
"""
import datetime
import numpy as np
import os
import pandas as pd
import re
import subprocess
import tempfile
import multiprocessing as mp
from functools import reduce

from ARC.mhc_G_domain import mhc_G_domain
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from Bio import SeqIO
from io import StringIO



class FastaClassifier:
    """Classifies input sequence/s into BCR, TCR or MHC chains.

    Attributes:
        seqfile (str): Input sequence file in FASTA format

        outfile (str): Name of output file

        threads (int): number of cpu threads to run on

        hmmer_path (str): optional argument containing the path to your HMMER installation

        blast_path (str): optional argument containing the path to your BLAST installation"""

    def __init__(
            self, seqfile=None, outfile=None, recalc_species=True, threads=1, hmmer_path=None, blast_path=None, hmm_path=None,
            speedy=True):
        # Relative paths and IO handling
        self.package_directory = '/mnt/BioAdHoc/Users/erichard/ARC2/ARC/' #os.path.dirname(os.path.abspath(__file__))
        self.recalc_species = recalc_species
        self.seqfile = seqfile
        self.outfile = outfile
        # HMM related scores and files
        self.hmm_score_threshold = 25
        self.hmm_path = hmm_path if hmm_path != None else os.path.join(self.package_directory,
                                                                       "data/HMMs/ALL_AND_C.hmm")
        self.mhc_I_hmm = os.path.join(
            self.package_directory, "data/MHC_HMMs/Pfam_MHC_I.hmm"
        )
        self.mhc_II_alpha_hmm = os.path.join(
            self.package_directory, "data/MHC_HMMs/Pfam_MHC_II_alpha.hmm"
        )
        self.mhc_II_beta_hmm = os.path.join(
            self.package_directory, "data/MHC_HMMs/Pfam_MHC_II_beta.hmm"
        )
        # G domain assignment files and IgNAR database
        self.mro_file = os.path.join(self.package_directory, "data/chain-sequence.tsv")
        self.mro_gdomain_file = os.path.join(
            self.package_directory, "data/MRO_Gdomain.csv"
        )
        self.mro_df = self.get_MRO_Gdomains(self.mro_file)
        self.ignar_db = os.path.join(self.package_directory, "data/IgNAR/IgNAR")
        self.b2m_db = os.path.join(self.package_directory, "data/blastdb/b2m.fasta")
        self.num_threads = threads
        if hmmer_path == None:
            self.hmmer_path = ""
        else:
            self.hmmer_path = hmmer_path
        if blast_path == None:
            self.blast_path = os.path.join(self.package_directory, 'bin/blastp')
        else:
            self.blast_path = blast_path
        self.blast_db_path = os.path.join(self.package_directory, 'data/imgt/blast_fasta')
        self.mmseqs_path = os.path.join(self.package_directory, 'bin/mmseqs')
        self.b2m_db_speedy = os.path.join(self.package_directory, 'data/mmseqsdb/b2m_mmseqs/b2m')
        self.ignar_db_speedy = os.path.join(self.package_directory, 'data/mmseqsdb/IgNAR_mmseqs/IgNAR')
        self.speedy = speedy

    def is_aa(self, seq_record):
        characters = set(str(seq_record.seq))
        if len(characters.intersection(set(['U', 'C', 'A', 'G', 'T']))) == len(characters):
            return False
        return True

    def check_seq(self, seq_record):
        """Checks validity of an amino acid sequence

        Args:
            seq_record: A biopython sequence record object
        """
        pattern = re.compile(
            r"[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]", re.IGNORECASE
        )
        if len(str(seq_record.seq)) > 0 and not pattern.findall(str(seq_record.seq)) and self.is_aa(seq_record):
            return True
        else:
            return False

    def run_cmd(self, cmd, input_string=""):
        """Runs a command using subprocess

        Args:
            cmd: The shell command to run

            input_string: String to pass in via stdin

        Raises:
            Exception: The command has failed to run. Stderr is printed
        """
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            close_fds=True,
            env=dict(os.environ, my_env_prop="value"),
            shell=True,
        )

        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception("Cmd {} failed: {}".format(cmd, stderr))

        return out

    def get_species(self, seq_record, locus = "IG"):

        blast_path = self.blast_db_path

        with tempfile.NamedTemporaryFile(mode="w") as temp_out:
            if not seq_record.seq:
                return False
            SeqIO.write(seq_record, temp_out.name, "fasta")
            path = os.path.join(blast_path, f'{locus}V.fasta')
            subprocess.call(f'{self.blast_path} -query {temp_out.name} -db {path} -evalue 1e-6 -num_threads 4 -out {temp_out.name}blast.txt -outfmt 6', shell =
                            True)
            if os.path.getsize(f'{temp_out.name}blast.txt') == 0:
                top_species = {'species':'none', 'bitscore': 0}
                return top_species
            output = pd.read_csv(f'{temp_out.name}blast.txt', sep = '\t', header = None)
            output.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send','evalue','bitscore']
            output['species'] = output['sseqid'].map(lambda x:x.split('|')[1])
            top_species = output.groupby('sseqid').apply(lambda x:x.loc[x['bitscore'].idxmax()], include_groups=False).sort_values('bitscore', ascending = False).iloc[0]
            return top_species

    def get_species_seqfile(self, seq_file, locus = "IG"):
        db_path = os.path.join(self.blast_db_path, locus + "V.fasta")
        with tempfile.NamedTemporaryFile(mode="w") as temp_out:
            subprocess.call(f'{self.blast_path} -query {seq_file} -db {db_path} -evalue 1e-6 -num_threads 4 -out {temp_out.name} -outfmt 6', shell =
                            True)
            if os.path.getsize(temp_out.name) == 0:
                top_species = {'species':'none', 'score': 0}
                return pd.DataFrame({p.description:top_species for p in SeqIO.parse(seq_file, format='fasta')}).T
            output = pd.read_csv(temp_out.name, sep = '\t', header = None)
            if output.shape[0] == 0:
                top_species = {'species':'none', 'score': 0}
                return pd.DataFrame({p.description:top_species for p in SeqIO.parse(seq_file, format='fasta')}).T
            output.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send','evalue','bitscore']
            output['species'] = output['sseqid'].map(lambda x:x.split('|')[1])
            top_species = output.groupby('qseqid').apply(lambda x: x.loc[x['bitscore'].idxmax()], include_groups=False).reset_index()
            return top_species

    def run_hmmscan(self, seq_records, hmm_out):
        """Runs hmmscan from the HMMER3 software suite

        Args:
            seq_record: A biopython sequence record object

            hmm_out: tempfile object for hmm output
        """
        with tempfile.NamedTemporaryFile(mode="w") as temp_out:
            SeqIO.write([p for p in seq_records if self.check_seq(p)], temp_out.name, "fasta")
            temp_out.flush()
            hmmer = self.hmmer_path + "hmmscan"
            args = [
                hmmer,
                "-o",
                hmm_out.name,
                self.hmm_path,
                temp_out.name,
            ]
            cmd = (" ").join(args)
            self.run_cmd(cmd)
            if not (os.path.exists(hmm_out.name) and os.access(hmm_out.name, os.R_OK)):
                return False
            if os.path.getsize(hmm_out.name) == 0:
                return False
            return True

    def domains_are_same(self, dom1, dom2):
        """Check to see if two domains are overlapping.

        Args:
            dom1:

            dom2:

        Returns:
            True (domains are same) or False (domains are not the same)
        """
        dom1, dom2 = sorted([dom1, dom2], key=lambda x: x.query_start)
        if dom2.query_start >= dom1.query_end:
            return False
        return True

    def parse_hmmer_query(self, query, bit_score_threshold=100):
        """Retrieves top hits of HMMER query

        The function will identify multiple domains if they have been found
        and provide the details for the best alignment for each domain.
        This allows the ability to identify single chain fvs and engineered
        antibody sequences as well as the capability to identify constant
        domains.

        Args:
            query: hmmer query object from Biopython

            bit_score_threshold: the threshold for which to consider a hit a hit.

        Returns:
            hit_table: format of HMMER output

            top_descriptions: which provides a list of the top domain's species and chain types
        """
        hit_table = [
            [
                "id",
                "description",
                "evalue",
                "bitscore",
                "bias",
                "query_start",
                "query_end",
            ]
        ]

        # Find the best hit for each domain in the sequence.

        top_descriptions, domains = [], []

        if query.hsps:  # We have some hits
            # Iterate over the matches of the domains in order of their e-value (most significant first)
            for hsp in sorted(query.hsps, key=lambda x: x.evalue):
                new = True
                # Only look at those with hits that are over the threshold bit-score.
                if hsp.bitscore >= bit_score_threshold:
                    # Check to see if we already have seen the domain
                    for i in range(len(domains)):
                        if self.domains_are_same(domains[i], hsp):
                            new = False
                            break
                    hit_table.append(
                        [
                            hsp.hit_id,
                            hsp.hit_description,
                            hsp.evalue,
                            hsp.bitscore,
                            hsp.bias,
                            hsp.query_start,
                            hsp.query_end,
                        ]
                    )
                    if (
                            new
                    ):  # It is a new domain and this is the best hit. Add it for further processing.
                        domains.append(hsp)
                        # Add the last added to the descriptions list.
                        top_descriptions.append(dict(zip(hit_table[0], hit_table[-1])))
                else:
                    hit_table.append(
                        [
                            hsp.hit_id,
                            hsp.hit_description,
                            hsp.evalue,
                            hsp.bitscore,
                            hsp.bias,
                            hsp.query_start,
                            hsp.query_end,
                        ]
                    )
                    return hit_table, top_descriptions
            # Reorder the domains according to the order they appear in the sequence.
            ordering = sorted(range(len(domains)), key=lambda x: domains[x].query_start)
            domains = [domains[_] for _ in ordering]
            top_descriptions = [top_descriptions[_] for _ in ordering]

        ndomains = len(domains)
        # If any significant hits were identified parse and align them to the reference state.
        for i in range(ndomains):
            domains[i].order = i
            if top_descriptions[i]["id"].count("_") == 1:
                species, chain = top_descriptions[i]["id"].split("_")
                top_descriptions[i]["species"] = species  # Reparse
                top_descriptions[i]["chain_type"] = chain

        return hit_table, top_descriptions



    def get_chain_type(self, top_hits):
        """Retrieves the chain type from the list of top hits

        Args:
            top_hits: the highest scoring hits per domain of an HMMER query

        Returns:
            The chain type as well as the specific domains present in a HMM hit
        """
        ndomains = len(top_hits)
        # set of the top domains - i.e. the unique types
        top_domains_set = {"".join(x["id"].split('_')[1:]) for x in top_hits}
        # set incl. species name
        top_domains = [x["id"].split('_')[0]+'_'+''.join(x['id'].split('_')[1:]) for x in top_hits]
        # species observed per unique
        chain_type_sp = {hit_type: [sp for sp, _ in (obj.split('_') for obj in top_domains) if _ == hit_type]
                         for hit_type in top_domains_set}

        # These sets simplify checking for various conditions
        bcr_constant = {
            "KCC": "Kappa C",
            "LCC": "Lambda C",
            "HCC": "Heavy C",
            "HC1": "Heavy C domain 1",
            "HC2": "Heavy C domain 2",
            "HC3": "Heavy C domain 3",
        }
        tcr_constant = {
            "TRAC": "Alpha C",
            "TRBC": "Beta C",
            "TRDC": "Delta C",
            "TRGC": "Gamma C",
        }
        tcr_var = {"A": "Alpha V", "B": "Beta V", "G": "Gamma V", "D": "Delta V"}
        bcr_var = {"H": "Heavy V", "K": "Kappa V", "L": "Lambda V"}

        var = {**bcr_var, **tcr_var}
        constant = {**bcr_constant, **tcr_constant}

        # We have no hits
        if ndomains == 0:
            return None, None, None

        v_domains = top_domains_set.intersection(var)
        # set of V domain types
        c_domains = top_domains_set.intersection(constant)
        # set of C domain types

        # if no V domain types recognised:
        if len(v_domains) < 1:
            return None, None, None

        # if only a single V domain type recognised,
        if len(v_domains) == 1:
            v_domain = next(iter(v_domains))
            v_type = var[v_domain]
            receptor_type = 'BCR' if v_domain in bcr_var else 'TCR'
            # if only one V domain hit
            if len(chain_type_sp[v_domain]) == 1:
                species = chain_type_sp[v_domain][0]
                if len(c_domains) != 0:
                    c_types = '/'.join(sorted([constant[c_domain] for c_domain in c_domains]))
                    return (receptor_type, v_type + ', ' + c_types, species)
                else:
                    return (receptor_type, v_type, species)
            # > one domain of the same type: a construct, not currently explicitly typed in IEDB schema.
            else:
                species = '/'.join(set(chain_type_sp[v_domain]))
                return (receptor_type, 'construct', species )

        # Possible scFv / tandem TCR.
        if len(v_domains) == 2:
            # heavy + light --> scFv.
            if v_domains == {"H", "L"}:
                # if one copy of each, traditional scFv
                if ( len(chain_type_sp['H']) == 1) & ( len(chain_type_sp['L']) == 1):
                    species = chain_type_sp['H'][0] + '/' + chain_type_sp['L'][0]
                    if len(c_domains) != 0:
                        c_types = '/'.join(sorted([constant[c_domain] for c_domain in c_domains]))
                        return ('BCR', 'scFv' + ', ' + c_types, species)
                    else:
                        return ('BCR', 'scFv', species)
                # else, a construct.
                else:
                    species = '/'.join(['/'.join(chain_type_sp[p]) for p in ['H', 'L']])
                    return ('BCR', 'construct', species)

            if (v_domains == {"G", "D"}) | (v_domains == {"A", "B"}):
                if v_domains == {"G", "D"}:
                    heavy = "G"
                    light = "D"
                else:
                    heavy = "A"
                    light = "B"
                if ( len(chain_type_sp[heavy]) == 1) & ( len(chain_type_sp[light]) == 1):
                    species = chain_type_sp[heavy][0] + '/' + chain_type_sp[light][0]
                    return ('TCR', 'TscFv', species)
                else:
                    species = '/'.join(['/'.join(chain_type_sp[p]) for p in [heavy, light]])
                    return ('TCR', 'construct', species)

            else:
                species = '/'.join([p.split('_')[0] for p in top_domains])
                receptor_types = '/'.join(set(['BCR' if v_domain in bcr_var else 'TCR' for v_domain in v_domains]))
                return (receptor_types, 'construct', species)

        if len(v_domains) > 2:
            species = '/'.join([p.split('_')[0] for p in top_domains])
            receptor_types = '/'.join(set(['BCR' if v_domain in bcr_var else 'TCR' for v_domain in v_domains]))
            return (receptor_types, 'construct', species)


    def assign_Gdomain(self, seq, seq_id=None):
        """Returns G domain of a MHC sequence.

        Args:
            seq: a biopython sequence record

            seq_id: the corresponding sequence ID

        Returns:
            The G domain results, or None if there are no hits
        """
        gd = mhc_G_domain(chain1_seq=seq, ch1_id=re.sub(r"[^\w|\.]", "", seq_id))
        res = gd.get_g_domain()
        if res:
            return res
        else:
            return None, None

    def get_MRO_Gdomains(self, mro_TSVfile):
        """
        Returns G domains of the MRO chain sequences.

        Args:
            mro_TSVfile: Contains MHC alleles and their Gdom seqs
                (Can be located in the MRO github repo if lost)
        """
        mro = pd.read_csv(mro_TSVfile, sep="\t", skiprows=[1])
        if (
                os.path.exists(self.mro_gdomain_file)
                and os.path.getsize(self.mro_gdomain_file) > 0
        ):
            mro_out = pd.read_csv(self.mro_gdomain_file)
            cnt = mro_out.Label.index[-1] + 1
        else:
            mro_out = pd.DataFrame(
                columns=["Label", "Sequence", "calc_mhc_class", "ch_g_dom"]
            )
            cnt = 0
        for i in list(range(mro.shape[0])):
            if pd.isnull(mro.loc[i, "Sequence"]):
                continue
            if mro.loc[i, "Label"] in list(mro_out["Label"]) and mro.loc[
                i, "Sequence"
            ] in list(mro_out["Sequence"]):
                continue
            mro_out.loc[cnt, "Label"] = mro.loc[i, "Label"]
            mro_out.loc[cnt, "Sequence"] = mro.loc[i, "Sequence"]
            (
                mro_out.loc[cnt, "calc_mhc_class"],
                mro_out.loc[cnt, "ch_g_dom"],
            ) = self.assign_Gdomain(mro.loc[i, "Sequence"], mro.loc[i, "Accession"])
            cnt += 1
        mro_out.to_csv(self.mro_gdomain_file, index=False)
        return mro_out

    def get_MRO_allele(self, mro_df, seq, seq_id=None):
        """Gets the allele corresponding to an MRO sequence

        Args:
            mro_df: pandas dataframe with MRO data

            seq: sequence to query against MRO

            seq_id: corresponding sequence ID

        Returns:
            The sequences allele if found in the MRO dataframe
        """
        mhc_class, pdb_g_dom = self.assign_Gdomain(seq, seq_id)
        if pdb_g_dom:
            mro_allele = (
                str(
                    list(
                        mro_df[
                            mro_df.fillna({"ch_g_dom": ""}).apply(
                                lambda r: r["ch_g_dom"] != ""
                                          and (
                                                  pdb_g_dom in r["ch_g_dom"]
                                                  and seq in r["Sequence"]
                                                  or r["ch_g_dom"] in pdb_g_dom
                                                  and r["Sequence"] in seq
                                          ),
                                axis=1,
                            )
                        ]["Label"]
                    )
                )
                .strip("[]")
                .replace(",", "#")
            )
            return mro_allele
        else:
            print("[INFO] Unable to assign G domain to the {} chain sequence".format(seq_id))
            return

    def is_MHC(self, sequence_records, hmm):
        """Checks if sequence is MHC using HMMER

        Args:
            sequence: sequence to query

            hmm: HMM used to query sequence using HMMER

        Returns:
            score: The bit score for sequence query against provided HMM
        """
        # create a temporary file
        with tempfile.NamedTemporaryFile(mode="w") as fp:
            SeqIO.write(sequence_records, fp, format='fasta')
            fp.flush()
            # Find MHC sequences
            hmmer = self.hmmer_path + "hmmscan"
            args = [hmmer, hmm, fp.name]
            cmd = " ".join(args)
            output = self.run_cmd(cmd)
            max_scores = {}
            for hmmer_query in SearchIO.parse(StringIO(output), 'hmmer3-text'):
                headers, values = self.parse_hmmer_query(hmmer_query)
                max_score = 0
                for hit in values:
                    if hit['bitscore'] >= max_score:
                        max_score = hit['bitscore']
                max_scores[hmmer_query.id] = max_score
        # close the file. When the file is closed it will be removed.
        fp.close()
        return max_scores

    def is_b2m(self, sequence_records):
        """Checks if sequence is b2m
        Uses BLAST as method rather than HMMER

        Args:
            sequence: Sequence to query

        Returns:
            True if sequence is b2m, False if sequence is not
        """
        print(f'[INFO] B2M check (BLASTP).')
        blast = self.blast_path
        hit_coverage = "75"
        hit_perc_id = 0.50
        b2m_set = set()
        with tempfile.NamedTemporaryFile(mode="w") as temp_in:
            with tempfile.NamedTemporaryFile(mode="r") as outfile:
                SeqIO.write(sequence_records, temp_in.name, "fasta")
                blast_cmd = [
                    blast,
                    "-db",
                    self.b2m_db,
                    "-query",
                    temp_in.name,
                    "-evalue",
                    "10e-4",
                    "-qcov_hsp_perc",
                    hit_coverage,
                    "-outfmt",
                    "5",
                    ">",
                    outfile.name,
                ]
                self.run_cmd((" ".join(blast_cmd)))
                rec = NCBIXML.parse(outfile)
                for res in rec:
                    for alignment in res.alignments:
                        for hsp in alignment.hsps:
        
                            if (
                                    float(hsp.identities) / float(hsp.align_length)
                                    > hit_perc_id
                            ):
                                b2m_set.add(res.query)
        return b2m_set


    def is_b2m_speedy(self, sequence_records):
        """Checks if sequence is b2m
        Uses MMseqs rather than BLAST

        Args:
            sequence: Sequence to query

        Returns:
            True if sequence is b2m, False if sequence is not
        """
        print(f'[INFO] B2M check (speedy).')
        mmseqs = self.mmseqs_path
        hit_coverage = "0.75"
        hit_perc_id = 0.50
        mmseqs_columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        with tempfile.NamedTemporaryFile(mode="w") as temp_in:
            with tempfile.NamedTemporaryFile(mode="r") as outfile:
                with tempfile.TemporaryDirectory() as temp_dir:
                    SeqIO.write(sequence_records, temp_in.name, "fasta")
                    mmseqs_cmd = [mmseqs,
                                'easy-search',
                                temp_in.name,
                                self.b2m_db_speedy,
                                outfile.name,
                                temp_dir,
                                '--min-seq-id',
                                str(hit_perc_id),
                                '-c',
                                hit_coverage,
                                '-e',
                                '10e-4'
                                ]
                    self.run_cmd((" ".join(mmseqs_cmd)))
                    res = pd.read_csv(outfile.name, sep='\t', names=mmseqs_columns)
                    return set(res[res['pident']>= hit_perc_id]['qseqid'].unique())

    def is_ignar_speedy(self, sequence):
        """Checks if sequence is shark antibody (IgNAR)

        Uses BLAST as method rather than HMMER due to lack of sequences

        Args:
            sequence: The sequence to query

        Returns:
            True if sequence is IgNAR, False if sequence is not IgNAR
        """
        print(f'[INFO] IgNAR check (speedy).')
        mmseqs = self.mmseqs_path
        hit_coverage = "0.75"
        hit_perc_id = 0.5
        mmseqs_columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        with tempfile.NamedTemporaryFile(mode="w") as temp_in:
            with tempfile.NamedTemporaryFile(mode="r") as outfile:
                with tempfile.TemporaryDirectory() as temp_dir:
                    SeqIO.write(sequence, temp_in.name, "fasta")
                    mmseqs_cmd = [mmseqs,
                                'easy-search',
                                temp_in.name,
                                self.ignar_db_speedy,
                                outfile.name,
                                temp_dir,
                                '--min-seq-id',
                                str(hit_perc_id),
                                '-c',
                                hit_coverage,
                                '-e',
                                '10e-4'
                                ]
                    self.run_cmd((" ".join(mmseqs_cmd)))
                    res = pd.read_csv(outfile.name, sep='\t', names=mmseqs_columns)
                    return set(res[res['pident']>= hit_perc_id]['qseqid'].unique())
            
    def is_ignar(self, sequence_records):
        """Checks if sequence is shark antibody (IgNAR)

        Uses BLAST as method rather than HMMER due to lack of sequences

        Args:
            sequence: The sequence to query

        Returns:
            True if sequence is IgNAR, False if sequence is not IgNAR
        """
        print(f'[INFO] IgNAR check (BLASTP).')
        blast = self.blast_path
        hit_coverage = "0.75"
        hit_perc_id = 0.5
        ignar_set = set()
        with tempfile.NamedTemporaryFile(mode="w") as temp_in:
            with tempfile.NamedTemporaryFile(mode="r") as outfile:
                SeqIO.write(sequence_records, temp_in.name, "fasta")
                blast_cmd = [
                    blast,
                    "-db",
                    self.ignar_db,
                    "-query",
                    temp_in.name,
                    "-evalue",
                    "10e-4",
                    "-qcov_hsp_perc",
                    hit_coverage,
                    "-outfmt",
                    "5",
                    ">",
                    outfile.name,
                ]
                self.run_cmd((" ".join(blast_cmd)))
                res = NCBIXML.parse(outfile)
                for record in res:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            if (
                                    float(hsp.identities) / float(hsp.align_length)
                                    > hit_perc_id
                            ):
                                ignar_set.add(record.query)
        return ignar_set

    def assign_class(self, seq_records, bit_score_threshold=100):
        """Classifies sequence as BCR, TCR, or MHC

        Args:
            seq_recored: A biopython sequence record object

        Returns:
            The receptor and chain type of input sequence, if available
        """
        seq_records = [p for p in seq_records if self.check_seq(p)]
        if len(seq_records) == 0:
            return {}
    
        print(f"[INFO] Running receptor scan.")
        with tempfile.NamedTemporaryFile(mode="w") as hmm_out:
            self.run_hmmscan(seq_records, hmm_out)
            hmmer_results = SearchIO.parse(hmm_out.name, "hmmer3-text")
            output_table = {}
            for hmmer_query in hmmer_results:
                hit_table, top_descriptions = self.parse_hmmer_query(hmmer_query, bit_score_threshold=bit_score_threshold)
                try:
                    score = int(hit_table[1][3] - 100)
                except:
                    score = int(0 - 100)
                receptor, chain_type, species = self.get_chain_type(top_descriptions)
                output_table[hmmer_query.id] = (receptor, chain_type, score, species)

        non_receptor = [p for p in seq_records if output_table[p.description][0] is None or output_table[p.description][1] is None]


        if len(non_receptor) != 0:
            if self.speedy:
                b2m_results = self.is_b2m_speedy(non_receptor)
                output_table.update({p.description: ("B2M", "-", 0, '') for p in non_receptor if p.description in b2m_results})
                ignar_results = self.is_ignar_speedy(non_receptor)
                output_table.update({p.description: ("BCR", "IgNAR", 0, '') for p in non_receptor if p.description in ignar_results})
            else:
                b2m_results = self.is_b2m(non_receptor)
                output_table.update({p.description: ("B2M", "-", 0, '') for p in non_receptor if p.description in b2m_results})
                ignar_results = self.is_ignar(non_receptor)
                output_table.update({p.description: ("BCR", "IgNAR", 0, '') for p in non_receptor if p.description in ignar_results})


        print(f"[INFO] Running MHC scan.")
        possible_mhc = [p for p in seq_records if output_table[p.description][0] is None or output_table[p.description][1] is None]
        if len(possible_mhc) == 0:
            return output_table
        
        mhc_I_score = None
        mhc_I_score = self.is_MHC(possible_mhc, self.mhc_I_hmm)
    
        output_table.update({p.description: ('MHC-I', 'alpha', mhc_I_score[p.description] - self.hmm_score_threshold, '') for p in possible_mhc if mhc_I_score[p.description] >= self.hmm_score_threshold})

        possible_mhc = [p for p in seq_records if output_table[p.description][0] is None or output_table[p.description][1] is None]
        if len(possible_mhc) == 0:
            return output_table

        mhc_II_alpha_score = self.is_MHC(possible_mhc, self.mhc_II_alpha_hmm)
        output_table.update({p.description: ('MHC-II', 'alpha', mhc_II_alpha_score[p.description] - self.hmm_score_threshold, '') for p in possible_mhc if mhc_II_alpha_score[p.description] >= self.hmm_score_threshold})

        possible_mhc = [p for p in seq_records if output_table[p.description][0] is None or output_table[p.description][1] is None]

        if len(possible_mhc) == 0:
            return output_table
        
        mhc_II_beta_score = self.is_MHC(possible_mhc, self.mhc_II_beta_hmm)
        output_table.update({p.description: ('MHC-II', 'beta', mhc_II_beta_score[p.description] - self.hmm_score_threshold, '') for p in possible_mhc if mhc_II_beta_score[p.description] >= self.hmm_score_threshold})

        possible_mhc = [p for p in seq_records if output_table[p.description][0] is None or output_table[p.description][1] is None]
        
        if len(possible_mhc) == 0:
            return output_table
        
        for p in possible_mhc:
            if mhc_II_alpha_score[p.description] == 0 and mhc_II_beta_score[p.description] == 0:
                output_table[p.description] = (None, None, score, '')
            elif mhc_II_alpha_score[p.description] >= mhc_II_beta_score[p.description]:
                output_table[p.description] =  (
                            None,
                            None,
                            int(mhc_II_alpha_score[p.description] - self.hmm_score_threshold),
                            ''
                        )
            else:
                output_table[p.description] = (
                            None,
                            None,
                            int(mhc_II_beta_score[p.description] - self.hmm_score_threshold),
                            ''
                        )
        return output_table

    def classify_seqfile(self, seq_file):
        """Classifies the sequences in a FASTA format file

        This method will write results of classificaiton to specified
        outfile in a tab separated format.

        Args:
            seq_file: the name of a FASTA file of sequences
        """
        print(f'[INFO] Reading {seq_file}.')
        seq_records = list(SeqIO.parse(seq_file, "fasta"))
        print(f'[INFO] Removing redundancy for speed.')
        inputs_unique = {}
        for p in seq_records:
            if p.seq not in inputs_unique:
                inputs_unique[p.seq] = []
            inputs_unique[p.seq].append(p)
        inputs_unique_records = [SeqRecord(description=f'seq{i}', seq=Seq(p), id=f'seq{i}') for i,p in enumerate(inputs_unique)]
        inputs_unique_records_mapping = {p.description: p for p in inputs_unique_records}
        rev_mapping = {p.description: inputs_unique[p.seq] for p in inputs_unique_records}
        output = pd.DataFrame(self.assign_class(inputs_unique_records)).T
        if output.shape[0] == 0:
            output = pd.DataFrame({p.description: {'class': None, 'chain_type': None, 'score': -100, 'species': np.nan, 'species_score':np.nan, 'calc_mhc_allele':np.nan}
                                 for p in seq_records}).T.reset_index().rename(columns={'index':'id'})
            output.to_csv(self.outfile, sep="\t", index=False)
            return output
        output.columns = ['class', 'chain_type', 'score', 'species']
        output = output.reset_index().rename(columns={'index':'id'})
        for k ,p in output.iterrows():
            if "MHC" in str(p['class']):
                seq_record = inputs_unique_records_mapping[p.id]
                calc_mhc_allele = self.get_MRO_allele(
                        self.mro_df, str(seq_record.seq), str(seq_record.description)
                    )
                output.at[k, 'calc_mhc_allele'] = calc_mhc_allele
        if self.recalc_species:
            print(f'[INFO] Starting species reassignment on BCRs and TCRs.')
            ig_tr = output[output['class'].isin(set(['BCR', 'TCR']))]
            if ig_tr.shape[0] == 0:
                print(f'[INFO] No BCRs or TCRs detected.')
                pass
            else:
                ig_tr_sp = []
                for locus, df in ig_tr.groupby('class'):
                    locus_name = 'IG' if locus == 'BCR' else 'TR'
                    ids = set(df['id'].unique())
                    rename_records = {p: f'seq{i}' for i,p in enumerate(ids)}
                    rename_records_rev = {rename_records[p]:p for p in rename_records}
                    with tempfile.NamedTemporaryFile(mode="w") as temp_out:
                        records = [inputs_unique_records_mapping[p] for p in ids]
                        new_records = []
                        for record in records:
                            record.id = rename_records[record.description]
                            new_records.append(record)
                        SeqIO.write(new_records, temp_out.name, "fasta")
                        species_reassignment = self.get_species_seqfile(seq_file = temp_out.name, locus = locus_name)
                        species_reassignment['qseqid'] = species_reassignment['qseqid'].map(str).map(rename_records_rev)
                        ig_tr_sp.append(species_reassignment)
                ig_tr_sp = pd.concat(ig_tr_sp)
                print(f'[INFO] Species reassignment complete.')
                output = pd.merge(left = output.drop(['species'], axis = 1), right = ig_tr_sp[['qseqid', 'species', 'bitscore']].rename(columns = {'bitscore':'species_score'}),
                               left_on = 'id', right_on = 'qseqid', how = 'left').drop(['qseqid'], axis = 1)
        if 'species_score' not in output.columns:
            output['species_score'] = np.nan
        if 'calc_mhc_allele' not in output.columns:
            output['calc_mhc_allele'] = np.nan
        order = ['id','class','chain_type','score','calc_mhc_allele','species','species_score']
        output_remapped = []
        for entry in rev_mapping:
            subset = output[output['id'] == entry][order]
            if subset.shape[0] == 0:
                for map in rev_mapping[entry]:
                    add = [map.description, None, None, None, None, None, None]
                    output_remapped.append(add)
            else:
                subset = subset.iloc[0].values
                for map in rev_mapping[entry]:
                    if subset.shape[0] != 0:
                        add =[ map.description] + list(subset[1:])
                    else:
                        add = [map.description, None, None, None, None, None, None, None]
                    output_remapped.append(add)
        output_remapped = pd.DataFrame(output_remapped, columns=['id','class','chain_type','score','calc_mhc_allele','species','species_score'])
        output_remapped.to_csv(self.outfile, sep="\t", index=False)
        return output_remapped