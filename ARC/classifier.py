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
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from Bio import SeqIO


class SeqClassifier:
    """Classifies input sequence/s into BCR, TCR or MHC chains.

    Attributes:
        seqfile (str): Input sequence file in FASTA format

        outfile (str): Name of output file

        threads (int): number of cpu threads to run on

        hmmer_path (str): optional argument containing the path to your HMMER installation

        blast_path (str): optional argument containing the path to your BLAST installation"""

    def __init__(
            self, seqfile=None, outfile=None, threads=1, hmmer_path=None, blast_path=None, hmm_path=None):
        # Relative paths and IO handling
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
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
            self.blast_path = ""
        else:
            self.blast_path = blast_path

    def check_seq(self, seq_record):
        """Checks validity of an amino acid sequence

        Args:
            seq_record: A biopython sequence record object
        """
        pattern = re.compile(
            r"[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]", re.IGNORECASE
        )
        if len(str(seq_record.seq)) > 0 and not pattern.findall(str(seq_record.seq)):
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

    def get_species(self, seq_record,
                    blast_path = 'test/data/imgt/blast_fasta/', locus = "IG"):
        with tempfile.NamedTemporaryFile(mode="w") as temp_out:
            if not seq_record.seq:
                return False
            SeqIO.write(seq_record, temp_out.name, "fasta")
            path = os.path.join(blast_path, f'{locus}V.fasta')
            subprocess.call(f'blastp -query {temp_out.name} -db {path} -evalue 1e-6 -num_threads 4 -out {temp_out.name}blast.txt -outfmt 6', shell =
                            True)
            if os.path.getsize(f'{temp_out.name}blast.txt') == 0:
                top_species = {'species':'none', 'score': 0}
                return top_species
            output = pd.read_csv(f'{temp_out.name}blast.txt', sep = '\t', header = None)
            output.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send','evalue','bitscore']
            output['species'] = output['sseqid'].map(lambda x:x.split('|')[1])
            top_species = output.groupby('sseqid').apply(lambda x:x.loc[x['bitscore'].idxmax()].drop('qseqid')).sort_values('bitscore', ascending = False).iloc[0]
            return top_species

    def get_species_seqfile(self, seq_file,
                    blast_path = 'test/data/imgt/blast_fasta/', locus = "IG"):
        db_path = os.path.join(blast_path, locus + "V.fasta")
        with tempfile.NamedTemporaryFile(mode="w") as temp_out:
            subprocess.call(f'blastp -query {seq_file} -db {db_path} -evalue 1e-6 -num_threads 4 -out {temp_out.name} -outfmt 6', shell =
                            True)
            if os.path.getsize(temp_out.name) == 0:
                top_species = {'species':'none', 'score': 0}
                return pd.DataFrame(top_species)
            output = pd.read_csv(temp_out.name, sep = '\t', header = None)
            output.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send','evalue','bitscore']
            output['species'] = output['sseqid'].map(lambda x:x.split('|')[1])
            top_species = output.groupby('qseqid').apply(lambda x: x.loc[x['bitscore'].idxmax()].drop(['qseqid'])).reset_index()
            return top_species

    def run_hmmscan(self, seq_record, hmm_out):
        """Runs hmmscan from the HMMER3 software suite

        Args:
            seq_record: A biopython sequence record object

            hmm_out: tempfile object for hmm output
        """
        with tempfile.NamedTemporaryFile(mode="w") as temp_out:
            if not seq_record.seq:
                return False
            SeqIO.write(seq_record, temp_out.name, "fasta")
            hmmer = self.hmmer_path + "hmmscan"
            args = [
                hmmer,
                "-o",
                hmm_out.name,
                self.hmm_path,
                temp_out.name,
            ]
            cmd = (" ").join(args)
            self.run_cmd(cmd, str(seq_record.seq))

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
                    print(hsp.hit_id)
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
            print("Unable to assign G domain to the {} chain sequence".format(seq_id))
            return

    def is_MHC(self, sequence, hmm):
        """Checks if sequence is MHC using HMMER

        Args:
            sequence: sequence to query

            hmm: HMM used to query sequence using HMMER

        Returns:
            score: The bit score for sequence query against provided HMM
        """
        score = 0
        # create a temporary file
        fp = tempfile.NamedTemporaryFile(mode="w")
        fp.write(">seq\n")
        fp.write("{}".format(sequence))
        fp.flush()

        # Find MHC sequences
        hmmer = self.hmmer_path + "hmmscan"
        args = [hmmer, hmm, fp.name]
        cmd = " ".join(args)
        output = self.run_cmd(cmd)
        aln = [line.split() for line in output.splitlines()]
        # Search for score to see if there is a match
        for i, line in enumerate(aln):
            if line[0:3] == ["E-value", "score", "bias"] and aln[i + 2]:
                try:
                    # E_value = float(aln[i+2][0])
                    score = float(aln[i + 2][1])
                    break
                except ValueError:
                    # E_value = float(aln[i+3][0])
                    score = float(aln[i + 3][1])
                    break

        # close the file. When the file is closed it will be removed.
        fp.close()
        return score

    def is_b2m(self, sequence):
        """Checks if sequence is b2m
        Uses BLAST as method rather than HMMER

        Args:
            sequence: Sequence to query

        Returns:
            True if sequence is b2m, False if sequence is not
        """
        blast = self.blast_path + "blastp"
        hit_coverage = "75"
        hit_perc_id = 0.50
        with tempfile.NamedTemporaryFile(mode="w") as temp_in:
            with tempfile.NamedTemporaryFile(mode="r") as outfile:
                SeqIO.write(sequence, temp_in.name, "fasta")
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
                res = NCBIXML.read(outfile).alignments
                for alignment in res:
                    for hsp in alignment.hsps:
                        if (
                                float(hsp.identities) / float(hsp.align_length)
                                > hit_perc_id
                        ):
                            return True
                return False

    def is_ignar(self, sequence):
        """Checks if sequence is shark antibody (IgNAR)

        Uses BLAST as method rather than HMMER due to lack of sequences

        Args:
            sequence: The sequence to query

        Returns:
            True if sequence is IgNAR, False if sequence is not IgNAR
        """
        blast = self.blast_path + "blastp"
        hit_coverage = "75"
        hit_perc_id = 0.50
        with tempfile.NamedTemporaryFile(mode="w") as temp_in:
            with tempfile.NamedTemporaryFile(mode="r") as outfile:
                SeqIO.write(sequence, temp_in.name, "fasta")
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
                res = NCBIXML.read(outfile).alignments
                for alignment in res:
                    for hsp in alignment.hsps:
                        if (
                                float(hsp.identities) / float(hsp.align_length)
                                > hit_perc_id
                        ):
                            return True
                return False

    def assign_class(self, seq_record, bit_score_threshold=100):
        """Classifies sequence as BCR, TCR, or MHC

        Args:
            seq_recored: A biopython sequence record object

        Returns:
            The receptor and chain type of input sequence, if available
        """
        with tempfile.NamedTemporaryFile(mode="w") as hmm_out:
            receptor, chain_type = None, None
            self.run_hmmscan(seq_record, hmm_out)
            hmmer_query = SearchIO.read(hmm_out.name, "hmmer3-text")
            hit_table, top_descriptions = self.parse_hmmer_query(hmmer_query, bit_score_threshold=bit_score_threshold)
            try:
                score = int(hit_table[1][3] - 100)
            except:
                score = int(0 - 100)
            receptor, chain_type, species = self.get_chain_type(top_descriptions)

            # We have no hits so now we check for MHC and IgNAR
            # This avoids excessive computations
            if not receptor or not chain_type:
                if self.is_b2m(seq_record):
                    return ("B2M", "-", 0, '')
                if self.is_ignar(seq_record):
                    return ("BCR", "IgNAR", 0, '')
                mhc_I_score = None
                mhc_I_score = self.is_MHC(str(seq_record.seq), self.mhc_I_hmm)
                if mhc_I_score >= self.hmm_score_threshold:
                    return (
                        "MHC-I",
                        "alpha",
                        int(mhc_I_score - self.hmm_score_threshold),
                        ''
                    )
                else:
                    mhc_II_alpha_score = None
                    mhc_II_alpha_score = self.is_MHC(
                        str(seq_record.seq), self.mhc_II_alpha_hmm
                    )
                    if (
                            mhc_II_alpha_score
                            and mhc_II_alpha_score >= self.hmm_score_threshold
                    ):
                        return (
                            "MHC-II",
                            "alpha",
                            mhc_II_alpha_score - self.hmm_score_threshold,
                            ''
                        )
                    else:
                        mhc_II_beta_score = None
                        mhc_II_beta_score = self.is_MHC(
                            str(seq_record.seq), self.mhc_II_beta_hmm
                        )
                        if (
                                mhc_II_beta_score
                                and mhc_II_beta_score >= self.hmm_score_threshold
                        ):
                            return (
                                "MHC-II",
                                "beta",
                                int(mhc_II_beta_score - self.hmm_score_threshold),
                                ''
                            )
                        else:
                            if mhc_II_alpha_score == 0 and mhc_II_beta_score == 0:
                                return (None, None, score, '')
                            if mhc_II_alpha_score >= mhc_II_beta_score:
                                return (
                                    None,
                                    None,
                                    int(mhc_II_alpha_score - self.hmm_score_threshold),
                                    ''
                                )
                            else:
                                return (
                                    None,
                                    None,
                                    int(mhc_II_beta_score - self.hmm_score_threshold),
                                    ''
                                )
            else:
                if score < 0:
                    score = 0

                return (receptor, chain_type, score, species)

    def gen_classify(self, seq, seq_id):
        """Returns BCR, TCR, or MHC class and chain type for input sequence without needing
        Biopython FASTA sequence object

        Args:
            seq: string containing protein sequence of interest

            seq_id: id of sequence (> in FASTA)

        Returns:
            Receptor, chain type, and calculated MHC allele if applicable
        """
        seq_record = SeqRecord(Seq(seq), id=seq_id)
        g_domain = ""
        calc_mhc_allele = ""
        receptor, chain_type, score, species = self.assign_class(seq_record)
        if receptor == "MHC-I" or receptor == "MHC-II":
            g_domain = self.assign_Gdomain(str(seq_record.seq), seq_record.id)
            calc_mhc_allele = self.get_MRO_allele(
                self.mro_df, str(seq_record.seq), str(seq_record.description)
            )
        return receptor, chain_type, calc_mhc_allele, species

    def classify(self, seq_record, bit_score_threshold=100, recalc_species=True):
        """Returns BCR, TCR or MHC class and chain type for an input sequence.

        If sequence is MHC, finds its g-domain and returns its corresponding
        allele

        Args:
            seq_record: a biopython sequence record object

        Returns:
            The receptor, chain type, and calculated MHC allele, if applicable
        """
        g_domain = ""
        calc_mhc_allele = ""
        receptor, chain_type, score, species = self.assign_class(seq_record, bit_score_threshold=bit_score_threshold)
        species_score = score
        if receptor == "MHC-I" or receptor == "MHC-II":
            g_domain = self.assign_Gdomain(str(seq_record.seq), seq_record.id)
            calc_mhc_allele = self.get_MRO_allele(
                self.mro_df, str(seq_record.seq), str(seq_record.description)
            )
        else:
            if recalc_species:
                locus = 'IG' if receptor == 'BCR' else 'TR'
                new_species = self.get_species(seq_record, locus = locus)
                species = new_species['species']
                species_score = new_species['bitscore']
        return receptor, chain_type, calc_mhc_allele, score, species, species_score

    def classify_multiproc(self, seq_list, recalc_species = False):
        out = pd.DataFrame(columns=["id", "class", "chain_type", "calc_mhc_allele", "species", "species_score"])
        cnt = 0
        for seq in seq_list:
            if seq.seq == "":
                print(f"{seq.description} has empty sequence. Skipping sequence.")
                continue
            if self.check_seq(seq):
                receptor, chain_type, calc_mhc_allele, score, species, species_score = self.classify(seq, recalc_species \
                = recalc_species)
            else:
                print(
                    f"{seq.description} contains invalid amino acid sequence. Skipping sequence."
                )
                continue

            out.loc[cnt, "id"] = seq.description
            out.loc[cnt, "class"] = receptor
            out.loc[cnt, "chain_type"] = chain_type
            out.loc[cnt, "calc_mhc_allele"] = calc_mhc_allele
            out.loc[cnt, "score"] = score
            out.loc[cnt, "species"] = species
            out.loc[cnt, "species_score"] = species_score
            cnt += 1

        return out

    def classify_seqfile(self, seq_file, recalc_species = True):
        """Classifies the sequences in a FASTA format file

        This method will write results of classificaiton to specified
        outfile in a tab separated format.

        Args:
            seq_file: the name of a FASTA file of sequences
        """
        seq_records = list(SeqIO.parse(seq_file, "fasta"))
        if len(seq_records) == 1:
            out = self.classify_multiproc(seq_records)
        else:
            num_records = len(seq_records)
            records_per_thread = max(num_records // self.num_threads, 1)  # ensure at least 1 record per thread

            # Create chunks of sequences to be processed by each thread
            chunks = [seq_records[i:i + records_per_thread] for i in range(0, num_records, records_per_thread)]
            with mp.Pool(processes=self.num_threads) as pool:
                # Use the map function to distribute the workload
                results = pool.map(self.classify_multiproc, chunks)
                out = pd.concat(results)
        if recalc_species:
            ig_tr = out[out['class'].isin(set(['BCR', 'TCR']))]
            if ig_tr.shape[0] == 0:
                pass
            else:
                ig_tr_sp = []
                for locus, df in ig_tr.groupby('class'):
                    locus_name = 'IG' if locus == 'BCR' else 'TR'
                    ids = set(df['id'].unique())
                    with tempfile.NamedTemporaryFile(mode="w") as temp_out:
                        records = [p for p in seq_records if p.description in ids]
                        SeqIO.write(records, temp_out.name, "fasta")
                        species_reassignment = self.get_species_seqfile(seq_file = temp_out.name, locus = locus_name)
                        ig_tr_sp.append(species_reassignment)
                ig_tr_sp = pd.concat(ig_tr_sp)
                ig_tr_sp['qseqid'] = ig_tr_sp['qseqid'].map(str)
                out = pd.merge(left = out.drop(['species', 'species_score'], axis = 1), right = ig_tr_sp[['qseqid', 'species', 'bitscore']].rename(columns = {'bitscore':'species_score'}),
                               left_on = 'id', right_on = 'qseqid', how = 'left').drop(['qseqid'], axis = 1)
        out.to_csv(self.outfile, sep="\t", index=False)