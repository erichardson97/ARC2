import requests
import os
import ssl
import certifi
import urllib
from bs4 import BeautifulSoup
import sys
import time
import os, sys
from subprocess import Popen, PIPE
import subprocess
from Bio import SeqIO
import yaml
from datetime import date
from importlib.resources import files


species_translations = {"Homo_sapiens": "human",
                    "Mus": "mouse",
                    "Rattus_norvegicus": "rat",
                    "Oryctolagus_cuniculus": "rabbit",
                    "Macaca_mulatta": "rhesus",
                    "Sus_scrofa": "pig",
                    "Vicugna_pacos": "alpaca",
                    "Bos_taurus": "cow"}


class DataDownloader():
    def __init__(self, overwrite = False):
        self.package_directory = files('ARC')
        self.overwrite_status = overwrite

    def validate_data(self, data_manifest):
        ## Needs to be finished.
        pass


    def download_MRO_tsv(self):
        outfile = os.path.join(self.package_directory, 'data/chain-sequence.tsv')
        if not self.overwrite_status:
            if os.path.exists(outfile):
                ## Replace with proper error handling
                return None
        mro_request = requests.get('https://raw.githubusercontent.com/IEDB/MRO/master/ontology/chain-sequence.tsv')
        if mro_request.status_code != 200:
            ## Replace with proper error handling:)
            return f'Unable to access IEDB MRO file.'
        else:
            with open(outfile, 'w') as k:
                k.writelines(mro_request.content.decode())

    def download_IG_TR_databases(self):
        db = IG_TR_Database()
        db.write_imgt_urls()
        db.download_imgt(species_list=['Homo+sapiens','Mus'], locus='TR')
        db.download_imgt(species_list=['Homo+sapiens', 'Mus', 'Macaca+mulatta', 'Oryctolagus+cuniculus', 'Vicugna+pacos', 'Sus+scrofa', 'Bos+taurus', 'Rattus+norvegicus'], locus='IG')
        db.build_BLAST_databases(species_list=['Homo+sapiens', 'Mus', 'Macaca+mulatta', 'Oryctolagus+cuniculus', 'Vicugna+pacos', 'Sus+scrofa', 'Bos+taurus', 'Rattus+norvegicus'], loci=['IG'])
        db.build_BLAST_databases(species_list=['Homo+sapiens','Mus'], loci=['TR'])
        db.run_ANARCI_processing()





class IG_TR_Database():

    def __init__(self):
        self.package_directory = files('ARC')
        if os.path.exists(os.path.join(self.package_directory, 'data/imgt')) is False:
            os.mkdir(os.path.join(self.package_directory, 'data/imgt'))
        self.source_fasta = os.path.join(self.package_directory, 'data/imgt/fasta')
        self.blast_path = os.path.join(self.package_directory, 'data/imgt/blast_fasta')
        if os.path.exists(self.source_fasta) is False:
            os.mkdir(self.source_fasta)
        if os.path.exists(self.blast_path) is False:
            os.mkdir(self.blast_path)
        self.context = ssl.create_default_context()
        self.context.load_verify_locations(certifi.where())

    def write_imgt_urls(self):
        ig_urls = {"HV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHV&species=%s",
                "HJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=%s",
                "KV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKV&species=%s",
                "KJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGKJ&species=%s",
                "LV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLV&species=%s",
                "LJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGLJ&species=%s",
                "HC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHC&species=%s",
                "KC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKC&species=%s",
                "LC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLC&species=%s",
                   }

        tr_urls = {
                "AV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRAV&species=%s",
                "AJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRAJ&species=%s",
                "BV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRBV&species=%s",
                "BJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRBJ&species=%s",
                "GV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRGV&species=%s",
                "GJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRGJ&species=%s",
                "DV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRDV&species=%s",
                "DJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRDJ&species=%s",
                "AC": "https://www.imgt.org/genedb/GENElect?query=7.3+TRAC&species=%s",
                "BC": "https://www.imgt.org/genedb/GENElect?query=7.3+TRBC&species=%s",
                "GC": "https://www.imgt.org/genedb/GENElect?query=7.3+TRGC&species=%s",
                "DC": "https://www.imgt.org/genedb/GENElect?query=7.3+TRGC&species=%s",
                }
        with open(os.path.join(self.package_directory, f'data', f'imgt_access.yaml'), 'w') as file:
            file.write(yaml.dump({'IG': ig_urls, 'TR': tr_urls, 'date': date.today()}))


    def download_imgt(self, species_list = ['Homo+sapiens','Mus+musculus'], locus = 'IG'):
        urls = yaml.safe_load(open(os.path.join(self.package_directory, f'data', 'imgt_access.yaml'), 'r'))[locus]
        for sp in species_list:
            for locus in urls:
                fasta_outfile = os.path.join(self.source_fasta, f'{sp}_{locus}.fasta').replace('+', '_')
                request_url = urls[locus] % sp
                try:
                    print(f"Working on {request_url}.")
                    with urllib.request.urlopen(request_url, context = self.context) as r:
                        if r.status != 200:
                            raise urllib.error.URLError(request_url)
                        content = r.read()
                        fasta = self.caveman_parse_imgt_html(content)
                        if fasta[0] == '<hr/><b>Number of results = 0</b><br /><pre>':
                            f"Requested species + locus n'existe pas."
                            continue
                        with open(fasta_outfile, 'w') as x:
                            x.writelines('\n'.join(fasta))
                        seqs = dict((p.description, str(p.seq)) for p in SeqIO.parse(fasta_outfile, format='fasta'))
                        with open(fasta_outfile, 'w') as x:
                            for p in seqs:
                                x.write('>'+p+'\n'+seqs[p]+'\n')
                except urllib.error.URLError as ex:
                    print("Failed to download from IMGT. %s\n" % ex)
                    print("Failed to download all files from IMGT. Exiting.")
                    sys.exit(1)
                time.sleep(1)

    def parse_imgt_html(self, html_str):
        soup = BeautifulSoup(html_str, 'html.parser')
        fasta_data = []
        for pre_tag in soup.find_all('pre'):
            fasta_data.append(pre_tag.text)

    def caveman_parse_imgt_html(self, html_str):
        html_str = html_str.decode().split('\n')
        fasta_block = [k for k, p in enumerate(html_str) if p == '</pre>']
        fasta_seqs = []
        seq = ''
        for row in html_str[fasta_block[0]+2: fasta_block[1]]:
            if row.startswith('>'):
                if seq != '':
                    fasta_seqs.append(seq)
                fasta_seqs.append(row)
                seq = ''
            else:
                if row == '':
                    continue
                if row[0] != '<':
                    seq += row
        fasta_seqs.append(seq)
        return fasta_seqs
    def run_ANARCI_processing(self):
        anarci_processing_module = ANARCI_Processing()
        anarci_processing_module.run_process()

    def remove_duplicates(self, fasta_file):
        sequence_dict = dict((p.description, str(p.seq)) for p in SeqIO.parse(fasta_file, 'fasta'))
        with open(fasta_file, 'w') as k:
            for seq in sequence_dict:
                k.write('>'+seq+'\n'+sequence_dict[seq]+'\n')

    def build_BLAST_databases(self, species_list = ['Homo+sapiens','Mus'], loci = ['IG']):
        urls = yaml.safe_load(open(os.path.join(self.package_directory, f'data', 'imgt_access.yaml'), 'r'))
        ## clear all files
        for sp in species_list:
            for locus in loci:
                for gene_name in urls[locus]:
                    originating_gene = gene_name[1]
                    species_name = species_translations[sp.replace('+','_')]
                    outpath = os.path.join(self.blast_path,
                                           f'imgt_{species_name}_{locus}_{originating_gene}_input.fasta')

                    if os.path.exists(outpath) is True:
                        os.remove(outpath)

        for sp in species_list:
            for locus in loci:
                for gene_name in urls[locus]:
                    originating_gene = gene_name[1]
                    species_name = species_translations[sp.replace('+','_')]
                    outpath = os.path.join(self.blast_path, f'imgt_{species_name}_{locus}_{originating_gene}_input.fasta')
                    originating_fasta = os.path.join(self.source_fasta, f'{sp}_{gene_name}.fasta').replace('+', '_')
                    subprocess.call(f'cat {originating_fasta} >> {outpath}', shell=True)
                input_file = os.path.join(self.blast_path, f'imgt_{species_name}_{locus}_V_input.fasta')
                output_file = os.path.join(self.blast_path, f'{species_name}_{locus}V.fasta')
                subprocess.call(f'edit_imgt_file.pl {input_file} > {output_file}', shell = True)
                self.remove_duplicates(output_file)
                subprocess.call(f'makeblastdb -parse_seqids -dbtype prot -in {output_file}', shell = True)

        for locus in loci:
            all_sp = os.path.join(self.blast_path, f'{locus}V.fasta')
            with open(all_sp, 'w') as k:
                k.writelines('')
            for species in species_list:
                species_name = species_translations[species.replace('+', '_')]
                output_file = os.path.join(self.blast_path, f'{species_name}_{locus}V.fasta')
                seqs = dict((p.id+'|'+species_name, str(p.seq)) for p in SeqIO.parse(output_file, format = 'fasta'))
                with open(all_sp, 'a') as k:
                    for x in seqs:
                        k.write('>'+x+'\n'+seqs[x]+'\n')
                subprocess.call(f'makeblastdb -parse_seqids -dbtype prot -in {all_sp}',
                    shell=True)


class ANARCI_Processing():
    """
    This code is taken directly from ANARCI: https://github.com/oxpig/ANARCI/blob/master/build_pipeline/FormatAlignments.py
    And is copied here in case the source files change.
    """

    def __init__(self):
        self.package_directory = files('ARC')
        self.fasta_path = os.path.join(self.package_directory, 'data', 'imgt', 'fasta')
        self.amino_acids = sorted(list("QWERTYIPASDFGHKLCVNM"))
        self.acid_set = set(self.amino_acids + ["."])
        self.alignment_path = os.path.join(self.package_directory, 'data', 'imgt', 'curated_alignments')
        if os.path.exists(self.alignment_path) is False:
            os.mkdir(self.alignment_path)
        self.muscle_alignments = os.path.join(self.package_directory, 'data', 'imgt', 'muscle_alignments')
        if os.path.exists(self.muscle_alignments) is False:
            os.mkdir(self.muscle_alignments)
        self.ig_species = ["Homo_sapiens",
                           "Mus",
                           "Rattus_norvegicus",
                           "Oryctolagus_cuniculus",
                           "Macaca_mulatta",
                           "Sus_scrofa",
                           "Vicugna_pacos",
                           "Bos_taurus"]
        self.tr_species = ["Homo_sapiens",
                           "Mus"]
        self.species_translations = {"Homo_sapiens": "human",
                    "Mus": "mouse",
                    "Rattus_norvegicus": "rat",
                    "Oryctolagus_cuniculus": "rabbit",
                    "Macaca_mulatta": "rhesus",
                    "Sus_scrofa": "pig",
                    "Vicugna_pacos": "alpaca",
                    "Bos_taurus": "cow"}

    def read_alignment(self, input_file, read_all=False, region_name=""):
        """
        """
        imgt_fields = ["accession_number",
                       "allele",
                       "species",
                       "functionality",
                       "region",
                       "start_and_end_positions_IMGT/LIGM-DB_accession_number",
                       "number_of_nucleotides",
                       "codon_start",
                       "number_nucleotides_added_in_5'_compared_IMGT/LIGM-DB",
                       "number_nucleotides_added_in_3'_compared_IMGT/LIGM-DB",
                       "number_nucleotides_to_correct_sequencing_errors",
                       "number_of_amino_acids",
                       "number_of_characters",
                       "partial",
                       "reverse"]

        records = {}
        try:
            handle = open(input_file, "r")
        except IOError:
            print('Warning file', input_file, 'could not be found')
            return records

        region = ""
        for record in SeqIO.parse(input_file, format = 'fasta'):

            fields = dict(list(zip(imgt_fields, record.description.split("|"))))
            sequence = str(record.seq)
            # These are the ones we care about and will be used
            try:
                if fields['accession_number'] == 'None': continue
                if fields["functionality"] == "F" and not fields["partial"].strip() and not fields["reverse"].strip():
                    if read_all:
                        pass
                    elif fields["allele"].split("*")[-1].strip() != "01":
                        continue
                    if set(list(sequence)) - self.acid_set:
                        #                print >> sys.stderr,"Unexpected character in sequence"
                        #                print >> sys.stderr,sequence
                        continue

                    if fields["region"] == region_name:
                        records[(fields["species"], fields["allele"])] = sequence
                    elif region_name.startswith("C"):
                        if len(sequence) < 100:
                            continue  # Filter out partial sequences that IMGT have not....
                    elif region:
                        assert fields["region"] == region, "The region for some the entries is different"

                    region = fields["region"]
                    records[(fields["species"], fields["allele"])] = sequence
            except KeyError:
                print("Something wrong with the file %s" % input_file)
                continue

        handle.close()
        return records

    def read_fasta(self, filename):
        """
        Read a sequence file and parse as description, string
        """
        records = [[s.description, str(s.seq).replace(" ", "")] for s in SeqIO.parse(filename, format = 'fasta')]
        return records

    def write_fasta(self, sequences):
        """
        Write a fasta file containing all sequences
        """
        filename = os.path.join(self.muscle_alignments, "all_js.fasta")
        with open(filename, "w") as outfile:
            for al in sequences:
                for s in sequences[al]:
                    print(">%s|%s|%s|%s" % tuple(list(al) + list(s)), file=outfile)
                    print(sequences[al][s], file=outfile)
        return filename

    def format_c_genes(self, calignments, gene_name=""):
        new_calignments = {}
        for entry in calignments:
            if len(calignments[entry]) == 0:
                continue
            new_calignments[entry] = {}
            for seq in calignments[entry]:

                cspecies, callele = seq
                sequence = calignments[entry][seq]

                # IMGT has a different alignment for C-genes than V and J genes.
                # This condition filters out the two sequeces that are not in the consistent format for C gene.
                # tttt = sequence[:132].ljust( 132 ).replace(" ",".")
                tttt = sequence[:149].ljust(149).replace(" ", ".")
                if tttt[104] != "C" or tttt[33] != "C":
                    print("Something wrong with ", entry, gene_name, sequence)
                    continue

                max_length = 149
                new_name = "%s_%s_%s" % (gene_name, cspecies, callele)
                new_calignments[entry][new_name] = sequence[:max_length].ljust(max_length).replace(" ", ".")

        return new_calignments

    def format_j_genes(self, jalignments):

        reference = ("WFAYWGQGTLVTVSA", 4, 19)
        #                 seq           start  end

        ffile = self.write_fasta(jalignments)
        al_filename = os.path.join(self.muscle_alignments, "all_js_aligned.fasta")
        pr = Popen(["muscle", "-in", ffile, "-gapopen", "-10", "-out", al_filename, ], stdout=PIPE, stderr=PIPE)
    #    pr = Popen(["muscle", "-align", ffile, "-gapopen", "-10", "-output", al_filename, ], stdout=PIPE, stderr=PIPE)
        o, e = pr.communicate()
        aligned = self.read_fasta(al_filename)
        new_jalignments = {}

        # Find the reference sequence and what we need to do to map
        for name, sequence in aligned:
            print(name)
            if name == "Mus|H|Mus musculus|IGHJ3*01":
                ref_aligned = sequence
                break
        print(ref_aligned, reference[0])
        start = ref_aligned.index(reference[0])
        if start > reference[1]:
            START = start + 1 - reference[1]
        else:
            START = 0
        END = start + 15

        for name, sequence in aligned:
            species, chain_type, id1, id2 = name.strip(">").split("|")

            if (species, chain_type) not in new_jalignments:
                new_jalignments[(species, chain_type)] = {}
            # We take the last 13 of the new alignment and pad into 20 long string
            new_jalignments[(species, chain_type)][(id1, id2)] = sequence[START: END][-14:].rjust(20).replace(" ", ".")
        return new_jalignments

    def format_v_genes(self, valignments):
        """
        Take upto and including imgt position 108 in the alignment. Pad with gaps on the right side
        """

        new_valignments = {}
        for entry in valignments:
            species, chain_type = entry
            new_valignments[entry] = {}
            for seq in valignments[entry]:
                sequence = valignments[entry][seq]
                if chain_type == "L" and self.species_translations[species] == "rhesus":
                    sequence = self.rhesus_lambda(sequence)
                elif chain_type == "A" and self.species_translations[species] == "mouse":
                    sequence = self.mouse_alpha(sequence)
                elif chain_type == "D" and self.species_translations[species] == "mouse":
                    sequence = self.mouse_delta(sequence)
                new_valignments[entry][seq] = sequence[:108].ljust(108).replace(" ", ".")
                if new_valignments[entry][seq][103] != "C" or new_valignments[entry][seq][22] != "C":
                    sys.stderr.write(
                        "Warning - this alignment doesn't feature CYS at position 23 and/or position 104.\n")
                    sys.stderr.write("%s,%s,%s\n" % (new_valignments[entry][seq], entry, seq))

        return new_valignments

    def mouse_delta(self, sequence):
        """
        Mouse delta chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.

        This is particularly bad because alignment is not even consistent within the chain and species!!!

        Remove and return
        """
        # Check in here because not all are bad...recheck again in the format v genes just to make sure.
        if sequence[103] != "C" or sequence[22] != "C":
            return sequence[: 8] + sequence[9:85] + sequence[86:]
        return sequence

    def rhesus_lambda(self, sequence):
        """
        Rhesus lambda chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
        Remove and return
        """
        return sequence[:20] + sequence[21:51] + sequence[53:]

    def mouse_alpha(self, sequence):
        """
        Mouse alpha chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
        Remove and return
        """
        return sequence[:8] + sequence[9:85] + sequence[86:]

    def combine_sequences(self, vsequences, jsequences):
        """
        Do a pairwise combination of the v and j sequences to get putative germline sequences for the species.

        """
        combined_sequences = {}
        for v in vsequences:
            vspecies, vallele = v
            for j in jsequences:
                _, jallele = j
                combined_sequences[("%s_%s_%s" % (vspecies, vallele, jallele)).replace(" ", "_")] = vsequences[v] + \
                                                                                                    jsequences[j]
        return combined_sequences

    def make_putative_alignments(self, vsequences, jsequences, calignments=None):
        all_sequences = {}
        for species, chain_type in vsequences:
            if (species, chain_type) not in jsequences or (species, chain_type) not in vsequences: continue
            combined_sequences = self.combine_sequences(vsequences[(species, chain_type)], jsequences[(species, chain_type)])
            all_sequences[(species, chain_type)] = combined_sequences
            self.output_stockholm(combined_sequences, "%s_%s" % (self.species_translations[species], chain_type))

        # Write just the V and J combinations
        self.output_stockholm_all(all_sequences)

        # Write the V and J combinations and the c-domains
        self.output_stockholm_all_and_C(all_sequences, calignments)

    def write_germlines(self, vsequences, jsequences):
        """
        Compile a dictionary containing all the v and j germline sequences.
        """

        all_gene_alignments = {"J": {}, "V": {}}

        for species, chain_type in vsequences:
            for ((_, gene), seq) in vsequences[(species, chain_type)].items():
                assert len(seq) == 108, species + _ + gene + chain_type + _ + seq + str(len(seq))
                try:
                    all_gene_alignments["V"][chain_type][self.species_translations[species]][gene] = seq.replace(".", "-") + "-" * 20
                except KeyError:
                    try:
                        all_gene_alignments["V"][chain_type][self.species_translations[species]] = {
                            gene: seq.replace(".", "-") + "-" * 20}
                    except KeyError:
                        try:
                            all_gene_alignments["V"][chain_type] = {
                                self.species_translations[species]: {gene: seq.replace(".", "-") + "-" * 20}}
                        except KeyError:
                            all_gene_alignments["V"] = {
                                chain_type: {self.species_translations[species]: {gene: seq.replace(".", "-") + "-" * 20}}}

            for ((_, gene), seq) in jsequences.get((species, chain_type), {}).items():
                assert len(seq) == 20
                try:
                    all_gene_alignments["J"][chain_type][self.species_translations[species]][gene] = "-" * 108 + seq.replace(".",
                                                                                                                "-")
                except KeyError:
                    try:
                        all_gene_alignments["J"][chain_type][self.species_translations[species]] = {
                            gene: "-" * 108 + seq.replace(".", "-")}
                    except KeyError:
                        try:
                            all_gene_alignments["J"][chain_type] = {
                                self.species_translations[species]: {gene: "-" * 108 + seq.replace(".", "-")}}
                        except KeyError:
                            all_gene_alignments["J"] = {
                                chain_type: {self.species_translations[species]: {gene: "-" * 108 + seq.replace(".", "-")}}}
            self.output_python_lookup(all_gene_alignments)

    def output_python_lookup(self, all_gene_alignments, path=None):
        """
        Format a lookup table for the germline sequences. This can then be used by the final program.
        """

        if path is None:
            path = self.alignment_path
        filename = os.path.join(path, "germlines.py")
        with open(filename, 'w') as outfile:
            print("all_germlines = " + repr(all_gene_alignments), file=outfile)

    def write_stockholm(self, sequences, ID, outfile):
        print("# STOCKHOLM 1.0", file=outfile)
        print("#=GF ID %s" % ID, file=outfile)

        pad_length = max(list(map(len, list(sequences.keys())))) + 1
        for s in sequences:
            print(s.replace(" ", "_").ljust(pad_length), sequences[s].replace(".", "-"), file=outfile)
        print("#=GC RF".ljust(pad_length), "x" * len(sequences[s]), file=outfile)
        print("//", file=outfile)

    def output_C_alignments(self, alignments, c_name):
        """
        Write a stockholm for all C domains.
        """
        for species, chain_type in alignments:
            self.output_stockholm(alignments[(species, chain_type)],
                             "%s_%s_%s" % (c_name, self.species_translations[species], chain_type))

    def output_stockholm_all_and_C(self, all_sequences, all_C_alignments, path=None):
        """
        Output a minimal stockholm alignment file for all sequences.
        """
        if path is None:
            path = self.alignment_path

        filename = os.path.join(path, "ALL_AND_C.stockholm")
        with open(filename, "w") as outfile:
            for species, chain_type in all_sequences:
                sequences = all_sequences[(species, chain_type)]
                l = len(list(sequences.values())[0])
                assert all([1 if l == len(sequences[s]) else 0 for s in
                            sequences]), "Not all sequences in alignment are the same length"
                self.write_stockholm(sequences, "%s_%s" % (self.species_translations[species], chain_type), outfile)

            for c_name in all_C_alignments:
                for species, chain_type in all_C_alignments[c_name]:
                    self.write_stockholm(all_C_alignments[c_name][(species, chain_type)],
                                    "%s_%s_%s" % (self.species_translations[species], chain_type, c_name), outfile)
        return filename

    def run_hmmbuild(self):
        output_hmm = os.path.join(self.package_directory, 'HMMs')
        if os.path.exists(output_hmm) is False:
            os.mkdir(output_hmm)
        o = os.path.join(output_hmm, 'ALL.hmm')
        i = os.path.join(self.alignment_path, 'ALL.stockholm')
        subprocess.call(f'hmmbuild --hand {o} {i}', shell = True)
        subprocess.call(f'hmmpress -f {o}', shell = True)
        o = os.path.join(output_hmm, 'ALL_AND_C.hmm')
        i = os.path.join(self.alignment_path, 'ALL_AND_C.stockholm')
        subprocess.call(f'hmmbuild --hand {o} {i}', shell=True)
        subprocess.call(f'hmmpress -f {o}', shell = True)
    def output_stockholm_all(self, all_sequences, path=None):
        """
        Output a minimal stockholm alignment file for all sequences.
        """
        if path is None:
            path = self.alignment_path

        filename = os.path.join(path, "ALL.stockholm")
        with open(filename, "w") as outfile:
            for species, chain_type in all_sequences:
                sequences = all_sequences[(species, chain_type)]
                l = len(list(sequences.values())[0])
                assert all([1 if l == len(sequences[s]) else 0 for s in
                            sequences]), "Not all sequences in alignment are the same length"
                self.write_stockholm(sequences, "%s_%s" % (self.species_translations[species], chain_type), outfile)

        return filename

    def output_stockholm(self, sequences, name, path=None):
        """
        Output a minimal stockholm alignment file.
        """
        if path is None:
            path = self.alignment_path

        filename = os.path.join(path, "%s.stockholm" % name)
        l = len(list(sequences.values())[0])

        assert all([1 if l == len(sequences[s]) else 0 for s in
                    sequences]), "Not all sequences in alignment are the same length"

        with open(filename, "w") as outfile:
            self.write_stockholm(sequences, name, outfile)

        return filename

    def run_process(self):
        """
        Read in the raw v and j alignments
        Format them and combine the sequences
        """
        print("\nFormatting alignments\n")
        valignments, jalignments = {}, {}
        all_valignments, all_jalignments = {}, {}
        ccalignments, c1alignments, c2alignments, c3alignments = {}, {}, {}, {}

        print("IGs")
        for species in self.ig_species:
            for chain_type in "HKL":
                if not os.path.isfile(os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type))):
                    print(f'missing', os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type)))
                    continue

                print(species, chain_type)
                valignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type)), region_name="V-REGION")
                jalignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sJ.fasta" % (species, chain_type)), region_name="J-REGION")

                ### Comment out if you want constant regions?
                if chain_type == "H":
                   c1alignments[ (species, chain_type) ]   = self.read_alignment( os.path.join( self.fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="CH1")
                   c2alignments[ (species, chain_type) ]   = self.read_alignment( os.path.join( self.fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="CH2")
                   c3alignments[ (species, chain_type) ]   = self.read_alignment( os.path.join( self.fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="CH3")
                else:
                   ccalignments[ (species, chain_type) ]   = self.read_alignment( os.path.join( self.fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="C-REGION")

                all_valignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type)), region_name="V-REGION",
                    read_all=True)
                all_jalignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sJ.fasta" % (species, chain_type)), region_name="J-REGION",
                    read_all=True)

        print("\nTRs")
        for species in self.tr_species:
            for chain_type in "ABGD":
                if not os.path.isfile(os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type))):
                    continue

                print(species, chain_type)
                valignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type)))
                jalignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sJ.fasta" % (species, chain_type)))
                all_valignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sV.fasta" % (species, chain_type)), read_all=True)
                all_jalignments[(species, chain_type)] = self.read_alignment(
                    os.path.join(self.fasta_path, "%s_%sJ.fasta" % (species, chain_type)), read_all=True)

        valignments = self.format_v_genes(valignments)
        self.j_alignments = jalignments
        jalignments = self.format_j_genes(jalignments)


        ccalignments = self.format_c_genes(ccalignments, 'CC')
        c1alignments = self.format_c_genes(c1alignments, 'C1')
        c2alignments = self.format_c_genes(c2alignments, 'C2')
        c3alignments = self.format_c_genes(c3alignments, 'C3')

        all_valignments = self.format_v_genes(all_valignments)
        all_jalignments = self.format_j_genes(all_jalignments)

        all_C_alignments = { "CC":ccalignments,"C1":c1alignments,"C2":c2alignments,"C2":c2alignments}

        self.write_germlines(all_valignments, all_jalignments)

        # Combine the alignments to make putative germline alignments (obviously no d gene in there for Hs)
        # Write them to a stockholm alignment file.
        combined_sequences = self.make_putative_alignments(valignments, jalignments, all_C_alignments)

        # Write the constant domains each to file.
        self.output_C_alignments(ccalignments, 'CC')
        self.output_C_alignments(c1alignments, 'C1')
        self.output_C_alignments(c2alignments, 'C2')
        self.output_C_alignments(c3alignments, 'C3')

        self.run_hmmbuild()




