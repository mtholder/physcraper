import sys
import re
import os
import io
import subprocess
import datetime
import glob
import json
import configparser
import pickle
import random
import contextlib
import time
import csv

from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel
from physcraper import ncbi_data_parser
from physcraper.opentree_helpers import get_mrca_ott
from physcraper.helpers import standardize_label, to_string, debug

_VERBOSE = 1
_DEBUG = 1

def generate_ATT_from_files(seqaln,
                            mattype,
                            workdir,
                            config_obj,
                            treefile,
                            otu_json,
                            schema_trf,
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem, use your own files instead.

    Spaces vs underscores kept being an issue, so all spaces are coerced to underscores when data are read in.

    Note: has test -> test_owndata.py

    :param seqaln: path to sequence alignment
    :param mattype: string containing format of sequence alignment
    :param workdir: path to working directory
    :param config_obj: config class including the settings
    :param treefile: path to phylogeny
    :param otu_json: path to json file containing the translation of tip names to taxon names, generated with OtuJsonDict()
    :param schema_trf: string defining the format of the input phylogeny
    :param ingroup_mrca: optional - OToL ID of the mrca of the clade of interest. If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca.

    :return: object of class ATT
    """

    # replace ? in seqaln with - : papara handles them as different characters

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    # use replaced aln as input
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
    assert aln.taxon_namespace
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore
    tre = Tree.get(path=treefile,
                   schema=schema_trf,
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    assert tre.taxon_namespace is aln.taxon_namespace, "tre and aln have not the same namespace."
    otu_newick = tre.as_string(schema=schema_trf)
    otu_dict = json.load(open(otu_json, "r"))
    if ingroup_mrca:
        mrca_ott = int(ingroup_mrca)
    else:
        ott_ids = [otu_dict[otu].get(u'^ot:ottId', ) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        mrca_ott = get_mrca_ott(ott_ids)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=mrca_ott, workdir=workdir,
                        config_obj=config_obj, schema=schema_trf)



class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.
        Hypothetically, all the keys in the  otu_dict should be clean.

        To build the class the following is needed:

          * **newick**: dendropy.tre.as_string(schema=schema_trf) object
          * **otu_dict**: json file including the otu_dict information generated earlier
          * **alignment**: dendropy :class:`DnaCharacterMatrix <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix>` object
          * **ingroup_mrca**: OToL identifier of the group of interest, either subclade as defined by user or of all tip labels in the phylogeny
          * **workdir**: the path to the corresponding working directory
          * **config_obj**: Config class
          * **schema**: optional argument to define tre file schema, if different from "newick"

        During the initializing process the following self objects are generated:

          * **self.aln**: contains the alignment and which will be updated during the run
          * **self.tre**: contains the phylogeny, which will be updated during the run
          * **self.otu_dict**: dictionary with taxon information and physcraper relevant stuff

               * key: a unique identifier (otu plus either "tiplabel of phylogeny" or for newly found sequences PS_number.
               * value: dictionary with the following key:values:

                    * '^ncbi:gi': GenBank identifier - deprecated by Genbank - only older sequences will have it
                    * '^ncbi:accession': Genbanks accession number
                    * '^ncbi:title': title of Genbank sequence submission
                    * '^ncbi:taxon': ncbi taxon identifier
                    * '^ot:ottId': OToL taxon identifier
                    * '^physcraper:status': contains information if it was 'original', 'queried', 'removed', 'added during filtering process'
                    * '^ot:ottTaxonName': OToL taxon name
                    * '^physcraper:last_blasted': contains the date when the sequence was blasted.
                    * '^user:TaxonName': optional, user given label from OtuJsonDict
                    * "^ot:originalLabel" optional, user given tip label of phylogeny
          * **self.ps_otu**: iterator for new otu IDs, is used as key for self.otu_dict
          * **self.workdir**: contains the path to the working directory, if folder does not exists it is generated.
          * **self.mrca_ott**: OToL taxon Id for the most recent common ancestor of the ingroup
          * **self.orig_seqlen**: list of the original sequence length of the input data
          * **self.gi_dict**: dictionary, that has all information from sequences found during the blasting.
            * key: GenBank sequence identifier
            * value: dictionary, content depends on blast option, differs between webquery and local blast queries
                * keys - value pairs for local blast:
                    * '^ncbi:gi': GenBank sequence identifier
                    * 'accession': GenBank accession number
                    * 'staxids': Taxon identifier
                    * 'sscinames': Taxon species name
                    * 'pident': Blast  percentage of identical matches
                    * 'evalue': Blast e-value
                    * 'bitscore': Blast bitscore, used for FilterBlast
                    * 'sseq': corresponding sequence
                    * 'title': title of Genbank sequence submission
                * key - values for web-query:
                    * 'accession':Genbank accession number
                    * 'length': length of sequence
                    * 'title': string combination of hit_id and hit_def
                    * 'hit_id': string combination of gi id and accession number
                    * 'hsps': Bio.Blast.Record.HSP object
                    * 'hit_def': title from GenBank sequence
                * optional key - value pairs for unpublished option:
                    * 'localID': local sequence identifier
          * **self._reconciled**: True/False,
          * **self.unpubl_otu_json**: optional, will contain the OTU-dict for unpublished data, if that option is used

        Following functions are called during the init-process:

          * **self._reconcile()**:
                removes taxa, that are not found in both, the phylogeny and the aln
          * **self._reconcile_names()**:
                is used for the own file stuff, it removes the character 'n' from tip names that start with a number

        The physcraper class is then updating:

          * self.aln, self.tre and self.otu_dict, self.ps_otu, self.gi_dict
    """

    def __init__(self, tree, otu_dict, alignment, ingroup_mrca, workdir, config_obj,
                 schema='newick', taxon_namespace=None):
        debug("build ATT class")
        self.aln = alignment
        assert isinstance(self.aln, datamodel.charmatrixmodel.DnaCharacterMatrix), \
                ("your aln '%s' is not a DnaCharacterMatrix" % alignment)
        self.tre = Tree.get(data=tree,
                            schema=schema,
                            preserve_underscores=True,
                            taxon_namespace=self.aln.taxon_namespace)
        assert (self.tre.taxon_namespace is self.aln.taxon_namespace), "tre and aln taxon_namespace are not identical"
        assert isinstance(otu_dict, dict), ("otu_dict '%s' is not of type dict" % otu_dict)
        self.otu_dict = otu_dict
        self.config = config_obj
        self.ps_otu = 1  # iterator for new otu IDs
        self._reconcile()
        self._reconcile_names()
        self.workdir = os.path.abspath(workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca), ("your ingroup_mrca '%s' is not an integer." % ingroup_mrca)
        self.mrca_ott = ingroup_mrca  # ott_ingroup mrca can be pulled directly from phylesystem
        self.orig_seqlen = []  # will get filled in later...
        self.gb_dict = {}  # has all info about new blast seq
        self._reconciled = False
        self.unpubl_otu_json = None

    def _reconcile(self):
        """Taxa that are only found in the tree, or only in the alignment are deleted.

        This checks that the tree "original labels" from phylesystem
        align with those found in the alignment.
        """
        debug("reconcile")
        treed_tax = set()
        for leaf in self.tre.leaf_nodes():
            treed_tax.add(leaf.taxon)
        aln_tax = set()
        for tax, seq in self.aln.items():
            aln_tax.add(tax)
        prune = treed_tax ^ aln_tax
        missing = [i.label for i in prune]
        if missing:
            errmf = 'NAME RECONCILIATION Some of the taxa in the tree are not in the alignment or vice versa' \
                    ' and will be pruned. Missing "{}"\n'
            errm = errmf.format('", "'.join(missing))
            sys.stderr.write(errm)
        del_aln = []
        del_tre = []
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_tax)
            if taxon in aln_tax:
                del_aln.append(taxon)
            if taxon in treed_tax:
                del_tre.append(taxon)
        self.aln.remove_sequences(del_aln)
        self.tre.prune_taxa(del_tre)
        for tax in prune:
            # potentially slow at large number of taxa and large numbers to be pruned
            found = 0
            for otu in self.otu_dict:
                if "^ot:originalLabel" in self.otu_dict[otu]:
                    if self.otu_dict[otu][u'^ot:originalLabel'] == tax.label:
                        self.otu_dict[otu]['^physcraper:status'] = "deleted in reconciliation"
                        found = 1
                elif otu == tax.label:
                    self.otu_dict[otu]['^physcraper:status'] = "deleted in reconciliation"
                    found = 1
            if found == 0:
                sys.stderr.write("lost taxon {} in reconcilliation \n".format(tax.label))
            self.aln.taxon_namespace.remove_taxon(tax)
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

    def _reconcile_names(self):
        """It rewrites some tip names, which kept being an issue when it starts with a number at the beginning.
        Then somehow a n was added to the tip names.

        :return: replaced tip names
        """
        for tax in self.aln.taxon_namespace:
            if tax.label in self.otu_dict.keys():
                pass
            else:
                found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ""
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for otu in self.otu_dict:
                    original = self.otu_dict[otu].get("^ot:originalLabel")
                    if original == tax.label or original == newname:
                        tax.label = otu
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))

    def prune_short(self):
        """Prunes sequences from alignment if they are shorter than specified in the config file,
         or if tip is only present in tre.

        Sometimes in the de-concatenating of the original alignment taxa with no sequence are generated
        or in general if certain sequences are really short. This removes those from both the tre and the alignment.

        has test: test_prune_short.py

        :return: prunes aln and tre
        """
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in
                            self.aln]
        avg_seqlen = sum(self.orig_seqlen) / len(self.orig_seqlen)
        #debug("average sequence length is {}".format(avg_seqlen))
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        prune = []
        aln_ids = set()
        for tax, seq in self.aln.items():
            aln_ids.add(tax.label)
            if len(seq.symbols_as_string().translate(None, "-?")) <= seq_len_cutoff:
                prune.append(tax)
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short "
                     "step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                self.remove_taxa_aln_tre(tax.label)
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]["^physcraper:status"] = "deleted in prune short"
        self.check_tre_in_aln()
        self._reconciled = 1

    def trim(self):
        """ It removes bases at the start and end of alignments, if they are represented by less than the value
        specified in the config file. E.g. 0.75 given in config means, that 75% of the sequences need to have a
        base present

        Ensures, that not whole chromosomes get dragged in. It's cutting the ends of long sequences.

        has test: test_trim.py
        """
        # debug('in trim')
        taxon_missingness = self.config.trim_perc
        i = 0
        seqlen = len(self.aln[i])
        while seqlen == 0:
            i = i + 1
            seqlen = len(self.aln[i])
        for tax in self.aln:
            if len(self.aln[tax]) != seqlen:
                sys.stderr.write("can't trim un-aligned inputs, moving on")
                return
        start = 0
        stop = seqlen
        cutoff = len(self.aln) * taxon_missingness
        for i in range(seqlen):
            counts = {"?": 0, "-": 0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ["?", "-"]:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:  # first ok column
                start = i
                break
        for i in range(seqlen, 0, -1):  # previously seqlen-1, that cuts off last character of aln, I changed it.
            counts = {'?': 0, '-': 0}
            for tax in self.aln:
                call = self.aln[tax][i - 1].label  # changing seqlen-1 to seqlen requires that we have here i-1
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                stop = i
                break
        # here alignment gets shortened to start:stop
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        # make sure that tre is presented in aln
        self.check_tre_in_aln()
        if _VERBOSE:
            sys.stdout.write("trimmed alignment ends to < {} missing taxa, "
                             "start {}, stop {}\n".format(taxon_missingness, start, stop))
        return


    def check_tre_in_aln(self):
        """Makes sure that everything which is in tre is also found in aln.

        Extracted method from trim. Not sure we actually need it there.
        """
        aln_ids = set()
        for taxon in self.aln:
            aln_ids.add(taxon.label)
        assert aln_ids.issubset(self.otu_dict.keys()), \
            ([x for x in aln_ids if x not in self.otu_dict.keys()], self.otu_dict.keys())

        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        for leaf in self.tre.leaf_nodes():
            if leaf.taxon not in aln_ids:
                self.tre.prune_taxa([leaf])
                self.tre.prune_taxa_with_labels([leaf.taxon])
                self.tre.prune_taxa_with_labels([leaf])
                treed_taxa.remove(leaf.taxon)
        assert treed_taxa.issubset(aln_ids), (treed_taxa, aln_ids)

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxa from aln and tre and updates otu_dict,
        takes a single taxon_label as input.

        note: has test, test_remove_taxa_aln_tre.py

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return: removes information/data from taxon_label
        """
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.discard_sequences([tax])
            self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found
            # the first prune does not remove it sometimes...
        if tax2:
            self.tre.prune_taxa([tax2])
            self.tre.prune_taxa_with_labels([taxon_label])
            self.tre.prune_taxa_with_labels([tax2])
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['^physcraper:status'] = "deleted, updated otu_dict but was never in tre or aln!"

    
    def get_otu_for_acc(self, gb_id):
        if gb_id in set([self.otu_dict[otu].get("^ncbi:accession",'UNK') for otu in self.otu_dict]):
            for otu in self.otu_dict:
                if self.otu_dict[otu].get("^ncbi:accession") == gb_id:
                    debug("tried to create OTU for {} but already had otu {}".format(gb_id, otu))
                    return otu
        else:
            return None

    def add_otu(self, gb_id, ids_obj):
        """ Generates an otu_id for new sequences and adds them into self.otu_dict.
        Needs to be passed an IdDict to do the mapping.

        :param gb_id: the Genbank identifier/ or local unpublished
        :param ids_obj: needs to IDs class to have access to the taxonomic information
        :return: the unique otu_id - the key from self.otu_dict of the corresponding sequence
        """
        # debug("add_otu function")
        otu_id = self.get_otu_for_acc(gb_id)
        if otu_id:
            return otu_id
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1
        ott_id = None
        #debug("trying to add an otu with accesion {}".format(gb_id))
        ncbi_id, tax_name = ncbi_data_parser.get_tax_info_from_acc(gb_id, self, ids_obj)
        if ncbi_id == None:
            debug("DID NOT ADD accession {} ncbi_id {}".format(gb_id, ncbi_id, tax_name))
            return None
        else:
            ncbi_id = int(ncbi_id)
        if ncbi_id in ids_obj.ncbi_to_ott.keys():
            #debug("ADDED OTU: accession {} ncbi_id {}".format(gb_id, ncbi_id, tax_name))
            ott_id = int(ids_obj.ncbi_to_ott[ncbi_id])
        else:
            debug("{} Ncbi id not found in ott_ncbi dictionaries\n".format(ncbi_id))
            ott_id = None
        if ott_id in ids_obj.ott_to_name:
            ott_name = ids_obj.ott_to_name[ott_id]
        else:
            ott_name = None
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]["^ncbi:title"] = self.gb_dict[gb_id]["title"]
        self.otu_dict[otu_id]["^ncbi:taxon"] = ncbi_id
        self.otu_dict[otu_id]["^ncbi:TaxonName"] = tax_name
        self.otu_dict[otu_id]["^ot:ottId"] = ott_id
        self.otu_dict[otu_id]["^physcraper:status"] = "query"
        self.otu_dict[otu_id]["^ot:ottTaxonName"] = ott_name
        self.otu_dict[otu_id]["^physcraper:last_blasted"] = None
        if gb_id[:6] == "unpubl":
            self.otu_dict[otu_id]["^physcraper:status"] = "local seq"
            self.otu_dict[otu_id]["^ot:originalLabel"] = self.gb_dict[gb_id]["localID"]
            self.otu_dict[otu_id]['^user:TaxonName'] = self.gb_dict[gb_id][u'^user:TaxonName']
        else:
            self.otu_dict[otu_id]["^ncbi:gi"] = self.gb_dict[gb_id]["^ncbi:gi"]
            self.otu_dict[otu_id]["^ncbi:accession"] = gb_id
        # get a name for the OTU, no matter from which source
        if tax_name is not None:
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = tax_name
        elif ott_name is not None:
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = ott_name
        elif self.otu_dict[otu_id].get('^user:TaxonName'):
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = self.otu_dict[otu_id]['^user:TaxonName']
        else:
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = "ACC_{}".format(gb_id)
        assert self.otu_dict[otu_id]["^physcraper:TaxonName"]  # is not None
        if _DEBUG >= 2:
            sys.stderr.write("acc:{} assigned new otu: {}\n".format(gb_id, otu_id))
        #debug("RETURNED OTU_ID {}".format(otu_id))
        return otu_id

    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """This writes out needed files for papara (except query sequences).
        Papara is finicky about trees and needs phylip format for the alignment.

        NOTE: names for tree and aln files should not be changed, as they are hardcoded in align_query_seqs().

        Is only used within func align_query_seqs.
        """
        #debug('write papara files')
        self.write_random_resolve_tre(treefilename)
        self.aln.write(path="{}/{}".format(self.workdir, alnfilename), schema="phylip")


    def write_random_resolve_tre(self, treefilename='random_resolve.tre'):
        self.tre.resolve_polytomies()
        self.tre.deroot()
        treepath= "{}/{}".format(self.workdir, treefilename)
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;", ";")  # Papara is diffffffficult about root
        tmptre = tmptre.replace("'", "_")
        fi = open("{}".format(treepath), "w")
        fi.write(tmptre)
        fi.close()
        return treepath

    def write_aln(self, alnname="physcraper.fas", alnschema="fasta"):
        alnpath = "{}/{}".format(self.workdir, alnname)
        self.aln.write(path=alnpath,
                       schema=alnschema)
        return os.path.abspath(alnpath)

    def write_files(self, treepath="physcraper.tre", treeschema="newick", alnpath="physcraper.fas", alnschema="fasta"):
        """Outputs both the streaming files, labeled with OTU ids.
        Can be mapped to original labels using otu_dict.json or otu_seq_info.csv"""
        #debug("write_files")
        self.tre.write(path="{}/{}".format(self.workdir, treepath),
                       schema=treeschema, unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.workdir, alnpath),
                       schema=alnschema)


    def write_labelled(self, label, filename = "labelled", direc='workdir', norepeats=True, add_gb_id=False):
        """output tree and alignment with human readable labels
        Jumps through a bunch of hoops to make labels unique.

        NOT MEMORY EFFICIENT AT ALL

        Has different options available for different desired outputs

        :param label: which information shall be displayed in labelled files: possible options:
                    '^ot:ottTaxonName', '^user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"
        :param treepath: optional: full file name (including path) for phylogeny
        :param alnpath:  optional: full file name (including path) for alignment
        :param norepeats: optional: if there shall be no duplicate names in the labelled output files
        :param add_gb_id: optional, to supplement tiplabel with corresponding GenBank sequence identifier
        :return: writes out labelled phylogeny and alignment to file
        """
        #debug("write labelled files")
        if direc == 'workdir':
            direc = self.workdir
        treepath = "{}/{}".format(direc, "{}.tre".format(filename))
        alnpath = "{}/{}".format(direc, '{}.fas'.format(filename))
        debug(treepath)
        assert label in ['^ot:ottTaxonName', '^user:TaxonName', '^physcraper:TaxonName',
                         "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        tmp_fasta = self.aln.as_string(schema="fasta")
        tmp_aln = DnaCharacterMatrix.get(data=tmp_fasta,
                                         schema="fasta",
                                         taxon_namespace=tmp_tre.taxon_namespace)
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:
            new_label = self.otu_dict[taxon.label].get(label, None)
            if new_label is None:
                if self.otu_dict[taxon.label].get("^ot:originalLabel"):
                    new_label = "orig_{}".format(self.otu_dict[taxon.label]["^ot:originalLabel"])
                else:
                    new_label = "ncbi_{}_ottname_{}".format(self.otu_dict[taxon.label].get("^ncbi:taxon", "unk"),
                                                            self.otu_dict[taxon.label].get('^physcraper:TaxonName', "unk"))
            new_label = str(new_label).replace(' ', '_')
            if add_gb_id:
                gb_id = self.otu_dict[taxon.label].get('^ncbi:accession')
                if gb_id is None:
                    gb_id = self.otu_dict[taxon.label].get("^ot:originalLabel")
                new_label = "_".join([new_label, str(gb_id)])
                sp_counter = 2
                if new_label in new_names and norepeats:
                    new_label = "_".join([new_label, str(sp_counter)])
                    sp_counter += 1
            else:
                if new_label in new_names and norepeats:
                    new_label = "_".join([new_label, taxon.label])
            taxon.label = new_label
            new_names.add(new_label)
        tmp_tre.write(path=treepath,
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
        tmp_aln.write(path=alnpath,
                      schema="fasta")

    def write_otus(self, filename = "otu_info", schema="table"):
        """Writes out OTU dict as json or table.

        :param filename: filename
        :param schema: either table or json format
        :return: writes out otu_dict to file
        """
        # TODO: schema is unused!
        assert schema in ["table", "json"]
        if schema == "json":
            with open("{}/{}.json".format(self.workdir, filename), "w") as outfile:
                json.dump(self.otu_dict, outfile)
        if schema == "table":
            all_keys =  set()
            for otu in self.otu_dict:
                all_keys.update(self.otu_dict[otu].keys())
            keys = list(all_keys) 
            header = ["otu_id"] + keys
            with open("{}/{}.csv".format(self.workdir, filename), "w") as outfile:
                outfile.write("\t".join(header)+"\n")
                for otu in self.otu_dict:
                    vals = [str(self.otu_dict[otu].get(key, "-")) for key in keys]
                    outfile.write("\t".join([to_string(otu)]+vals)+"\n")


    def dump(self, filename=None):
        """writes pickled files from att class"""
        if filename:
            if not os.path.exists(os.path.dirname(filename)):
                os.makedirs(os.path.dirname(filename))
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/att_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self, ofi)

