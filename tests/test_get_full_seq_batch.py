import sys
import os
#from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#
import physcraper
import physcraper.filter_by_local_blast as local_blast
from copy import deepcopy
from physcraper import get_acc_from_blast, get_gi_from_blast


# tests select_seq_local_blast_test
# tests if the building of select_seq_from local blast is selecting the right amount of species
workdir = "tests/output/test_get_full_seq_batch"
configfi = "tests/data/test.config"
threshold = 2
selectby = "blast"
downtorank = None

absworkdir = os.path.abspath(workdir)

from pytest import mark
localblast = mark.localblast

@localblast
def test_get_full_seq_batch():
    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape =  physcraper.FilterBlast(data_obj, ids)
    filteredScrape.add_setting_to_self(downtorank, threshold)

    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
   

    # simplified copy from read_blast_wrapper
    if blast_dir:
        filteredScrape.blast_subdir = os.path.abspath(blast_dir)
    else:
        if not os.path.exists(filteredScrape.blast_subdir):
            os.mkdir(filteredScrape.blast_subdir)
    
    if not filteredScrape._blasted:
        filteredScrape.run_blast_wrapper()
    assert os.path.exists(filteredScrape.blast_subdir)
    for taxon in filteredScrape.data.aln:
        fn = None
        if filteredScrape.config.blast_loc == "local":
            file_ending = "txt"
        if filteredScrape.config.gb_id_filename is True: #TODO what is this doing?
            fn = filteredScrape.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label) 
            if fn is None:
                fn = filteredScrape.data.otu_dict[taxon.label].get('^user:TaxonName', taxon.label)
            fn_path = "{}/{}.{}".format(filteredScrape.blast_subdir, fn, file_ending)
        else:
            fn_path = "{}/{}.{}".format(filteredScrape.blast_subdir, taxon.label, file_ending)
        if os.path.isfile(fn_path):
            if filteredScrape.config.blast_loc == 'local':  # new method to read in txt format
                # simplified copy from read_local_blast_query
                query_dict = {}
                with open(fn_path, mode="r") as infile:
                    for lin in infile:
                        sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, salltitles, sallseqid = lin.strip().split('\t')
                        gb_acc = get_acc_from_blast(sseqid)
                        gi_id = get_gi_from_blast(sseqid)                
                        sseq = sseq.replace("-", "") #TODO here is where we want to grab the full sequence MK: I wrote a batch query for the seqs we are interested. Makes it faster.
                        sscinames = sscinames.replace(" ", "_").replace("/", "_")
                        pident = float(pident)
                        evalue = float(evalue)
                        bitscore = float(bitscore)
                        stitle = salltitles

                        # get additional info only for seq that pass the eval
                        if evalue < float(filteredScrape.config.e_value_thresh):
                            if gb_acc not in filteredScrape.acc_ncbiid.keys(): # do not do it for gb_ids we already considered
                                # NOTE: sometimes there are seq which are identical & are combined in the local blast db...
                                # Get all of them! (they can be of a different taxon ids = get redundant seq info)
                                if len(sallseqid.split(";")) > 1:
                                    var_list = [gb_acc, gi_id, sseq, staxids, sscinames, pident, evalue, bitscore, stitle, sallseqid]
                                    query_dict = filteredScrape.get_new_seqs_for_mergedseq(var_list, query_dict)
                                    taxids_l = staxids.split(";")
                                else:  # if there are no non-redundant data
                                    staxids = int(staxids)
                                    filteredScrape.ids.spn_to_ncbiid[sscinames] = staxids
                                    if gb_acc not in filteredScrape.ids.acc_ncbi_dict:  # fill up dict with more information.
                                        filteredScrape.ids.acc_ncbi_dict[gb_acc] = staxids
                                    if gb_acc not in query_dict and gb_acc not in filteredScrape.newseqs_acc:
                                        query_dict[gb_acc] = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                                              'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                                              'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
                before = deepcopy(query_dict)

                # add data which was not added before and that passes the evalue threshhold
                gb_acc_d = {}
                for key in query_dict.keys():
                    if float(query_dict[key]["evalue"]) < float(filteredScrape.config.e_value_thresh):
                        gb_acc = query_dict[key]["accession"]
                        if gb_acc not in filteredScrape.data.gb_dict.keys()  or filteredScrape.config.add_lower_taxa is True:
                            # make dict with queries for full seq batch
                            gb_acc_d[gb_acc] = query_dict[key]["sseq"]
                if gb_acc_d != {}:
                    seq_d = filteredScrape.get_full_seq_batch(gb_acc_d)  # currently we get full seqs of all seqs found in a single blast search and which were not discarded until this point.
                    assert len(seq_d.keys()) == len(gb_acc_d.keys()), (len(seq_d.keys()), len(gb_acc_d.keys()))
                    for key in seq_d:
                        query_dict[key]["sseq"] = seq_d[key]  
                        filteredScrape.new_seqs[key] = query_dict[key]["sseq"]
                        filteredScrape.data.gb_dict[key] = query_dict[key]

                after = deepcopy(seq_d)

                different = False
                for item in before.keys():
                    seq_before = before[item]
                    seq_after = after[item]
                    if seq_before != seq_after:
                        different = True
                assert different is True

                assert before != after