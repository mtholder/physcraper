import sys
import os
#from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#
import physcraper
import physcraper.filter_by_local_blast as local_blast


workdir = "tests/output/test_select_seq_randomly"
configfi = "tests/data/test.config"
threshold = 2
selectby = "random"
downtorank = None

absworkdir = os.path.abspath(workdir)

from pytest import mark
localblast = mark.localblast

def test_select_seq_randomly():
    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))


    filteredScrape =  physcraper.FilterBlast(data_obj, ids)
    filteredScrape.add_setting_to_self(downtorank, threshold)

    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict()

    ##this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
    count = 0
    for tax_id in filteredScrape.sp_d:
        count_dict = filteredScrape.count_num_seq(tax_id)
        # print(tax_id)
        # print(count_dict)
        if count_dict["new_taxon"]:
            if count_dict["query_count"] <= threshold:
                count += count_dict["query_count"]
            if count_dict["query_count"] > threshold:
                count += threshold
        else:
            if count_dict["query_count"] >= 1:
                if count_dict["seq_present"] < threshold:
                    count += threshold-count_dict["seq_present"]
                if count_dict["seq_present"] > threshold:
                    count += 0
    
    filteredScrape.how_many_sp_to_keep(selectby)
    print(count, len(filteredScrape.filtered_seq))
    assert count == len(filteredScrape.filtered_seq) and count>0
 
