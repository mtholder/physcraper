import sys
import os
#from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#
import physcraper
import physcraper.filter_by_local_blast as local_blast


# tests select_seq_local_blast_test
# tests if the building of select_seq_from local blast is selecting the right amount of species
workdir = "tests/output/test_get_taxid_from_acc"
configfi = "tests/data/test.config"
threshold = 2
selectby = "random"
downtorank = None
absworkdir = os.path.abspath(workdir)

gb_acc = "FJ980341.1"
tax_id = 189247

gb_acc_2 = "JF284826.1"
tax_id2 = 292617


# # used to generate the file, used below

# conf = physcraper.ConfigObj(configfi, interactive=False)
# data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
# data_obj.workdir = absworkdir
# ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
# ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
# filteredScrape =  physcraper.FilterBlast(data_obj, ids)

# if not os.path.exists("{}/tmp".format(filteredScrape.workdir)):
#     os.mkdir("{}/tmp".format(filteredScrape.workdir))
# filteredScrape.config.blastdb == "/home/blubb/local_blast_db/"
# file_to_open = "./tests/data/precooked/testing_localdb/tax_id_{}.csv".format(gb_acc)
# # if not os.path.exists(file_to_open) or os.stat(file_to_open).st_size == 0:
# fn = "{}/tmp/tmp_search.csv".format(filteredScrape.workdir)
# fn_open = open(fn, "w+")                                        
# fn_open.write("{}\n".format(gb_acc))
# fn_open.close()
# cmd1 = "blastdbcmd -db /home/blubb/local_blast_db/nt  -entry_batch {} -outfmt %T -out ./tests/data/precooked/testing_localdb/tax_id_{}.csv".format(fn, gb_acc)
# os.system(cmd1)


# if not os.path.exists("{}/tmp".format(filteredScrape.workdir)):
#     os.mkdir("{}/tmp".format(filteredScrape.workdir))
# filteredScrape.config.blastdb == "/home/blubb/local_blast_db/"
# file_to_open = "./tests/data/precooked/testing_localdb/tax_id_{}.csv".format(gb_acc)
# # if not os.path.exists(file_to_open) or os.stat(file_to_open).st_size == 0:
# fn = "{}/tmp/tmp_search.csv".format(filteredScrape.workdir)
# fn_open = open(fn, "w+")                                        
# fn_open.write("{}\n".format(gb_acc_2))
# fn_open.close()
# cmd1 = "blastdbcmd -db /home/blubb/local_blast_db/nt  -entry_batch {} -outfmt %T -out ./tests/data/precooked/testing_localdb/tax_id_{}.csv".format(fn, gb_acc_2)
# os.system(cmd1)


def test_get_taxid_from_acc():
   

    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape =  physcraper.FilterBlast(data_obj, ids)

    tax_id_l = filteredScrape.get_taxid_from_acc(gb_acc)
    assert tax_id == tax_id_l[0], (tax_id, tax_id_l)

    tax_id_l = filteredScrape.get_taxid_from_acc(gb_acc_2)
    assert tax_id2 == tax_id_l[0], (tax_id2, tax_id_l)
