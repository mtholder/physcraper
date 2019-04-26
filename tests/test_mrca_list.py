# package import
import os
import json
import pickle
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts, FilterBlast

# define here your files
def test_mrca_list():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = [723076, 710505, 187044, 4727685, 4728090, 711399 ]

    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    ids = IdDicts(conf, workdir=workdir)

    # print(ids.mrca_ott, ids.mrca_ncbi)


    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = FilterBlast(data_obj, ids, ingroup_mrca)
    assert len(filteredScrape.mrca_ncbi_list) >= 2
    assert filteredScrape.mrca_ott_list == ingroup_mrca
    assert filteredScrape.mrca_ott_list != filteredScrape.mrca_ncbi_list
    assert filteredScrape.mrca_ott_list != filteredScrape.data.ott_mrca
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) == 63
    #fixed length to reflect new length after taxon_id as integer fix in line 2190


def test_no_mrca():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = None
    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    ids = IdDicts(conf, workdir=workdir)

    # print(ids.mrca_ott, ids.mrca_ncbi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = FilterBlast(data_obj, ids, ingroup_mrca)
    assert len(filteredScrape.mrca_ncbi_list) == 1
    assert filteredScrape.mrca_ott_list == ingroup_mrca
    assert filteredScrape.mrca_ott_list != filteredScrape.mrca_ncbi_list
    assert filteredScrape.mrca_ott_list != filteredScrape.data.ott_mrca

    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) <= 63

def test_higher_mrca():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = 557768
    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    ids = IdDicts(conf, workdir=workdir)

    # print(ids.mrca_ott, ids.mrca_ncbi)


    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = FilterBlast(data_obj, ids, ingroup_mrca)
    assert filteredScrape.mrca_ott_list == ingroup_mrca
    assert filteredScrape.mrca_ott_list != filteredScrape.mrca_ncbi_list
    assert filteredScrape.mrca_ott_list != filteredScrape.data.ott_mrca

    assert len(filteredScrape.mrca_ncbi_list) == 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) > 63

