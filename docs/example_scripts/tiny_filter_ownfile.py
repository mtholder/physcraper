# package import
import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts

# define here your files
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
workdir="tests/output/tiny_filter_own"
configfi = "tests/data/test.config"

otu_jsonfi = "{}/otu_dict.json".format(workdir)

# change to your filtering criteria
threshold = 2
selectby = "blast"
downtorank = "species"
shared_blast_folder = "/home/blubb/Documents/gitdata/physcraper/shared_runs/"
ingroup_mrca = 723076


# setup the run
if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi, interactive=False)
ids = IdDicts(conf, workdir=workdir, mrca=ingroup_mrca)


if os.path.exists(otu_jsonfi):
    print("load json")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, ids)
    json.dump(otu_json, open(otu_jsonfi, "w"))


# select a wrapper function, depending on what you want to do, see short tutorial:
wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         otu_jsonfi,
                         configfi,
                         selectby=selectby,
                         downtorank=downtorank,
                         ingroup_mrca=ingroup_mrca,
                         shared_blast_folder=shared_blast_folder)
