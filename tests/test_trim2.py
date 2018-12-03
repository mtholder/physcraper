from dendropy import Tree, \
      DnaCharacterMatrix, \
      DataSet, \
      datamodel
from physcraper import generate_ATT_from_files, AlignTreeTax, OtuJsonDict, IdDicts, ConfigObj
import os
import json



seqaln= "tests/data/Senecioneae_ets_update/Senecioneae_ets2.fas"
mattype="fasta"
treefile= "tests/data/Senecioneae_ets_update/2taxon.tre"
schema_trf = "newick"
workdir="tests/output/test_trim2"
configfi = "tests/data/localblast.config"
id_to_spn = r"tests/data/Senecioneae_ets_update/nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)



if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi, interactive=False)
ids = IdDicts(conf, workdir=workdir)

otu_json = OtuJsonDict(id_to_spn, ids)
json.dump(otu_json, open(otu_jsonfi,"w"))


data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                 mattype=mattype, 
                                 workdir=workdir,
                                 treefile=treefile,
                                 schema_trf = schema_trf,
                                 otu_json=otu_jsonfi,
                                 ingroup_mrca=None)


for tax, seq in data_obj.aln.items():
  len_start = len(seq)


data_obj.trim()

for tax, seq in data_obj.aln.items():
  len_end = len(seq)

assert len_start ==  len_end


for tax, seq in data_obj.aln.items():
	len_start = len(seq)


data_obj.trim(taxon_missingness=0.5)

for tax, seq in data_obj.aln.items():
	len_end = len(seq)

assert len_start >  len_end
