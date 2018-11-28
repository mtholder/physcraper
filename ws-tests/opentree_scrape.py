import physcraper
import sys
from dendropy import DnaCharacterMatrix


# Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir = "tests/output/opentree"
configfi = "tests/data/remotencbi.config"

sys.stdout.write("\nTesting 'opentree scrape (1 round)'\n")
conf = physcraper.ConfigObj(configfi, interactive=False)
print "1. {}".format(conf.email)
      
    
aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                                    workdir=workdir,
                                                    study_id=study_id,
                                                    tree_id=tree_id,
                                                    phylesystem_loc=conf.phylesystem_loc)

ids = physcraper.IdDicts(conf, workdir=workdir)

print "3. {}".format(ids.config.email)

data_obj.prune_short()
assert len(data_obj.aln) == 9
data_obj.write_files()
try:
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    scraper.run_blast_wrapper()
    scraper.read_blast_wrapper()
    scraper.remove_identical_seqs()
    scraper.generate_streamed_alignment()
    sys.stdout.write("\nTest opentree_scrape.py (round 1) passed\n")
except:
    sys.stdout.write("\nTest opentree_scrape.py (round 1) FAILED'\n")
    