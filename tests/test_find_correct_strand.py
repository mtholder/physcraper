import sys
import os
#from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#
import physcraper
import shutil
import physcraper.filter_by_local_blast as local_blast
from Bio.Seq import Seq
from copy import deepcopy
from Bio.Alphabet import generic_dna


# tests select_seq_local_blast_test
# tests if the building of select_seq_from local blast is selecting the right amount of species
workdir = "tests/output/test_find_correct_strand"
configfi = "tests/data/test.config"
threshold = 2
selectby = "random"
downtorank = None

absworkdir = os.path.abspath(workdir)


def test_find_correct_strand():
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



    full_seqs = {}


    orig = Seq("TTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGA", generic_dna)
    
    reverse = Seq("TTTGCCAATATAAAAACTTCTTTGTATCCTTATGAATCGGATAATACTATGTTATTTCCAATACTTATATTAATTCTATTTACTTTGTTCGTTGGATTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGTTAACCCCTTCTATAAATCTTTTACATAAAAATTCATACAATTCAATAGATTGGTATGAATTTTGTAAAGATGCAGTTTTTTCGGTTAGTATAGCCTCTTTCGGAATATTTATAGCATTTTTTTTATATCAACCTGTTTATTCATCTTTTCAAAATTTGGACTTAATTAATTCATTTGTTAAAATGGGTCCTAAGAGAAATTTTTCTGACAAAATAAAAAATGCTATATATGATTGGTCATATAATCGAGGTTACATAGATGCCTTTTATGGAACATTCTTAATTGCGGGGATGAGA", generic_dna)
    reverse = reverse[::-1]

    comp = Seq("TTTACTTTGTTTGTTGGATTCTTAGGAATTCCTTTCAATGAAGACATGGATATATTATCCAAATGGTTAACCCCTTCTATAAATCTTTTACATAAAAATTCATACAATTCAATAGATTGGTATGAATTTTGTAAAGATGCAGTTTTTTCGGTTAGTATAGCCTCTTTCGGAATATTTATAGCATTTTTTTTATATCAACCTGTTTATTCATCTTTTCAAAATTTGGACTTAATTAATTCATTTGTTAAAATGGGTCCTAAGAGGAATTTTTCTGACAAAATAAAAAATGCTATATA", generic_dna)
    comp = comp.complement()

    rcomp = Seq("AAAAATCTTATATCTGCCCTTACCGAGAAAGATCAAAAATTTTTTTTTTTTTTTGTAAAGGTCGATGGGGAAAAAAAGACTCCCCACCACACCACCAAAATCTGAGTGAAAATATACTAAAAAAAGAAAAAATGATTTTTTTATAATTCAGAAAACGGAAGAATTCTTCTTTTAGACATCTTATAATCTTATAATTAAGAGTCTATTTCATATTGATTAGAATAGGGA", generic_dna)
    rcomp = rcomp.reverse_complement()

    gb_acc_dict = {"EF538330.1": str(orig),
                    "JN790046.1": str(reverse),
                    "JN790051.1": str(comp),
                    "LM999594.1": str(rcomp),
                    }


    # copy from get_full_seq_batch
    if not os.path.exists("./tests/data/precooked/long_seqs/"):
        os.mkdir("./tests/data/precooked/long_seqs/")
    fn = "./tests/data/precooked/long_seqs/tmp_search.csv"
    fn_open = open(fn, "w+")
    for gb_acc in gb_acc_dict:
        fn_open.write("{}\n".format(gb_acc))
    fn_open.close()
    
    if not os.path.exists("./tests/data/precooked/long_seqs/all_full_seqs.fasta"):
        db_path = "{}/nt".format(filteredScrape.config.blastdb)
        if len(gb_acc_dict.keys()) >= 2:
            cmd1 = "blastdbcmd -db {}  -entry_batch {} -outfmt %f -out ./tests/data/precooked/long_seqs/tmp_full_seqs_perfile.fasta".format(db_path, fn)
        else:
            cmd1 = "blastdbcmd -db {}  -entry {} -outfmt %f -out ./tests/data/precooked/long_seqs/tmp_full_seqs_perfile.fasta".format(db_path, gb_acc_dict.keys()[0])
        print(cmd1)
        os.system(cmd1)
        # this is code used to copy needed file for test: cmd1 output is named: tmp_full_seqs_perfile.fasta, 
        f = open("./tests/data/precooked/long_seqs/all_full_seqs.fasta", "a+")
        internal_file = open("./tests/data/precooked/long_seqs/tmp_full_seqs_perfile.fasta")
        f.write(internal_file.read())

    file_name = "./tests/data/precooked/long_seqs/all_full_seqs.fasta"
    if (os.path.isfile(file_name)):
        dest = "{}/tmp/all_full_seqs.fasta".format(filteredScrape.workdir)
        shutil.copy(file_name, dest)

    # read in file to get full seq        
    fn = "{}/tmp/all_full_seqs.fasta".format(filteredScrape.workdir)

    # # assert that every gb_acc is found in file, sometimes only one of the methods seem to work, even though acc is in file.
    for item in gb_acc_dict:
        found = False
        with open(fn) as myfile:
            if item.split(".")[0] in myfile.read():
                found = True
        found2 = False
        df = file(fn)
        for line in df:
            if item.split(".")[0] in line:
                found2 = True
        assert found == True or found2 == True, (item, fn, found, found2)


    # full_seqs_func = filteredScrape.get_full_seq_batch(gb_acc_dict)
    f = open(fn)
    seq = ""
    full_seqs = {}  # dictionary to be filled with full seqs that are in the correct direction
    gb_acc_l_intern = set()
    count = 0

    for line in iter(f):
        line = line.rstrip().lstrip()
        if line[0]  != ">":  # '>' delimits identifier in those files
            seq += line
        elif line[0]  == ">":
            count += 1
            if seq != "" and len(gb_acc_l_intern) != 0:
                # seq = str(seq)
                before = deepcopy(full_seqs)
                full_seqs = filteredScrape.find_correct_strand(full_seqs, gb_acc_l_intern, gb_acc_dict, seq)
                assert before != full_seqs, (before, full_seqs)
    
                x = iter(gb_acc_l_intern).next()
                assert gb_acc_dict[x] != full_seqs[x], (gb_acc_dict[x], full_seqs[x])

            seq = ""
    
            # now get new gb_acc for this line
            splitline = line.split(">")
            del splitline[0]  # remove first element of list, which is an empty string bc of split ">"

            gb_acc_l_intern = set()  # reassign for new list
            # print(splitline)
            for item in splitline:
                if "emb" in item or "gb" in item or "dbj" in item:
                    gb_acc = physcraper.get_acc_from_blast(item)
                    if gb_acc in gb_acc_dict:
                        gb_acc_l_intern.add(gb_acc)
            seq = ""  # make new empty string for next seq
    f.close()
    if seq != "" and len(gb_acc_l_intern) != 0:
        full_seqs = filteredScrape.find_correct_strand(full_seqs, gb_acc_l_intern, gb_acc_dict, seq)

    assert set(full_seqs.keys()) == set(gb_acc_dict.keys()), ("missing in query_acc:", [x for x in set(full_seqs.keys()) if x not in set(gb_acc_dict.keys())], "missing in full seqs:", [x for x in set(gb_acc_dict.keys()) if x not in set(full_seqs.keys())])
    assert set(gb_acc_dict.keys()) == set(full_seqs.keys()), ([x for x in  gb_acc_dict if x not in full_seqs.keys()])

