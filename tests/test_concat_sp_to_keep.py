from physcraper.concat import Concat
from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 
localblast = mark.localblast



#
workdir_its = "tests/data/precooked/concat_pre"
workdir_ets = "tests/data/precooked/concat_pre"
email = "martha.kandziora@yahoo.com"
pickle_fn = "final_ATT_checkpoint.p"

workdir_comb = "tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": "its_{}".format(pickle_fn)}, 
            "ets": {"workdir": workdir_ets, "pickle": "ets_{}".format(pickle_fn)}}

@localblast
def test():
    # get to test status

    concat = Concat(workdir_comb, email)
    for item in genelist.keys():
        concat.load_single_genes(genelist[item]['workdir'], genelist[item]["pickle"], item)

    concat.combine()
    concat.sp_seq_counter()
    sp_to_keep = concat.sp_to_keep()

    # print(sp_to_keep.keys())
    #
    # print("tests sp_to_keep")

    counter = 0
    sp_keep = []
    for sp in concat.sp_counter:
        for gene in concat.sp_counter[sp]:
            if concat.sp_counter[sp][gene] == 0:
                sp_keep.append(sp)

    assert set(sp_to_keep.keys()) == set(sp_keep)
       