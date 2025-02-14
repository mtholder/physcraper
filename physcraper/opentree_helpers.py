import requests
import json
import sys
import os
import physcraper
from peyotl.api.phylesystem_api import PhylesystemAPI, APIWrapper
from peyotl.sugar import tree_of_life, taxomachine, treemachine, oti
from peyotl.nexson_syntax import (
    extract_tree,
    get_subtree_otus,
    extract_otu_nexson,
    PhyloSchema
)



from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel

from physcraper.helpers import cd, standardize_label, to_string




_VERBOSE = 1
_DEBUG = 1
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)



if sys.version_info < (3,):
    from urllib2 import HTTPError
#    from urllib2 import ConnectionError
else:
    from urllib.error import HTTPError

phylesystemref = "McTavish EJ, Hinchliff CE, Allman JF, Brown JW, Cranston KA, Holder MT,  Phylesystem: a gitbased data store for community curated phylogenetic estimates. Bioinformatics. 2015 31 2794-800. doi: 10.1093/bioinformatics/btv276\n"
synthref = "Redelings BD, Holder MT. A supertree pipeline for summarizing phylogenetic and taxonomic information for millions of species. PeerJ. 2017;5:e3058. https://doi.org/10.7717/peerj.3058 \n"


def bulk_tnrs_load(filename, ids_obj = None):
    otu_dict = {}
    with open(filename) as data_file:
        input_dict = json.load(data_file)
    for name in input_dict["names"]:
        i = 1
        otu = "Otu" + name['id']
        while otu in otu_dict.keys():
            otu = "{}_{}".format(otu, i)
            i+=1
        otu_dict[otu]={"^ot:originalLabel":name["originalLabel"]}
        if name.get("ottTaxonName"):
            otu_dict[otu]["^ot:ottTaxonName"] = name["ottTaxonName"]
        if name.get("ottId"):
            otu_dict[otu]["^ot:ottId"] = name["ottId"]
        for source in name.get("taxonomicSources"):
            taxsrc = source.split(":")
            otu_dict[otu]["^{}:taxon".format(taxsrc[0])] = source.strip(taxsrc[1])
    return otu_dict


def get_tree_from_synth(ott_ids, label_format="name", citation="cites.txt",):
    assert label_format in ['id', 'name', 'name_and_id']
    resp = treemachine.get_synth_tree_pruned(ott_ids=ott_ids, label_format=label_format)
    cites = ''
    sys.stdout.write("gathering citations")
    for study in resp['supporting_studies']:
        sys.stdout.write('.')
        study = study.split('@')[0]
        query = {"ot:studyId":study}
        new_cite = oti.find_studies(query_dict = query, verbose=True)
        #print new_cite[0].keys()
        cites = cites + '\n' + to_string(new_cite[0]['ot:studyPublicationReference']) + new_cite[0]['ot:studyPublication']
  #  cites = cites + '\n' +phylesystemref + synthref
    with open(citation,'w') as citfile:
        citfile.write(cites)
    sys.stdout.write("citations printed to {}\n".format(citation))
    tre = Tree.get(data=resp['newick'],
                   schema="newick",
                   suppress_internal_node_taxa=True)
    tre.suppress_unifurcations()
    return tre




# ATT is a dumb acronym for Alignment Tree Taxa object
def generate_ATT_from_phylesystem(aln,
                                  workdir,
                                  config_obj,
                                  study_id,
                                  tree_id,
                                  phylesystem_loc='api',
                                  ingroup_mrca=None):
    """gathers together tree, alignment, and study info - forces names to otu_ids.

    Study and tree ID's can be obtained by using python ./scripts/find_trees.py LINEAGE_NAME

    Spaces vs underscores kept being an issue, so all spaces are coerced to underscores when data are read in.

    :param aln: dendropy :class:`DnaCharacterMatrix <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix>` alignment object
    :param workdir: path to working directory
    :param config_obj: config class containing the settings
    :param study_id: OToL study id of the corresponding phylogeny which shall be updated
    :param tree_id: OToL corresponding tree ID as some studies have several phylogenies
    :param phylesystem_loc: access the github version of the OpenTree data store, or a local clone
    :param ingroup_mrca: optional.  OToL identifier of the mrca of the clade that shall be updated (can be subset of the phylogeny)
    :return: object of class ATT
    """
    assert isinstance(aln, datamodel.charmatrixmodel.DnaCharacterMatrix), \
            "your alignment `%s` ist not of type DnaCharacterMatrix" % aln
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore
    nexson = get_nexson(study_id, phylesystem_loc)
    newick = extract_tree(nexson,
                          tree_id,
                          PhyloSchema('newick',
                                      output_nexml2json='1.2.1',
                                      content="tree",
                                      tip_label="ot:originalLabel"))
    newick = newick.replace(" ", "_")  # UGH Very heavy handed, need to make sure happens on alignment side as well.
    tre = Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    # this gets the taxa that are in the subtree with all of their info - ott_id, original name,
    otus = get_subtree_otus(nexson, tree_id=tree_id)
    otu_dict = {}
    orig_lab_to_otu = {}
    treed_taxa = {}
    for otu_id in otus:
        otu_dict[otu_id] = extract_otu_nexson(nexson, otu_id)[otu_id]
        otu_dict[otu_id]["^physcraper:status"] = "original"
        otu_dict[otu_id]["^physcraper:last_blasted"] = None
        orig = otu_dict[otu_id].get(u"^ot:originalLabel").replace(" ", "_")
        orig_lab_to_otu[orig] = otu_id
        treed_taxa[orig] = otu_dict[otu_id].get(u"^ot:ottId")
    for tax in aln.taxon_namespace:
        if tax .label in otu_dict:
            sys.stdout.write("{} aligned\n".format(tax.label))
        else:
            try:
                tax.label = orig_lab_to_otu[tax.label].encode("ascii")
            except KeyError:
                sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignment. "
                                 "This may indicate a mismatch between tree and alignment\n".format(tax.label))
    # need to prune tree to seqs and seqs to tree...
    otu_newick = tre.as_string(schema="newick")
    ott_ids = get_subtree_otus(nexson,
                               tree_id=tree_id,
                               subtree_id="ingroup",
                               return_format="ottid")
    if ingroup_mrca:
        if type(ingroup_mrca) == list:
            ott_ids = set(ingroup_mrca)
            ott_mrca = get_mrca_ott(ott_ids)
        else:
            ott_mrca = int(ingroup_mrca)
    elif ott_ids:  # if no ingroup is specified, ott_ids will be none
        ott_mrca = get_mrca_ott(ott_ids)
    else:  # just get the mrca for teh whole tree
        ott_mrca = get_mrca_ott([otu_dict[otu_id].get(u"^ot:ottId") for otu_id in otu_dict])
    workdir = os.path.abspath(workdir)
    return physcraper.aligntreetax.AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, config_obj=config_obj)
    # newick should be bare, but alignment should be DNACharacterMatrix



def get_dataset_from_treebase(study_id, phylesystem_loc="api"):
    """Function is used to get the aln from treebase, for a tree that OpenTree has the mapped tree.
    """
    try:
        nexson = get_nexson(study_id, phylesystem_loc)
    except HTTPError as err:
        sys.stderr.write(err)
        sys.stderr.write("couldn't find study id {} in phylesystem location {}\n".format(study_id, phylesystem_loc))
    treebase_url = nexson['nexml'][u'^ot:dataDeposit'][u'@href']
    if 'treebase' not in nexson['nexml'][u'^ot:dataDeposit'][u'@href']:
        sys.stderr.write("No treebase record associated with study ")
        sys.exit(-2)
    else:
        tb_id = treebase_url.split(':S')[1]
        url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexml".format(tb_id)
        if _DEBUG:
            sys.stderr.write(url + "\n")
        dna = DataSet.get(url=url, schema="nexml")
        return dna



def OtuJsonDict(id_to_spn, id_dict):
    """Makes otu json dict, which is also produced within the openTreeLife-query.

     This function is used, if files that shall be updated are not part of the OpenTreeofLife project.
    It reads in the file that contains the tip names and the corresponding species names.
    It then tries to get the different identifier from the OToL project or if not from ncbi.

    Reads input file into the var sp_info_dict, translates using an IdDict object
    using web to call Open tree, then ncbi if not found.

    :param id_to_spn: user file, that contains tip name and corresponding sp name for input files.
    :param id_dict: uses the id_dict generated earlier
    :return: dictionary with key: "otu_tiplabel" and value is another dict with the keys '^ncbi:taxon',
                                                    '^ot:ottTaxonName', '^ot:ottId', '^ot:originalLabel',
                                                    '^user:TaxonName', '^physcraper:status', '^physcraper:last_blasted'
    """
    sys.stdout.write("Set up OtuJsonDict \n")
    sp_info_dict = {}
    nosp = []
    with open(id_to_spn, mode="r") as infile:
        for lin in infile:
            ottid, ottname, ncbiid = None, None, None
            tipname, species = lin.strip().split(",")
            clean_lab = standardize_label(tipname)
            assert clean_lab not in sp_info_dict, ("standardized label ('{}') of `{}` already exists".format(clean_lab, tipname))
            otu_id = "otu{}".format(clean_lab)
            spn = species.replace("_", " ")
            info = get_ott_taxon_info(spn)
            if info:
                ottid, ottname, ncbiid = info
            if not info:
                ncbi = NCBITaxa()
                name2taxid = ncbi.get_name_translator([spn])
                if len(name2taxid.items()) >= 1:
                    ncbiid = name2taxid.items()[0][1][0]
                else:
                    sys.stderr.write("match to taxon {} not found in open tree taxonomy or NCBI. "
                                     "Proceeding without taxon info\n".format(spn))
                    nosp.append(spn)
            ncbi_spn = None
            if ncbiid is not None:
                ncbi_spn = spn
            else:
                ncbi_spn = id_dict.ott_to_ncbi[ottid]
            sp_info_dict[otu_id] = {
                "^ncbi:taxon": ncbiid,

                "^ot:ottTaxonName": ottname,
                "^ot:ottId": ottid,
                "^ot:originalLabel": tipname,
                "^user:TaxonName": species,
                "^physcraper:status": "original",
                "^physcraper:last_blasted": None,
                }
            if ncbi_spn is not None:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = ncbi_spn
                sp_info_dict[otu_id]["^ncbi:TaxonName"] = ncbi_spn
            elif ottname is not None:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = ottname
            elif sp_info_dict[otu_id]['^user:TaxonName']:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = sp_info_dict[otu_id]['^user:TaxonName']
            assert sp_info_dict[otu_id]["^physcraper:TaxonName"]  # is not None
    return sp_info_dict


#####################################
def get_nexson(study_id, phylesystem_loc):
    """Grabs nexson from phylesystem"""
    phy = PhylesystemAPI(get_from=phylesystem_loc)
    nexson = phy.get_study(study_id)["data"]
    return nexson


def get_mrca_ott(ott_ids):
    """finds the mrca of the taxa in the ingroup of the original
    tree. The blast search later is limited to descendants of this
    mrca according to the ncbi taxonomy

    Only used in the functions that generate the ATT object.

    :param ott_ids: list of all OToL identifiers for tip labels in phylogeny
    :return: OToL identifier of most recent common ancestor or ott_ids
    """
    debug("get_mrca_ott")
    if None in ott_ids:
        ott_ids.remove(None)
    synth_tree_ott_ids = []
    ott_ids_not_in_synth = []
    for ott in ott_ids:
        r = check_if_ottid_in_synth(ott)
        if r == 1:
            synth_tree_ott_ids.append(ott)
        else:
            ott_ids_not_in_synth.append(ott)
    if len(synth_tree_ott_ids) == 0:
        sys.stderr.write('No sampled taxa were found in the current synthetic tree. '
                         'Please find and input and appropriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit(-3)
    mrca_node = tree_of_life.mrca(ott_ids=synth_tree_ott_ids, wrap_response=False)  # need to fix wrap eventually
    if u'nearest_taxon' in mrca_node.keys():
        tax_id = mrca_node[u'nearest_taxon'].get(u'ott_id')
        if _VERBOSE:
            sys.stdout.write('(v3) MRCA of sampled taxa is {}\n'.format(mrca_node[u'nearest_taxon'][u'name']))
    elif u'taxon' in mrca_node['mrca'].keys():
        tax_id = mrca_node['mrca'][u'taxon'][u'ott_id']
        if _VERBOSE:
            sys.stdout.write('(v3) MRCA of sampled taxa is {}\n'.format(mrca_node['mrca'][u'taxon'][u'name']))
    else:
        # debug(mrca_node.keys())
        sys.stderr.write('(v3) MRCA of sampled taxa not found. Please find and input an '
                         'appropriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit(-4)
    return tax_id

def get_ott_taxon_info(spp_name):
    """get ottid, taxon name, and ncbid (if present) from Open Tree Taxonomy.
    ONLY works with version 3 of Open tree APIs

    :param spp_name: species name
    :return:
    """
    #This is only used to write out the opentree info file. Could use NCBI id's instead of name, and likely be quicker.
    # debug(spp_name)
    try:
        res = taxomachine.TNRS(spp_name)["results"][0]
    except IndexError:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
        return 0
    if res['matches'][0]['is_approximate_match'] == 1:
        sys.stderr.write("""exact match to taxon {} not found in open tree taxonomy.
                          Check spelling. Maybe {}?""".format(spp_name, res['matches'][0][u'ot:ottTaxonName']))
        return 0
    if res["matches"][0]["is_approximate_match"] == 0:
        ottid = res["matches"][0]["taxon"][u"ott_id"]
        ottname = res["matches"][0]["taxon"][u"unique_name"]
        ncbi_id = None
        for source in res["matches"][0]["taxon"][u"tax_sources"]:
            if source.startswith("ncbi"):
                ncbi_id = source.split(":")[1]
        return ottid, ottname, ncbi_id
    else:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
        return 0

def check_if_ottid_in_synth(ottid):
    url = 'https://api.opentreeoflife.org/v3/tree_of_life/node_info'
    payload = json.dumps({"ott_id":int(ottid)})
    headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
    try:
        r = requests.post(url, data=payload, headers=headers)
        if r.status_code == 200:
            return 1
        elif r.status_code == 400:
            return 0
        elif r.status_code == 502:
            sys.stderr.write("Bad OpenTree taxon ID: {}".format(ottid))
            return 0
        else:
            sys.stderr.write("unexpected status code from node_info call: {}".format(r.status_code))
            return 0
    except requests.ConnectionError:
        sys.stderr.write("Connection Error - coud not get taxon information from OpenTree\n")
