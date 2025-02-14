"""functions to write/run/read files using the BLAST+ and a local database.
"""

import os
import sys
import subprocess
import numpy
import shutil
from Bio.Blast import NCBIXML

_DEBUG_MK = 1

"""Functions are used within the FilterBlast class to select sequences based on a local blast filtering step."""


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)
    

debug("Current local_blast version number: 10252018.0")


def del_blastfiles(workdir):
    """Deletes all files in the local blast folder.
    """
    try:
        shutil.rmtree(os.path.join(workdir, "blast"))
    except: 
        sys.stderr.write("Blast folder was not removed. Maybe it was not present?")


def run_filter_blast(workdir, blast_seq, blast_db, output=None):
    """Runs  a local blast to get measurement of differentiation between available sequences for the same taxon concept.

    The blast run will only be run if there are more sequences found than specified by the threshold value.
    When several sequences are found per taxon, blasts each seq against all other ones found for that taxon.
    The measure of differentiation will then be used to select a random representative from the taxon concept,
    but allows to exclude potential mis-identifications.
    In a later step (select_seq_by_local_blast) sequences will be selected based on the blast results generated here.

    # Note: has test, runs -> test_run_local_blast.py

    :param workdir: working directory
    :param blast_seq: Name of the file which contains the query seq
    :param blast_db:  Name of the file which contains the seq in fasta format that shall be blasted against
    :param output: Optional outpu filename
    """

    debug("run_filter_blast")
    general_wd = os.getcwd()
    os.chdir(os.path.join(workdir, "blast"))
    out_fn = "{}_tobeblasted".format(str(blast_seq))
    assert os.path.isfile(out_fn)
    cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(blast_db)
    os.system(cmd1)
    if output is None:
        cmd2 = "blastn -query {} -db {}_db -out output_{}.xml -outfmt 5".format(out_fn, blast_db, out_fn)
    else:
        cmd2 = "blastn -query {} -db {}_db -out {} -outfmt 5".format(out_fn, blast_db, output)
    os.system(cmd2)
    os.chdir(general_wd)


def calculate_mean_sd(hsp_scores):
    """Calculates standard deviation and mean of scores which are used as a measure of sequence differentiation
    for a given taxon.

    This is being used to select a random representative of a taxon later.

    Note: has test, runs: test_calculate_mean_sd.py

    :param hsp_scores: is a dict generated in read_filter_blast:
                        hsp_scores[gi_id] = {'hsp.bits': hsp.bits, 'hsp.score': hsp.score,
                                            'alignment.length': alignment.length, 'hsp.expect': hsp.expect}
    """
    total_seq = 0
    bit_sum = 0
    bit_l = []
    for gi_id in hsp_scores:
        total_seq += 1
        bit_sum += hsp_scores[gi_id]["hsp.bits"]
        bit_l.append(hsp_scores[gi_id]["hsp.bits"])
    bit_sd = float(numpy.std(bit_l))
    mean_hsp_bits = float(bit_sum / total_seq)
    mean_sd = {"mean": mean_hsp_bits, "sd": bit_sd}
    return mean_sd


def read_filter_blast(workdir, seq_d, fn):
    """Reads the files of the local blast run and returns sequences below a value
    (within the sd of the mean scores of the hsp.bit scores at the moment).

    (this is to make sure seqs chosen are representative of the taxon)

    Note: has test, runs: test_read_local_blast.py

    :param workdir: working directory
    :param seq_d: sequence dictionary from FilterBlast
    :param fn: filename
    :return seq_blast_score, dictionary with the blast values
    """
    general_wd = os.getcwd()
    os.chdir(os.path.join(workdir, "blast"))
    output_blast = "output_{}_tobeblasted.xml".format(fn)
    xml_file = open(output_blast)
    os.chdir(general_wd)
    blast_out = NCBIXML.parse(xml_file)
    hsp_scores = {}
    tries = 5
    for i in range(tries):
        try:
            for record in blast_out:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        gi_id = alignment.title.split(" ")[1]
                        if gi_id.isdigit():
                            gi_id = int(gi_id)
                        hsp_scores[gi_id] = {'hsp.bits': hsp.bits, 'hsp.score': hsp.score,
                                             'alignment.length': alignment.length, 'hsp.expect': hsp.expect}
        except ValueError:
            debug("rebuild the local blast db and try again")
            sys.stderr.write("{} blast file has a problem. Redo running it".format(fn))
            general_wd = os.getcwd()
            os.chdir(os.path.join(workdir, "blast"))
            subprocess.call(["rm", "{}_db.*".format(fn)])
            cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(fn)
            os.system(cmd1)
            cmd2 = "blastn -query {}_tobeblasted -db {}_db -out output_{}tobeblasted.xml -outfmt 5".format(fn, fn, fn)
            os.system(cmd2)
            os.chdir(general_wd)
            if i < tries - 1:  # i is zero indexed
                continue
            else:
                # debug("im going to raise")
                raise
        break
    # make values to select for blast search, calculate standard deviation,mean
    seq_blast_score = {}

    if hsp_scores != {}:
        mean_sd = calculate_mean_sd(hsp_scores)
        # select which sequences to use
        for gi_id in hsp_scores:  # use only seq that are similar to mean plus minus sd
            if (hsp_scores[gi_id]['hsp.bits'] >= mean_sd['mean'] - mean_sd['sd']) & \
                    (hsp_scores[gi_id]['hsp.bits'] <= mean_sd['mean'] + mean_sd['sd']):
                if gi_id in seq_d:
                    seq_blast_score[gi_id] = seq_d[gi_id]

    return seq_blast_score




def write_filterblast_query(workdir, seq_name, seq, fn):
    """Writes local query files which will be read by run_filter_blast.

    :param workdir: working directory
    :param seq_name: sequence identifier
    :param seq: sequence to write
    :param fn: file name
    :return: files with sequences written to it in fasta format
    """
    """
    """
    if not os.path.exists("{}/blast".format(workdir)):
        os.makedirs("{}/blast/".format(workdir))
    fnw = "{}/blast/{}_tobeblasted".format(workdir, fn)
    fi_o = open(fnw, "w")
    fi_o.write(">{}\n".format(seq_name))
    fi_o.write("{}\n".format(str(seq).replace("-", "")))
    fi_o.close()
    debug("writing query file for {} to {}\n".format(seq_name, fnw))

