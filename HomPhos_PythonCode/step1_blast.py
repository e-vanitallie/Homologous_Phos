#!/usr/bin/env python3

"""
    Code to generate the files for matching,
    call BLASTP for the paired Xenopus Human matches via the command line,
    and parse the file with Xenopus phospho-residues as dictionary
"""

import csv
import fastaparser
from pathlib import Path
import pathlib
import os
import subprocess
import time
from multiprocessing import Queue, Process
import multiprocessing as mp
import sys
import queue


# Create and submit the blast command
def call_blast(i_f, name_part, ref, eval, out_string, o_f):
    """
    Inputs:
    ip - input folder
    name_part - the Xenopus reference name to create input and output file names
    ref - the Human reference
    eval - the e-value to use a for the blast cut off
    out_string - the output format for the blast result
    op - output folder
    """

    # Create the input file object
    in_file = pathlib.Path(i_f, name_part + ".fa")
    # Create the output file object
    out_file = pathlib.Path(o_f, name_part + "-out.txt")

    # Create the blast command
    command = ['blastp', '-query', in_file, '-subject', ref, \
      '-evalue',str(eval),'-outfmt',out_string, '-out', out_file]

    # Submit the command using subprocess
    #subprocess.call(command)

    try:
        subprocess.check_call(command, stderr=subprocess.STDOUT)
    except CalledProcessError as err:
        print("Error ocurred for command:" + command)

# Function that manages the blast command queue
def worker_2(queue_in):

     print('A worker has started!')
     # conditional so that nothing is executed until the queue has been filled
     while True:
        #try:
        # get the variables from the next item in the queue
        i_f, name_part, ref, eval, out_string, o_f = queue_in.get()

        if i_f is None:
            print('Waiting for last blastp to finish.')
            time.sleep(30)
            print('A worker has ended!')
            break
        # submit the inputs to the blast function
        call_blast(i_f, name_part, ref, eval, out_string, o_f)
        # not sure, remove the last object from the queue
        queue_in.task_done()

def main_2(argumentList):

    # First create the dictionary of files that need to be blasted

    dict_for_blast = read_phosresidues_as_dict(argumentList[0], argumentList[1], \
    argumentList[2], argumentList[3])

    number_of_processes = int(argumentList[9])

    num_cores_computer = mp.cpu_count()
    print('This computer has {} cores, and {} were requested for mulitprocessing blast execution.'.\
    format(num_cores_computer, number_of_processes))

    tasks_to_accomplish_1 = mp.JoinableQueue()
    processes = []

    for qr_ref in dict_for_blast.keys():

        tasks_to_accomplish_1.put((argumentList[4], qr_ref, argumentList[5], argumentList[7], \
                argumentList[8], argumentList[6]),)

    # add the SENTINELS at the end of the queue, one each per worker, so that
    #   the process workers can end
    for add_sent in range(number_of_processes):
        tasks_to_accomplish_1.put((None, None, None, None, None, None,))

    # create the list of workers one at a time
    for w in range(number_of_processes):
         p = Process(target=worker_2, args=(tasks_to_accomplish_1,))
         processes.append(p)
         p.start()

    #
    for m in processes:
        m.join()

    print('Blasting complete!!')

    return True

# read the file with xenopus references and phosphrylated residues as a dictionary
#   where the keys are the xenopus references and the values are the
#   phosphorylated residues
def read_phosresidues_as_dict(filename, key_ind, val_in, val_sep):
    """
    Inputs:
    filename  - name of file
    separator - character that separates fields
    quote     - character used to optionally quote fields
    key_ind   - index to the key in the rows of the file
    val_ind   - index to the values in the rows of the file
    val_sep   - the seperator to separate the string of residues
    Output:
    Returns a list of lists where each item in the list
    corresponds to a row in the file.
    """
    separator = "\t"
    quote = '"'
    table = {}
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=separator, quotechar=quote)
        for row in csvreader:
            phos_vals = row[int(val_in)].split(val_sep)
            table[row[int(key_ind)]] = phos_vals[1:]
    return table

if __name__ == "__main__":
    # capture the inputs from the command line
    argumentList = sys.argv[1:]

    # call main_2 with the arguments
    return_val = main_2(argumentList)
