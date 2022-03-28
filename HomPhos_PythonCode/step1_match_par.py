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
import multiprocessing
import time
from concurrent.futures import ThreadPoolExecutor
#from multiprocessing import Queue
import multiprocessing as mp
from threading import Thread
import threading, queue

# Initialize the Queue
q = queue.Queue()

# signal(SIGPIPE,SIG_DFL)

# read a file as a list of lists
def read_csv_as_list(filename, separator, quote):
    """
    Inputs:
    filename  - name of file
    separator - character that separates fields
    quote     - character used to optionally quote fields
    Output:
    Returns a list of lists where each item in the list
    corresponds to a row in the file.
    """
    table = []
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=separator, quotechar=quote)
        for row in csvreader:
            table.append(row)
    return table


# generate the single reference fasta files for the matched xenopus references
def generate_files4blast_xen(dict_in, fasta_in, path_for_files):
    """
    Inputs:
        dict_in - the dictionary where the keys are the xenopus reference names that have phosphrylated residue
        fasta_in - the fasta file that has the xenopus reference information
        path_for_files
    Outputs:
        Nothing is returned
    """
    # read the fasta file into a dictionary
    fasta_dict = {}
    with open(fasta_in) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
        # seq is a FastaSequence object
            fasta_dict[seq.id] = seq.sequence_as_string()

    for key in dict_in.keys():
        if key in fasta_dict.keys():
        # generate the file
        #file_name_1 = row[col_ref] + ".fa"
            file_name = pathlib.Path(path_for_files, key + ".fa")
        # use fastaparser to write the single entry
            with open(file_name, 'w') as fasta_file:
                writer = fastaparser.Writer(fasta_file)
                writer.writefasta((key,fasta_dict[key]))

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
    out_file = pathlib.Path(o_f, name_part + ".fa")

    # Create the blast command
    command = ['blastp', '-query', in_file, '-subject', ref, \
      '-evalue',str(eval),'-outfmt',out_string, '-out', out_file]

    # Submit the command using subprocess
    subprocess.call(command)

# Function that manages the blast command queue
def worker():
     i = 0

     # conditional so that nothing is executed until the queue has been filled
     while True:

        # store initial time for estimates of duration and progress
        if i == 0:
            t = time.time()

        # get the variables from the next item in the queue
        i_f, name_part, ref, eval, out_string, o_f = q.get()
        # submit the inputs to the blast function
        call_blast(i_f, name_part, ref, eval, out_string, o_f)

        # not sure, remove the last object from the queue
        q.task_done()

        # make an estimate about the amount it will take for all of the blast files
        if i == 0:
            time_elapsed = time.time() - t
            q_size = q.qsize()
            expect_dur_hours = round((q_size * time_elapsed)/(60*60),2)

            print('Expected duration of xenopus sequence blasting is {} hours.'\
                .format(expect_dur_hours))

        if (i%200 == 0):

            t_elapsed_min = round((time.time() - t)/60,2);
            print('{} minutes have elapsed since blasting started.'.format(t_elapsed_min))

        i = i+1


def blast_all(dict_in, input_folder, sub_fasta_file, output_folder, e_val, ofmt_string):

    num_worker_threads = 2

    if __name__ == "step1_match_par":

        # turn-on the worker thread
        threading.Thread(target=worker, daemon=True).start()

        for qr_ref in dict_in.keys():
            q.put((input_folder, qr_ref, sub_fasta_file, e_val, ofmt_string, output_folder,))

        print('All task requests sent\n', end='')

        # block until all tasks are done
        q.join()
        print('All work completed')

            # q.put(ip = input_folder, name_part = qr_ref, \
            # ref = sub_fasta_file, eval = e_val, out_string =  ofmt_string, \
            # op = output_folder)



    # if __name__ == "step1_match_par":
    #
    #     with ThreadPoolExecutor(max_workers=4) as executor:
    #         for qr_ref in dict_in.keys():
    #
    #             blast_command = ['blastp', '-query', pathlib.Path(input_folder, qr_ref + ".fa"), \
    #              '-subject', sub_fasta_file, \
    #             '-evalue',str(e_val),'-outfmt',ofmt_string,\
    #             '-out',pathlib.Path(output_folder, qr_ref + "-out.txt")]
    #
    #             executor.submit(call_blast, ip = input_folder, name_part = qr_ref, \
    #             ref = sub_fasta_file, eval = e_val, out_string =  ofmt_string, \
    #             op = output_folder)
    #
    #             time.sleep( 30 )

    # NUMBER_OF_TASKS = 10; #len(dict_in.keys())
    #
    # if __name__ == '__main__':
    #
    # pool = mp.Pool(NUMBER_OF_TASKS)
    #
    # for qr_ref in dict_in.keys():
    #     pool.apply_async(call_blast, (qr_ref, input_folder, e_val, ofmt_string, output_folder))
    #     print(qr_ref)
    #
    # pool.close()
    # pool.join()

    #child_processes = []
    # for qr_ref in dict_in.keys():
    #
    #     q_file_name = pathlib.Path(input_folder, qr_ref + ".fa")
    #     o_file_name = pathlib.Path(output_folder, qr_ref + "-out.txt")
    #     print(q_file_name)
    #
    #     command_line = ['blastp', '-query', q_file_name, '-subject', sub_fasta_file, \
    #     '-evalue',str(e_val),'-outfmt',ofmt_string,'-out',o_file_name]
    #     subprocess.call(command_line)
        #p = subprocess.Popen(command_line, shell=True, stdout = PIPE, stderr = PIPE)
        #child_processes.append(p)
        #time.sleep( 30 )

    #for cp in child_processes:
    #    cp.wait()

# read the file with xenopus references and phosphrylated residues as a dictionary
#   where the keys are the xenopus references and the values are the
#   phosphorylated residues
def read_phosresidues_as_dict(filename, separator, quote, key_ind, val_in, val_sep):
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
    table = {}
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=separator, quotechar=quote)
        for row in csvreader:
            phos_vals = row[val_in].split(val_sep)
            table[row[key_ind]] = phos_vals[1:]
    return table
