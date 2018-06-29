import os
import sys
import re
from os import listdir
from os.path import isfile, join
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import nermal
import bisect
import random
import itertools
from collections import defaultdict
import argparse
import csv
import gzip
import string
from datetime import datetime
import sqlite3

#--------------------------------------------------------------------------------------------------
#Command line input parameters
parser = argparse.ArgumentParser(description='Spike-ins Tag Reader')
parser.add_argument('--m',
                    dest='m',
                    action='store_true',
                    help='Turns on one mismatch per seed length.')
parser.add_argument('--d',
                    dest='d',
                    action='store_true',
                    help='Turns on the distribution of reads for multiple alignments. (not recommended)')
parser.add_argument('--l',
                    metavar='-len',
                    required=True,
                    help='Seed length. Default: 12 (recommended)')

parser.add_argument('--i',
                    metavar='-input',
                    required=True,
                    help="Input directory, contains reads in fastq.gz " + \
                    "format. Do not include I1/I2 index file" + \
                    "s. The software will automatically filter the R1 and/" + \
                    "or R2 files with respect to paired-ended-ness. For ex" + \
                    "ample, if only R1s are present in the input directory" + \
                    ", there will not be any respect of index when filterin" + \
                    "g tagged reads. However, if R1 and R2 for a library a" + \
                    "re both present, the reads that are tagged in either " + \
                    "will be removed in both.")
parser.add_argument('--o',
                    metavar='-output',
                    required=True,
                    help="Name of output report.\nAll output reports will " + \
                    "contain the out_suffix specified. Option available fo" + \
                    "r labelling purposes.")
parser.add_argument('--f',
                    dest='f',
                    action='store_true',
                    help="Print filtered .fastq files (Optional). This wil" + \
                    "l produce fastqs (uncompressed) with spike-ins filter" +\
                    "ed out.")
parser.set_defaults(f=False)
parser.set_defaults(m=False)
parser.set_defaults(d=False)
#--------------------------------------------------------------------------------------------------

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def main():
    start_time = datetime.now()
    # print "calling ref gen"
    print "Parsing arguments"
    options = parser.parse_args()
    seed_length = int(options.l)
    distribute_multiples = options.d
    SEED_LENGTHS = [2,3,4,6,9,12]
    tag_bin_bin = {}
    total_seqs = {}

    if seed_length not in SEED_LENGTHS:
        seed_length = seedFix(seed_length)
    mismatch = 0
    if options.m:
        mismatch = 1

    db_name = "spikedb/noah_spike" + str(seed_length) + "." + str(mismatch) + ".db"

    print "Accessing " + str(db_name)
    
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    filt_ref_pos = {}
    cur.execute("SELECT * FROM filter_ref")
    for row in cur.fetchall():
        #print row
        if row[0] not in filt_ref_pos:
            filt_ref_pos[row[0]] = []
        filt_ref_pos[row[0]].append(row[2])
    cur.execute("SELECT * FROM match_dict")
    match_ref_pos = {}
    ids = {}
    for row in cur.fetchall():
        if row[0] not in match_ref_pos:
            match_ref_pos[row[0]] = []
        if row[5] not in ids:
            ids[row[5]] = 0
        match_ref_pos[row[0]].append((row[2],row[3],row[4],row[5]))

    cur.execute("SELECT id FROM ids")
    ref_names = cur.fetchall()
    ref_names = [ r[0] for r in ref_names ]
    collapsed_ref_names = list(set([ int(''.join(re.split('[a-z]+', str(r), flags=re.IGNORECASE))) if isInt(''.join(re.split('[a-z]+', str(r), flags=re.IGNORECASE))) else r for r in ref_names ]))
    collapsed_ref_names.sort()
    collapsed_ref_names = [ str(c) for c in collapsed_ref_names ]
    

    print("Counting tags")
    only_files = [f for f in listdir(options.i) if isfile(join(options.i, f))]
    only_files.sort()
    filter_reads = {}
    for o in only_files:
        print("Parsing " + str(o))

        seqs = []
        titles = []
        qualities = []

        
        
        handle = str(options.i) + "/" +str(o)

        name_fields = o.split('_')
        paired_end_reads_half_finished = False
        root_name = ""
        for e in xrange(len(name_fields)-2):
            root_name += str(name_fields[e])
            root_name += "_"
        root_name = root_name[:-1]
        print "Root:\t" + str(root_name)
        if root_name not in filter_reads:
            filter_reads[root_name] = []
        else:
            paired_end_reads_half_finished = True
        


        for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
            seqs.append(sequence)
            titles.append(title)
            qualities.append(quality)
            if paired_end_reads_half_finished == False:    
                filter_reads[root_name].append(False)

        wow_count = 0
        jeez_count = 0
        little_guys = 0
        big_boys = 0
        processing = []
        score = float(1.0)
        total_reads = len(seqs)
        total_tags = 0
        tag_bin = {"multiple" : 0, "not_aligned" : 0}

        
        for s in range(len(seqs)):
            if len(seqs[s]) <= seed_length:
                little_guys += 1
                
            else:
                big_boys += 1
                primer_scan = True
                wow = False
                z = 0
                screened_bases = []
                amplicon_lengths = []
                score = float(1.0)
                candidates = {}
                new_candidates = {}
                #print "***********" + str(seqs[s])
                cand_path = []
                meta_path = []
                iters = 1


                while primer_scan:
                    #print z
                    seq = seqs[s][z:z+seed_length]

                            
                    temp = []
                    temp_seqs = seqs[s]

                    new_candidates = {}
                    if z == 0:
                        if seq in filt_ref_pos:
                            for y in filt_ref_pos[seq]:
                                for m in match_ref_pos[y]:
                                    if m[3] == "1F" or m[3] == "1R":
                                        temp.append((z,seq,m[3],m[2],m[0]))
                                    #print "old"
                                    if m[3] not in tag_bin:
                                        tag_bin[m[3]] = 0
                                    if m[3] not in candidates:
                                        candidates[m[3]] = {}
                                    if m[0] > 0: 
                                        #if m[2]+seed_length
                                        candidates[m[3]][m[2]+seed_length] = m[0]
                                    else:
                                        candidates[m[3]][m[2]-seed_length] = m[0]
                                    #temp.append("")
                                    #print "eeeeeee: " + str(m[2]+seed_length)
                                    
                                                

                    else:
                        if seq in filt_ref_pos:
                            temp = []
                            temp_seqs = seqs[s]
                            
                            for y in filt_ref_pos[seq]:
                                for m in match_ref_pos[y]:
                                    if m[3] in candidates:
                                        temp.append((z,seq,m[3],m[2],m[0]))
                                        if m[2] in candidates[m[3]]:
                                            if candidates[m[3]][m[2]] == m[0]:
                                                if m[3] not in new_candidates:
                                                    new_candidates[m[3]] = {}
                                                
                                                if m[0] > 0: 
                                                    new_candidates[m[3]][m[2]+seed_length] = m[0]
                                                else:
                                                    new_candidates[m[3]][m[2]-seed_length] = m[0]
                                                


                    if z != 0:
                        candidates = new_candidates
                    cand_path.append(candidates)
                    meta_path.append(temp)
                    #print candidates
                    if z + seed_length >= 32 or len(candidates) == 0:
                        primer_scan = False
                        if len(candidates) != 0:
                            wow = True
                    z = seed_length * iters
                    iters += 1
                        
                if wow:
                    total_tags += 1
                    wow_count += 1

                    c_keys = candidates.keys()  
                    c_keys.sort()

                    if len(c_keys) == 1:
                        tag_bin[c_keys[0]] += 1.0
                    elif len(c_keys) > 1:
                        if distribute_multiples:
                            d_mult = 1.0 / float(len(c_keys))
                            for c in c_keys:
                                tag_bin[c] += d_mult
                        tag_bin["multiple"] += 1.0
                    filter_reads[root_name][s] = True
                else:
                    jeez_count += 1
                    tag_bin["not_aligned"] += 1.0
                    #print "cand path" 
                    f_true = False
                    r_true = False
                    hey_look_at_this = False

        tag_keys = tag_bin.keys()
        #tag_keys = [ str(f) if (f == "multiple" or f == "not_aligned") else int(f) for f in tag_keys ]
        tag_keys.sort()
        #tag_keys = [ str(f) for f in tag_keys ]
        for t in tag_keys:
            if tag_bin[t] != 0:
                print "\t* " + str(t) + ":\t" +str(tag_bin[t]) + "\t% total: " + str(format(float(tag_bin[t])/float(len(seqs))*100.0,'.2f'))
            

        print "\t* total reads: " + str(len(seqs))

        #print("small reads: " + str(float(little_guys) / float(little_guys+big_boys)))
        #print("primers: " + str(float(wow_count) / float(wow_count+jeez_count)))
        #print "read processing: " + str(processing)
        
       

        tag_bin_bin[o] = tag_bin
        total_seqs[o] = len(seqs)

    filt_final_csv = [ [f] for f in collapsed_ref_names ]
    filt_final_csv = [["libs"]] + filt_final_csv + [["multiple"],["not_aligned"],["total_reads"],["percent_aligned"],["top_tag"],["purity"]]
    #print filt_final_csv
    filt_percent_final_csv = [ [f] for f in collapsed_ref_names ]
    filt_percent_final_csv = [["libs"]] + filt_percent_final_csv + [["multiple"],["not_aligned"],["total_reads"],["percent_aligned"],["top_tag"],["purity"]]
    f_keys = [ f[0] for f in filt_final_csv ]
    for p in only_files:
        most_abundant = 0
        most_abundant_hits = 0
        tag_total = 0
        for q in xrange(len(f_keys)):
            if q == 0:
                filt_final_csv[q].append(p)
                filt_percent_final_csv[q].append(p)
            elif q == len(f_keys)-4:
                filt_final_csv[q].append(total_seqs[p])
                filt_percent_final_csv[q].append(total_seqs[p])
            elif q == len(f_keys)-3:
                if total_seqs[p] == 0:
                    filt_final_csv[q].append(0.0)
                    filt_percent_final_csv[q].append(0.0)
                else:
                    filt_final_csv[q].append(float(tag_total)/float(total_seqs[p])*100.0)
                    filt_percent_final_csv[q].append(float(tag_total)/float(total_seqs[p])*100.0)
            elif q == len(f_keys)-2:
                filt_final_csv[q].append(most_abundant)
                filt_percent_final_csv[q].append(most_abundant)
            elif q == len(f_keys)-1:
                if tag_total == 0:
                    filt_final_csv[q].append(0.0)
                    filt_percent_final_csv[q].append(0.0)
                else:
                    filt_final_csv[q].append(float(most_abundant_hits)/float(tag_total)*100.0)
                    filt_percent_final_csv[q].append(float(most_abundant_hits)/float(tag_total)*100.0)
            else:
                ff = unicode(str(f_keys[q]) + "F", 'utf-8')
                rr = unicode(str(f_keys[q]) + "R", 'utf-8')
                temp_hits = 0
                if rr in tag_bin_bin[p] and ff in tag_bin_bin[p]:
                    filt_final_csv[q].append(tag_bin_bin[p][ff] + tag_bin_bin[p][rr])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][ff] + tag_bin_bin[p][rr])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][ff] + tag_bin_bin[p][rr]
                elif rr in tag_bin_bin[p]:
                    filt_final_csv[q].append(tag_bin_bin[p][rr])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][rr])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][rr]     
                elif ff in tag_bin_bin[p]:
                    filt_final_csv[q].append(tag_bin_bin[p][ff])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][ff])/float(total_seqs[p])*100.0)
                    temp_hits = tag_bin_bin[p][ff]
                elif f_keys[q] == "multiple" or f_keys[q] == "not_aligned":
                    filt_final_csv[q].append(tag_bin_bin[p][f_keys[q]])
                    filt_percent_final_csv[q].append(float(tag_bin_bin[p][f_keys[q]])/float(total_seqs[p])*100.0)
                else:
                    filt_final_csv[q].append(0.0)
                    filt_percent_final_csv[q].append(0.0)
                tag_total += temp_hits
                if temp_hits > most_abundant_hits:
                    most_abundant = q
                    most_abundant_hits = temp_hits




    count_out = str(options.o) + "_counts.csv"
    with open(count_out, 'w') as of:
        a = csv.writer(of, delimiter=',')
        a.writerows(filt_final_csv)
        of.close()
    percent_out = str(options.o) + "_percent.csv"
    with open(percent_out, 'w') as of:
        a = csv.writer(of, delimiter=',')
        a.writerows(filt_percent_final_csv)
        of.close()

    print "\n\nPrinting Trimmed fastqs"
    if options.f:
        for o in only_files:
            name_fields = o.split('_')
            root_name = ""
            for e in xrange(len(name_fields)-2):
                root_name += str(name_fields[e])
                root_name += "_"
            root_name = root_name[:-1]
            p = o.split("/")
            out_name = str(options.o) + "_" + str(p[-1])
            print "Printing:\t" + str(out_name)
            f = 0
            handle = str(options.i) + "/" +str(o)
            with open(out_name, 'w') as nfq:
                for (title, sequence, quality) in FastqGeneralIterator(gzip.open(handle)):
                    if filter_reads[root_name][f] == False:
                        nfq.write("@"+str(title)+"\n")
                        nfq.write(str(sequence)+"\n")
                        nfq.write("+\n")
                        nfq.write(str(quality)+"\n")
                    f += 1
                nfq.close()

    print "Run time: " + str(datetime.now() - start_time)

def seedFix(seed_length):
    ret_length = 2
    print seed_length
    if seed_length < 2:
        print "Seed length is too small; defaulting to 1. Advising seed_length of 12."
    if seed_length == 5:
        print "Seed length not a factor of 36; defaulting to 6. Advising seed_length of 12."
        ret_length = 6
    if seed_length > 6 and seed_length < 9:
        print "Seed length not a factor of 36; defaulting to 8. Advising seed_length of 12."
        ret_length = 9
    if seed_length < 12 and seed_length > 9:
        print "Seed length not a factor of 36; defaulting to 12."
        ret_length = 12
    if seed_length > 12:
        print "Seed length capped at 12; defaulting to 12."
        ret_length = 12
    return ret_length





#--------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()