
from pylab import *
import numpy as np
import os
import sys
import re
import time
import gzip
import shutil

from multiprocessing import Pool
from multiprocessing import Process
import subprocess

##
# \brief main
# \details
class merge_BMIQ:
    def __init__(self, settings):
        start = time.time()

        probename = settings.data_container.probes
        # for line in open( settings.datamatdir+'focal_probe.tsv', 'r'):
        #     probename.append(line.split(",")[0].replace('"',''))
        # for k in probename:
        #     print k
        # exit()

        files = []
        args_set = []
        i = 0
        print (settings.sampleinfo)
        for k in settings.data_container.sample_info:
            i = i + 1
            print(k)
            # if re.match("#", k):
            #     print (k[:-1]+" is skipped")
            #     continue
            # if i > 10:
            #     break
            filename = k[0].replace('"','')
            files.append(filename)
            filename = settings.datadirloc+filename+".tsv"

            args_set.append([probename, filename, settings.datamatdir, start])

        if os.path.isdir(settings.datamatdir+"infi_data/") is False:
            os.makedirs(settings.datamatdir+"infi_data/")

        ## most important process in this script
        p = Pool(settings.usecpu)
        p.map(pick_bata_value, args_set)
        p.close()

        print("beta_file was completed")

        beta_sum = {}
        for k in probename:
            beta_sum[k] = k
        for k in files:
            for line in open( settings.datamatdir+"infi_data/"+k+".tsv", 'r'):
                beta_sum[line.split("\t")[0]] += "\t" + line[:-1].split("\t")[1]

        f = open(settings.merged_data.replace(".gz",""), 'w')
        f.write("ID_REF")
        for k in files:
            f.write("\t"+k)
        f.write("\n")
        i = 0
        for k in probename:
            # print beta_sum[k]
            m = re.search("NA", beta_sum[k])
            if "NA" in beta_sum[k].split("\t"):
                i += 1
                print (beta_sum[k].split("\t")[0], " is skipped", i)
            # if m == True:
            #     print beta_sum[k]
            # f.write(k+beta_sum[k]+"\n")
                continue
            f.write(beta_sum[k]+"\n")
        f.close()

        if settings.merged_data[-3:] == ".gz":
            subprocess.call("gzip "+settings.merged_data.replace(".gz",""), shell=True)


        os.remove(settings.datamatdir+"focal_probe.tsv")
        shutil.rmtree(settings.datamatdir+"infi_data/")

def pick_bata_value(args):
    [probename, filename, datadir, start] = args
    start = float(start)

    i = 1
    j = 0
    log_bin = 100000
    sep = " "

    print ("BMIQ filename", filename)
    if os.path.exists(filename):
        # print filename[-2:]
        probeset = {}
        if filename[-2:] == "gz":
            for line in gzip.open(filename, 'rt'):
                j += 1
                if j % log_bin == 0:
                    elapsed_time = time.time() - start
                    print ("\t", j, elapsed_time, " in ", filename)
                if line.split(sep)[0] in probename:
                    probeset[line.split(sep)[0]] = line.split(sep)[0] + "\t" + line[:-1].split(sep)[1]
        else:
            for line in open(filename, 'r'):
                j += 1
                if j % log_bin == 0:
                    elapsed_time = time.time() - start
                    print ("\t", j, elapsed_time, " in ", filename)
                if line.split(sep)[0] in probename:
                    probeset[line.split(sep)[0]] = line.split(sep)[0] + "\t" + line[:-1].split(sep)[1]

        outname = filename.split("/")[-1].split(".")[0]
        f = open(datadir+"infi_data/"+outname+".tsv", 'w')
        for k in probename:
            if k not in probeset.keys():
            # if len(probeset[k].split("\t")) < i+1:
                # print "error", k
                probeset[k] = k + "\tNA"
            f.write(probeset[k]+"\n")
        f.close()

    elif os.path.exists(filename+".gz"):
        probeset = {}
        for line in gzip.open(filename+".gz", 'rt'):
            j += 1
            if j % log_bin == 0:
                elapsed_time = time.time() - start
                print ("\t", j, elapsed_time, " in ", filename)
            if line.split(sep)[0] == "":
                continue
            elif line.split(sep)[0] in probename:
                probeset[line.split(sep)[0]] = line.split(sep)[0] + "\t" + line[:-1].split(sep)[1]


        ## output each data file for debug
        outname = filename.split("/")[-1].split(".")[0]
        f = open(datadir+"infi_data/"+outname+".tsv", 'w')
        for k in probename:
            if k not in probeset.keys():
            # if len(probeset[k].split("\t")) < i+1:
                probeset[k] = k + "\tNA"
            f.write(probeset[k]+"\n")
        f.close()
    else:
        print (filename + " is absent")
        sys.exit()
