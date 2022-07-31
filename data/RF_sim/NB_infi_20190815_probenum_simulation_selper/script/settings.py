# -*- coding: utf-8 -*-
import re

import os
import shutil
import sys
import gzip

import methylome_data

## \brief infiniumの結果をまとめるための情報をsetting fileから読み込む。
# 設定を読み込むだけ。
class settings:

    def __init__(self, setting_file):
        self.setting_file = setting_file
        self.initialize_variable()
        self.data_container = methylome_data.methylome_data()
        self.import_settings(setting_file)
        self.dir_check()
        self.get_sample_info()

    ## \brief 設定の初期化を行う。
    # setting fileで指定されるべき設定を予め設定しておく。ほとんどの情報はsetting fileの情報で上書きする。
    def initialize_variable(self):
        print("settings.initialize_variable:")

        ## infinium dataの一時保管場所
        self.datamatdir = "result/tmp/"
        ## 使うcpuの数。全部使いたいときはmultiprocessing.cpu_count()。
        self.usecpu = 5
        ## 出力するディレクトリを初期化するかどうか。データを使いまわしたいときは０にしておくこと。
        self.force_initialize_dir = 0
        ## サンプルの情報を含んだファイル。
        self.sampleinfo = "data/training.tsv"
        ## サンプルごとのbeta valueを含んだファイルのディレクトリ。
        self.datadirloc = ""

    ## setting fileの情報を読み込む。
    # setting fileは " = " で区切られていることを想定している。左側にパラメーターの名前、右側にパラメーター。
    def import_settings(self, setting_file):
        print ("settings.import_settings:")
        # if self.datamatdir == "":
        #     print ("\t\tSetting file isn't indicated.")
        #     print ("\t\tThis process is quited.")
        #     exit()
        if os.path.isfile(setting_file) is False:
            print ("\t\t", setting_file, "cannot be found.")
            print ("\t\tPlease make setting file. This process is quited.")
            exit()

        for line in open(setting_file, 'r'):
            if line == "\n":
                continue
            if line[0] == "#":
                continue

            # print(line)

            if line.split(" = ")[0] == "datamatdir":
                self.datamatdir = line[:-1].split(" = ")[1]
                if self.datamatdir[-1] != "/":
                    self.datamatdir += "/"
            elif line.split(" = ")[0] == "usecpu":
                self.usecpu = int(line[:-1].split(" = ")[1])
            elif line.split(" = ")[0] == "manifest":
                self.data_container.manifest = line[:-1].split(" = ")[1]
            elif line.split(" = ")[0] == "sampleinfo":
                self.sampleinfo = line[:-1].split(" = ")[1]
                print(line)
            elif line.split(" = ")[0] == "datadirloc":
                self.datadirloc = line[:-1].split(" = ")[1]
            elif line[:-1].split(" = ")[0] == "Design_Type" or line[:-1].split(" = ")[0] == "Infinium_Design_Type":
                if line[:-1].split(" = ")[1] == "\n":
                    continue
                elif line[:-1].split(" = ")[1] == "I":
                    self.data_container.Design_Type = "I"
                elif line[:-1].split(" = ")[1] == "II":
                    self.data_container.Design_Type = "II"
                else:
                    print ("\tWrong option in Design_Type", line[:-1].split(" = ")[1])
                    exit()
            elif line.split(" = ")[0] == "chromosome":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    tmp = line[:-1].split(" = ")[1].split(",")
                    for k in tmp:
                        self.data_container.chromosome.append(str(k))
            elif line.split(" = ")[0] == "position_gene":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    tmp = line[:-1].split(" = ")[1].split(",")
                    for k in tmp:
                        self.data_container.position_gene.append(str(k))
            elif line.split(" = ")[0] == "CG_island":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    tmp = line[:-1].split(" = ")[1].split(",")
                    for k in tmp:
                        self.data_container.CG_island.append(str(k))
            elif line.split(" = ")[0] == "Phantom5_Enhancers":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.Phantom5_Enhancers = "+"
            elif line.split(" = ")[0] == "DMR":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.DMR = "+"
            elif line.split(" = ")[0] == "450k_Enhancer":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.k_Enhancer = "+"
            elif line.split(" = ")[0] == "DNase_Hypersensitivity_NAME":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.DNase_Hypersensitivity_NAME = "+"
            elif line.split(" = ")[0] == "OpenChromatin_NAME":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.OpenChromatin_NAME = "+"
            elif line.split(" = ")[0] == "TFBS_NAME":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.TFBS_NAME = "+"
            elif line.split(" = ")[0] == "SNP_ID":
                if line.split(" = ")[1] == "\n":
                    continue
                elif line.split(" = ")[1] == '""':
                    continue
                else:
                    self.data_container.SNP_ID = "+"
            elif line.split(" = ")[0] == "merged_data":
                self.merged_data = line[:-1].split(" = ")[1]
            elif line.split(" = ")[0] == "result_dir":
                self.result_dir = line[:-1].split(" = ")[1]

    def dir_check(self):
        if os.path.isdir(self.datamatdir) is False:
            # shutil.rmtree(self.datamatdir)
            os.makedirs(self.datamatdir)
            os.makedirs(self.datamatdir+"infi_data/")


    ## \brief サンプルの情報をimport
    ## サンプルの情報をimportする。separatorはここで指定する。
    def get_sample_info(self):
        if os.path.isfile(self.sampleinfo) is False:
            print ("\t\t", self.sampleinfo, "is absent.")
            print ("\t\tPlease make sample info file. The process is quited.")
            exit()

        sep = "\t"
        if self.sampleinfo[-3:] == "csv":
            sep = ","
        # print(self.sampleinfo[-3:])
        for line in open(self.sampleinfo, 'r'):
            # if len(line.split(sep)) != 6:
            #     print ("\t\tWrong format at sample info.")
            #     print (line)
            #     exit()

            l = line[:-1].split(sep)
            if l[0] == "ID":
                continue
            if l[1] not in self.data_container.sample_category_name:
                # if l[1] not in self.data_container.sample_group_name:
                self.data_container.sample_category_name.append(l[1])
                self.data_container.sample_loc_bycategory[l[1]] = []
            if l[2] not in self.data_container.sample_group_name:
                # if l[1] not in self.data_container.sample_group_name:
                self.data_container.sample_group_name.append(l[2])
                self.data_container.sample_loc_bygroup[l[2]] = []
            self.data_container.sample_info.append(l)
            self.data_container.sample_name.append(l[0])
            self.data_container.sample_category.append(l[1])
            self.data_container.sample_group.append(l[2])
            self.data_container.sample_data[l[0]]=line[:-1].split(sep)
            self.data_container.sample_tested[l[0]]=0
            self.data_container.sample_correct[l[0]]=0
            self.data_container.sample_assigned[l[0]]={}

        i = 0
        while i < len(self.data_container.sample_name):
            self.data_container.sample_loc_bycategory[self.data_container.sample_category[i]].append(i)
            self.data_container.sample_loc_bygroup[self.data_container.sample_group[i]].append(i)
            # self.data_container.sample_loc_bygroup[self.data_container.sample_category[i]].append(i)
            i += 1

        for k in self.data_container.sample_loc_bygroup:
            print(self.data_container.sample_loc_bygroup[k])


    ## infiniumのdata matrixが読み込まれているかを調べ、読み込まれていなければ読み込む。
    # check_infinium_probes_data -> import_data_matrix
    def import_data_matrix(self):
        print ("\tsettings.data_check:")
        # print("pgrpgr")

        if len(self.data_container.probes_data) > 0:
            print ("\t\tInfinium data looks present")
            return

        print ("\t\tInfinium data is null, then import from infinium_data.txt")
        print(self.merged_data)
        if os.path.isfile(self.merged_data) is False:
            print ("\t\tmerged_data isn't annotated.")
            print ("\t\tPlease make setting file. This process is quited.")
            exit()


        # pgr = 1
        print(self.merged_data[-3:], self.merged_data)
        if self.merged_data[-3:] == ".gz":
            for line in gzip.open(self.merged_data, 'rt'):
                # pgr += 1
                # print(pgr, line.split("\t")[0])
                if line.split("\t")[0] == "ID_REF":
                    self.check_sample_order(line[:-1].split("\t")[1:])
                    self.data_container.probes = []
                    tp1=[]
                    tp1=line[:-1].split("\t")
                    for i1 in range(1,len(tp1)):
                        self.data_container.samples.append(tp1[i1])
                        self.data_container.probes_bysample[tp1[i1]]=''
                    continue

                if line == "\n":
                    continue

                self.data_container.probes.append(line.split("\t")[0])
                self.data_container.probes_data[line.split("\t")[0]] = list(map(float, line[:-1].split("\t")[1:]))

        else:
            for line in open(self.merged_data, 'r'):
                # pgr += 1
                # print(pgr, line.split("\t")[0])
                if line.split("\t")[0] == "ID_REF":
                    self.check_sample_order(line[:-1].split("\t")[1:])
                    self.data_container.probes = []
                    tp1=[]
                    tp1=line[:-1].split("\t")
                    for i1 in range(1,len(tp1)):
                        self.data_container.samples.append(tp1[i1])
                        self.data_container.probes_bysample[tp1[i1]]=''
                    continue

                if line == "\n":
                    continue

                self.data_container.probes.append(line.split("\t")[0])
                self.data_container.probes_data[line.split("\t")[0]] = list(map(float, line[:-1].split("\t")[1:]))
            # settings.data_container.probes_data[line.split("\t")[0]] = line[:-1].split("\t")[1:]

        #     tp1=line[:-1].split("\t")
        #     for i1 in range(1,len(tp1)):
        #         # print(i1, tp1[i1], tp1)
        #         # print(self.data_container.probes_bysample)
        #         # print(self.data_container.probes_bysample[self.data_container.samples[i1-1]])
        #         if(tp1[i1]==''):
        #             print("null data in beta matrix",i1)
        #         else:
        #             # print(self.data_container.probes_bysample[self.data_container.samples[i1-1]])
        #             # print(tp1[i1])
        #             self.data_container.probes_bysample[self.data_container.samples[i1-1]]+="\t"+tp1[i1]
        #
        # #print (probes_bysample)
        # for k1 in sorted(self.data_container.probes_bysample):
        #     tp1=[]
        #     tp1=self.data_container.probes_bysample[k1].split("\t")



    def check_sample_order(self, header):
        print ("\tcal_stat.check_sample_order:")

        if len(self.data_container.sample_name) != len(header):
            print ("\t\tWrong sample name order. Please check the order of sample_info.tsv and infinium_data.txt")
            print ("\t\tself.data_container.sample_name:",len(self.data_container.sample_name))
            print ("\t\theader:",len(header))
            exit()

        i = 0
        while i < len(header):
            if header[i] != self.data_container.sample_name[i]:
                print ("\t\tWrong sample name order at position", i)
                print (header[i], "and" ,self.data_container.sample_name[i])
            i += 1
