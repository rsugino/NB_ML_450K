# -*- coding: utf-8 -*-
import re

import os
import shutil
import sys

import numpy as np
import random


## infiniumの情報を管理する。
# This class includes infinium probe condition. This class is carried for later process.
class methylome_data:
    def __init__(self):
        print ("settings.methylome_data:")
        self.manifest = "/media/sugino/data_storage/ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/MethylationEPIC_v-1-0_B4.csv" ## illuminaが提供しているmanifest file。EPIC。
        # self.manifest = "GPL13534_HumanMethylation450_15017482_v.1.1.csv" ## illuminaが提供しているmanifest file。

        self.focal_gene = ""
        self.Design_Type = "" ## I or II, Infinium_Design_Type -> Design_Type
        self.chromosome = [] ## always character not number
        self.position_gene = [] ## TSS1500, TSS200, 5'UTR, 3'UTR, 1stExon, Body
        self.CG_island = [] ## N_Shelf, N_Shore, Island, S_Shore, S_Shelf. If set as null character "", it returns sea.
        # self.position_gene = ["TSS1500", "TSS200"] ##< TSS1500, TSS200, 5'UTR, 3'UTR, 1stExon, Body
        # self.CG_island = ["Island"] ##< N_Shelf, N_Shore, Island, S_Shore, S_Shelf

        # for EPIC manifest
        self.Phantom5_Enhancers = ""
        self.DMR = ""
        self.k_Enhancer = ""
        self.DNase_Hypersensitivity_NAME = ""
        self.OpenChromatin_NAME = ""
        self.TFBS_NAME = ""
        self.SNP_ID = ""

        # infinium_の部分を除いた。
        self.probes = [] ## List of focal probes Ordered. select_probe, import_data_matrix
        self.probes_anno = {}
        self.probes_data = {}
        self.probes_bysample={} # 使うサンプルのデータだけを抜き出してタブ区切りにする。
        self.samples = []
        self.sample_info = []
        self.sample_name = []
        self.sample_category = [] # 推定したいもの。INSS stage。
        self.sample_loc_bycategory = {} # 推定したいもの。INSS stage。
        self.sample_category_name = [] # 推定したいものの名前をkey,場所をリストで。INSS stage。

        self.sample_group = [] # data source (eg. target)を想定
        self.sample_loc_bygroup = {} # sample_groupごとのリスト。
        self.sample_group_name = [] ## グループのリスト。
        self.sample_data={} # sampleの情報をIDをkeyとするstringで読み込む
        self.sample_tested={}
        self.sample_correct={}
        self.sample_assigned={}

        ## key is the group name, item is plot information [num, color, marker]
        self.sample_mark = {}

    ## データを3つのグループにわける。
    #  分類の比率はこの関数の中で指定する。
    # データのグループは、コラムの2番めで指定する。ステージかつデータソースで分類したいときはそこのコラムで指定する。
    def split3(self):
        ## training用にデータを揃える。

        self.train_ratio = 0.7
        self.tune_ratio = 0
        self.test_ratio = 1 - self.train_ratio - self.tune_ratio

        self.train_samples_pos = []
        self.tune_samples_pos = []
        self.test_samples_pos = []

        balance_sample = "category" ##
        print("balance_sample option", balance_sample)
        if balance_sample == "category":
            self.sample_balance_category()
        elif balance_sample == "same":
            self.same_train_size()
        elif balance_sample == "group":
            self.sample_balance_group()
        elif balance_sample == "no_balance":
            self.sample_no_balance()
        elif balance_sample == "twoclass":
            self.sample_no_balance()
        elif balance_sample == "nonNA":
            self.nonNA()

        return self.train_samples_pos,self.tune_samples_pos,self.test_samples_pos

    # カテゴリー虫で分ける
    def sample_no_balance(self):

        train_datnum = int(len(self.sample_name)*0.7)
        tune_datnum = int(len(self.sample_name)*0)
        test_datnum = len(self.sample_name) - train_datnum - tune_datnum

        print("train_datnum", train_datnum)
        print("tune_datnum", tune_datnum)
        print("test_datnum", test_datnum)

        self.train_samples_pos.extend(random.sample(range(len(self.sample_name)),train_datnum))

        tmp1 = list(set(range(len(self.sample_name))).difference(self.train_samples_pos))
        self.tune_samples_pos.extend(random.sample(tmp1,tune_datnum))

        self.test_samples_pos.extend(list(set(tmp1).difference(self.tune_samples_pos)))

    def same_train_size(self):
        train_datnum={"MYCNA":5,"4":5,"4s":5,"other":5}
        tune_datnum={"MYCNA":5,"4":5,"4s":5,"other":5}
        test_datnum={"MYCNA":5,"4":5,"4s":5,"other":5}

        # train_datnum={"4Amp":5,"4Noamp":5,"3":5,"4s":5,"other":5}
        # tune_datnum={"4Amp":5,"4Noamp":5,"3":5,"4s":5,"other":5}
        # test_datnum={"4Amp":5,"4Noamp":5,"3":5,"4s":5,"other":5}

        for k in train_datnum:
            train_datnum[k] = 40
            # train_datnum[k] = int(len(self.sample_loc_bycategory[k])*0.7)
            tune_datnum[k] = int(len(self.sample_loc_bycategory[k])*0)
            test_datnum[k] = len(self.sample_loc_bycategory[k]) - train_datnum[k] - tune_datnum[k]

        print("train_datnum", train_datnum)
        print("tune_datnum", tune_datnum)
        print("test_datnum", test_datnum)

        for k1 in train_datnum:
            # print("train_datnum", k1,train_datnum[k1])
            # train_datnumで指定した数だけ抜き出す
            self.train_samples_pos.extend(random.sample(self.sample_loc_bycategory[k1],train_datnum[k1]))

            tmp1 = list(set(self.sample_loc_bycategory[k1]).difference(self.train_samples_pos))
            self.tune_samples_pos.extend(random.sample(tmp1,tune_datnum[k1]))

            self.test_samples_pos.extend(list(set(tmp1).difference(self.tune_samples_pos)))


    def sample_balance_category(self):
        train_datnum={"MYCNA":5,"4":5,"4s":5,"other":5}
        tune_datnum={"MYCNA":5,"4":5,"4s":5,"other":5}
        test_datnum={"MYCNA":5,"4":5,"4s":5,"other":5}

        # train_datnum={"4Amp":5,"4Noamp":5,"3":5,"4s":5,"other":5}
        # tune_datnum={"4Amp":5,"4Noamp":5,"3":5,"4s":5,"other":5}
        # test_datnum={"4Amp":5,"4Noamp":5,"3":5,"4s":5,"other":5}

        for k in train_datnum:
            train_datnum[k] = 40
            train_datnum[k] = int(len(self.sample_loc_bycategory[k])*0.7)
            tune_datnum[k] = int(len(self.sample_loc_bycategory[k])*0)
            test_datnum[k] = len(self.sample_loc_bycategory[k]) - train_datnum[k] - tune_datnum[k]

        print("train_datnum", train_datnum)
        print("tune_datnum", tune_datnum)
        print("test_datnum", test_datnum)

        for k1 in train_datnum:
            # print("train_datnum", k1,train_datnum[k1])
            # train_datnumで指定した数だけ抜き出す
            self.train_samples_pos.extend(random.sample(self.sample_loc_bycategory[k1],train_datnum[k1]))

            tmp1 = list(set(self.sample_loc_bycategory[k1]).difference(self.train_samples_pos))
            self.tune_samples_pos.extend(random.sample(tmp1,tune_datnum[k1]))

            self.test_samples_pos.extend(list(set(tmp1).difference(self.tune_samples_pos)))


    def sample_balance_group(self):

        train_datnum={"Target":5,"Henrich":5,"JNB":5,"Ackerman":5}
        tune_datnum={"Target":5,"Henrich":5,"JNB":5,"Ackerman":5}
        test_datnum={"Target":5,"Henrich":5,"JNB":5,"Ackerman":5}

        for k in train_datnum:
            train_datnum[k] = int(len(self.sample_loc_bygroup[k])*0.7)
            tune_datnum[k] = int(len(self.sample_loc_bygroup[k])*0)
            test_datnum[k] = len(self.sample_loc_bygroup[k]) - train_datnum[k] - tune_datnum[k]

        print("train_datnum", train_datnum)
        print("tune_datnum", tune_datnum)
        print("test_datnum", test_datnum)

        for k1 in train_datnum:
            print("train_datnum", k1,train_datnum[k1])
            # train_datnumで指定した数だけ抜き出す
            self.train_samples_pos.extend(random.sample(self.sample_loc_bygroup[k1],train_datnum[k1]))

            tmp1 = list(set(self.sample_loc_bygroup[k1]).difference(self.train_samples_pos))
            self.tune_samples_pos.extend(random.sample(tmp1,tune_datnum[k1]))

            self.test_samples_pos.extend(list(set(tmp1).difference(self.tune_samples_pos)))

    def twoclass(self):

        category = list(set(self.sample_category))
        print(category)

        train_datnum={}
        tune_datnum={}
        test_datnum={}
        for k in category:
            train_datnum[k]=0
            tune_datnum[k]=0
            test_datnum[k]=0

        for k in train_datnum:
            train_datnum[k] = int(len(self.sample_loc_bygroup[k])*0.7)
            tune_datnum[k] = int(len(self.sample_loc_bygroup[k])*0)
            test_datnum[k] = len(self.sample_loc_bygroup[k]) - train_datnum[k] - tune_datnum[k]

        print("train_datnum", train_datnum)
        print("tune_datnum", tune_datnum)
        print("test_datnum", test_datnum)

        for k1 in train_datnum:
            print("train_datnum", k1,train_datnum[k1])
            # train_datnumで指定した数だけ抜き出す
            self.train_samples_pos.extend(random.sample(self.sample_loc_bygroup[k1],train_datnum[k1]))

            tmp1 = list(set(self.sample_loc_bygroup[k1]).difference(self.train_samples_pos))
            self.tune_samples_pos.extend(random.sample(tmp1,tune_datnum[k1]))

            self.test_samples_pos.extend(list(set(tmp1).difference(self.tune_samples_pos)))

    def nonNA(self):

        category = list(set(self.sample_category))
        print("category", category)

        train_datnum={}
        tune_datnum={}
        test_datnum={}
        for k in category:
            if k == "NA":
                continue
            train_datnum[k]=0
            tune_datnum[k]=0
            test_datnum[k]=0

        for k in train_datnum:
            train_datnum[k] = int(len(self.sample_loc_bycategory[k])*0.7)
            tune_datnum[k] = int(len(self.sample_loc_bycategory[k])*0)
            test_datnum[k] = len(self.sample_loc_bycategory[k]) - train_datnum[k] - tune_datnum[k]

        print("train_datnum", train_datnum)
        print("tune_datnum", tune_datnum)
        print("test_datnum", test_datnum)

        for k1 in train_datnum:

            print("train_datnum", k1,train_datnum[k1])
            # train_datnumで指定した数だけ抜き出す
            self.train_samples_pos.extend(random.sample(self.sample_loc_bycategory[k1],train_datnum[k1]))

            tmp1 = list(set(self.sample_loc_bycategory[k1]).difference(self.train_samples_pos))
            self.tune_samples_pos.extend(random.sample(tmp1,tune_datnum[k1]))

            self.test_samples_pos.extend(list(set(tmp1).difference(self.tune_samples_pos)))


    ## データをとってくる。
    def data_pickup(self, train_samples_pos, tune_samples_pos, test_samples_pos):
        X_train_tmp = []
        X_tuning_tmp = []
        X_test_tmp = []
        ## data_matrixのサブセットを作る
        for k in self.probes_data.keys():
            X_train_tmp.append([self.probes_data[k][i] for i in train_samples_pos])
            X_tuning_tmp.append([self.probes_data[k][i] for i in tune_samples_pos])
            X_test_tmp.append([self.probes_data[k][i] for i in test_samples_pos])

        X_train = np.array(X_train_tmp).transpose()
        X_tuning = np.array(X_tuning_tmp).transpose()
        X_test = np.array(X_test_tmp).transpose()

        y_train = [self.sample_data[self.sample_name[i]][1] for i in train_samples_pos]
        y_tuning = [self.sample_data[self.sample_name[i]][1] for i in tune_samples_pos]
        y_test = [self.sample_data[self.sample_name[i]][1] for i in test_samples_pos]


        return X_train,y_train,X_tuning,y_tuning,X_test,y_test

    ## データをランダムにとってくる。
    def data_pickup_sim(self, train_samples_pos, tune_samples_pos, test_samples_pos, probe_num):

        probes_index = np.random.choice(len(self.probes_data), probe_num, replace=False)

        # print("num probes", len(self.probes_data))
        # print(probes_index)

        X_train_tmp = []
        X_tuning_tmp = []
        X_test_tmp = []
        ## data_matrixのサブセットを作る
        for k in self.probes_data.keys():
            X_train_tmp.append([self.probes_data[k][i] for i in train_samples_pos])
            X_tuning_tmp.append([self.probes_data[k][i] for i in tune_samples_pos])
            X_test_tmp.append([self.probes_data[k][i] for i in test_samples_pos])

        X_train = np.array(X_train_tmp).transpose()[:, probes_index]
        X_tuning = np.array(X_tuning_tmp).transpose()[:, probes_index]
        X_test = np.array(X_test_tmp).transpose()[:, probes_index]

        y_train = [self.sample_data[self.sample_name[i]][1] for i in train_samples_pos]
        y_tuning = [self.sample_data[self.sample_name[i]][1] for i in tune_samples_pos]
        y_test = [self.sample_data[self.sample_name[i]][1] for i in test_samples_pos]


        # print(X_test)
        # print(X_test.shape)
        # np.savetxt('test.csv', X_test, delimiter=',')

        return X_train,y_train,X_tuning,y_tuning,X_test,y_test

    ## データをランダムにとってくる。
    def data_pickup_sim_selper(self, train_samples_pos, tune_samples_pos, test_samples_pos, probe_num):

        probes_index = np.random.choice(len(self.probes_data), probe_num, replace=False)

        # print("num probes", len(self.probes_data))
        # print(probes_index)

        X_train_tmp = []
        X_tuning_tmp = []
        X_test_tmp = []
        ## data_matrixのサブセットを作る
        for k in self.probes_data.keys():
            X_train_tmp.append([self.probes_data[k][i] for i in train_samples_pos])
            X_tuning_tmp.append([self.probes_data[k][i] for i in tune_samples_pos])
            X_test_tmp.append([self.probes_data[k][i] for i in test_samples_pos])

        X_train = np.array(X_train_tmp).transpose()[:, probes_index]
        X_tuning = np.array(X_tuning_tmp).transpose()[:, probes_index]
        X_test = np.array(X_test_tmp).transpose()[:, probes_index]

        y_train = [self.sample_data[self.sample_name[i]][1] for i in train_samples_pos]
        y_tuning = [self.sample_data[self.sample_name[i]][1] for i in tune_samples_pos]
        y_test = [self.sample_data[self.sample_name[i]][1] for i in test_samples_pos]


        # print(X_test)
        # print(X_test.shape)
        # np.savetxt('test.csv', X_test, delimiter=',')

        return X_train,y_train,X_tuning,y_tuning,X_test,y_test


    def print_condition(self):
        print ("Design_Type", self.Design_Type)
        print ("chromosome", self.chromosome)
        print ("position_gene", self.position_gene)
        print ("CG_island", self.CG_island)
