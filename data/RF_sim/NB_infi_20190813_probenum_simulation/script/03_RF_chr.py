# -*- coding: utf-8 -*-
## ＠mainpage GEOから落としてきたデータをBMIQでnormalizeする、ためのファイルを作る。ID\tDesign_type\tbeta
# \version 1
# merge infinium data for PCA analysis

import time
import os
import shutil

import settings
import random_forest_probe_num_sim

from multiprocessing import Pool
from multiprocessing import Process

##
# \brief main
# \details
def main():
    start = time.time()

    # datadir = './settings/'
    datadir = '/media/sugino/HDD2/project/machine_learning/NB_RF/NB_infi_20190311_probenum_simulation/settings/'
    setting_files = os.listdir(datadir)
    setting_files = ['4data_all_20190813.txt']
    # setting_files = ['4data_enha.txt']
    # setting_files = ['4data_SNP.txt']
    print(setting_files)
    args_set = []
    # for setting_file in setting_files:
        # if setting_file[-1] == "~" or setting_file[0] == ".":
        #     continue
        # if setting_file == "settings_event.txt":
        #     continue

    setting_file = setting_files[0]

    print ("settings:", setting_file)
    a = settings.settings(datadir+setting_file)
    a.import_data_matrix()
    print("check done")
    # for probe_num in [50000,39719,2070,16376,2259]:
    # for probe_num in [100]:
    # for probe_num in [50000,100000,150000,200000,250000,300000,350000,400000,36313,26416,54338,34643,138386,78798,95986,23349,95986,275166,14017,58921,18534,369310,58724,78798,7269,58151,452453,39719,2070,163726,2259]:
    for probe_num in [36313,26416,54338,34643,138386,78798,95986,23349,95986,275166,14017,58921,18534,369310,58724,78798,7269,58151,452453,39719,2070,163726,2259]:
        # a.data_container.probe_num = probe_num
        args_set.append([a, probe_num])
        # random_forest_score_sample.run_random_forest(a)


    print ("\tcheck")
    # for k in args_set:

    p = Pool(5)
    p.map(random_forest_probe_num_sim.run_random_forest, args_set)
    p.close()


    elapsed_time = time.time() - start
    print (("\telapsed_time:{0}".format(elapsed_time)) + "[sec]")


if __name__ == '__main__':
    main()

#plt.show()
