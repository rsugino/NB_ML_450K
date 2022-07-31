# -*- coding: utf-8 -*-
## ＠mainpage GEOから落としてきたデータをBMIQでnormalizeする、ためのファイルを作る。ID\tDesign_type\tbeta
# \version 1
# merge infinium data for PCA analysis
# testtest


import time
import os

def main():
    start = time.time()

    params = {}
    params[chromosome] = ["01","10","21"]
    params[Infinium_Design_Type] = [""]
    params[CG_island] = [""]
    params[Phantom5_Enhancers] = ["","+"]
    params[DMR] = [""]
    params[450k_Enhancer] = [""]
    params[DNase_Hypersensitivity_NAME] = [""]
    params[OpenChromatin_NAME] = [""]
    params[TFBS_NAME] = [""]
    params[SNP_ID] = [""]

    for k in params.keys():


        f = open("/media/sugino/HDD2/project/machine_learning/NB_infi_IIBMP/01_dataset/settings/settings_"+str(k)+".tsv", "w")
        f.write("projectdir = /media/sugino/HDD2/project/machine_learning/NB_infi_IIBMP/01_dataset/each_chr/\n\n")

        f.write("usecpu = 30\n")
        f.write("manifest = /media/sugino/data_storage/ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv\n")
        f.write("sampleinfo = /media/sugino/HDD2/project/machine_learning/NB_infi_IIBMP/01_dataset/3dat_set.csv\n\n")

        f.write("datadirloc = /media/sugino/HDD2/project/machine_learning/NB_infi_IIBMP/preprocessIllumina/\n")
        f.write("datafullpath = 0\n")
        f.write("datagz = 1\n")
        f.write("merged_data = /media/sugino/HDD2/project/machine_learning/NB_infi_IIBMP/01_dataset/each_chr/chr"+str(k)+".tsv\n\n")

        f.write("Infinium_Design_Type =\n")
        f.write("chromosome = "+str(k)+"\n")
        f.write("CG_island = \n")

        f.close()


if __name__ == '__main__':
    main()

#plt.show()
