# -*- coding: utf-8 -*-
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import font_manager
import re
import pylab
# import Image, ImageDraw
from PIL import Image

##
# @package infinium
# @brief handling infinium data of annotation and data.

## @brief 図の設定情報をもったクラス
class conditions:
    """UCSC refseq annotation"""
    def __init__(self):
        self.focal_gene = ""
        self.Infinium_Design_Type = "" ##< I or II
        self.chromosome = [] ##< always character not number
        self.position_gene = [] ##< TSS1500, TSS200, 5'UTR, 3'UTR, 1stExon, Body
        self.CG_island = [] ##< N_Shelf, N_Shore, Island, S_Shore, S_Shelf. If set as null character "", it returns sea.
        # self.position_gene = ["TSS1500", "TSS200"] ##< TSS1500, TSS200, 5'UTR, 3'UTR, 1stExon, Body
        # self.CG_island = ["Island"] ##< N_Shelf, N_Shore, Island, S_Shore, S_Shelf

        self.infinium_probes = [] ##< 領域に含まれるinfinium probeのリスト. Ordered.
        self.infinium_probes_anno = {}
        self.infinium_probes_data = {}
        self.sample_name = []
        self.sample_category = []


##
# @brief infiniumのアノテーションデータを取り込む。450K。
# @details アノテーションデータはNCBIからダウンロードしてきたものを使う（/ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv）。
class infinium_annotation:
    """UCSC refseq annotation"""
    def __init__(self, line):
        itemList = line[:-1].split(',')

        if len(itemList) < 12:
            self.Name = "false"
            return
        if itemList[10] == "":
            self.Name = "false"
            return
        if itemList[0] == "IlmnID":
            self.Name = "false"
            return

        self.IlmnID = itemList[0]
        self.Name = itemList[1]
        self.AddressA_ID = itemList[2]
        self.AlleleA_ProbeSeq = itemList[3]
        self.AddressB_ID = itemList[4]
        self.AlleleB_ProbeSeq = itemList[5]
        self.Infinium_Design_Type = itemList[6]
        self.Next_Base = itemList[7]
        self.Color_Channel = itemList[8]
        self.Forward_Sequence = itemList[9]
        self.Genome_Build = itemList[10]
        self.CHR = itemList[11]
        self.MAPINFO = int(itemList[12])
        self.SourceSeq = itemList[13]
        self.Chromosome_36 = itemList[14]
        self.Coordinate_36 = itemList[15]
        self.Strand = itemList[16]
        self.Probe_SNPs = itemList[17]
        self.Probe_SNPs_10 = itemList[18]
        self.Random_Loci = itemList[19]
        self.Methyl27_Loci = itemList[20]
        self.UCSC_RefGene_Name = [ item for item in itemList[21].split(';') ]
        self.UCSC_RefGene_Accession = [ item for item in itemList[22].split(';') ]
        self.UCSC_RefGene_Group = [ item for item in itemList[23].split(';') ]
        self.UCSC_CpG_Islands_Name = itemList[24]
        self.Relation_to_UCSC_CpG_Island = itemList[25]
        self.Phantom = itemList[26]
        self.DMR = itemList[27]
        self.Enhancer = itemList[28]
        self.HMM_Island = itemList[29]
        self.Regulatory_Feature_Name = itemList[30]
        self.Regulatory_Feature_Group = itemList[31]
        self.DHS = itemList[32]

# @brief infiniumのアノテーションデータを取り込む。EPIC。
# @details アノテーションデータはNCBIからダウンロードしてきたものを使う（/ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv）。
class infinium_annotation_EPIC:
    """UCSC refseq annotation"""
    def __init__(self, line):
        itemList = line[:-1].split(',')

        if len(itemList) < 12:
            self.Name = "false"
            return
        if itemList[10] == "":
            self.Name = "false"
            return
        if itemList[0] == "IlmnID":
            self.Name = "false"
            return

        self.IlmnID = itemList[0]
        self.Name = itemList[1]
        self.AddressA_ID = itemList[2]
        self.AlleleA_ProbeSeq = itemList[3]
        self.AddressB_ID = itemList[4]
        self.AlleleB_ProbeSeq = itemList[5]
        self.Infinium_Design_Type = itemList[6]
        self.Next_Base = itemList[7]
        self.Color_Channel = itemList[8]
        self.Forward_Sequence = itemList[9]
        self.Genome_Build = itemList[10]
        self.CHR = itemList[11]
        self.MAPINFO = int(itemList[12])
        self.SourceSeq = itemList[13]
        self.Strand = itemList[14]
        self.UCSC_RefGene_Name = [ item for item in itemList[15].split(';') ]
        self.UCSC_RefGene_Accession = [ item for item in itemList[16].split(';') ]
        self.UCSC_RefGene_Group = [ item for item in itemList[17].split(';') ]
        self.UCSC_CpG_Islands_Name = itemList[18]
        self.Relation_to_UCSC_CpG_Island = itemList[19]
        self.Phantom4_Enhancers = itemList[20]
        self.Phantom5_Enhancers = itemList[21]
        self.DMR = itemList[22]
        self.k_Enhancer = itemList[23] #450k_enhancer
        self.HMM_Island = itemList[24]
        self.Regulatory_Feature_Name = itemList[25]
        self.Regulatory_Feature_Group = itemList[26]
        self.GencodeBasicV12_NAME = [ item for item in itemList[27].split(';') ]
        self.GencodeBasicV12_Accession = [ item for item in itemList[28].split(';') ]
        self.GencodeBasicV12_Group = [ item for item in itemList[29].split(';') ]
        self.GencodeCompV12_NAME = [ item for item in itemList[30].split(';') ]
        self.GencodeCompV12_Accession = [ item for item in itemList[31].split(';') ]
        self.GencodeCompV12_Group = [ item for item in itemList[32].split(';') ]
        self.DNase_Hypersensitivity_NAME = itemList[33]
        self.DNase_Hypersensitivity_Evidence_Count = itemList[34]
        self.OpenChromatin_NAME = itemList[35]
        self.OpenChromatin_Evidence_Count = itemList[36]
        self.TFBS_NAME = itemList[37]
        self.TFBS_Evidence_Count = itemList[38]
        self.Methyl27_Loci = itemList[39]
        self.Methyl450_Loci = itemList[40]
        self.Chromosome_36 = itemList[41]
        self.Coordinate_36 = itemList[42]
        self.SNP_ID = itemList[43]
        self.SNP_DISTANCE = itemList[44]
        self.SNP_MinorAlleleFrequency = itemList[45]


##
# \brief ゲノムブラウザにinfiniumの情報を加える
def get_infinium_anno():
    print ("get_infinium_anno:")

    infi_anno = {}
    for line in open('/media/sugino/data_storage/ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv', 'r'):
        itemList = line.replace('"','').split(",")
        if len(itemList) < 12:
            continue

        infi_anno[itemList[0]] = "\t"+itemList[21]+"\t"+itemList[22]+"\t\""+itemList[23].replace("'","_")+"\"\t"+itemList[6]+"\t"+itemList[25]

    return infi_anno


class CpG_class:
    def __init__(self):
        # self.gene_pos = ["TSS1500","TSS200","5_UTR_1stExon","Body_1stExon","5_UTR","Body","3_UTR",""]
        self.gene_pos = ["TSS1500","TSS200","1stExon","5_UTR","Body","3_UTR","","other"]
        self.CGI = ["N_Shelf","N_Shore","Island","S_Shore","S_Shelf",""]

        self.count_class = {}

        # print self.gene_pos
        self.initialize()

    def initialize(self):
        # print self.gene_pos
        for k in self.gene_pos:
            for l in self.CGI:
                # print k,l
                self.count_class[(k,l)] = 0

    def check_class(self, gene_info, CGI_info, ID):
        gene_info_list = gene_info.replace("'","_").split(";")
        uni_list = list(set(gene_info_list))


        # print type(gene_info)
        # print type(CGI_info)
        # print gene_info, CGI_info
        if gene_info == "":
            self.count_class[(gene_info,CGI_info)] += 1
        # elif len(uni_list) == 1 and uni_list[0] != "1stExon":
        elif len(uni_list) == 1:
            # print gene_info, CGI_info
            self.count_class[(gene_info_list[0],CGI_info)] += 1
        elif len(uni_list) == 2 and "1stExon" in uni_list and ("5_UTR" or "Body") in uni_list:
            # print ID, gene_info, CGI_info
            self.count_class[("1stExon",CGI_info)] += 1
        else :
            # print "skipped", gene_info, CGI_info
            # print "skipped", gene_info, CGI_info
            self.count_class[("other",CGI_info)] += 1
