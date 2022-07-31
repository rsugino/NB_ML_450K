
import os
import sys
import re

import settings
import infinium

##
# @brief select probe
# @details Probes are selected by some conditions, like chromosome, probe type so on.
class select_probe:
    def __init__(self, settings):

        # settings = infinium.conditions()
        self.choose_probe(settings.data_container)

        f = open(settings.datamatdir+"focal_probe.tsv", 'w')
        j = 0
        for k in settings.data_container.probes:
            f.write(settings.data_container.probes_anno[k][0])
            i = 1
            while i < len(settings.data_container.probes_anno[k]):
                f.write("\t"+settings.data_container.probes_anno[k][i])
                i += 1
            f.write("\n")
            i += 1
            # f.write(settings.data_container.probes_anno[k])
        f.close()
        print("check select_probe")

    ##
    # \brief choose probes
    def choose_probe(self,settings):
        print ("choose_probe: ")

        if os.path.isfile(settings.manifest) is False:
            print ("\t\t", settings.manifest, "is absent.")
            print ("\t\tPlease indicate manifest file. This process is quited.")
            exit()

        # print settings.chromosome
        # exit()

        self.header = []
        for line in open(settings.manifest, 'r'):
            itemList = []

            if line.split(",")[0] == "IlmnID":
                self.header = line[:-1].split(",")
                continue

            probe_anno = infinium.infinium_annotation_EPIC(line)

            #450K以外を除く
            if probe_anno.Name == "false":
                continue
            if probe_anno.Methyl450_Loci != "TRUE":
                continue

            if settings.Design_Type != "":
                if probe_anno.Design_Type != settings.Design_Type:
                    continue
                    print (probe_anno.IlmnID, probe_anno.CHR, settings.chromosome)

            if len(settings.chromosome) > 0:
                check = 0
                for j in settings.chromosome:
                    if probe_anno.CHR == j:
                        check = 1
                        break
                    if re.match("!",j):
                        if probe_anno.CHR == j.replace("!",""):
                            print (j, j.replace("!",""))
                            break
                        else:
                            check = 1
                            break

                if check == 0:
                    continue

            if len(settings.position_gene) > 0:
                query = set(probe_anno.UCSC_RefGene_Group)
                subject = set(settings.position_gene)
                if len(query & subject) == 0:
                    continue

            if len(settings.CG_island) > 0:
                if probe_anno.Relation_to_UCSC_CpG_Island not in settings.CG_island:
                    continue

            # EPIC
            if settings.Phantom5_Enhancers == "+":
                if probe_anno.Phantom5_Enhancers == "":
                    continue
            if settings.DMR == "+":
                if probe_anno.DMR == "":
                    continue
            if settings.k_Enhancer == "+":
                if probe_anno.k_Enhancer == "":
                    continue
            if settings.DNase_Hypersensitivity_NAME == "+":
                if probe_anno.DNase_Hypersensitivity_NAME == "":
                    continue
            if settings.OpenChromatin_NAME == "+":
                if probe_anno.OpenChromatin_NAME == "":
                    continue
            if settings.TFBS_NAME == "+":
                if probe_anno.TFBS_NAME == "":
                    continue
            if settings.SNP_ID == "+":
                if probe_anno.SNP_ID == "":
                    continue


            settings.probes.append(probe_anno.IlmnID)
            settings.probes_anno[probe_anno.IlmnID] = line[:-1].split(",")
