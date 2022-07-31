# -*- coding: utf-8 -*-

import time
import os

import settings
import select_probe
import merge_BMIQ


from multiprocessing import Pool
from multiprocessing import Process

##
# @package pick_data
# \version 2
# EPICのアノテーションを使う。
# @brief サンプルの情報を読み込み、machine learning用のdata matrixを作る。
# \details
# 解析に用いるinfiniumのデータを（サンプル、プローブの両方）を選んできて一つのmatrixにまとめる。
# 必要なファイル
# - 抽出情報など、設定をまとめたsetting file。
# - サンプルごとに分割されたinfiniumのデータファイル。beta valueでnomarizeされたものを想定している。tsvやcsv、またはgzなどの形式はファイル名から判断する。一つのディレクトリにまとめておくこと。
# - サンプルの情報を含んだスプレッドシート形式のファイル。一列目にサンプルの情報を含んでいること。この名前でデータを取ってくることになっている。各列には以下の情報が含まれていること。
#   -# サンプルの名前。拡張子を除いたファイル名と同じであること。
#   -# 分類するためのグループ。
#   -# サンプル自体の属性。日本の患者、アメリカの患者など。
#
# プログラムの流れ。
# - 抜き出したいprobeの情報を読み込む。
# - 抜き出したprobeのデータを連結する。\link select_probe \endlink
def main():
    start = time.time()

    ## \var settingdir setting fileをおいている場所。setting fileは複数あっても対応する。
    settingdir = '/media/sugino/HDD2/project/machine_learning/NB_infi_data_310219/settings/'
    setting_files = os.listdir(settingdir)
    print(setting_files)

    for setting_file in setting_files:
    # for setting_file in setting_files:
        if setting_file[-1] == "~" or setting_file[0] == ".": ## ubuntuのファイルシステムのせいか、一時ファイルが作られることがあるので、それらのファイルはスキップする。
            continue

        print ("settings:")
        a = settings.settings(settingdir+setting_file)

        print("output data file", a.merged_data)
        if os.path.isfile(a.merged_data):
            print(a.merged_data," is present")
            continue

        ## select_probe
        print ("select_probe:")
        select_probe.select_probe(a)
        elapsed_time = time.time() - start
        print (("\telapsed_time:{0}".format(elapsed_time)) + "[sec]")

        ## merge_BMIQ
        print ("merge_BMIQ:")
        merge_BMIQ.merge_BMIQ(a)
        elapsed_time = time.time() - start
        print ("\telapsed_time:{0}".format(elapsed_time) + "[sec]")
        # exit()

        a.import_data_matrix()


        elapsed_time = time.time() - start
        print (("\telapsed_time:{0}".format(elapsed_time)) + "[sec]")



if __name__ == '__main__':
    main()

#plt.show()
