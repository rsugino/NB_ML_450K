#from pylab import *
import numpy as np
import time
import os
import shutil

import random
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import f1_score
from sklearn.metrics import log_loss

from sklearn.ensemble import RandomForestClassifier

from sklearn.tree import export_graphviz
import subprocess

from multiprocessing import Pool
from multiprocessing import Process
import settings
import methylome_data

import copy

## This script separate train-test dataset by own script without using build-in classifier.
## This script prepare same number train samples from each dataset, each sampleclass.

class run_random_forest:
	def __init__(self, tmp):
		settings, probe_num = tmp
		start = time.time()

		self.settings = settings
		self.data_container = self.settings.data_container

		self.settings.result_dir = "result/RF/"
		self.initialize_data()

		self.out_dir = settings.result_dir+self.data_name
		if os.path.isdir(self.out_dir) is False:
			os.makedirs(self.out_dir)

		print(self.settings.usecpu)
		shutil.copyfile(self.settings.setting_file, self.settings.result_dir+self.data_name+"/settings.txt")

		# best_result = []
		# importances = []
		# scores = []
		#f1= open(settings.result_dir+self.data_name+"/false_list.tsv","w")
		#f2=open(settings.result_dir+self.data_name+"/true_rate_list.tsv","w")

		# self.arr = np.array([self.data_container.probes_data[k] for k in self.data_container.probes_data])

		print("simulation start: ",probe_num)

		data_tmp = self.data_name
		self.scores = []
		self.data_name = data_tmp+str(probe_num)
		if os.path.isdir(self.settings.result_dir+self.data_name) is False:
			os.makedirs(self.settings.result_dir+self.data_name)
		for self.rep in range(0,100):
			#X_train, X_test, y_train, y_test = train_test_split(self.arr.transpose(),  self.data_container.sample_category)
			print (self.rep,"th itaration", probe_num)

			## それぞれに割り振ったサンプルのlistの場所を取ってくる。
			train_samples_pos,tune_samples_pos,test_samples_pos = self.data_container.split3()
			## とってきた場所をもとに、データとカテゴリーをとってくる
			# X_train,y_train,X_tuning,y_tuning,X_test,y_test = self.data_container.data_pickup(train_samples_pos,tune_samples_pos,test_samples_pos)
			X_train,y_train,X_tuning,y_tuning,X_test,y_test = self.data_container.data_pickup_sim(train_samples_pos,tune_samples_pos,test_samples_pos,probe_num)

			# X_train, X_test, y_train, y_test = train_test_split(self.arr.transpose(),  self.data_container.sample_category)

			##### Random forest: start training #####
			a = RandomForestClassifier(n_estimators=10000, n_jobs=-1, class_weight='balanced')
			# a = RandomForestClassifier(n_estimators=1000, n_jobs=-1, max_depth=4, oob_score=True, class_weight='balanced')
			a.fit(X_train, y_train)
			train_score = a.score(X_train, y_train)

			if self.rep == 0:
				f = open(self.settings.result_dir+self.data_name+"/params.txt","w")
				f.write(str(a.get_params()))
				f.close()

			y_tuning_name = [self.data_container.sample_data[self.data_container.sample_name[i]][0] for i in tune_samples_pos]
			threshold = self.determine_threshold(a, X_tuning, y_tuning, y_tuning_name)

			y_test_name = [self.data_container.sample_data[self.data_container.sample_name[i]][0] for i in test_samples_pos]
			self.output_test(a, X_test, y_test, y_test_name, threshold)



	# predict_probにつけるthresholdをきめる。
	def determine_threshold(self, a, X_tuning, y_tuning, y_tuning_name):

		threshold = 0

		if len(y_tuning) == 0:
			return 0

		pred_list = a.predict(X_tuning) #予測した結果
		pred_proba = a.predict_proba(X_tuning) #各分類群に対しての確率

		max_probs = np.array([pred_proba[i].max() for i in range(len(pred_proba))])
		print(max_probs.argsort())

		f = open(self.settings.result_dir+self.data_name+"/pred_prob.tsv","w")
		for j in max_probs.argsort():
			f.write(y_tuning_name[j]+"\t"+y_tuning[j]+"\t"+pred_list[j]+"\t"+str(pred_proba[j].max())+"\t"+str(pred_proba[j])+"\n")
		f.close()



		i = 0
		for i in range(len(X_tuning)):
			print(y_tuning[i],pred_list[i],max_probs.argsort()[i],pred_proba[i].max(),pred_proba[i])

		f = open(self.settings.result_dir+self.data_name+"/evaluate.tsv","w")
		f.write("num\ttest_samples\ttest_samples_group\tlowerst_prob\taccuracy\tf1\tf1s\tlogloss\n")
		for i in range(len(max_probs)):
			tmp = max_probs.argsort()[i:]
			f.write(str(i)+"\t"+str(len(y_tuning)-i))

			labelname= list(self.data_container.sample_loc_bycategory.keys())

			label_out = ""
			for l in range(len(labelname)):
				if l != 0:
					label_out += ","
				label_out += labelname[l]+":"+str(len([j for j in tmp if y_tuning[j] == labelname[l]]))
			f.write("\t"+label_out)

			f.write("\t"+str(max_probs[tmp[0]]))

			score = accuracy_score([y_tuning[j] for j in tmp], [pred_list[j] for j in tmp])
			f.write("\t"+str(score))
			if score < 0.9:
				threshold = max_probs[tmp[0]]

			score = f1_score([y_tuning[j] for j in tmp], [pred_list[j] for j in tmp], average="macro")
			f.write("\t"+str(score))

			score = f1_score([y_tuning[j] for j in tmp], [pred_list[j] for j in tmp], average=None)
			f.write("\t"+str(score))

			# for k in labelname:
			# 	score = f1_score([y_tuning[j] for j in tmp], [pred_list[j] for j in tmp], labels=labelname)
			# 	f.write("\t"+str(score))

			check_sample_num = 0
			for k in labelname:
				if len([j for j in tmp if y_tuning[j] == k]) == 0:
					check_sample_num += 1
			if check_sample_num != 0:
				score = "-"
			else:
				score = log_loss([y_tuning[j] for j in tmp], np.array([list(pred_proba[j]) for j in tmp]))
			f.write("\t"+str(score))

			f.write("\n")

		return threshold


	def output_test(self, a, X_test, y_test, y_test_name, threshold):

		self.label_num = {"MYCNA":5,"4":5,"4s":5,"other":5}
		self.label_true = {"MYCNA":5,"4":5,"4s":5,"other":5}
		for k in list(set(self.data_container.sample_category)):
			self.label_num[k] = 0
			self.label_true[k] = 0

		print("threshold: ", threshold)
		tmp_score = a.score(X_test, y_test)
		pred_list = a.predict(X_test)
		pred_proba = a.predict_proba(X_test)
		print("score: ", tmp_score)

		print("Confusion matrix:\n{}".format(confusion_matrix(y_test, pred_list)))

		max_probs = np.array([pred_proba[i].max() for i in range(len(pred_proba))])

		f = open(self.settings.result_dir+self.data_name+"/pred_prob_test.tsv","w")
		for j in max_probs.argsort():
			f.write(y_test_name[j]+"\t"+y_test[j]+"\t"+pred_list[j]+"\t"+str(pred_proba[j].max())+"\t"+str(pred_proba[j])+"\n")
			if pred_proba[j].max() < threshold:
				continue
			self.label_num[y_test[j]] += 1
			self.data_container.sample_tested[y_test_name[j]] += 1
			if y_test[j] == pred_list[j]:
				self.label_true[y_test[j]] += 1
				self.data_container.sample_correct[y_test_name[j]] += 1
		f.close()

		true_num = 0
		total_num = 0
		for k in self.label_num:
			print(k,self.label_true[k],self.label_num[k])
			true_num += self.label_true[k]
			total_num += self.label_num[k]

		self.scores.append(true_num/total_num)
		# print("score: ", str(true_num/total_num), " for ", self.rep)
		#
		print(classification_report(y_test,pred_list))
		# print(y_test)
		#
		# print(precision_recall_fscore_support(y_test, pred_list, average=None))
		prec_tmp,reca_tmp,f1_tmp,support_tmp = precision_recall_fscore_support(y_test, pred_list, average=None, labels=list(self.label_num.keys()))
		prec_cat = {}
		reca_cat = {}
		f1_cat = {}
		support_cat = {}


		print("label_num", self.label_num)
		i = 0
		for k in self.label_num:
			prec_cat[k] = prec_tmp[i]
			reca_cat[k] = reca_tmp[i]
			f1_cat[k] = f1_tmp[i]
			support_cat[k] = support_tmp[i]
			i += 1
		# print("prec_cat", prec_cat)
		# print("reca_cat", reca_cat)
		# print("f1_cat", f1_cat)
		# print("support_cat", support_cat)

		if self.rep == 0:
			f = open(self.settings.result_dir+self.data_name+"/scores.tsv","w")
			f.write("rep\tmean_score")
			for k in self.label_num:
				f.write("\t"+str(k)+"_score")
			for k in self.label_num:
				f.write("\t"+str(k)+"_precision")
			for k in self.label_num:
				f.write("\t"+str(k)+"_recall")
			for k in self.label_num:
				f.write("\t"+str(k)+"_f1")
			for k in self.label_num:
				f.write("\t"+str(k)+"_support")
			f.write("\n")
			f.close()

		f = open(self.settings.result_dir+self.data_name+"/scores.tsv","a")
		f.write(str(self.rep)+"\t"+str(true_num/total_num))
		for k in self.label_num:
			f.write("\t"+str(self.label_true[k]/self.label_num[k]))
		for k in self.label_num:
			f.write("\t"+str(prec_cat[k]))
		for k in self.label_num:
			f.write("\t"+str(reca_cat[k]))
		for k in self.label_num:
			f.write("\t"+str(f1_cat[k]))
		for k in self.label_num:
			f.write("\t"+str(support_cat[k]))
		f.write("\n")
		f.close()

	## draw tree
	def draw_tree(self):
		print(self.arr.shape,len(self.data_container.sample_category))
		X_train, X_test, y_train, y_test = train_test_split(self.arr.transpose(),  self.data_container.sample_category)
		forest = RandomForestClassifier(n_estimators=self.num_tree, max_depth=self.depth,n_jobs=-1, oob_score=True)
		# forest = RandomForestClassifier(n_estimators=self.num_tree, max_depth=self.depth,n_jobs=-1, oob_score=True, class_weight = "balanced")
		forest.fit(X_train, y_train)

		for i in range(1,3):
			filename = self.tree_dir+str(self.num_tree)+"_"+str(self.depth)+"_v"+str(i)+".dot"
			# print(filename)
			# print("forest.estimators_[0]",i,forest.estimators_)
			export_graphviz(forest.estimators_[0], out_file=filename,  class_names=["4Amp", "4Noamp", "4s", "other"], feature_names=self.probename, impurity=False, filled=True)
#			export_graphviz(forest.estimators_[0], out_file=filename,  class_names=["4Amp", "4Noamp", "4s", "other"], feature_names=self.probename, impurity=False, filled=True)
			command = "dot -T pdf "+filename+" -o "+filename.replace(str(i)+".dot",str(i)+".pdf")
			subprocess.call(command, shell=True)


	## Initialize variables
	def initialize_data(self):
		self.betalist = []
		self.data_name = self.settings.merged_data.split("/")[-1].split(".")[0]
		# self.data_name = self.settings.merged_data.split("/")[-1].replace(".tsv","")
		# self.data_name = self.settings.merged_data.split("/")[-1].replace(".txt","")
		# self.data_name = self.settings.merged_data.split("/")[-1].replace(".csv","")

		# random_forest_dir -> result_dir
		# if os.path.isdir(self.settings.result_dir+self.data_name) is False:
		# 	os.makedirs(self.settings.result_dir+self.data_name)

		group_num = {}
		for k in self.data_container.sample_category:
			if k not in group_num:
				group_num[k] = 1
			else:
				group_num[k] += 1

		sample_num = len(self.data_container.sample_category)
		# print(group_num)
		# self.base_Gini = cal_Gini(group_num)

		tmp_key = ""
		self.probename = self.data_container.probes
		args_set = []
		for k in self.data_container.probes_data.keys():
			if tmp_key == "":
				tmp_key = k
			self.betalist.append(self.data_container.probes_data[k])


			# if self.args_set == []:
			# args_set.append([k,self.base_Gini,group_num,self.data_container.sample_name,self.data_container.probes_data[k],self.data_container.sample_category])

		self.make_category()

	#	self.select_top_probe()
		print(group_num)

		self.arr = np.array(self.betalist)
		# self.arr = np.array(self.selected)

	## 上位いくつかのprobeのみで解析を行う。
	def select_top_probe(self):
		p = Pool(30)
		self.sort_by_Gini = p.map(Gini_decrement, args_set)
		p.close()

		self.sort_by_Gini.sort(key = lambda x:x[2])

		top_Gini = []
		f = open(self.settings.result_dir+self.data_name+"Gini.tsv","w")
		i = 0
		while i < len(self.sort_by_Gini):
			k = self.sort_by_Gini[i]
		# for k in self.sort_by_Gini:
			f.write(k[0]+"\t"+str(k[1])+"\t"+str(k[2])+"\n")
			if i < 1000:
				top_Gini.append(k[0])
			i += 1
		f.close()

		## top 1000だけを用いる。
		self.selected = []
		for k in top_Gini:
			if tmp_key == "":
				tmp_key = k
			self.selected.append(self.data_container.probes_data[k])

	def make_category(self):
		for k in self.data_container.sample_mark.keys():
			self.data_container.sample_mark[k][1] = -0

		for k in self.data_container.sample_info:
			if k[1] not in self.data_container.sample_mark.keys():
				self.data_container.sample_mark[k[1]] = [k[1],1,k[2]]
				continue

			self.data_container.sample_mark[k[1]][1] += 1

	def plot_feature_importances_cancer(self, model):
		# n_features = cancer.data.shape[1]
		# plt.barh(range(n_features), model.feature_importances_, align='center')
		n_features = self.arr.transpose().shape[1]
		plt.barh(range(n_features), model.feature_importances_, align='center')
		plt.yticks(np.arange(n_features), self.probename)
		plt.xlabel("Feature importance")
		plt.ylabel("Feature")
		plt.ylim(-1, n_features)
