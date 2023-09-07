
#!/usr/bin/env python3

# import libraries 
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import datetime
import xgboost
import sklearn
import sklearn.ensemble as ensemble
#from imblearn.ensemble import BalancedRandomForestClassifier


from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from xgboost import XGBRFClassifier

from sklearn.feature_selection import RFE
from platform import python_version

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--model", "-m", type=str, required=True, choices=["TME_included", "TME_excluded", "Simple"], help="Model selection (TME_included, TME_excluded, Simple)")

# Warnings for input files
if not args.model:
    warnings.warn("Model not provided. Please specify the --model or -m argument.")

#aign parsed values 
model = args.model

  
print('The scikit-learn version is {}.'.format(sklearn.__version__))

cols_to_include = ['cohort','Patient','HLA_allele','Mut_peptide','response','Partition',
                   'Aro','mw','pI', 'Inst', 'CysRed','RankEL','RankBA','NetMHCExp',
                        'Expression','SelfSim','Prime','PropHydroAro','HydroCore',
                        'PropSmall','PropAro','PropBasic','PropAcidic','DAI','Stability','Foreigness',
                        'CelPrev','PrioScore','CYT','HLAexp','Monocytes']



if model == "TME_included": 
    cols_to_drop = ['cohort','Patient','HLA_allele','Mut_peptide','response','Partition']
    Model = "TME_included"
if model == "TME_excluded": 
    cols_to_drop = ['cohort','Patient','HLA_allele','Mut_peptide','response','Partition',
                   'CYT','HLAexp','Monocytes']
    Model = "TME_excluded"
if model == "Simple": 
    cols_to_drop = ['cohort','Patient','HLA_allele','Mut_peptide','response','Partition',
                   'CelPrev','PrioScore','CYT','HLAexp','Monocytes']
    Model = "Simple"

# load in data 
data = pd.read_csv('../../data/03_data_for_CV/03_final_peptide_features_Partition.txt', sep='\t')

Var_importance_filename = '../../results/5_fold_CV/'+Model+'/Feature_importance_' + Model +'.txt'



# run 5 CV random forrest 


outfile = open(Var_importance_filename ,'w')
n_sample = 50
partitioning_list = data.Partition.unique()
#partitioning_list.sort()
print(partitioning_list)
pred_df = pd.DataFrame()
for i in partitioning_list:
    av_pred_rf = 0
    print(i)
    for n in range(n_sample):  
        test = data[data.Partition == i]
        train = data[data.Partition != i]
        info = test[['cohort',"Mut_peptide", "HLA_allele","Patient","Partition"]]
        train_pos = train[train.response==1]
        train_neg = train[train.response==0].sample(n=500, random_state=n) 
        train = shuffle(pd.concat([train_pos, train_neg], axis=0), random_state = 42).reset_index(drop=True)
        train = train[~train.Mut_peptide.isin(test["Mut_peptide"])]
        X_train = train.drop(cols_to_drop, axis=1).reset_index(drop=True)
        
        y_train = train["response"]
        feature_list = list(X_train.columns)
        


        # Instantiate the RF and the MLP # 
        rf = RandomForestClassifier(random_state = 42, max_depth = 6, min_samples_leaf = 6,
                                    n_estimators = 2000, n_jobs=-1) #class_weight={0:1,1:3}

                # Train the models on training data
        rf.fit(X_train, y_train)
        
        # pickle dumb model 
        pickle_model = '../../models/'+ Model +'/rf' + str(i) + '_' + Model +'_' + str(n) + '.pkl'
        pickle.dump(rf, open(pickle_model, 'wb'))
        
        X_test = test.drop(cols_to_drop, axis=1).reset_index(drop=True)
        y_test = test["response"]  
        len(X_test)

        #Get fetures importance per partition
        importances = list(rf.feature_importances_)
        feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]
        # Sort
        feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
     #   [print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances]
      #  gather all feature importance 
        for pair in feature_importances:
            outfile.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(i) + "\n")

        prediction_rf = rf.predict_proba(X_test)

        
        av_pred_rf += prediction_rf
    
    pred_val_rf = av_pred_rf/n_sample
        
    pred_val_rf = pd.DataFrame(pred_val_rf[:,1])
    pred_val_rf.reset_index(drop=True, inplace=True)
    mapping = {pred_val_rf.columns[0]: "prediction_rf"}
    pred_val_rf = pred_val_rf.rename(columns=mapping)

    X_test.reset_index(drop=True, inplace=True)
    y_test.reset_index(drop=True, inplace=True)
    info.reset_index(drop=True, inplace=True)
    print(len(X_test))
    ped_df = pd.concat([info,X_test,y_test,pred_val_rf],axis = 1)
    ped_df["Partition"] = i
    pred_df = pred_df.append(ped_df)
    print(X_train.columns)

        
#print("RF AUC:",round(roc_auc_score(pred_df.response, pred_df.prediction_rf),4))    
print("n sample:",n_sample,"RF AUC:",round(roc_auc_score(pred_df.response, pred_df.prediction_rf),4))  
outfile.close()  


# save output to file 
pred_df.to_csv(r'../../results/5_fold_CV/'+Model+'/pred_df_'+ Model +'.txt', index=None, sep=' ', mode='w')


