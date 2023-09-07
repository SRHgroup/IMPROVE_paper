
#!/usr/bin/env python3

# run script
#python3 bin/python_script/feature_calculations.py --file data/01_data/01_Validated_neoepitopes.txt --dataset "Validated_neoepitopes" --ProgramDir "/Users/annieborch/Documents/programs/" --outfile "data/02_feature_data/02_Validated_neoepitopes_calculated_features.tsv" --TmpDir "/Users/annieborch/Documents/programs/"

# ----------
# Load global modules 
import random
import os
import sys
import importlib
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('bin/python_script/src')

import multimerPatientTools
import partitionTools
import procTools
import physiochemical_properties
import kernelSim

# append directories 

# arg parse 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--file", "-i", type=str, required=True)
parser.add_argument("--outfile", "-o", type=str, required=True)
parser.add_argument("--dataset", "-d", type=str, required=True)
parser.add_argument("--ProgramDir", "-p", type=str, required=True)
parser.add_argument("--TmpDir", "-t", type=str, required=True)
args = parser.parse_args()


#aign parsed values 
infile = args.file
dataset_name = args.dataset
output_file = args.outfile
ProgramDirectory = args.ProgramDir
TemptDirectory = args.TmpDir

homeDir = '.'

pd.set_option("display.max_columns",999)
pd.set_option("display.max_rows",None)
sns.set_context(context='talk',rc={"lines.linewidth": 2})

#df = pd.read_csv(os.path.join(homeDir,'benchmark_data.txt'),sep='\t')
df = pd.read_csv(os.path.join(homeDir,infile),sep='\t')
# down sample for testrun 
#df = df.sample(n = 2000)

#correct if Norm peptide is NA
df['PeptNorm'] = df.apply(lambda row: row['Mut_peptide'] if pd.isna(row['Norm_peptide']) else row['Norm_peptide'], axis=1)


#df.to_csv(os.path.join('.',"df.txt"),sep='\t',index=False)
#print(df.columns)

#def procMultimer2Uniq(df):
df = df.rename(columns={'HLA_allele':'MHC',
                            'Mut_peptide':'PeptMut'})
cols = ['MHC', 'PeptMut','PeptNorm']
dfdup = df[cols]
dfdup = dfdup[dfdup. duplicated(keep=False)]
dfSel = df[cols].drop_duplicates()
dfSel['PeptLen'] = dfSel['PeptMut'].apply(len)
dfSel['HLA'] = dfSel['MHC']
dfSel['Patient'] = 1

#dfSel = procMultimer2Uniq(df)
#print("Hola")
#print(dfSel.columns)
#print(kernelSimOut)
#dfdup.to_csv(os.path.join('.',"dfdup.txt"),sep='\t',index=False)


importlib.reload(multimerPatientTools)
dfMerge = multimerPatientTools.runFeatureGeneration_NetMHCpan_stab_prime_molecular(dfSel,os.path.join(homeDir,'results','predictions'),dataSet=dataset_name,utilsDir = ProgramDirectory,tmpDir = TemptDirectory,runPred=True,plot=False,clean=True)






#dfMerge_with_input = pd.concat([df,dfMerge],axis = 1)


# make nice names 
dfMerge = dfMerge.rename(columns={'MHC':'HLA_allele','PeptMut':'Mut_peptide',
  'aro':'Aro', 'inst':'Inst', 'cys_red':'CysRed','%Rank_EL_mut':'RankEL_4.1','%Rank_EL_wt':'RankEL_wt_4.1','%Rank_BA_mut':'RankBA_4.1',
                        'Expression_Level':'Expression','Self_Similarity':'SelfSim','Score_PRIME':'Prime',
                        'helix':'PropHydroAro','MeanHydroph_coreNoAnc':'HydroCore', 'MeanHydroph':'HydroAll',
                        'Prop_Small':'PropSmall','Prop_Aromatic':'PropAro','Prop_Basic':'PropBasic',
                        'Prop_Acidic':'PropAcidic','Agrotopicity':'DAI_4.1','%Rank_Stab':'Stability'})


cols_to_include =['HLA_allele', 'Mut_peptide', 'PeptNorm', 'PeptLen', 
       'Core', 'Of', 'Gp', 'Gl', 'Ip',
       'Il',  'RankEL_4.1', 'RankBA_4.1', 'RankEL_wt_4.1', 'Stability',
       'Prime', 'DAI_4.1','CoreNonAnchor', 'Loci', 'HydroAll', 'HydroCore',
       'PropSmall', 'PropAro', 'PropBasic', 'PropAcidic', 'SelfSim', 'mw',
       'Aro', 'Inst', 'PropHydroAro', 'CysRed', 'pI']
dfMerge = dfMerge[cols_to_include]

df = df.rename(columns={'MHC':'HLA_allele','PeptMut':'Mut_peptide'})

#make idenifier and merge input and output to keep all columns 
#df['id'] = df['HLA_allele'] + '_' + df['Mut_peptide']+ '_' + df['Norm_peptide']
#dfMerge['id'] = dfMerge['HLA_allele'] + '_' + dfMerge['Mut_peptide']+ '_' + dfMerge['Norm_peptide']

#print(len(df))
#print(len(dfMerge))

#print(dfMerge.columns)
#print(df.columns)
#dfMerge.to_csv(os.path.join('.',"dfMerge.txt"),sep='\t',index=False)

#dfMerge_with_input = df.merge(dfMerge, on='id', how='left')
dfMerge_with_input = df.merge(dfMerge, on=['HLA_allele','Mut_peptide','PeptNorm'], how='left')
#dfMerge_with_input = dfMerge_with_input.drop_duplicates(subset =['Patient','HLA_allele','Mut_peptide','PeptNorm'])
# correct NA values when MutPeptide is NA 
dfMerge_with_input['RankEL_wt_4.1'] = dfMerge_with_input.apply(lambda row: pd.NA if pd.isna(row['Norm_peptide']) else row['RankEL_wt_4.1'], axis=1)
dfMerge_with_input['DAI_4.1'] = dfMerge_with_input.apply(lambda row: pd.NA if pd.isna(row['Norm_peptide']) else row['DAI_4.1'], axis=1)

#fMerge_with_input = dfMerge_with_input.drop(columns=['PeptNorm'])

#dfMerge_with_input = pd.merge(df,dfMerge, on="id")

print(len(dfMerge_with_input))
print(dfMerge_with_input.columns)



# write output file 

dfMerge_with_input.to_csv(os.path.join('.',output_file),sep='\t',index=False)
#dfMerge.to_csv(os.path.join('.','results',output_file),sep='\t',index=False)


