

''' Set the environment  '''
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem,Draw
from rdkit.Chem import MCS
import pandas as pd
import numpy as np
import copy,re
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol,MakeScaffoldGeneric,MurckoScaffoldSmiles,MurckoScaffoldSmilesFromSmiles
from my_toolset.my_utils import get_mol,compute_FP,canonic_smiles,mapper,weight
from my_toolset.drawing_utils import show_mols
import os,sys
import argparse


from IPython.display import display, SVG, display_svg
import seaborn as sns
from matplotlib import pyplot
from pathlib import Path
import glob
from functools import partial
from pandarallel import pandarallel
n_jobs=50
pandarallel.initialize(nb_workers=n_jobs)

import math
import itertools
import matplotlib.pyplot as plt
import plotly.express as px
import json

from utils.common_utils import mapper,csvToExcel,get_mol,compute_FP,canonic_smiles,float_row,heavy_atm_smiles,kekulize_smi,sort_mol_sim_df
from utils.sarm_utils import get_core,fragmentize,create_SARM,create_SARM_MLP
from utils.dual_target_utils import kekulize_smi,remove_dummy,get_RGrps_dummy,get_frag,combine_CombTab_Fround1,replace_nH,add_singleBond_dummy,match_frag,get_activity_info,add_singleBond_dummy,get_full_core,is_number,is_parent,is_child,remove_dummy,real_sonNode,remove_ionizationForm,find_children_single,find_children,find_parents_single,find_parents,if_root,find_root,smiles2shapeSmarts,scaffoldFamily,scaffold_info,get_full_core,get_atomRingInfo,match_core,get_atomPos
from my_toolset.my_utils import canonic_smiles,df_valid

jupyterMode=0


class argNamespace:
    ''' for test in jupyter notebook  '''
    def __init__(self):  
        self.protein = '/mnt/public-bg6/jay.zhang/Codes-public/Vina/Test/8jzx.pdb'  # 'self'引用的是类实例自身  
        self.ligand = '/mnt/public-bg6/jay.zhang/Codes-public/Vina/Test/c5.pdb' 
        self.addH = 1


def main(args):
    ''' change the work directory  '''
    act_csv_file=args.csvFile
    value_cols=args.column
    log=args.log

    dfAcat=pd.read_csv(act_csv_file)
    # print(dfAcat)
    dfAcat=float_row(dfAcat,cols=value_cols)
    for icol in value_cols:
        dfAcat[icol]=dfAcat[icol].apply(lambda x: math.log(x) if log else x)
    df_sele=df_valid(dfAcat,row_smi='smiles')

    # df_sele[value_col]=df_sele[value_col].parallel_apply(is_number)
    # df_sele=pd.DataFrame(df_sele[df_sele[value_col]!=''])

    if 'Cano_SMILES' not in df_sele.columns:
        df_sele['Cano_SMILES']=df_sele['smiles'].parallel_apply(canonic_smiles)
    df_sele = df_sele.drop_duplicates(subset=['Cano_SMILES'])
    
    if args.type=='scaffold':
        df_sele['Scaffold']=[MurckoScaffoldSmiles(ismi) for ismi in df_sele["Cano_SMILES"]]
        act_CPDs=df_sele['Scaffold']
        smi_col='Scaffold'

    if args.type=='smiles':
        act_CPDs=df_sele['Cano_SMILES']
        smi_col='Cano_SMILES'

    df_round1, df_round2=fragmentize(act_CPDs, n_jobs=args.n_jobs, drop_duplicate=False, pos_args={'RR':True, 'nRnR':True})
    print(df_round1)
    print("Creating SARM!")
    df_table_info_left, df_table_info_right, df_table_info_combine=create_SARM(df_round1, df_round2, df_sele, save_folder=args.save_folder, smi_col=smi_col, value_col=value_cols, minimum_count_site1=args.minimumSite1, minimum_count_site2=args.minimumSite2, csv2excel=args.csv2excel, cal_table_stats=True, n_jobs=args.n_jobs)



def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csvFile", help="the path of csv file with smiles",
                        default='')
    parser.add_argument("--type", help="the type of molecule, [smiles, scaffold] default is smiles", default='smiles')
    parser.add_argument("--column",  nargs='+',  help="the column to create SAR table", required=True, default=[])
    parser.add_argument("--log", help="logrithm of the value of the selected column",default=0,type=int)
    parser.add_argument("--minimumSite1", help="the minmum count of fragments for keeping",default=3,type=float)
    parser.add_argument("--minimumSite2", help="the minmum count of fragments for keeping",default=3,type=float)
    parser.add_argument("--n_jobs", help="the number of cores will be used",default=8,type=int)
    parser.add_argument("--save_folder", help="the folder to save the results and SAR tables",default='SAR_Tables')
    parser.add_argument("--csv2excel", help="if transform csv file into excel file", type=int, default=0)

    args = parser.parse_args()
    return args


if jupyterMode:
    args=argNamespace()
    main(args)

# %%

if __name__ == "__main__":
    args = get_parser()
    main(args)
