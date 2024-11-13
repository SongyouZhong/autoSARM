import pandas as pd
from .tree_utils import match_frag
from .common_utils import mapper,csvToExcel,get_mol,compute_FP,canonic_smiles,float_row,heavy_atm_smiles,kekulize_smi,sort_mol_sim_df



def combine_CombTab_Fround1(combTable,FRound1):
    df_combTable=pd.read_csv(combTable)
    df_combTable.columns=['Table','SMILES','Size','Count']
    df_FRound1=pd.read_csv(FRound1)
    df_FRound1.columns=['SMILES','Count']
    df_FRound1['Size']=[f'({i})' for i in df_FRound1['Count']]
    
    df_FRound1=df_FRound1[['SMILES','Size']]
    df_combTable=df_combTable[['SMILES','Size']]
    df_comb=pd.concat([df_FRound1,df_combTable])
    return df_comb

def get_activity_info(frag, file_act, cols):
    df_act=pd.read_csv(file_act)
    if 'Cano_SMILES' not in df_act.columns:
        df_act['Cano_SMILES']=df_act['SMILES'].apply(canonic_smiles)
    df_act['matched']=df_act.apply(lambda x: match_frag(x['Cano_SMILES'],ismarts=frag),axis=1)
    df_match=df_act[df_act['matched']==1]
    # print(df_match)
    if len(df_match)==0:
        res=[]
        for icol in cols:
            res+=['','',''] 
            return res
    if len(df_match)==1:
        means=df_match.iloc[0]
        stds=df_match.iloc[0]
        medians=df_match.iloc[0]
    if len(df_match)>1:
        means=df_match.mean().round(2)
        stds=df_match.std().round(2)
        medians=df_match.median().round(2)
    res=[]
    for icol in cols:
        res+=[means[icol],stds[icol],medians[icol]] 
    return res

def relative_score(frag1,frag2):
    frag1_FP=compute_FP(frag1)
    frag2_FP=compute_FP(frag2)
    