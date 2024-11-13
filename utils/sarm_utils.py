from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import copy,re
from pathlib import Path
from .common_utils import mapper,get_mol,compute_FP,mol_with_atom_index
from functools import partial
from pandarallel import pandarallel

''' Common fragmentation functions '''
def load_csv_df(csv_or_df):
    if isinstance(csv_or_df, str):
        df=pd.read_csv(csv_or_df)
    elif isinstance(csv_or_df, pd.DataFrame):
        df=csv_or_df
    else:
        raise Exception("The file should be String or Dataframe.")
    return df

def frags_count(file,minimum_count=5):
    ''' A function to count the fragments and filter the less frequent ones '''
    df=load_csv_df(file)
    counts=df['Key'].value_counts()
    df_count=pd.DataFrame(counts)
    df_count=df_count[df_count['Key']>minimum_count-1]
    df_count['Count']=df_count['Key']
    df_count['Key']=df_count.index
    df_count.columns=['Key','Count']
    # op_file=file.replace('.csv','')+'_count.csv'
    # df_count.to_csv(op_file, index=None)
    return df_count

def attach_left_right(frag, smi):
    ''' Detect if the fragments attached to right or left. '''
    frag_mol=Chem.MolFromSmarts(frag)
    matchDummy = frag_mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    if len(matchDummy)!=2:
        print(f'frag= {frag}; less than 2 dummy atoms were detected in the fragment!')
        return None
    leftDummy=matchDummy[0][0]
    rightDummy=matchDummy[1][0]
    mol=Chem.MolFromSmiles(smi)
    match = mol.GetSubstructMatches(frag_mol)
    Atoms=mol.GetAtoms()
    leftAtmIdx=match[0][leftDummy]
    rightAtmIdx=match[0][rightDummy]
    if len(Atoms[leftAtmIdx].GetNeighbors())>1:
        return 0
    if len(Atoms[rightAtmIdx].GetNeighbors())>1:
        return 1

def sort_mol_sim(smi_list):
    ''' Sort the compounds based on the similarity! '''
    df_sorted=pd.DataFrame.from_dict({'Smi':[],'Mol':[],'Fp':[]})
    for idx,ismi in enumerate(smi_list):
        mol=get_mol(ismi)
        Fp=compute_FP(mol)
        if len(df_sorted)==0:
            df_sorted.loc[0]=[ismi, mol, Fp]
            continue
        Sim_np=np.array([DataStructs.TanimotoSimilarity(
                Fp,iFp) for iFp in df_sorted['Fp']])
        SimArgmax=np.argmax(Sim_np)
        df_sorted.loc[SimArgmax+0.5]=[ismi, mol, Fp]
        df_sorted=df_sorted.sort_index()
        df_sorted=df_sorted.reset_index(drop=True)
    return pd.DataFrame(df_sorted['Smi'])
    
def sort_mol_sim_df(df, col='OrgSmi'):
    ''' Sort the dataframe based on the similarity of CPDs!  '''
    df['index']=df[col]
    df=df.set_index('index')
    df_sorted=pd.DataFrame.from_dict({'Key2':[],'Value2':[],'Smi':[],'Mol':[],'Fp':[]})
    for idx,row in df.iterrows():
        #for idx,ismi in enumerate(smi_list):
        ismi=row[col]
        key2=row['Key']
        value2=row['Value']
        mol=get_mol(ismi)
        Fp=compute_FP(mol)
        if len(df_sorted)==0:
            df_sorted.loc[0]=[key2,value2,ismi, mol, Fp]
            continue
        Sim_np=np.array([DataStructs.TanimotoSimilarity(
                Fp,iFp) for iFp in df_sorted['Fp']])
        SimArgmax=np.argmax(Sim_np)
        df_sorted.loc[SimArgmax+0.5]=[key2,value2,ismi, mol, Fp]
        df_sorted=df_sorted.sort_index()
        df_sorted=df_sorted.reset_index(drop=True)
    return pd.DataFrame(df_sorted[['Key2','Value2','Smi']]) 

def frag_mol_near_ring(smi, pos_args={'RR':True, 'nRnR':False}):
    '''A function to fragmentize molecules via delete the non-ring bond near the ring'''  
    try: 
        # smi='CC(C1=C(C2=CN3C(C(C)=C2)=NC=N3)NC4=CC=C(C5CCN(C(CCO)CCO)CC5)C=C14)C'
        mol=get_mol(smi)
        if pos_args["nRnR"]:  ## All single bonds except those in ring
            bs=[]
            for bond in mol.GetBonds():
                idx = bond.GetIdx()
                if bond.IsInRing():
                    continue
                btype = bond.GetBondType()
                if btype == Chem.BondType.SINGLE:
                    bs.append(idx)
        else:
            bis = mol.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
            bs = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bis]
            if pos_args['RR']:  #single bond between two ring
                bis = mol.GetSubstructMatches(Chem.MolFromSmarts('[R]!@;-[R]'))
                bs.extend([mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bis])
            
        FragPairs=[]
        for ibond in bs:
            nm = Chem.FragmentOnBonds(mol,[ibond], dummyLabels=[(0, 0)])
            # smis = Chem.MolToSmiles(nm,True)
            Mols=Chem.GetMolFrags(nm,asMols=True)
            dummpy_pair=0
            for imol in Mols:
                if imol.GetSubstructMatches(Chem.MolFromSmarts('[#0][#0]')):
                    dummpy_pair=1
            if dummpy_pair:
                continue
            smi0=Chem.MolToSmiles(Mols[0])
            smi1=Chem.MolToSmiles(Mols[1])
            if Mols[0].GetNumAtoms() < Mols[1].GetNumAtoms() and len(Mols[0].GetSubstructMatches(Chem.MolFromSmarts('[#0]')))<2:
                FragPairs.append([smi1,smi0,smi])
            else:
                FragPairs.append([smi0,smi1,smi])
            # ms = BRICS.BRICSBuild(Mols)
            # print(Chem.MolToSmiles(ms))
        return FragPairs
    except Exception as e:
        print(e)
        return None
    
def fragmentize(act_CPDs, n_jobs=20, drop_duplicate=True, pos_args={'RR':True, 'nRnR':False}):
    ''' A function to fragmentize the active compounds twice '''
    ## first round fragmentation  -->  Key1,Value1
    frag_mol_near_ring_p=partial(frag_mol_near_ring, pos_args=pos_args)
    fragPairs=mapper(n_jobs)(frag_mol_near_ring_p, act_CPDs)
    fragPairs=[ifrag_pair for ifrag_pair in fragPairs if ifrag_pair!=None]
    fragPairs=[jfrag for ifrag in fragPairs for jfrag in ifrag]
    df_round1=pd.DataFrame(fragPairs, columns=['Key','Value','OrgSmi'])
    # df_round1.to_csv('Results/Fragment_round1.csv', index=None)
    ## Second round fragmentation  -->  Key2,Value2; Key1 will be fragmentized in this round.
    if drop_duplicate:
        df_round1=df_round1.drop_duplicates('Key')
    fragPairs=mapper(n_jobs)(frag_mol_near_ring_p, df_round1['Key'])
    fragPairs=[ifrag_pair for ifrag_pair in fragPairs if ifrag_pair!=None]
    fragPairs=[jfrag for ifrag in fragPairs for jfrag in ifrag]
    df_round2=pd.DataFrame(fragPairs, columns=['Key','Value','OrgSmi'])
    return df_round1, df_round2
    # df.to_csv('Results/Fragment_round2.csv', index=None)

def set_isotope(mol, int_isotope):
    ''' reset the dummy atom isotope; Only the first one! '''
    match_list=mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    atoms=mol.GetAtoms()
    dummy_idx=match_list[0][0]  ### first idx of first match
    atoms[dummy_idx].SetIsotope(int_isotope)
    return mol

def get_dummy_negb(atom):
    ''' Get the neighbor index of the dummy atom '''
    negb=atom.GetNeighbors()[0]
    return int(negb.GetIdx())

def connect_R1(R, core='', left_or_right=0, return_type='smiles'):
    ''' Connect single R group to the core
    '''
    # print("R=", R)
    # print("core=", core)
    core_mol=get_mol(core)
    core_mol = copy.deepcopy(core_mol)  ### To protect the inpur molecule object
    matchDummy = core_mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    if len(matchDummy)!=2:
        print(f'frag= {core}; less than 2 dummy atoms were detected in the fragment!')
        return None
    leftDummy=matchDummy[0][0]
    rightDummy=matchDummy[1][0]
    core_mol_atoms=core_mol.GetAtoms()
    core_mol_atoms[leftDummy].SetIsotope(0)
    core_mol_atoms[rightDummy].SetIsotope(1)
    
    r_mol=Chem.MolFromSmiles(R)  ### the isotope of dummy atom is zero
    r_mol=set_isotope(r_mol, 2)
    combo = Chem.CombineMols(core_mol, r_mol)
    match = combo.GetSubstructMatches(Chem.MolFromSmarts('[#0]')) ### detect the dummy atoms
    combo_atoms=combo.GetAtoms()
    dummy_pair=[]
    dummy_negb=[]            ### store the idx of connect dummy atoms
    for imatch in match:     ### look through all the dummy atoms
        atm_idx=imatch[0]
        isotope=combo_atoms[atm_idx].GetIsotope()
        if isotope in [2, left_or_right]:
            dummy_pair.append(atm_idx)
            dummy_negb.append(get_dummy_negb(combo_atoms[atm_idx]))
        else:
            combo_atoms[atm_idx].SetAtomicNum(50)  ## protect this dummy atom
            
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(dummy_negb[0],dummy_negb[1],order=Chem.rdchem.BondType.SINGLE)   
    dummy_pair.sort(reverse=True) 
    for idummy in dummy_pair:
        edcombo.RemoveAtom(idummy)
    combo = edcombo.GetMol()
    ''' Replace dummy atom with hydrogen '''
    products = Chem.ReplaceSubstructs(combo,Chem.MolFromSmarts('[#0]'),Chem.MolFromSmarts('[#1]'),replaceAll=True)
    combo=products[0]
    products = Chem.ReplaceSubstructs(combo,Chem.MolFromSmarts('[#50]'),Chem.MolFromSmarts('[#0]'),replaceAll=True)
    combo=products[0]
    combo_smi=Chem.MolToSmiles(combo)  ### move the hydrogen
    combo=Chem.MolFromSmiles(combo_smi) 
    combo=Chem.RemoveHs(combo)
    if return_type=='mol':
        return combo
    if return_type=='smiles':
        combo_smi=Chem.MolToSmiles(combo)
        return combo_smi

def key_filter(key, ring=True, num_atoms=True, aromatic=False):
    mol=get_mol(key)
    ringinfo = mol.GetRingInfo()
    if mol.GetNumAtoms() < num_atoms:
        return False
    if ring:
        if len(ringinfo.AtomRings())<1:
            return False
        if aromatic:
            for iring in ringinfo.AtomRings():
                aromatic_flag=0
                for iatm in iring:
                    if mol.GetAtomWithIdx(iatm).GetIsAromatic():
                        print(f'Aromatic atom found! {iatm}')
                        aromatic_flag=1
                        break
                if aromatic_flag==1:
                    break
                    # return True
            if aromatic_flag==0:
                    return False
    return True
 
def create_SARM(df_FragsRound1, df_FragsRound2, df_active, save_folder, smi_col="SMILES", value_col="Value", minimum_count=5, csv2excel=True, cal_table_stats=True, n_jobs=50):
    if not csv2excel:
        def csvToExcel(*args1,**args2):
            pass
    else:
        from .common_utils import csvToExcel
    pandarallel.initialize(nb_workers=n_jobs)
    ''' Create Path '''
    root_path=Path(save_folder).absolute()
    root_path.joinpath('Left_Table').mkdir(exist_ok=True, parents=True)
    root_path.joinpath('Right_Table').mkdir(exist_ok=True, parents=True)
    root_path.joinpath('Combine_Table').mkdir(exist_ok=True, parents=True)
    ''' Loading active set '''
    # df_active=pd.read_csv(act_file)
    # df_active=act_file
    df_active=df_active.set_index(smi_col)
    
    ''' Loading round1 fragments '''
    print("Loading round1 fragments")
    # Frag_round1=pd.read_csv(df_FragsRound1)
    Frag_round1=load_csv_df(df_FragsRound1)
    Frag_round1['index']=Frag_round1['Key']
    Frag_round1=Frag_round1.set_index('index')
    Frag_round1_count=frags_count(df_FragsRound1, minimum_count=minimum_count)
    Frag_round1_count.to_csv(root_path.joinpath(f'Frag_round1_count.csv'),index=None)
    csvToExcel(root_path.joinpath(f'Frag_round1_count.csv'), imgCols=['Key'],save_file=root_path.joinpath(f'Frag_round1_count.xlsx'))
    
    ''' Loading round2 fragments '''
    print("Loading round2 fragments")
    Frag_round2=load_csv_df(df_FragsRound2)
    Frag_round2_count=frags_count(df_FragsRound2, minimum_count=minimum_count)
    index_list=sort_mol_sim(Frag_round2_count['Key'])
    Frag_round2_count['index']=Frag_round2_count['Key']
    Frag_round2_count=Frag_round2_count.set_index('index')
    Frag_round2_count=Frag_round2_count.loc[index_list['Smi']]
    Frag_round2_count=Frag_round2_count.reset_index(drop=True)
    Frag_round2_count.to_csv(root_path.joinpath(f'Frag_round2_count.csv'),index=None)
    csvToExcel(root_path.joinpath(f'Frag_round2_count.csv'), imgCols=['Key'],save_file=root_path.joinpath(f'Frag_round2_count.xlsx'))
    
    ''' Build SARMs ''' 
    df_table_info_left=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    df_table_info_right=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    df_table_info_combine=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    
    for idx,row in Frag_round2_count.iterrows():
        print(f'Working on row: {idx}')
        ''' To enumerate the key2s, each will generate a table! '''
        Key2=row['Key']
        if not key_filter(Key2, ring=True, num_atoms=6, aromatic=True):
            continue  ### Filter some smaller, non-ring, non-aromatic keys
        df_key2=Frag_round2[Frag_round2['Key']==Key2]  ## the CPDs with the same Frag B
        df_key2=pd.DataFrame(df_key2)
        # print(f"Stuck after {idx} 1")
        BreakBondPos=[attach_left_right(row['Key'], row['OrgSmi']) for idx,row in df_key2.iterrows()]
        df_key2['BreakBondPos']=BreakBondPos
        # print(f"Stuck after {idx} 2")
        df_sarm_list=[]
        for attach_id in [0,1]:
            df_left=df_key2[df_key2['BreakBondPos']==attach_id] ## attach to the left
            df_key1_sort=sort_mol_sim_df(df_left, col='OrgSmi')
            df_sarm=df_key1_sort.copy()
            df_sarm['index']=df_sarm['Smi']
            df_sarm=df_sarm.set_index('index')
            ####  SARMs table cycle
            def cal_stats(irow):
                # irow=df_key1_sort.loc[iidx]
                df_key1_sel = Frag_round1[Frag_round1['Key']==irow['Smi']]
                for jdx, jrow in df_key1_sel.iterrows():
                    # print(irow['Smi'],str(jrow['Value']),jrow['OrgSmi'],value_col)
                    imax=round(df_active.loc[jrow['OrgSmi'],value_col].max(),1)
                    try:
                        imedian=round(df_active.loc[jrow['OrgSmi'],value_col].median(),1)
                        icount=df_active.loc[jrow['OrgSmi'],value_col].count()
                    except:
                        imedian=''
                        icount=''
                return [irow['Smi'],str(jrow['Value']),f"max:{imax}  median:{imedian}  count:{icount}"]
                
            if cal_table_stats and len(df_key1_sort)>0:                
                table_stats_list=df_key1_sort.parallel_apply(lambda x: cal_stats(x), axis=1)
                for istat in table_stats_list:
                    df_sarm.loc[istat[0],istat[1]]=istat[2]
            
            # for iidx,irow in df_key1_sort.iterrows():
            #     df_key1_sel = Frag_round1[Frag_round1['Key']==irow['Smi']]
            #     for jdx, jrow in df_key1_sel.iterrows():
            #         # print(irow['Smi'],str(jrow['Value']),jrow['OrgSmi'],value_col)
            #         imax=round(df_active.loc[jrow['OrgSmi'],value_col].max(),1)
            #         try:
            #             imedian=round(df_active.loc[jrow['OrgSmi'],value_col].median(),1)
            #             icount=df_active.loc[jrow['OrgSmi'],value_col].count()
            #         except:
            #             imedian=''
            #             icount=''
            #         df_sarm.loc[irow['Smi'],str(jrow['Value'])]=f"max:{imax}  median:{imedian}  count:{icount}"
            # print(list(df_sarm.count(axis=1)))
            
            df_sarm=df_sarm.reset_index(drop=True)
            df_sarm.loc[-0.5]=list(df_sarm.count(axis=0))
            df_sarm=df_sarm.sort_index()
            df_sarm.insert(1,'Count',list(df_sarm.count(axis=1)-3))
            df_tmp = pd.DataFrame({'Table':[f'Table_{idx}_right.csv'],'Key2':[Key2],'Size':[df_sarm.shape],'Items_count':[df_sarm.iloc[0].astype(int)[4:].sum()]}) 
            df_sarm=df_sarm.drop_duplicates('Value2')
            if attach_id==0:
                df_sarm.to_csv(root_path.joinpath(f'Left_Table/Table_{idx}_left.csv'),index=None)
                csvToExcel(root_path.joinpath(f'Left_Table/Table_{idx}_left.csv'), imgCols=['Value2','Smi'],save_file=root_path.joinpath(f'Left_Table/Table_{idx}_left.xlsx'))
                df_table_info_left=df_table_info_left.append(df_tmp)
            if attach_id==1:
                df_sarm.to_csv(root_path.joinpath(f'Right_Table/Table_{idx}_right.csv'),index=None)
                csvToExcel(root_path.joinpath(f'Right_Table/Table_{idx}_right.csv'), imgCols=['Value2','Smi'],save_file=root_path.joinpath(f'Right_Table/Table_{idx}_right.xlsx'))
                df_table_info_right=df_table_info_right.append(df_tmp)
            print(f'Processing Table_{idx}, size: {df_sarm.shape}, count: {df_sarm.size}')
            df_sarm_list.append(df_sarm.copy()) ## store the left and right table
        ###  combine two dataframe
        if df_sarm_list[0].size > df_sarm_list[1].size:
            df_sarm_main=df_sarm_list[0]
            df_sarm_sub=df_sarm_list[1]
            connect_site=0
        else:
            df_sarm_main=df_sarm_list[1]
            df_sarm_sub=df_sarm_list[0]
            connect_site=1
        frag_cols=df_sarm_sub.columns[4:]
        add_frag_cols=[ifrag for ifrag in frag_cols if ifrag not in list(df_sarm_main['Value2'])]
        print("frag_cols:",frag_cols)
        df_tmp = pd.DataFrame({'Value2':add_frag_cols})
        df_sarm_main=df_sarm_main.append(df_tmp)
        df_sarm_main['index']=df_sarm_main['Value2']
        df_sarm_main=df_sarm_main.set_index('index')
        df_sarm_sub['index']=df_sarm_sub['Value2']
        df_sarm_sub=df_sarm_sub.set_index('index')
        for iidx,row in df_sarm_sub.iterrows():
            if isinstance(iidx,str):  ## skip first row
                for icol in frag_cols:
                    iKey2 = df_sarm_sub.loc[row['Value2'],'Key2']
                    df_sarm_main.loc[icol,'Key2'] = iKey2
                    iValue2 = df_sarm_sub.loc[row['Value2'],icol]
                    df_sarm_main.loc[icol,row['Value2']]= iValue2
                    df_sarm_main.loc[icol,'Smi']= connect_R1(R=icol, core=iKey2, left_or_right=connect_site, return_type='smiles')
        df_sarm_main=df_sarm_main.reset_index(drop=True)
        df_sarm_main.loc[0]=list(df_sarm_main.count(axis=0))
        df_sarm_main=df_sarm_main.sort_index()
        # df_sarm_main['Count']=list(df_sarm_main.count(axis=1)-3)
        df_sarm_main=df_sarm_main.drop('Count',axis=1)
        df_sarm_main.insert(1,'Count',list(df_sarm_main.count(axis=1)-3))
        df_sarm_main.to_csv(root_path.joinpath(f'Combine_Table/Table_{idx}_combine.csv'),index=None)
        csvToExcel(root_path.joinpath(f'Combine_Table/Table_{idx}_combine.csv'), imgCols=['Value2','Smi'],save_file=root_path.joinpath(f'Combine_Table/Table_{idx}_combine.xlsx'))
        df_tmp = pd.DataFrame({'Table':[f'Table_{idx}_combine.csv'],'Key2':[Key2],'Size':[df_sarm_main.shape],'Items_count':[df_sarm_main.iloc[0].astype(int)[4:].sum()]})
        df_table_info_combine=df_table_info_combine.append(df_tmp)
        
    df_table_info_left=df_table_info_left.sort_values('Items_count', ascending=False)
    df_table_info_right=df_table_info_right.sort_values('Items_count', ascending=False)
    df_table_info_combine=df_table_info_combine.sort_values('Items_count', ascending=False)
    df_table_info_left.to_csv(root_path.joinpath('Left_Table_info.csv'),index=None)   
    csvToExcel(root_path.joinpath('Left_Table_info.csv'), imgCols=['Key2'],save_file=root_path.joinpath(f'Left_Table_info.xlsx'))
    df_table_info_right.to_csv(root_path.joinpath('Right_Table_info.csv'),index=None) 
    csvToExcel(root_path.joinpath('Right_Table_info.csv'), imgCols=['Key2'],save_file=root_path.joinpath(f'Right_Table_info.xlsx'))
    df_table_info_combine.to_csv(root_path.joinpath('Combine_Table_info.csv'),index=None)
    csvToExcel(root_path.joinpath('Combine_Table_info.csv'), imgCols=['Key2'],save_file=root_path.joinpath(f'Combine_Table_info.xlsx'))
    return df_table_info_left,df_table_info_right,df_table_info_combine


def process_row(row,Frag_round1, Frag_round2,df_active,root_path, value_col="Value",cal_table_stats=True, csv2excel=True):
    # pandarallel.initialize(nb_workers=10)
    if not csv2excel:
        def csvToExcel(*args1,**args2):
            pass
    else:
        from .common_utils import csvToExcel
    # global Frag_round2
    # global df_active
    # global df_table_info_left
    # global df_table_info_right
    # global df_table_info_combine
    df_table_info_left=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    df_table_info_right=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    df_table_info_combine=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    idx=row['index']
    print(f'Working on row: {idx}')
    ''' To enumerate the key2s, each will generate a table! '''
    Key2=row['Key']
    if not key_filter(Key2, ring=True, num_atoms=6, aromatic=True):
        return ''  ### Filter some smaller, non-ring, non-aromatic keys
    df_key2=Frag_round2[Frag_round2['Key']==Key2]  ## the CPDs with the same Frag B
    df_key2=pd.DataFrame(df_key2)
    # print(f"Stuck after {idx} 1")
    BreakBondPos=[attach_left_right(row['Key'], row['OrgSmi']) for idx,row in df_key2.iterrows()]
    df_key2['BreakBondPos']=BreakBondPos
    # print(f"Stuck after {idx} 2")
    df_sarm_list=[]
    for attach_id in [0,1]:
        df_left=df_key2[df_key2['BreakBondPos']==attach_id] ## attach to the left
        df_key1_sort=sort_mol_sim_df(df_left, col='OrgSmi')
        df_sarm=df_key1_sort.copy()
        df_sarm['index']=df_sarm['Smi']
        df_sarm=df_sarm.set_index('index')
        ####  SARMs table cycle
        def cal_stats(irow):
            # irow=df_key1_sort.loc[iidx]
            df_key1_sel = Frag_round1[Frag_round1['Key']==irow['Smi']]
            for jdx, jrow in df_key1_sel.iterrows():
                # print(irow['Smi'],str(jrow['Value']),jrow['OrgSmi'],value_col)
                imax=round(df_active.loc[jrow['OrgSmi'],value_col].max(),1)
                try:
                    imedian=round(df_active.loc[jrow['OrgSmi'],value_col].median(),1)
                    icount=df_active.loc[jrow['OrgSmi'],value_col].count()
                except:
                    imedian=''
                    icount=''
            return [irow['Smi'],str(jrow['Value']),f"max:{imax}  median:{imedian}  count:{icount}"]
            
        if cal_table_stats and len(df_key1_sort)>0:  
            table_stats_list=df_key1_sort.apply(lambda x: cal_stats(x), axis=1)    
            if len(df_key1_sort)>100:        
                table_stats_list=df_key1_sort.parallel_apply(lambda x: cal_stats(x), axis=1)
            for istat in table_stats_list:
                df_sarm.loc[istat[0],istat[1]]=istat[2]
        
        # for iidx,irow in df_key1_sort.iterrows():
        #     df_key1_sel = Frag_round1[Frag_round1['Key']==irow['Smi']]
        #     for jdx, jrow in df_key1_sel.iterrows():
        #         # print(irow['Smi'],str(jrow['Value']),jrow['OrgSmi'],value_col)
        #         imax=round(df_active.loc[jrow['OrgSmi'],value_col].max(),1)
        #         try:
        #             imedian=round(df_active.loc[jrow['OrgSmi'],value_col].median(),1)
        #             icount=df_active.loc[jrow['OrgSmi'],value_col].count()
        #         except:
        #             imedian=''
        #             icount=''
        #         df_sarm.loc[irow['Smi'],str(jrow['Value'])]=f"max:{imax}  median:{imedian}  count:{icount}"
        # print(list(df_sarm.count(axis=1)))
        
        df_sarm=df_sarm.reset_index(drop=True)
        df_sarm.loc[-0.5]=list(df_sarm.count(axis=0))
        df_sarm=df_sarm.sort_index()
        df_sarm.insert(1,'Count',list(df_sarm.count(axis=1)-3))
        df_tmp = pd.DataFrame({'Table':[f'Table_{idx}_right.csv'],'Key2':[Key2],'Size':[df_sarm.shape],'Items_count':[df_sarm.iloc[0].astype(int)[4:].sum()]}) 
        df_sarm=df_sarm.drop_duplicates('Value2')
        if attach_id==0:
            df_sarm.to_csv(root_path.joinpath(f'Left_Table/Table_{idx}_left.csv'),index=None)
            csvToExcel(root_path.joinpath(f'Left_Table/Table_{idx}_left.csv'), imgCols=['Value2','Smi'],save_file=root_path.joinpath(f'Left_Table/Table_{idx}_left.xlsx'))
            df_table_info_left=df_table_info_left.append(df_tmp)
        if attach_id==1:
            df_sarm.to_csv(root_path.joinpath(f'Right_Table/Table_{idx}_right.csv'),index=None)
            csvToExcel(root_path.joinpath(f'Right_Table/Table_{idx}_right.csv'), imgCols=['Value2','Smi'],save_file=root_path.joinpath(f'Right_Table/Table_{idx}_right.xlsx'))
            df_table_info_right=df_table_info_right.append(df_tmp)
        print(f'Processing Table_{idx}, size: {df_sarm.shape}, count: {df_sarm.size}')
        df_sarm_list.append(df_sarm.copy()) ## store the left and right table
    ###  combine two dataframe
    if df_sarm_list[0].size > df_sarm_list[1].size:
        df_sarm_main=df_sarm_list[0]
        df_sarm_sub=df_sarm_list[1]
        connect_site=0
    else:
        df_sarm_main=df_sarm_list[1]
        df_sarm_sub=df_sarm_list[0]
        connect_site=1
    frag_cols=df_sarm_sub.columns[4:]
    add_frag_cols=[ifrag for ifrag in frag_cols if ifrag not in list(df_sarm_main['Value2'])]
    print("frag_cols:",frag_cols)
    df_tmp = pd.DataFrame({'Value2':add_frag_cols})
    df_sarm_main=df_sarm_main.append(df_tmp)
    df_sarm_main['index']=df_sarm_main['Value2']
    df_sarm_main=df_sarm_main.set_index('index')
    df_sarm_sub['index']=df_sarm_sub['Value2']
    df_sarm_sub=df_sarm_sub.set_index('index')
    for iidx,row in df_sarm_sub.iterrows():
        if isinstance(iidx,str):  ## skip first row
            for icol in frag_cols:
                iKey2 = df_sarm_sub.loc[row['Value2'],'Key2']
                df_sarm_main.loc[icol,'Key2'] = iKey2
                iValue2 = df_sarm_sub.loc[row['Value2'],icol]
                df_sarm_main.loc[icol,row['Value2']]= iValue2
                df_sarm_main.loc[icol,'Smi']= connect_R1(R=icol, core=iKey2, left_or_right=connect_site, return_type='smiles')
    df_sarm_main=df_sarm_main.reset_index(drop=True)
    df_sarm_main.loc[0]=list(df_sarm_main.count(axis=0))
    df_sarm_main=df_sarm_main.sort_index()
    # df_sarm_main['Count']=list(df_sarm_main.count(axis=1)-3)
    df_sarm_main=df_sarm_main.drop('Count',axis=1)
    df_sarm_main.insert(1,'Count',list(df_sarm_main.count(axis=1)-3))
    df_sarm_main.to_csv(root_path.joinpath(f'Combine_Table/Table_{idx}_combine.csv'),index=None)
    csvToExcel(root_path.joinpath(f'Combine_Table/Table_{idx}_combine.csv'), imgCols=['Value2','Smi'],save_file=root_path.joinpath(f'Combine_Table/Table_{idx}_combine.xlsx'))
    df_tmp = pd.DataFrame({'Table':[f'Table_{idx}_combine.csv'],'Key2':[Key2],'Size':[df_sarm_main.shape],'Items_count':[df_sarm_main.iloc[0].astype(int)[4:].sum()]})
    df_table_info_combine=df_table_info_combine.append(df_tmp)
    # return 1
    df_tables=[df_table_info_left,df_table_info_right,df_table_info_combine]
    return df_tables

def create_SARM_MLP(df_FragsRound1, df_FragsRound2, df_active, save_folder, smi_col="SMILES", value_col="Value", minimum_count=5, csv2excel=True, cal_table_stats=True, n_jobs=50):
    if not csv2excel:
        def csvToExcel(*args1,**args2):
            pass
    else:
        from .common_utils import csvToExcel
    pandarallel.initialize(nb_workers=n_jobs)
    ''' Create Path '''
    root_path=Path(save_folder).absolute()
    root_path.joinpath('Left_Table').mkdir(exist_ok=True, parents=True)
    root_path.joinpath('Right_Table').mkdir(exist_ok=True, parents=True)
    root_path.joinpath('Combine_Table').mkdir(exist_ok=True, parents=True)
    ''' Loading active set '''
    # df_active=pd.read_csv(act_file)
    # df_active=act_file
    df_active=df_active.set_index(smi_col)
    
    ''' Loading round1 fragments '''
    print("Loading round1 fragments")
    # Frag_round1=pd.read_csv(df_FragsRound1)
    Frag_round1=load_csv_df(df_FragsRound1)
    Frag_round1['index']=Frag_round1['Key']
    Frag_round1=Frag_round1.set_index('index')
    Frag_round1_count=frags_count(df_FragsRound1, minimum_count=minimum_count)
    Frag_round1_count.to_csv(root_path.joinpath(f'Frag_round1_count.csv'),index=None)
    csvToExcel(root_path.joinpath(f'Frag_round1_count.csv'), imgCols=['Key'],save_file=root_path.joinpath(f'Frag_round1_count.xlsx'))
    
    ''' Loading round2 fragments '''
    print("Loading round2 fragments")
    Frag_round2=load_csv_df(df_FragsRound2)
    Frag_round2_count=frags_count(df_FragsRound2, minimum_count=minimum_count)
    index_list=sort_mol_sim(Frag_round2_count['Key'])
    Frag_round2_count['index']=Frag_round2_count['Key']
    Frag_round2_count=Frag_round2_count.set_index('index')
    Frag_round2_count=Frag_round2_count.loc[index_list['Smi']]
    Frag_round2_count=Frag_round2_count.reset_index(drop=True)
    Frag_round2_count.to_csv(root_path.joinpath(f'Frag_round2_count.csv'),index=None)
    csvToExcel(root_path.joinpath(f'Frag_round2_count.csv'), imgCols=['Key'],save_file=root_path.joinpath(f'Frag_round2_count.xlsx'))
    
    ''' Build SARMs ''' 
    print(" Build SARMs ")
    df_table_info_left=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    df_table_info_right=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    df_table_info_combine=pd.DataFrame.from_dict({'Table':[],'Key2':[],'Size':[],'Items_count':[]})
    
    # for idx,row in Frag_round2_count.iterrows():
    Frag_round2_count['index']=[irow for irow in range(len(Frag_round2_count))]
    
        
    process_row_p=partial(process_row,Frag_round1=Frag_round1,Frag_round2=Frag_round2,df_active=df_active,root_path=root_path,value_col=value_col,cal_table_stats=cal_table_stats,csv2excel=csv2excel)
    df_tables=Frag_round2_count.apply(lambda x: process_row_p(x), axis=1)
    # process_row(row,Frag_round1, Frag_round2,df_active,root_path, value_col="Value",cal_table_stats=True, csv2excel=True)
    df_rows=[row.to_dict() for irow,row in Frag_round2_count.iterrows()]
    df_tables=mapper(n_jobs)(process_row_p,df_rows)
    '''  Combine the tables '''
    for itable in df_tables:
        df_table_info_left.append(itable[0])
        df_table_info_right.append(itable[1])
        df_table_info_combine.append(itable[2])
        
    df_table_info_left=df_table_info_left.sort_values('Items_count', ascending=False)
    df_table_info_right=df_table_info_right.sort_values('Items_count', ascending=False)
    df_table_info_combine=df_table_info_combine.sort_values('Items_count', ascending=False)
    df_table_info_left.to_csv(root_path.joinpath('Left_Table_info.csv'),index=None)   
    csvToExcel(root_path.joinpath('Left_Table_info.csv'), imgCols=['Key2'],save_file=root_path.joinpath(f'Left_Table_info.xlsx'))
    df_table_info_right.to_csv(root_path.joinpath('Right_Table_info.csv'),index=None) 
    csvToExcel(root_path.joinpath('Right_Table_info.csv'), imgCols=['Key2'],save_file=root_path.joinpath(f'Right_Table_info.xlsx'))
    df_table_info_combine.to_csv(root_path.joinpath('Combine_Table_info.csv'),index=None)
    csvToExcel(root_path.joinpath('Combine_Table_info.csv'), imgCols=['Key2'],save_file=root_path.joinpath(f'Combine_Table_info.xlsx'))
    return df_table_info_left,df_table_info_right,df_table_info_combine

def detect_bond(frag, smarts=False):
    ''' Find the dummy bonds of the fragment  '''
    if smarts:
        mol=Chem.MolFromSmarts(frag)
    else:
        mol=get_mol(frag)
    mol_atoms=mol.GetAtoms()
    matched=mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    bond_atms=[]
    for imatch in matched:
        iatom=imatch[0]
        ibond_atm=get_dummy_negb(mol_atoms[iatom])
        bond_atms.append([iatom,ibond_atm])
    return bond_atms

def get_core(mol, frags):
    ''' remove fragments and keep the core  '''
    mol=get_mol(mol)
    mol=mol_with_atom_index(mol, Idx=None, start=0)
    # display(mol)
    bond_atms_list=[]  ##  detected bond will be break
    for ifrag in frags:
        # print("ifrag=", ifrag)
        frag_bond_atms=detect_bond(ifrag) 
        matched=mol.GetSubstructMatches(Chem.MolFromSmarts(ifrag))
        # print("matched=", matched)
        for imatch in matched:
            matched_atms=np.array(imatch)
            bond_atms=[list(matched_atms[ibd_atm]) for ibd_atm in frag_bond_atms]
            [bond_atms_list.append(ibd_atm) for ibd_atm in bond_atms if ibd_atm not in bond_atms_list]
    # bond_atms_list=list(set(bond_atms_list))
    bonds=[mol.GetBondBetweenAtoms(int(x),int(y)).GetIdx() for x,y in bond_atms_list]
    dummyLabels=[(0, 0) for ibond in bonds]
    Frag_mol = Chem.FragmentOnBonds(mol, bonds, dummyLabels=dummyLabels)
    Frag_mols=Chem.GetMolFrags(Frag_mol, asMols=True)
    Frag_smis=[Chem.MolToSmiles(imol) for imol in Frag_mols]
    frag_connect_atoms=set([ibond_atm[1] for ibond_atm in bond_atms_list])
    Frag_smis_mem=copy.deepcopy(Frag_smis)
    for ifrag in Frag_smis:
        # print("ifrag=", ifrag)
        re_p=re.compile(r'\(\*\)|\*|\[\*\]')
        ifrag_nodummy = re.sub(re_p, '', ifrag)  ## remove dummy atoms
        imatched=mol.GetSubstructMatches(Chem.MolFromSmarts(ifrag_nodummy))
        # print("imatched=",imatched)
        if len(imatched)>0:
            imatched = set([i for ii in imatched for i in ii])
            if len(frag_connect_atoms & imatched)>0:
                Frag_smis_mem.remove(ifrag)
    if len(Frag_smis_mem)==1:
        return Frag_smis_mem[0]
    return ''

def connect2Frags(R, core='', Rsite=0, return_type='smiles'):
    ''' Connect single R group to the core, only one bond will be formed.
    Note: the isotopes cannot be the same except only two zeros
    '''
    core_mol=get_mol(core)
    core_mol = copy.deepcopy(core_mol)  ## To protect the inpur molecule object
    r_mol=Chem.MolFromSmiles(R)  # the isotope of dummy atom is zero
    combo = Chem.CombineMols(core_mol, r_mol)
    match = combo.GetSubstructMatches(Chem.MolFromSmarts('[#0]')) ## detect the dummy atoms
    combo_atoms=combo.GetAtoms()
    dummy_pair=[]
    dummy_negb=[]  # store the idx of connect dummy atoms
    for imatch in match: # look through all the dummy atoms
        atm_idx=imatch[0]
        isotope=combo_atoms[atm_idx].GetIsotope()
        if isotope in [0, Rsite]:
            dummy_pair.append(atm_idx)
            dummy_negb.append(get_dummy_negb(combo_atoms[atm_idx]))
            
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(dummy_negb[0],dummy_negb[1],order=Chem.rdchem.BondType.SINGLE)   
    dummy_pair.sort(reverse=True) 
    for idummy in dummy_pair:
        edcombo.RemoveAtom(idummy)
    combo = edcombo.GetMol()
    ''' Replace dummy atom with hydrogen '''
    products = Chem.ReplaceSubstructs(combo,Chem.MolFromSmarts('[#0]'),Chem.MolFromSmarts('[#1]'),replaceAll=True)
    combo=products[0]
    combo_smi=Chem.MolToSmiles(combo)  ## To remove the hydrogen
    combo=Chem.MolFromSmiles(combo_smi) 
    combo=Chem.RemoveHs(combo)
    if return_type=='mol':
        return combo
    if return_type=='smiles':
        combo_smi=Chem.MolToSmiles(combo)
        return combo_smi #[combo_smi,Rsite,R]

def frag_pos_double(mol, frags, bonds):
    ''' Detect [Left,Core,Right]fragment position '''
    def bond_pos(matched, bonds):
        left=0;right=0
        for imatch in matched:
            if bonds[0][0] in imatch:
                left=1
                return 0 ## means "Left"
            if bonds[1][0] in imatch:
                right=1  ## means "Right"
                return 1
            
    left_match=0;right_match=0;match_count=0
    frag_ordered=['','','']
    frags_mem=copy.deepcopy(frags)
    for ifrag in frags:
        re_p=re.compile(r'\(\*\)|\*|\[\*\]')
        ifrag_nodummy = re.sub(re_p, '', ifrag)  ## remove dummy atoms
        matched=mol.GetSubstructMatches(Chem.MolFromSmarts(ifrag_nodummy))
        if len(matched)==1:
            if bond_pos(matched, bonds)==0:
                frag_ordered[0]=ifrag
                frags_mem.remove(ifrag)
            if bond_pos(matched, bonds)==1:
                frag_ordered[2]=ifrag
                frags_mem.remove(ifrag)
        if len(matched)==2:
            frag_ordered[0]=ifrag
            frag_ordered[2]=ifrag
            frags_mem.remove(ifrag)
    # print("frags_mem",frags_mem)
    # print("frags_mem0",frags_mem[0])
    if len(frags_mem)==1:
        frag_ordered[1]=frags_mem[0]
    return frag_ordered       

def frag_pos_singe(mol, Frag_smis, bond_atms):
    ''' Detect [Core, Rgrp] fragment position '''
    frag_ordered=['','']
    frags_mem=copy.deepcopy(Frag_smis)
    for ifrag in Frag_smis:
        re_p=re.compile(r'\(\*\)|\*|\[\*\]')
        ifrag_nodummy = re.sub(re_p, '', ifrag)  ## remove dummy atoms
        matched=mol.GetSubstructMatches(Chem.MolFromSmarts(ifrag_nodummy))
        matched = set([i for ii in matched for i in ii])
        if bond_atms[0][0] in matched:
            frags_mem.remove(ifrag)
            frag_ordered[1]=ifrag
    if len(frags_mem)==1:
        frag_ordered[0]=frags_mem[0]
    return frag_ordered  
    
def fragmize_on_bonds(mol, smarts, frag_smi):
    ''' Fragmentize the molecule based on the fragment.  Comparable with 1-2 cut sites.'''
    '''  Analyze the fragment '''
    re_p=re.compile(r'\[([^\*HC]*?)\]')
    smarts_readable = re.sub(re_p, 'c', frag_smi)   ### frag_smi is a readable version of smarts
    frag_bond_atms=detect_bond(smarts_readable) 

    ''' bonded atoms '''
    mol=get_mol(mol)
    mol=mol_with_atom_index(mol, Idx=None, start=0)
    matched=mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    matched_atms=np.array(matched[0])
    bond_atms=[list(matched_atms[ibd_atm]) for ibd_atm in frag_bond_atms]
    bond_atms=list(bond_atms)
    bonds=[mol.GetBondBetweenAtoms(int(x),int(y)).GetIdx() for x,y in bond_atms]
    dummyLabels=[(0, 0) for ibd in bonds]
    Frag_mol = Chem.FragmentOnBonds(mol, bonds, dummyLabels=dummyLabels)
    Frag_mols=Chem.GetMolFrags(Frag_mol, asMols=True)
    Frag_smis=[Chem.MolToSmiles(imol) for imol in Frag_mols]
    # print("Frag_smis=", Frag_smis)
    if len(Frag_smis)==3:
        Frag_smis_ordered=frag_pos_double(mol, Frag_smis, bond_atms)
    elif len(Frag_smis)==2:
        Frag_smis_ordered=frag_pos_singe(mol, Frag_smis, bond_atms)
    else:
        return None
    # print(Frag_smis_ordered)
    return Frag_smis_ordered

def get_complete_frag(mol, frag_smi):
    '''  '''
    '''  Analyze the fragment '''
    mol=get_mol(mol)
    mol=mol_with_atom_index(mol, Idx=None, start=0)
    # display(mol)
    frag_bond_atms=detect_bond(frag_smi) 
    # print("frag_bond_atms=",frag_bond_atms)
    matched=mol.GetSubstructMatches(Chem.MolFromSmarts(frag_smi))
    # print("matched=",matched)
    res_frags=[]
    for imatch in matched:
        matched_atms=np.array(imatch)
        bond_atms=[list(matched_atms[ibd_atm]) for ibd_atm in frag_bond_atms]
        bond_atms=bond_atms[0]
        # print("bond_atms=",bond_atms)
        bond=mol.GetBondBetweenAtoms(int(bond_atms[0]),int(bond_atms[1])).GetIdx()
        Frag_mol = Chem.FragmentOnBonds(mol, [bond], dummyLabels=[(0, 0)])
        Frag_mols=Chem.GetMolFrags(Frag_mol, asMols=True)
        Frag_smis=[Chem.MolToSmiles(imol) for imol in Frag_mols]
        for ifrag in Frag_smis:
            # print("ifrag=", ifrag)
            re_p=re.compile(r'\(\*\)|\*|\[\*\]')
            ifrag_nodummy = re.sub(re_p, '', ifrag)  ## remove dummy atoms
            imatched=mol.GetSubstructMatches(Chem.MolFromSmarts(ifrag_nodummy))
            # print("imatched=",imatched)
            if len(imatched)>0:
                imatched = [i for ii in imatched for i in ii]
                if bond_atms[1] in imatched:
                    res_frags.append(ifrag)
    return list(set(res_frags))

def get_single_frag(smi, smarts):
    print(smi,'   ',smarts)
    ''' single dummy atom and the dummy atom in not included in the R group. For example, "[cH]1c(c2cc(cccc3)c3[nH]2)[cH]n4n[cH]nc4[*;!#1]1" '''
    re_p=re.compile(r'\[(\*((?!nH).)*?)\]') 
    re_p_other=re.compile(r'\[((?!\*).)*?\]') 
    smarts_readable = re.sub(re_p, '[1*]', smarts)  
    smarts_readable = re.sub(re_p_other, '*', smarts_readable) 
    print(f"smarts_readable={smarts_readable}")
    smarts_mol=Chem.MolFromSmarts(smarts_readable)
    # smarts_mol=mol_with_atom_index(smarts_mol, start=0)
    # display(smarts_mol)
    smart_mol_atoms=smarts_mol.GetAtoms()
    smart_dummies=smarts_mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    # print(smart_dummies)
    for idummy in smart_dummies:
        isotope=smart_mol_atoms[idummy[0]].GetIsotope()
        # print(f"isotope= {isotope}")
        if isotope==1:
            dummy_index=idummy[0]
    # print(Chem.MolToSmiles(smarts_mol))
    ###  molecule match section
    mol=get_mol(smi)
    mol=mol_with_atom_index(mol, start=0)
    mol_atoms=mol.GetAtoms()
    # display(mol)
    Smarts_tmp=re.sub(re_p, '[*]', smarts) 
    print(f"Smarts_tmp={Smarts_tmp}")
    mol_match=mol.GetSubstructMatches(Chem.MolFromSmarts(Smarts_tmp))[0]
    # print(f"mol_match={mol_match}")
    dummyIdx_mol=mol_match[dummy_index]
    idummy_neighs=mol_atoms[dummyIdx_mol].GetNeighbors()
    idummy_neighs=[iatom.GetIdx() for iatom in idummy_neighs]
    # print(f"idummy_neighs={idummy_neighs}")
    for idummy_neighb in idummy_neighs:
        if idummy_neighb in mol_match:
            R_atom=idummy_neighb   ## the atom in R group bonded to the core
    bond=mol.GetBondBetweenAtoms(int(dummyIdx_mol),int(R_atom)).GetIdx()
    Frag_mol = Chem.FragmentOnBonds(mol, [bond], dummyLabels=[(0, 0)])
    Frag_mols=Chem.GetMolFrags(Frag_mol, asMols=True)
    Frag_smis=[Chem.MolToSmiles(imol) for imol in Frag_mols]
    R_smi=''
    for ifrag in Frag_smis:
        ifrag_mol=Chem.MolFromSmiles(ifrag)
        imatched=ifrag_mol.GetSubstructMatches(Chem.MolFromSmarts(Smarts_tmp))
        # print(f"ifrag={ifrag}")
        # print("imatched=",imatched)
        # print(f"len(imatched)={len(imatched)}")
        if len(imatched)==0:
            R_smi=ifrag
    return R_smi


