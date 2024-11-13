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
sys.path.append("../")
from utils.common_utils import mapper,csvToExcel,get_mol,compute_FP,canonic_smiles,float_row,heavy_atm_smiles,kekulize_smi,sort_mol_sim_df
from utils.sarm_utils import get_core

from IPython.display import display, SVG, display_svg
import seaborn as sns
from matplotlib import pyplot
from pathlib import Path
import glob
from functools import partial
from pandarallel import pandarallel
# n_jobs=40
# pandarallel.initialize(nb_workers=n_jobs)

def kekulize_smi(smi):
    try:
        mol=get_mol(smi)
        Chem.Kekulize(mol)
        smi=Chem.MolToSmiles(mol, kekuleSmiles=True)
        return smi
    except Exception as e:
        print(e)
        return smi

def remove_dummy(smi):
    FragMol=Chem.MolFromSmiles(smi)
    matched=FragMol.GetSubstructMatches(Chem.MolFromSmarts("[#7][#0]"))
    atoms=FragMol.GetAtoms()
    for imatch in matched:
        atoms[imatch[1]].SetAtomicNum(1)
    # imol.UpdatePropertyCache()
    Chem.SanitizeMol(FragMol)
    FragMol=Chem.RemoveHs(FragMol)
    smi=Chem.MolToSmiles(FragMol)
    re_p=re.compile(r'\(\*\)|\*|\[\*\]|\-')
    ifrag_nodummy = re.sub(re_p, '', smi)  ## remove dummy atoms
    return ifrag_nodummy
    
def get_RGrps_dummy(mol,Core_smi):
    mol=get_mol(mol)
    # FragMol=get_mol(Frag_smi)
    CoreMol=Chem.MolFromSmarts(Core_smi)  ### Must be read as SMARTS
    matchedFrag=mol.GetSubstructMatches(CoreMol)[0]  ## Only the first match will be considered!
    print(matchedFrag)
    matchedDummy=CoreMol.GetSubstructMatches(Chem.MolFromSmarts('[*][#0]'))  ##[CoreAtom, DummyAtom]
    if len(matchedDummy)==2:
        if matchedDummy[0][1]>matchedDummy[1][1]:
            matchedDummy=[matchedDummy[1],matchedDummy[0]]
    print(matchedDummy)
    BondList=[]
    AtomPairList=[]
    '''  Put the dummy atoms in order  '''
    
    for imatDummy in matchedDummy:
        bond_coreAtom=matchedFrag[imatDummy[0]]
        bond_RAtom=matchedFrag[imatDummy[1]]
        AtomPairList.append([bond_coreAtom, bond_RAtom])
        bond=mol.GetBondBetweenAtoms(bond_coreAtom, bond_RAtom).GetIdx()
        BondList.append(bond)
    dummyLabels=[(0, 0) for ibond in BondList]    
    FragedMol= Chem.FragmentOnBonds(mol, BondList, dummyLabels=dummyLabels)
    '''  Get the Fragments  '''
    fragSmis=[get_frag(FragedMol, iatomPair[1]) for iatomPair in AtomPairList]
    return fragSmis

def get_frag(fraged_mol, atomIdx):
    Atoms=fraged_mol.GetAtoms()
    count=0
    newAtoms=[Atoms[atomIdx]]
    RAtomIdx=[atomIdx]
    while len(newAtoms) > 0:
        count+=1
        tmpAtoms=[]
        for iAtom in newAtoms:
            negbs=iAtom.GetNeighbors()
            for inegb in negbs:
                inegbIdx=inegb.GetIdx()
                if inegbIdx not in RAtomIdx:
                    RAtomIdx.append(inegbIdx)
                    tmpAtoms.append(inegb)
        newAtoms=tmpAtoms
        if count>100:  ## get avoid of dead circle
            break
        # frag_smi = Chem.MolFragmentToSmiles(mol, RAtomIdx, canonical=True,isomericSmiles=False,kekuleSmiles=True) 
        frag_smi = Chem.MolFragmentToSmiles(fraged_mol, RAtomIdx, canonical=True,isomericSmiles=False) 
    return kekulize_smi(frag_smi)

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

def replace_nH(ismarts):
    ismarts=ismarts.replace('[nH]','n')
    ismarts=ismarts.replace('-','')
    return ismarts

def add_singleBond_dummy(smarts):
    newSmarts=''
    for idx,ilet in enumerate(smarts):
        if idx==0 and ilet=='*':
            newSmarts+='*-'
        elif idx==len(smarts)-1 and ilet=='*':
            newSmarts+='-*'
        elif ilet=='*':
            newSmarts+='-*'
        else:
            newSmarts+=ilet
    return newSmarts
        
def match_frag(ismi,ismarts):
    # print("ismi",ismi)
    # print("ismarts",ismarts)
    mol = Chem.MolFromSmiles(ismi)
    ismarts=replace_nH(ismarts)
    smartsMol = Chem.MolFromSmarts(ismarts) #,sanitize=False
    matched=mol.GetSubstructMatches(smartsMol)
    if len(matched)>0:
        # print(matched)
        return 1
    else:
        return 0
    
def get_activity_info(frag, file_act, cols):
    df_act=pd.read_csv(file_act)
    df_act['matched']=df_act.apply(lambda x: match_frag(x['Cano_SMILES'],ismarts=frag),axis=1)
    df_match=df_act[df_act['matched']==1]
    # print(df_match)
    if len(df_match)==0:
        return [['','',''] for icol in cols]
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

def add_singleBond_dummy(smarts):
    newSmarts=''
    for idx,ilet in enumerate(smarts):
        if idx==0 and ilet=='*':
            newSmarts+='*!@'
        elif idx==len(smarts)-1 and ilet=='*':
            newSmarts+='!@*'
        elif ilet=='*':
            newSmarts+='!@*'
        else:
            newSmarts+=ilet
    return newSmarts


    
def get_full_core(smi, Core_smi):
    mol=get_mol(smi)
    # FragMol=get_mol(Frag_smi)
    CoreMol=Chem.MolFromSmarts(Core_smi)  ### Must be read as SMARTS
    
    # matchedFrag=mol.GetSubstructMatches(CoreMol)[0]  ## Only the first match will be considered!
    for matchedFrag in mol.GetSubstructMatches(CoreMol):
        # print(matchedFrag)
        matchedDummy=CoreMol.GetSubstructMatches(Chem.MolFromSmarts('[*][#0]'))  ##[CoreAtom, DummyAtom]
        if len(matchedDummy)==2:
            if matchedDummy[0][1]>matchedDummy[1][1]:
                matchedDummy=[matchedDummy[1],matchedDummy[0]]
        # print(matchedDummy)
        BondList=[]
        AtomPairList=[]
        '''  Put the dummy atoms in order  '''
        
        for imatDummy in matchedDummy:
            bond_coreAtom=matchedFrag[imatDummy[0]]
            bond_RAtom=matchedFrag[imatDummy[1]]
            AtomPairList.append([bond_coreAtom, bond_RAtom])
            bond=mol.GetBondBetweenAtoms(bond_coreAtom, bond_RAtom)
            if not bond.IsInRing():
                BondList.append(bond.GetIdx())
        
    if len(BondList)>0:
        dummyLabels=[(0, 0) for ibond in BondList]    
        FragedMol= Chem.FragmentOnBonds(mol, BondList, dummyLabels=dummyLabels)
        '''  Get the Fragments  '''
        coreSmi=get_frag(FragedMol, AtomPairList[0][0])
        return coreSmi
    else:
        print(smi,Core_smi)
        return ''

def is_number(InString):
    try:
        string_num=float(InString)
        return string_num
    except Exception as e:
        return ''
    
'''  Create Tree Utils  '''
def is_parent(parent,child):
    '''  If the compound is a substructure of another and don't have a molecule as its own substurcture, it's a parent.  '''
    parent=parent.replace('[nH]','n')  ## Correct the influence of ionization
    mol_parent=Chem.rdmolfiles.MolFromSmarts(parent)
    mol_child=Chem.rdmolfiles.MolFromSmiles(child,sanitize=False)
    mol_child.UpdatePropertyCache()
    matched=mol_child.GetSubstructMatches(mol_parent)
    if len(matched)>0:
        # print(matched)
        return True
    else:
        return False
    
def is_child(child,parent):
    '''  If the compound is a substructure of another and don't have a molecule as its own substurcture, it's a parent.  '''
    parent=parent.replace('[nH]','n')  ## Correct the influence of ionization
    mol_parent=Chem.rdmolfiles.MolFromSmarts(parent)
    mol_child=Chem.rdmolfiles.MolFromSmiles(child,sanitize=False)
    mol_child.UpdatePropertyCache()
    matched=mol_child.GetSubstructMatches(mol_parent)
    if len(matched)>0:
        # print(matched)
        return True
    else:
        return False

def remove_dummy(smi):
    FragMol=Chem.MolFromSmiles(smi)
    matched=FragMol.GetSubstructMatches(Chem.MolFromSmarts("[#7][#0]"))
    atoms=FragMol.GetAtoms()
    for imatch in matched:
        atoms[imatch[1]].SetAtomicNum(1)
    # imol.UpdatePropertyCache()
    Chem.SanitizeMol(FragMol)
    FragMol=Chem.RemoveHs(FragMol)
    smi=Chem.MolToSmiles(FragMol)
    re_p=re.compile(r'\(\*\)|\*|\[\*\]|\-')
    ifrag_nodummy = re.sub(re_p, '', smi)  ## remove dummy atoms
    return ifrag_nodummy

def real_sonNode(smi,children_smis,parent_dict):
    has_parent=set(children_smis) & set(parent_dict[smi])
    # print(has_parent)
    # sys.exit()
    if len(has_parent)>0:
        return False
    else:
        return True

def remove_ionizationForm(smiList):
    '''  Remove the same compounds with different ionization state '''
    lenSmiList=len(smiList)
    rmIndex=[]
    for i,ismi in enumerate(smiList):
        for j in range(i+1,lenSmiList,1):
            jsmi=smiList[j]
            if is_child(ismi,jsmi) and is_child(jsmi,ismi):
                rmIndex.append(j)
    smiList=np.array(smiList)
    smiList=np.delete(smiList, list(set(rmIndex)))
    return smiList

def find_children_single(ismi,smi_list):
    children_list=[]
    for jsmi in smi_list:
        if ismi==jsmi: 
            continue
        if is_child(ismi,jsmi):
            children_list.append([ismi, jsmi]) 
    return children_list

def find_children(smi_list):
    find_children_single_p=partial(find_children_single,smi_list=smi_list)
    key_values=mapper(n_jobs)(find_children_single_p,smi_list)
    res_dict={}
    for ismi in smi_list:
        res_dict[ismi]=[]
    for ikey_value_list in key_values:
        if len(ikey_value_list)==0:
            continue
        for ikey_value in ikey_value_list:
            res_dict[ikey_value[0]].append(ikey_value[1])
    return res_dict

def find_parents_single(ismi,smi_list):
    # res_dict[ismi]=[]
    parents_list=[]
    for jsmi in smi_list:
        if ismi==jsmi: 
            continue
        if is_child(jsmi,ismi):
            parents_list.append([ismi,jsmi]) 
    return parents_list

def find_parents(smi_list):
    find_parents_single_p=partial(find_parents_single, smi_list=smi_list)
    key_values=mapper(n_jobs)(find_parents_single_p, smi_list)
    res_dict={}
    for ismi in smi_list:
        res_dict[ismi]=[]
    for ikey_value_list in key_values:
        if len(ikey_value_list)==0:
            continue
        for ikey_value in ikey_value_list:
            res_dict[ikey_value[0]].append(ikey_value[1])
    return res_dict

def if_root(ismi,smi_list):
    has_parent=False
    has_child=False 
    for jsmi in smi_list:
        if ismi==jsmi:
            continue
        # print(ismi,jsmi)
        if is_child(ismi,jsmi):
            has_child=True
        if is_child(jsmi,ismi):
            has_parent=True
        if has_child and has_parent:
            return None
    if has_child and not has_parent:
        return ismi    

def find_root(smi_list): 
    if_root_p=partial(if_root,smi_list=smi_list)
    roots=mapper(n_jobs)(if_root_p,smi_list)      
    roots=list(set(roots))
    roots=[iroot for iroot in roots if iroot!=None]
    return roots

def smiles2shapeSmarts(smi):
    ''' The SMILES should be aromatic '''
    for ichar in ['[nH]','c','n']:
        smi=smi.replace(ichar,'$')
    smi=smi.replace('$',"[c,n]")
    return smi    

def scaffoldFamily(scalffod, scaffolds):
    family=[]
    for iscaf in scaffolds:
        if is_parent(iscaf,scalffod) :#or is_parent(scalffod,iscaf)
            family.append(iscaf)
        elif abs(weight(scalffod)-weight(iscaf))<24:
            if is_parent(smiles2shapeSmarts(iscaf),scalffod) or is_parent(smiles2shapeSmarts(scalffod),iscaf):
                family.append(iscaf)
    return family

def scaffold_info(scaffolds):
    scaff_NoDumScaff={}
    NoDumScaff_scaff={}
    NoDumScaffs=[remove_dummy(ismi) for ismi in scaffolds]
    for idx,iscaf in enumerate(scaffolds):
        scaff_NoDumScaff[scaffolds[idx]]=NoDumScaffs[idx]
        NoDumScaff_scaff[NoDumScaffs[idx]]=scaffolds[idx]
    return scaff_NoDumScaff,NoDumScaff_scaff

'''  Scaffold manipulation  '''
def get_full_core(smi, Core_smi):
    '''  get all the complete fragments with R group in the dataset  '''
    mol=get_mol(smi)
    # FragMol=get_mol(Frag_smi)
    CoreMol=Chem.MolFromSmarts(Core_smi)  ### Must be read as SMARTS
    BondList=[]
    # matchedFrag=mol.GetSubstructMatches(CoreMol)[0]  ## Only the first match will be considered!
    for matchedFrag in mol.GetSubstructMatches(CoreMol):
        # print(matchedFrag)
        matchedDummy=CoreMol.GetSubstructMatches(Chem.MolFromSmarts('[*][#0]'))  ##[CoreAtom, DummyAtom]
        if len(matchedDummy)==2:
            if matchedDummy[0][1]>matchedDummy[1][1]:
                matchedDummy=[matchedDummy[1],matchedDummy[0]]
        # print(matchedDummy)
        BondList=[]
        AtomPairList=[]
        '''  Put the dummy atoms in order  '''
        
        for imatDummy in matchedDummy:
            bond_coreAtom=matchedFrag[imatDummy[0]]
            bond_RAtom=matchedFrag[imatDummy[1]]
            AtomPairList.append([bond_coreAtom, bond_RAtom])
            bond=mol.GetBondBetweenAtoms(bond_coreAtom, bond_RAtom)
            if not bond.IsInRing():
                BondList.append(bond.GetIdx())
        
    if len(BondList)>0:
        dummyLabels=[(0, 0) for ibond in BondList]    
        FragedMol= Chem.FragmentOnBonds(mol, BondList, dummyLabels=dummyLabels)
        '''  Get the Fragments  '''
        coreSmi=get_frag(FragedMol, AtomPairList[0][0])
        return coreSmi
    else:
        print(smi,Core_smi)
        return ''

def get_atomRingInfo(smi):
    mol=get_mol(smi)
    rInfo=mol.GetRingInfo()
    atomRingInfo=np.zeros(mol.GetNumAtoms())
    for iring in rInfo.AtomRings():
        iring=list(iring)
        atomRingInfo[iring]+=1
    return atomRingInfo
       
def match_core(ismi,ismarts):
    # print("ismi",ismi)
    # print("ismarts",ismarts)
    # mol = Chem.MolFromSmiles(ismi)
    mol=get_mol(ismi)
    # ismarts=replace_nH(ismarts)
    if '*' in ismarts:
        ismartsNodummy=remove_dummy(ismarts)
    else:
        ismartsNodummy=ismarts
    smartsMol = Chem.MolFromSmarts(ismartsNodummy) #,sanitize=False
    matched=mol.GetSubstructMatches(smartsMol)
    if len(matched)>0:
        smartRI=get_atomRingInfo(ismartsNodummy)
        # print(f"smartRI= {smartRI}")
        molRI=get_atomRingInfo(mol)
        for imatch in matched:
            submolRI=molRI[list(imatch)]
            # print(f"submolRI= {submolRI}")
            diff=np.absolute(smartRI-submolRI).sum()
            # print(matched)
            if diff==0:
                return 1
        return 0
    else:
        return 0
    
def get_atomPos(sdf):
    suppl = Chem.SDMolSupplier(sdf)
    mol=suppl[0]
    molNumAtoms=mol.GetNumAtoms()
    molAtoms=mol.GetAtoms()
    atomProps=[]
    if molNumAtoms>150:
        return atomProps
    conf = mol.GetConformer()
    # print(f"Molecule name: {mol.GetProp('_Name')}")
    # print(f"Number of atoms: {mol.GetNumAtoms()}")
    for i in range(molNumAtoms):
        molName=mol.GetProp('_Name')
        ipos=conf.GetAtomPosition(i)
        iAtomicNum=molAtoms[i].GetAtomicNum()
        iFormalCharge=molAtoms[i].GetFormalCharge()
        iIsAromatic=int(molAtoms[i].GetIsAromatic())
        iIsInRing=int(molAtoms[i].IsInRing())
        atomProps.append([molName,i,ipos.x,ipos.y,ipos.z,iAtomicNum,iFormalCharge,iIsAromatic,iIsInRing])
    # print(f"{i}:  {atomProps}")
    return atomProps