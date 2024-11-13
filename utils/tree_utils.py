import pandas as pd
import numpy as np
from my_toolset.my_utils import get_mol,mapper
from my_toolset.drawing_utils import show_mols
import re
from rdkit import Chem
import os,sys
import logging
from functools import reduce,partial
# sys.path.append("../")
from .common_utils import mapper,csvToExcel,get_mol,compute_FP,canonic_smiles,float_row,heavy_atm_smiles,kekulize_smi,remove_dummy
from IPython.display import display
from rdkit.Chem import AllChem,Draw,rdFMCS
import graphviz
from pathlib import Path
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import MolDraw2DCairo
from pandarallel import pandarallel
# n_jobs=40
# pandarallel.initialize(nb_workers=n_jobs)

# def fetch_table_info(df_table,smi):
#     info=df_table.loc[smi,['Size','Items_count']]
#     # res_str=''
#     # print(info)
#     if type(info['Items_count'])==list:
#         res_str=''+str(info['Size'].to_list())+str(info['Items_count'].to_list())
#     else:
#         res_str=''+str(info['Size'])+str(info['Items_count'])
#     return res_str

class node():
    def __init__(self,smi) -> None:
        self.root=True
        self.SMILES=smi
        self.children=[]
        self.parents=[]
        
    def add_child(self,node):
        node.root=False
        self.children.append(node)
        
    def print_tree(self,df_table,opfile):
        empty_node=node('')
        tree_list=[]
        level=1
        print(f"Level {level}:  {self.SMILES} {fetch_table_info(df_table,self.SMILES)}")
        tree_list.append([[self.SMILES]])
        ilevel_smiles=[]
       
        for idx,inode in enumerate(self.children):
            print(f'node_{idx}',end=' ')
            opfile.write(f"node_{idx}  ")
            print(f"{inode.SMILES} {fetch_table_info(df_table, inode.SMILES)}", end=' ')
            opfile.write(f"{inode.SMILES} {fetch_table_info(df_table, inode.SMILES)} ")
            ilevel_smiles.append(inode.SMILES)
            # tmp_next_level_inodes.append(inode)
            print(' |',end=' ')
            opfile.write(' |')
        tree_list.append([ilevel_smiles])
        opfile.write(f"Level {level}:  {self.SMILES} {fetch_table_info(df_table,self.SMILES)}\n")
        next_level_nodes=self.children
        beforeLevel_with_child=[idx for idx,inode in enumerate(next_level_nodes)]
        smis_next_level_node=1
        count=0
        while smis_next_level_node>0:
            level+=1
            print(f"Level {level}:", end='  ')
            opfile.write(f"Level {level}:  ")
            tmp_next_level_inodes=[]
            ilevel_smiles=[]
            ilevel_with_child=[]
            inode_with_child_count=-1
            for idx,inode in enumerate(next_level_nodes):
                # print(inode.SMILES, end='  ')
                print(f'node_{idx}',end=' ')
                opfile.write(f"node_{idx}  ")
                if len(inode.children)>0:
                    
                    inode_smiles=[]
                    for inode_child in inode.children:
                        inode_with_child_count+=1
                        ilevel_with_child.append(inode_with_child_count)
                        print(f"{inode_child.SMILES} {fetch_table_info(df_table, inode_child.SMILES)}", end=' ')
                        opfile.write(f"{inode_child.SMILES} {fetch_table_info(df_table, inode_child.SMILES)} ")
                        inode_smiles.append(inode_child.SMILES)
                        tmp_next_level_inodes.append(inode_child)
                    ilevel_smiles.append(inode_smiles)
                    
                else:
                    ilevel_smiles.append([])
                    ilevel_with_child.append(-1)
                    # tmp_next_level_inodes.append(empty_node)
                print(' |',end=' ')
                opfile.write(' |')
                
            '''  Put current child SMILES in right position  '''
            print('beforeLevel_with_child', beforeLevel_with_child)
            print('ilevel_smiles', ilevel_smiles)
            ilevel_smiles_withPos=[[] for i in beforeLevel_with_child]
            for idxSmi, smiList in enumerate(ilevel_smiles):
                ilevel_smiles_withPos[beforeLevel_with_child.index(idxSmi)]=smiList
            beforeLevel_with_child=ilevel_with_child
            
                
                
            tree_list.append(ilevel_smiles)
            print(' ')
            opfile.write('\n')
            next_level_nodes=tmp_next_level_inodes
            # for i_next_level_node in next_level_nodes:
            smis_next_level_node=[i_next_level_node.SMILES for i_next_level_node in tmp_next_level_inodes if i_next_level_node.SMILES!='' ]
            smis_next_level_node=len(smis_next_level_node)
            
            
            count+=1
            if count>500:
                break
        return tree_list

def replace_nH(ismarts,nH=1,DH=1):
    ''' fix the ionization state of nitrogen '''
    ismarts=ismarts.replace('-','')
    ismarts=ismarts.replace('\\','')
    ismarts=ismarts.replace('\/','')  ## remove isometric 
    if DH:
        re_p=re.compile(r'\(\[2H\]\)|\[2H\]|\(\)')
        ismarts = re.sub(re_p, '', ismarts)  ## remove hydrogen isotope 'Deuterium'
    if nH:
        ismarts=ismarts.replace('[nH]','n')
    return ismarts


def get_atomRingInfo(smi):
    mol=get_mol(smi)
    rInfo=mol.GetRingInfo()
    atomRingInfo=np.zeros(mol.GetNumAtoms())
    for iring in rInfo.AtomRings():
        iring=list(iring)
        atomRingInfo[iring]+=1
    return atomRingInfo
       
def match_core(ismi,ismarts):
    '''  exactly math the ring system  '''
    # print("ismi",ismi)
    # print("ismarts",ismarts)
    # mol = Chem.MolFromSmiles(ismi)
    mol=get_mol(ismi)
    
    if '*' in ismarts:    ## dummy atom has ring info which is not necessary
        ismartsNodummy=remove_dummy(ismarts)
        # ismartsNodummy=replace_nH(ismartsNodummy)
    else:
        ismartsNodummy=ismarts
    # print(f"ismartsNodummy:  {ismartsNodummy}")
    # ismartsNodummy=ismarts
    smartsMol = Chem.MolFromSmarts(replace_nH(ismartsNodummy)) #math different ionization state
    matched=mol.GetSubstructMatches(smartsMol)
    if len(matched)>0:
        smartRI=get_atomRingInfo(ismartsNodummy)  ## ionization state is necessary for ring ifno retrieval
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

def is_parent(parent,child,ring=0):
    '''  If the compound is a substructure of another and don't have a molecule as its own substurcture, it's a parent.  '''
    # parent=parent.replace('[nH]','n')  ## Correct the influence of ionization， adding this will lead to ring system failed!
    mol_parent=Chem.rdmolfiles.MolFromSmarts(parent)
    mol_child=Chem.rdmolfiles.MolFromSmiles(child,sanitize=False)
    mol_child.UpdatePropertyCache()
    matched=mol_child.GetSubstructMatches(mol_parent)
    if len(matched)>0:
        if ring:
            try:
                if match_core(child,parent):
                    return True
            except Exception as e:
                print(f"child: {child}")
                print(f"parent: {parent}")
                print(e)
                return False
        else:
            return True
    return False

# def find_root(smi_list): 
#     roots=[]      
#     for ismi in smi_list:
#         has_parent=False
#         has_child=False 
#         for jsmi in smi_list:
#             if ismi==jsmi:
#                 continue
#             # print(ismi,jsmi)
#             if is_parent(ismi,jsmi):
#                 has_child=True
#             if is_parent(jsmi,ismi):
#                 has_parent=True
#             if has_child and has_parent:
#                 break
#         if has_child and not has_parent:
#             roots.append(ismi)
#     roots=list(set(roots))
#     return roots

# def remove_dummy(smi):  ## moved to common utils
#     FragMol=Chem.MolFromSmiles(smi)
#     matched=FragMol.GetSubstructMatches(Chem.MolFromSmarts("[#7][#0]"))
#     atoms=FragMol.GetAtoms()
#     for imatch in matched:
#         atoms[imatch[1]].SetAtomicNum(1)
#     # imol.UpdatePropertyCache()
#     Chem.SanitizeMol(FragMol)
#     FragMol=Chem.RemoveHs(FragMol)
#     smi=Chem.MolToSmiles(FragMol)
#     re_p=re.compile(r'\(\*\)|\*|\[\*\]|\-')
#     ifrag_nodummy = re.sub(re_p, '', smi)  ## remove dummy atoms
#     return ifrag_nodummy

# def find_children(smi_list):
#     res_dict={}
#     for ismi in smi_list:
#         res_dict[ismi]=[]
#         for jsmi in smi_list:
#             if ismi==jsmi: 
#                 continue
#             if is_parent(ismi,jsmi):
#                 res_dict[ismi].append(jsmi)
#     return res_dict

# def find_parents(smi_list):
#     res_dict={}
#     for ismi in smi_list:
#         res_dict[ismi]=[]
#         for jsmi in smi_list:
#             if ismi==jsmi: 
#                 continue
#             if is_parent(jsmi,ismi):
#                 res_dict[ismi].append(jsmi)
#     return res_dict

def real_sonNode(smi,children_smis,parent_dict):
    has_parent=set(children_smis) & set(parent_dict[smi])
    # print(has_parent)
    # sys.exit()
    if len(has_parent)>0:
        return False
    else:
        return True
    
''' Functions for bulding a complete Trees  '''


def fetch_table_SMILES(df_table,smi):
    info=df_table.loc[smi,'SMILES']
    # print(info)
    try:
        # print(type(info))
        res_str=info.values[0]
        return res_str
    except:
        res_str=info
        return res_str

# def fetch_table_info(df_table,smi):
#     info=df_table.loc[smi,['Size','Items_count']]
#     # res_str=''
#     # print(info)
#     try:
#         res_str=''+"%10s"%str(info['Size'].values)+"%10s"%str(info['Items_count'].values)
#     except:
#         res_str=''+"%10s"%str(info['Size'])+"%10s"%str(info['Items_count'])
#     return res_str

def show_mols(smiles_mols,legends=[],subImgSize=(600,600)):
    '''Display multiple mols with legends'''
    mols=[Chem.rdmolfiles.MolFromSmiles(ismi,sanitize=False) for ismi in smiles_mols]
    mol_cls,legends_cls=[],[]
    for i in range(len(mols)):
        if mols[i]==None:
            continue
        mol_cls.append(mols[i])
        if len(legends)>i:
            legends_cls.append(legends[i])
        else:
            legends_cls.append('')
        Draw.MolToFile(mols[i], 'Images/test.png',size=(100,100))  ## strange thing withtout this below not working
    # svg=Draw.MolsToGridImage(mol_cls,subImgSize=subImgSize, molsPerRow=3,useSVG=False,legends=legends_cls,sanitize=False)
    # png=Draw.MolsToGridImage(mol_cls,subImgSize=subImgSize,useSVG=False,molsPerRow=3,legends=legends_cls,sanitize=False)
    if len(mols)>3:
        molsPerRow=3
    else:
        molsPerRow=len(mols)
    png=Draw.MolsToGridImage(mols,molsPerRow=molsPerRow,subImgSize=subImgSize,legends=legends_cls,returnPNG=False) 
    return png

def show_mols_cairo(ismi,legend):
    mol = Chem.MolFromSmiles(ismi)
    drawer = MolDraw2DCairo(400, 400)
    drawer.drawOptions().legendFontSize = 160
    drawer.DrawMolecule(mol,legend=legend)
    drawer.FinishDrawing()
    img = drawer.GetDrawingText()
    # with open('mol.png', 'wb') as f:
    #     f.write(img)
    return img

# ''' Create tree functions  '''
# def replace_nH(ismarts):
#     ismarts=ismarts.replace('[nH]','n')
#     ismarts=ismarts.replace('-','')
#     return ismarts

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
    
def remove_ionizationForm(smiList):
    '''  Remove the same compounds with different ionization state '''
    lenSmiList=len(smiList)
    rmIndex=[]
    for i,ismi in enumerate(smiList):
        for j in range(i+1,lenSmiList,1):
            jsmi=smiList[j]
            if is_parent(ismi,jsmi) and is_parent(jsmi,ismi):
                rmIndex.append(j)
    smiList=np.array(smiList)
    smiList=np.delete(smiList, list(set(rmIndex)))
    return smiList
    
def fetch_table_info(df_table,smi):
    info=df_table.loc[smi,['Size','Items_count']]
    # res_str=''
    # print(info)
    try:
        res_str=''+'{0: >15}'.format(str(info['Size'].values)+str(info['Items_count'].values))+'\n'
    except:
        res_str=''+'{0: >15}'.format(str(info['Size'])+str(info['Items_count']))+'\n'
    return res_str

def get_activity_info(df_act,frag, actCols=[]):
    df_act=df_act.copy()
    df_act['matched']=df_act.apply(lambda x: match_frag(x['Cano_SMILES'],ismarts=frag),axis=1)
    df_match=df_act[df_act['matched']==1]
    # print(df_match)
    if len(df_match)==0:
        return ''
    if len(df_match)==1:
        means=df_match.iloc[0]
        stds=df_match.iloc[0]
        medians=df_match.iloc[0]
    if len(df_match)>1:
        means=df_match.mean()
        stds=df_match.std()
        medians=df_match.median()
    actInfo=''
    for itarget in actCols:#["JAK1ToJAK2","JAK1","JAK2"]:
        iInfo='{0: <15}'.format(f"{itarget}:")
        for imetric in [means,stds,medians]:
            iInfo+='{0: <10}'.format(f"{round(imetric[itarget],1)}")
        iInfo+="|\n"
        actInfo+=iInfo
    return actInfo

def find_children_single(ismi,smi_list,ring=0):
    children_list=[]
    for jsmi in smi_list:
        if ismi==jsmi: 
            continue
        if is_parent(ismi,jsmi,ring):
            children_list.append([ismi,jsmi]) 
    return children_list

def find_children(smi_list,ring=0,n_jobs=20):
    find_children_single_p=partial(find_children_single,smi_list=smi_list,ring=ring)
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

def find_parents_single(ismi,smi_list,ring=0):
    # res_dict[ismi]=[]
    parents_list=[]
    for jsmi in smi_list:
        if ismi==jsmi: 
            continue
        if is_parent(jsmi,ismi,ring):
            parents_list.append([ismi,jsmi]) 
    return parents_list

def find_parents(smi_list,n_jobs=20,ring=0):
    find_parents_single_p=partial(find_parents_single,smi_list=smi_list,ring=ring)
    key_values=mapper(n_jobs)(find_parents_single_p,smi_list)
    res_dict={}
    for ismi in smi_list:
        res_dict[ismi]=[]
    for ikey_value_list in key_values:
        if len(ikey_value_list)==0:
            continue
        for ikey_value in ikey_value_list:
            res_dict[ikey_value[0]].append(ikey_value[1])
    return res_dict

def if_root(ismi,smi_list,ring=1):
    has_parent=False
    has_child=False 
    for jsmi in smi_list:
        if ismi==jsmi:
            continue
        # print(ismi,jsmi)
        if is_parent(ismi,jsmi,ring=ring):
            has_child=True
        if is_parent(jsmi,ismi,ring=ring):
            has_parent=True
        if has_child and has_parent:
            return None
    if has_child and not has_parent:
        return ismi    

def find_root(smi_list,n_jobs=20,ring=1): 
    if_root_p=partial(if_root,smi_list=smi_list,ring=ring)
    roots=mapper(n_jobs)(if_root_p,smi_list)      
    roots=list(set(roots))
    roots=[iroot for iroot in roots if iroot!=None]
    return roots

def show_cpd(img_icpd,df_table,df_act,actCols):
    imgPath=img_icpd[0]
    icpd=img_icpd[1]
    if icpd=='':
        return None
    icpd_dummy=fetch_table_SMILES(df_table,icpd)
    img=show_mols_cairo(icpd_dummy,fetch_table_info(df_table,icpd)+get_activity_info(df_act,icpd,actCols))
    with open(imgPath, 'wb') as f:
        f.write(img)
    return [imgPath,icpd_dummy]
    