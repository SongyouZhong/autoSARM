import os,sys
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem,Draw
from rdkit.Chem import rdFMCS as MCS
import pandas as pd
import numpy as np
import glob
import math
import itertools
from .tree_utils import get_atomRingInfo,replace_nH
from .common_utils import remove_dummy
import json
import re



schrodinger_root="/public/software/schrodinger/2021-2"
extract_sdf_fromMaegz="/public/home/zhangjie/Projects/autosarm/scripts/extract_sdf_fromMaegz.py"

def extract_sdf(maegz,saveDir,scoreLimit):
    Path(saveDir).mkdir(exist_ok=True,parents=True)
    os.system(f"{schrodinger_root}/run  {extract_sdf_fromMaegz} --maegz {maegz} --saveDir {saveDir} --scoreLimit  {scoreLimit}")

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

def grid_atom(sdfFolder):
    ''' grid the atom position of sdfs in a folder  '''
    # savePath=Path(sdfFolder)
    sdfFiles=glob.glob(f"{sdfFolder}/*.sdf")
    atomPos_list=[]     
    for isdf in sdfFiles:
        atomPos=get_atomPos(isdf)
        for iatmPos in atomPos:
            atomPos_list.append(iatmPos)
    df_ActAtomPos=pd.DataFrame(atomPos_list,columns=['molName','atomIndex','X','Y','Z','AtomicNum','FormalCharge','IsAromatic','IsInRing'])
    df_ActAtomPos.to_csv(f"{sdfFolder}/AMPK_actSet_atomPos.csv", index=None)
    print(f"Total atoms in the active set:  {len(atomPos_list)}.")
    grid_parameters={}  ## to save the grid parameters
    df_ActAtomPos=pd.read_csv(f"{sdfFolder}/AMPK_actSet_atomPos.csv")
    df_ActAtomPos[['X','Y','Z']]=df_ActAtomPos[['X','Y','Z']].round(1)
    xyzMax=df_ActAtomPos[['X','Y','Z']].max()
    xyzMin=df_ActAtomPos[['X','Y','Z']].min()
    grid_parameters['xyzMax']=xyzMax.to_dict()
    grid_parameters['xyzMin']=xyzMin.to_dict()
    print(f"xyzMax= {xyzMax}")
    print(f"xyzMin= {xyzMin}")
    cellWidth=1
    grid_parameters['cellWidth']=cellWidth
    xyzRange=(xyzMax-xyzMin)/cellWidth
    for iax in ['X','Y','Z']:
        xyzRange[iax]=int(math.ceil(xyzRange[iax]))
    xList=[xyzMin['X']+cellWidth*i for i in range(int(xyzRange['X']))]
    yList=[xyzMin['Y']+cellWidth*i for i in range(int(xyzRange['Y']))]
    zList=[xyzMin['Z']+cellWidth*i for i in range(int(xyzRange['Z']))]
    xyzIndex=list(itertools.product(xList,yList,zList))
    df_atmCount=pd.DataFrame(xyzIndex,columns=['X','Y','Z'])
    df_atmCount[['X','Y','Z']]=df_atmCount[['X','Y','Z']].astype(str)
    df_atmCount.set_index(['X','Y','Z'], drop=False, inplace=True)
    df_atmCount['Count_Carbon']=0  ## All Carbon
    df_atmCount['Count_ac']=0  ## All Carbon
    df_atmCount['Count_C']=0  ## All Carbon
    df_atmCount['Count_Oxygen']=0  ## All Carbon
    df_atmCount['Count_ao']=0  ## All Carbon
    df_atmCount['Count_O']=0  ## All Carbon
    df_atmCount['Count_Nitrogen']=0  ## All Carbon
    df_atmCount['Count_an']=0  ## All Carbon
    df_atmCount['Count_N']=0  ## All Carbon
    df_atmCount['Count_All']=0   ## All Atom
    for idx,irow in df_ActAtomPos.iterrows():
        ixPos=int((irow['X']-xyzMin['X'])/cellWidth)*cellWidth+xyzMin['X']
        ixPos=str(ixPos)
        iyPos=int((irow['Y']-xyzMin['Y'])/cellWidth)*cellWidth+xyzMin['Y']
        iyPos=str(iyPos)
        izPos=int((irow['Z']-xyzMin['Z'])/cellWidth)*cellWidth+xyzMin['Z']
        izPos=str(izPos)
        try:
            tmp=df_atmCount.loc[(ixPos,iyPos,izPos),'Count_All']
        except Exception as e:
            print(e)
            continue
        df_atmCount.loc[(ixPos,iyPos,izPos),'Count_All']+=1
        if irow['AtomicNum']==6:
            df_atmCount.loc[(ixPos,iyPos,izPos),'Count_Carbon']+=1
            if irow['IsAromatic']==1:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_ac']+=1
            else:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_C']+=1
        if irow['AtomicNum']==7:
            df_atmCount.loc[(ixPos,iyPos,izPos),'Count_Nitrogen']+=1
            if irow['IsAromatic']==1:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_an']+=1
            else:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_N']+=1
        if irow['AtomicNum']==8:
            df_atmCount.loc[(ixPos,iyPos,izPos),'Count_Oxygen']+=1
            if irow['IsAromatic']==1:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_ao']+=1
            else:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_O']+=1
                
    df_atmCount=df_atmCount[df_atmCount['Count_All']>0]
    return df_atmCount,grid_parameters

def create_grid_pse(gridCountFile,savePseFile,group='atomDensity'):
    from importlib import reload
    import pymol
    import pandas as pd
    from pymol import cmd, stored
    reload(pymol)
    
    # pymol.finish_launching() 
    # pl = pymol2.PyMOL()
    # pl.start()
    # cmd=pl.cmd
    
    df_ActAtomPos=pd.read_csv(gridCountFile)
    countCols=['Count_All','Count_Carbon','Count_ac','Count_C','Count_Nitrogen','Count_an','Count_N','Count_Oxygen','Count_ao','Count_O']
    
    countColMax=df_ActAtomPos[countCols].max().max()
    for icountCol in countCols:
        # countColMax=df_ActAtomPos[icountCol].max()
        df_ActAtomPos[icountCol]=df_ActAtomPos[icountCol]/countColMax
        
    greyColors=['tin','black','zinc','Grey50']
    orangeColors=['lightorange','brightorange','olive','deepolive','orange']
    blueColors=['lightblue','slate','tv_blue','skyblue','blue']
    redColors=['salmon','deepsalmon','tv_red','raspberry','red']
    # pymol.finish_launching()
    for idx,irow in df_ActAtomPos.iterrows():
        iX=irow['X']
        iY=irow['Y']
        iZ=irow['Z']
        
        if irow['Count_All']>0:
            objName=f"Count_All_Site_{group}"
            cmd.pseudoatom(objName, name=f"DUM{idx}", pos=[iX, iY, iZ], vdw=irow['Count_All'], color="gray80", label="")
            cmd.show("spheres", objName)
            cmd.color("gray80", objName)
        
        for idxCol,icol in enumerate(['Count_Carbon','Count_C','Count_ac']):
            if irow[icol]>0:
                objName=f"{icol}_Site_{group}"
                cmd.pseudoatom(objName, name=f"DUM{idx}", pos=[iX, iY, iZ], vdw=irow[icol], color="", label="")
                cmd.show("spheres", objName)
                cmd.color(greyColors[-idxCol-1], objName)
        for idxCol,icol in enumerate(['Count_Nitrogen','Count_an','Count_N']):
            if irow[icol]>0:
                objName=f"{icol}_Site_{group}"
                cmd.pseudoatom(objName, name=f"DUM{idx}", pos=[iX, iY, iZ], vdw=irow[icol], color="", label="")
                cmd.show("spheres", objName)
                cmd.color(blueColors[-idxCol-1], objName)
        for idxCol,icol in enumerate(['Count_Oxygen','Count_ao','Count_O']):
            if irow[icol]>0:
                objName=f"{icol}_Site_{group}"
                cmd.pseudoatom(objName, name=f"DUM{idx}", pos=[iX, iY, iZ], vdw=irow[icol], color="", label="")
                cmd.show("spheres", objName)
                cmd.color(redColors[-idxCol-1], objName)
    cmd.bg_color("white")            
    cmd.set("sphere_transparency", 0.3, "all") 
    objs=cmd.get_names("all")
    if f"Count_All_Site_{group}" in objs:
        cmd.center(f"Count_All_Site_{group}")
    print(objs)
    ObjStr=''
    for iObj in ["Count_All",'Count_Carbon','Count_C','Count_ac','Count_Nitrogen','Count_an','Count_N','Count_Oxygen','Count_ao','Count_O']:
        if f"{iObj}_Site_{group}" in objs:
            if ObjStr=='':
                ObjStr=f"{iObj}_Site_{group}"
            else:
                ObjStr+=f" or {iObj}_Site_{group}"
    # cmd.create(group,ObjStr)
    cmd.group(group,ObjStr,action='add')
    # try:       
    #     cmd.center("Count_All_Site")
    #     cmd.create(group,f"Count_All_Site or Count_Carbon or Count_C or Count_ac or Count_Nitrogen or Count_an or Count_N or Count_Oxygen or Count_ao or Count_O")
    # except Exception as e:
    #     print(e)
    
    cmd.save(savePseFile)   
    pymol.finish_launching()
    
def delete_dummy_direct(ifrag):
    ''' remove dummy function will change ionization state of fragment '''
    re_p=re.compile(r'\(\*\)|\*|\[\*\]')
    ifrag_nodummy = re.sub(re_p, '', ifrag)  ## remove dummy atoms
    return ifrag_nodummy
    
def get_fragPos(frag,sdf):
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
    if '*' in frag:
        frag=remove_dummy(frag)
        # frag=delete_dummy_direct(frag)
    frag=replace_nH(frag, nH=0, DH=1) ## Hydrogen isotope needs to be removed 
    smartMol=Chem.MolFromSmarts(replace_nH(frag))
    if mol==None:
        return atomProps
    matched=mol.GetSubstructMatches(smartMol)
    if len(matched)<1:
        return atomProps
    smartRI=get_atomRingInfo(frag)
    # print(f"smartRI= {smartRI}")
    molRI=get_atomRingInfo(mol)
    for imatch in matched:
        submolRI=molRI[list(imatch)]
        # print(f"submolRI= {submolRI}")
        diff=np.absolute(smartRI-submolRI).sum()
        # print(matched)
        if diff==0:
            for i in imatch:
                molName=mol.GetProp('_Name')
                ipos=conf.GetAtomPosition(i)
                iAtomicNum=molAtoms[i].GetAtomicNum()
                iFormalCharge=molAtoms[i].GetFormalCharge()
                iIsAromatic=int(molAtoms[i].GetIsAromatic())
                iIsInRing=int(molAtoms[i].IsInRing())
                atomProps.append([molName,i,ipos.x,ipos.y,ipos.z,iAtomicNum,iFormalCharge,iIsAromatic,iIsInRing])
    # print(f"{i}:  {atomProps}")
    return atomProps

def get_frag_gridPos(frag,sdfFolder,atmCount_file,grid_parameters_json):
    sdfFiles=glob.glob(f"{sdfFolder}/*.sdf")
    with open(grid_parameters_json, "r") as outfile:
        grid_parameters=json.load(outfile)   ## Load grid parameters
        xyzMin=grid_parameters['xyzMin']
        xyzMax=grid_parameters['xyzMax']
        cellWidth=grid_parameters['cellWidth']
    df_atmCount=pd.read_csv(atmCount_file)  ## Keep all the dataframe size are the same
    countCols=['Count_All','Count_Carbon','Count_ac','Count_C','Count_Nitrogen','Count_an','Count_N','Count_Oxygen','Count_ao','Count_O']  
    df_atmCount[countCols]=0
    df_atmCount[['X','Y','Z']]=df_atmCount[['X','Y','Z']].astype(str)
    df_atmCount.set_index(['X','Y','Z'], drop=False, inplace=True)
    atomPos_list=[]
    fragCount=0
    for isdf in sdfFiles:
        atomPos=get_fragPos(frag,isdf)
        if len(atomPos)>0:
            fragCount+=1
        for iatmPos in atomPos:
            atomPos_list.append(iatmPos)
    fragPercent=round(fragCount/len(sdfFiles),2)
    df_ActAtomPos=pd.DataFrame(atomPos_list,columns=['molName','atomIndex','X','Y','Z','AtomicNum','FormalCharge','IsAromatic','IsInRing'])
    ''' Calculate mean XYZ '''
    meanXYZ=df_ActAtomPos[['X','Y','Z']].mean()
    
    
    for idx,irow in df_ActAtomPos.iterrows():
        ixPos=int((irow['X']-xyzMin['X'])/cellWidth)*cellWidth+xyzMin['X']
        ixPos=str(ixPos)
        iyPos=int((irow['Y']-xyzMin['Y'])/cellWidth)*cellWidth+xyzMin['Y']
        iyPos=str(iyPos)
        izPos=int((irow['Z']-xyzMin['Z'])/cellWidth)*cellWidth+xyzMin['Z']
        izPos=str(izPos)
        try:
            tmp=df_atmCount.loc[(ixPos,iyPos,izPos),'Count_All']
        except Exception as e:
            print(e)
            continue
        df_atmCount.loc[(ixPos,iyPos,izPos),'Count_All']+=1
        if irow['AtomicNum']==6:
            df_atmCount.loc[(ixPos,iyPos,izPos),'Count_Carbon']+=1
            if irow['IsAromatic']==1:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_ac']+=1
            else:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_C']+=1
        if irow['AtomicNum']==7:
            df_atmCount.loc[(ixPos,iyPos,izPos),'Count_Nitrogen']+=1
            if irow['IsAromatic']==1:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_an']+=1
            else:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_N']+=1
        if irow['AtomicNum']==8:
            df_atmCount.loc[(ixPos,iyPos,izPos),'Count_Oxygen']+=1
            if irow['IsAromatic']==1:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_ao']+=1
            else:
                df_atmCount.loc[(ixPos,iyPos,izPos),'Count_O']+=1
                
    countCols=['Count_All','Count_Carbon','Count_ac','Count_C','Count_Nitrogen','Count_an','Count_N','Count_Oxygen','Count_ao','Count_O']
    countColMax=df_atmCount[countCols].max().max()
    for icountCol in countCols:
        # countColMax=df_ActAtomPos[icountCol].max()
        df_atmCount[icountCol]=df_atmCount[icountCol]/countColMax
    return  df_atmCount,meanXYZ,fragPercent


def frag_position_score(frag,motif,sdfFolder,atmCount_file,grid_parameters_json):
    # sdfFiles=glob.glob(f"{sdfFolder}/*.sdf")
    '''  Fragment position score relative to ref motif  '''
    countCols=['Count_All','Count_Carbon','Count_ac','Count_C','Count_Nitrogen','Count_an','Count_N','Count_Oxygen','Count_ao','Count_O']
    df_motif,_,_=get_frag_gridPos(frag=motif,sdfFolder=sdfFolder,atmCount_file=atmCount_file,grid_parameters_json=grid_parameters_json) #(frag,sdfFolder,atmCount_file,grid_parameters_json)
    df_frag,_,_=get_frag_gridPos(frag=frag,sdfFolder=sdfFolder,atmCount_file=atmCount_file,grid_parameters_json=grid_parameters_json) 
    score=(df_motif[countCols].astype(float)*df_frag[countCols].astype(float)).sum().sum()
    return score