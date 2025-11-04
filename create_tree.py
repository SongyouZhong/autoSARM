

''' Set the environment  '''
from my_toolset.drawing_utils import show_mols as display_mols
import pandas as pd
import numpy as np
import os,sys
from pathlib import Path
import glob
from my_toolset.my_utils import canonic_smiles,df_valid
import graphviz
from functools import partial
from pandarallel import pandarallel
import argparse
import time
n_jobs=30
pandarallel.initialize(nb_workers=n_jobs)

sys.path.append("../")
from utils.tree_utils import replace_nH,find_root,find_children,find_parents,node,real_sonNode,show_cpd,match_frag,remove_ionizationForm
from utils.grid_pos_utils import extract_sdf,get_atomPos,grid_atom,create_grid_pse,get_frag_gridPos,frag_position_score
from utils.common_utils import mapper,csvToExcel,get_mol,float_row,remove_dummy # ,df_match



jupyterMode=0


fragment_core = "*CN1CCC(c2ccc3[nH]c(-c4cc(CO*)c5ncnn5c4)c(C(C)C)c3c2)CC1"
rootTitle = "Table_775_combine"
workFolder = "/home/jay/Vue-molOpt/API/RunFolder/SARM/TLR78-selectivity_2025-01-13-18-20-13"
highlightDict = "[{'col':'TLR7IC50', 'type':'means', 'relation':'<', 'value':10}]"  ### type: means,stds,medians; relation: < = >; 

class argNamespace:
    ''' for test in jupyter notebook  '''
    def __init__(self):  
        self.protein = '/mnt/public-bg6/jay.zhang/Codes-public/Vina/Test/8jzx.pdb'  # 'self'引用的是类实例自身  
        self.ligand = '/mnt/public-bg6/jay.zhang/Codes-public/Vina/Test/c5.pdb' 
        self.addH = 1


def main(args):
    fragment_core = args.fragment_core
    rootTitle = args.rootTitle   ## 'Table_775_combine' table name
    workFolder = os.path.abspath(args.workFolder)
    highlightDict = eval(args.highlightDict)
    max_level = args.maxLevel
    treeContent = eval(args.treeContent)

    actCols = []
    for icol in highlightDict:
        actCols.append(icol['col'])

    treePath=Path(f"{workFolder}/Trees/FragTree_{rootTitle}").absolute()
    treePath.mkdir(exist_ok=True, parents=True)

    if 1:
        roots=[fragment_core]
        dfList = []
        if 'double-cut' in treeContent:
            df_tmp=pd.read_csv(f"{workFolder}/SAR_Tables/Combine_Table_info.csv")
            df_tmp["SMILES"] = df_tmp["Key2"]
            dfList.append(df_tmp)
        if 'single-cut' in treeContent:
            df_tmp=pd.read_csv(f"{workFolder}/SAR_Tables/singleCut_Table_info.csv")
            df_tmp["SMILES"] = df_tmp["Key2"]
            dfList.append(df_tmp)
        if 'whole-compound' in treeContent:
            df_tmp=pd.read_csv(f"{workFolder}/input.csv")
            df_tmp["SMILES"] = df_tmp["smiles"].parallel_apply(canonic_smiles)
            df_tmp["Items_count"] = 1
            dfList.append(df_tmp)
        df_table=pd.concat(dfList)

        opfile=open(f'{treePath.as_posix()}/Combine_Table_info_Tree_{rootTitle}.txt','w')
        ## active dataset
        df_act=pd.read_csv(f"{workFolder}/input.csv")
        # df_act["actValue"]=9-df_act["actValue"].apply(np.log10)
        df_act=float_row(df_act, cols=actCols,dropna=False)
        actCols=actCols
        if 'SMILES' not in df_act.columns and 'smiles' in df_act.columns:
            df_act['SMILES'] = df_act['smiles']
        df_act['Cano_SMILES'] =df_act['SMILES'].parallel_apply(canonic_smiles)
        df_act = df_act.dropna(subset=['Cano_SMILES'])

        print('#'*5+' Working on create the tree!') 
        # df_table["SMILESnoDummy"]=df_table["SMILES"].parallel_apply(remove_dummy)
        df_table['Matched']=df_table.parallel_apply(lambda x:match_frag(x['SMILES'],ismarts=fragment_core),axis=1)
        df_table=df_table[df_table['Matched']==1]
        # df_table["SMILESnoDummy"]=df_table["SMILESnoDummy"].apply(remove_ionizationForm)


        SMILESnoDummy=df_table["SMILES"].to_list()
        df_table['index']=df_table["SMILES"]
        df_table.set_index('SMILES',inplace=True,drop=False)
        df_table['Items_count']=df_table['Items_count'].astype(int)
        # SMILESnoDummy=list(set(SMILESnoDummy))
        # print('#'*5 + "remove_ionizationForm")
        # SMILESnoDummy=remove_ionizationForm(SMILESnoDummy)
        print('#'*5 + "find children")
        children_dict=find_children(SMILESnoDummy)
        print('#'*5 + "find parents")
        parent_dict=find_parents(SMILESnoDummy)

        print('#'*5+' Working on building the tree!')
        tree_list=[]
        tree_smile_list=[]
        for itree,ismi in enumerate(roots):
            opfile.write('\n\n\n'+'-'*8+f" Tree {itree}: {ismi} "+'-'*8+'\n')
            inode=node(ismi)
            # inode.print_tree()
            next_level_inodes=[inode]
            while len(next_level_inodes)>0:
                for jdx,jnode in enumerate(next_level_inodes):
                    children_smis=children_dict[jnode.SMILES]
                    if len(children_smis)>0:
                        for ichild_smi in children_smis:
                            if real_sonNode(ichild_smi,children_smis,parent_dict):
                                next_level_inodes[jdx].add_child(node(ichild_smi))
                tmp_next_level_inodes=[]
                for jnode in next_level_inodes:
                    tmp_next_level_inodes.extend(jnode.children)
                next_level_inodes=tmp_next_level_inodes    
            itree_smiles=inode.print_tree(df_table=df_table,opfile=opfile)
            tree_smile_list.append(itree_smiles)
            tree_list.append(inode)
        opfile.close()


    ''' View the tree '''
    if 1:
        os.chdir(treePath)
        Path(f'Images').mkdir(exist_ok=True, parents=True)
        
        print('#'*5+' Working on plot the tree!')
        for itree,itree_smiles in enumerate(tree_smile_list):  ## enumerate root
            d = graphviz.Digraph(filename=f'{rootTitle}')
            d.node_attr["shape"] = "plaintext" # "box"
            d.node_attr["fixedsize"] = 'true'
            d.node_attr["height"] = '1'
            d.node_attr["width"] = '1'
            d.node_attr["label"] = ''

            # 全局设置（整个图的默认字体大小）
            d.attr(fontsize='20')  # 图的整体字体大小（影响标题等）
            # 设置所有节点的默认属性（包括字体大小）
            d.node_attr.update(fontsize='18', fontname='SimHei')  # 节点默认字体14pt，支持中文
            # 设置所有边的默认属性
            d.edge_attr.update(fontsize='18')  # 边标签默认字体12pt


            levelDummyNodeList = []   ## the list of dummy nodes all graph level
            for ilevel_idx,ilevel_smiles in enumerate(itree_smiles):     ###     遍历每一层 
                # ilevel_smiles=[fetch_table_info_local(df_table,ismi) for ismi in ilevel_smiles]
                # print(ilevel_smiles)
                if ilevel_idx > max_level:
                    break
                with d.subgraph() as s:
                    s.attr(rank='same')
                    count=-1
                    levelNodeDict = {}   ## for saving the cpds have been displayed
                    
                    for inode,cpds_node in enumerate(ilevel_smiles):
                        cpdsNode_countList=[]
                        count_tmp=count
                        for icpd in cpds_node:
                            count_tmp+=1
                            cpdsNode_countList.append([f"Images/L{ilevel_idx}_{count_tmp}.png",icpd])
                        # show_cpd_p=partial(show_cpd, df_table=df_table, df_act=df_act, actCols=actCols, highlightDictList=highlightDict)
                        # cpdsNode_countList_res=mapper(40)(show_cpd_p, cpdsNode_countList)

                        '''  For debugging purpose '''
                        cpdsNode_countList_res =[]
                        for icpd_countList in cpdsNode_countList:
                            print(f"Processing icpd: {icpd_countList}")
                            print(f"actCols: {actCols}")
                            icpdsNode_countList_res = show_cpd(icpd_countList, df_table=df_table, df_act=df_act, actCols=actCols, highlightDictList=highlightDict)
                            cpdsNode_countList_res.append(icpdsNode_countList_res)

                        

                        for img_icpd in cpdsNode_countList_res:
                            # print(f"icpd={icpd}")
                            count+=1
                            if img_icpd==None:
                                continue
                            # icpd_dummy=fetch_table_SMILES(df_table,icpd)
                            # print(icpd)
                            # img=show_mols([icpd_dummy],legends=[fetch_table_info(df_table,icpd)])
                            # img.save(f"Images/L{ilevel_idx}_{count}.png")
                            # img=show_mols_cairo(icpd_dummy,fetch_table_info(df_table,icpd)+get_activity_info(df_act,icpd))
                            # with open(f"Images/L{ilevel_idx}_{count}.png", 'wb') as f:
                            #     f.write(img)
                            
                            cpdSmi =img_icpd[1]
                            if cpdSmi not in levelNodeDict.keys():
                                node_label = f"L{ilevel_idx}_{count}"
                                s.node(node_label, image=img_icpd[0])
                                
                                if ilevel_idx >0:
                                    d.edge(f'L{ilevel_idx-1}_{inode}', node_label, penwidth='0.2', arrowsize='0.2')
                                levelNodeDict[cpdSmi] = node_label
                            else:
                                node_label = f"L{ilevel_idx}_{count}"
                                s.node(node_label)   ## add dummy node to keep the node order
                                levelDummyNodeList.append(node_label)
                                node_label = levelNodeDict[cpdSmi]

                                if ilevel_idx>0 and f'L{ilevel_idx-1}_{inode}' not in levelDummyNodeList:
                                    d.edge(f'L{ilevel_idx-1}_{inode}', node_label, penwidth='0.2', arrowsize='0.2')
            # d.view()
            d.render(view=False)
            # time.sleep(1)      
            # os.system(f"mv *.gv *.pdf *.txt Images FragTree_{rootTitle}")  



def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragment_core", help="the fragment core file",
                        default='')
    parser.add_argument("--rootTitle", help="the title of pdf as rootTile", default='')
    parser.add_argument("--maxLevel", help="max level of the SAR tree", default=5, type=int)
    parser.add_argument("--workFolder", help="the workfolder of sar project", default='')
    parser.add_argument("--treeContent", help="the content in the tree: ['double-cut','','']", type=str, default="['double-cut']")

    parser.add_argument("--highlightDict", help="the dict to control the highlight", default='', type=str)

    args = parser.parse_args()
    return args


if jupyterMode:   ### test without input of args run as juputer notebook
    args=argNamespace()
    main(args)

# %%

if __name__ == "__main__":
    args = get_parser()
    main(args)
