from schrodinger import structure
from schrodinger.structutils import measure
from schrodinger.structutils import analyze
from schrodinger.test import mmshare_data_file
from pathlib import Path
import argparse
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser()
    ''' The section of docking parameters '''
    parser.add_argument("--saveDir", help="root folder of mol2 files", type=str, default='')
    parser.add_argument("--format", help="file formate to export: ['sd','mol2']", type=str, default='sd')
    parser.add_argument("--maegz", help="maegz file to unzip", type=str, default="")
    parser.add_argument("--scoreLimit", help="the minimum docking score", type=float, default=6)
    args = parser.parse_args()
    return args
args = get_parser()

saveDir=args.saveDir
savePath=Path(saveDir)
resList=[]

with structure.StructureReader(args.maegz) as reader:
    stCount_Dict={}
    for st in reader:
        title=st.title
        energy=st.property.get('r_i_docking_score')
        if not energy:
            energy=st.property.get('r_sd_docking_score')
        if energy != None:
            energy=round(-energy/1.36357,2)
            # print(energy)
            if energy<args.scoreLimit:
                continue
            print(f"Processing {title}")
            if title not in stCount_Dict.keys():   
                stCount_Dict[title]=1
            else:
                stCount_Dict[title]+=1
            if args.format=='sd':
                mol2_file=f"{title}-{stCount_Dict[title]}.sdf"
            if args.format=='mol2':
                mol2_file=f"{title}-{stCount_Dict[title]}.mol2"
            # ligPath=savePath.joinpath(title)
            # if not ligPath.exists():
            #     ligPath.mkdir(exist_ok=True, parents=True)
            st.write(savePath.joinpath(mol2_file), format=args.format)
            resList.append([f"{title}-{stCount_Dict[title]}",energy])
dfRes=pd.DataFrame(resList, columns=['fullName','Energy'])
dfRes.to_csv(savePath.joinpath('dock_score.csv'), index=None)
        
        
        