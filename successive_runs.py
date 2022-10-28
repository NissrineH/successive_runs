#Packages
import PySimpleGUI as sg                        
from pymol import cmd
import sys, os, random
from openbabel import pybel
from rdkit import Chem
import py3Dmol
import warnings
warnings.filterwarnings("ignore")
import subprocess
import pandas as pd
import csv

#Working folder
os.chdir('/home/nissrine/Jupyter_Dock/GUI/successive_runs/')
#Functions

#Get pdb format of the receptor 
def load_receptor(receptor):
    pdb_lists = [receptor] 
    for x in pdb_lists:
        cmd.fetch(code=x,type='pdb')
        cmd.select(name='Prot',selection=select)
        cmd.save(filename='rec.pdb',format='pdb',selection='Prot')
        cmd.delete('all')
        
# Receptor preparation: addH and charges
def receptor_preparation():    
    os.system("../bin/lepro_linux_x86 rec.pdb")
    os.rename('pro.pdb','rec_clean_H.pdb')

#Get ligands
def load_ligands(ligands, smiles):
    df = pd.read_excel(ligands)
    df_= df.sample(int(round((len(df)*ratio)/100)))
    
    smiles = df_['smiles'].values.tolist()
    id_smiles = df_['id'].values.tolist()
    return(df_, smiles, id_smiles)

# Ligands preparation: addH and charges
def ligands_preparation(smiles):
    out=pybel.Outputfile(filename='InputMols.mol2',format='mol2',overwrite=True)
    for index,smi in enumerate(smiles):
        mol=pybel.readstring(string=smi,format='smiles')
        mol.title='mol_'+str(index)
        mol.make3D('mmff94s')
        mol.localopt(forcefield='mmff94s', steps=500)
        out.write(mol)
    out.close()

# Docking box definition  
def getbox(selection='sele', extending = 6.0, software='vina'):
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)
    minX = minX - float(extending)
    minY = minY - float(extending)
    minZ = minZ - float(extending)
    maxX = maxX + float(extending)
    maxY = maxY + float(extending)
    maxZ = maxZ + float(extending)
    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX =  (maxX + minX)/2
    CenterY =  (maxY + minY)/2
    CenterZ =  (maxZ + minZ)/2
    cmd.delete('all')
    if software == 'vina':
        return {'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}
    elif software == 'ledock':
        return {'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ}
    elif software == 'both':
        return ({'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}),({'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ})
    else:
        print('software options must be "vina", "ledock" or "both"')

def box():    
    cmd.load(filename='rec_clean_H.pdb',format='pdb',object='prot') 
    center, size=getbox(selection='prot',extending=6.0,software='vina')
    cmd.delete('all')
    return (center, size)

    
def dock(center, size, num_exec):      
    random_seed=[]
    for i in range(num_exec):
        rs=random.randint(0, 1000000)
        commande="../bin/smina -r rec_clean_H.pdb -l InputMols.mol2 -o {6}.sdf --center_x {0} --center_y {1} --center_z {2} --size_x {3} --size_y {4} --size_z {5} --exhaustiveness 8 --num_modes 1 --seed {6}".format(center['center_x'], center['center_y'], center['center_z'], size['size_x'], size['size_y'], size['size_z'], rs)
        os.system(commande)
        random_seed.append(str(rs))
    return (random_seed)


def get_scores(random_seed):
    docking_scores=[]
    for i in range(len(random_seed)):
        poses=Chem.SDMolSupplier(str(random_seed[i])+'.sdf',True)
        sc=[]
        for p in list(poses)[::]:
            pose_1=Chem.MolToMolBlock(p)
            sc.append(p.GetProp('minimizedAffinity'))
        print(sc)
        docking_scores.append(sc)
    return (docking_scores)


def df_scores(docking_scores, random_seed, id_smiles):
    df=pd.DataFrame(docking_scores).T
    df.columns=random_seed
    df = df.astype(float)
    df['ids']=id_smiles
    first_column = df.pop('ids')
    df.insert(0, 'ids', first_column)
    best_scores = df.min(axis=1)
    df['best_scores']=best_scores
    best_value = df['best_scores'].min()
    ids = df.index[df['best_scores'] == best_value].tolist()
    df.to_csv("docking_results.csv")
    best_conf = df.loc[:, df.columns!='best_scores'].apply(lambda row: row[row==best_value].index.values, axis=1)
    for i in range(len(best_conf)):
        if ((best_conf).values[i]!='[]'):
            rs_best_conf=int(best_conf[i])
    
    return(df, ids[0], best_value, rs_best_conf)



def docking(receptor, ligands, ratio, num_exec):
    smiles=[]
    load_receptor(receptor)
    rec_clean=receptor_preparation() 
    df_ratio, sm, id_sm=load_ligands(ligands, smiles)
    ligands_preparation(sm)
    center, size=box()
    rs=dock(center, size, num_exec)
    docking_scores=get_scores(rs)
    scores, ids, best_value, rs_best_conf=df_scores(docking_scores, rs, id_sm)
    return(scores, ids, best_value, rs_best_conf)



    
# User's GUI    
    
sg.theme('SystemDefaultForReal')
    
# Define the window's contents
layout = [  [sg.Text("Molecular Docking")],    
            [sg.Text('Receptor ID'), sg.Input(key='receptor')],
            [sg.Text('Receptor\'s selection'), sg.Input(key='selection')],
            [sg.Text('Path to ligands file'), sg.In(key='ligands'), sg.FileBrowse()],
            [sg.Text('Number of runs'), sg.Input(key='num_exec')],
            [sg.Text('Ratio'), sg.Slider(range=(0,100.00),orientation='h', resolution=0.001, default_value=0.01,key='ratio')],
            [sg.Button('Get Docking Results')] ]

# Create the window
window = sg.Window('Molecular Docking', layout)    

# Display and interact with the Window
event, values = window.read()                   
args = values
receptor=args['receptor']
select=args['selection']
ligands=args['ligands']
ratio=args['ratio']
num_exec=args['num_exec']

scores, ids, best_value, rs_best_conf = docking(receptor, ligands, ratio, int(num_exec))
poses=Chem.SDMolSupplier('{0}.sdf'.format(rs_best_conf),True)
#for p in list(poses)[ids[0]]:
p=poses[ids]
pose_1=Chem.MolToMolBlock(p)
sg.popup(' Molecule id: {}'.format(scores['ids'][ids]),' Random Seed: {}'.format(rs_best_conf),' Score      : {}'.format(p.GetProp('minimizedAffinity')))
