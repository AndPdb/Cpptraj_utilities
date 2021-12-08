#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import re

# FUNCTIONS FOR INPUT CREATION

def input_intrahb(parm7, TRAJ_PATH, len_res):
    """Create input for cpptraj hbond analyses. $parm7 is the parm path, $TRAJ_PATH is the path for the trajectories, $len_res is the total number of protein residues"""
        
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond Bridges :1-{len_res} solventdonor :WAT solventacceptor :WAT@O \
solvout solvent_avg.dat bridgeout bridge.dat
hbond Protein :1-{len_res} out hb.dat avgout avghb.dat series uuseries hb_series.dat
run
lifetime Protein[solutehb] out Protein.lifetime.dat
runanalysis
quit"""
    return txt
    
def input_intrasb(parm7, TRAJ_PATH, len_res):
    """Create input for cpptraj saltbridges analyses. $parm7 is the parm path, $TRAJ_PATH is the path for the trajectories, $len_res is the total number of protein residues"""
    
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond donormask :ARG,LYS@N*&!@N acceptormask :ASP,GLU@O*&!@O dist 4 angle -1 out saltbridges.dat avgout avg_sb.dat series uuseries sb_series.dat
run
quit"""
    return txt
    
def input_intrahp(parm7, TRAJ_PATH, len_res):
    """Create input for cpptraj hydropobics contact analyses. $parm7 is the parm path, $TRAJ_PATH is the path for the trajectories, $len_res is the total number of protein residues"""
    
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
nativecontacts :1-{len_res}@/C&!@CA,C distance 4 resoffset 1 writecontacts contacts.dat out numcontacts.dat first savenonnative series seriesout series.dat seriesnnout nnseries.dat
run
quit"""
    return txt
    
def input_lighb(parm7, TRAJ_PATH, len_res, ligname):
    """Create input for cpptraj hbond analyses. $parm7 is the parm path, $TRAJ_PATH is the path for the trajectories, $len_res is the total number of protein residues, $ligname is the ligand name in the PDB"""
    
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond Bridges :1-{len_res},:{ligname} solventdonor :WAT solventacceptor :WAT@O \
solvout solvent_ligavg.dat bridgeout ligbridge.dat
hbond Substrate :1-{len_res},:{ligname} out hb.dat avgout avg_hb.dat series uuseries hb_series.dat nointramol
run
lifetime Substrate[solutehb] out Substrate.lifetime.dat
runanalysis
quit
"""
    return txt
    
def input_ligsb(parm7, TRAJ_PATH, len_res, ligA=None, ligD=None):
    """Create input for cpptraj salt bridges analyses. $parm7 is the parm path, $TRAJ_PATH is the path for the trajectories, $len_res is the total number of protein residues, $ligname is the ligand name in the PDB. $dlig is the mask for the ligand group positively charged; $alig is the mask for the ligand group negatively charged."""
    if ligA == None:
        txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond donormask {ligD} acceptormask :ASP,GLU@O* dist 4 angle -1 out num_sb1.dat avgout avg_sb1.dat series uuseries uuhbonds_1.dat
run
quit
"""
    elif ligD == None:
        txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond donormask :ARG,LYS@N* acceptormask {ligA} dist 4 angle -1 out num_sb2.dat avgout avg_sb2.dat series uuseries uuhbonds_2.dat
run
quit
"""
    else:  
        txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond donormask {ligD} acceptormask :ASP,GLU@O* dist 4 angle -1 out num_sb1.dat avgout avg_sb1.dat series uuseries uuhbonds_1.dat
hbond donormask :ARG,LYS@N* acceptormask {ligA} dist 4 angle -1 out num_sb2.dat avgout avg_sb2.dat series uuseries uuhbonds_2.dat
run
quit
"""
    return txt
    
def input_lighp(parm7, TRAJ_PATH, len_res, ligname):
    """Create input for cpptraj hydrophobics contact analyses. $parm7 is the parm path, $TRAJ_PATH is the path for the trajectories, $len_res is the total number of protein residues, $ligname is the ligand name in the PDB"""
    
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
nativecontacts :{ligname}@/C :1-{len_res}@/C&!@CA,C distance 4 resoffset 1 writecontacts contacts.dat out numcontacts.dat first savenonnative series seriesout series.dat seriesnnout nnseries.dat
run
quit"""
    return txt


# FUNCTIONS FOR OUTPUT ANALYSES



def intra_hbonds(avghb, frac=0.4):
    """Return a DataFrame with HBonds with avg duration higher than a percentage (default 0.4)"""
    df = pd.read_csv(avghb, delim_whitespace=True)
    hb = df[df.Frac >= frac]
    hb = hb.rename(columns = {'#Acceptor':'Acceptor'})
    hb = hb[(hb['Acceptor'].str.endswith('@O') == False) & (hb['Donor'].str.endswith('@N') == False)]
    return hb.sort_values('Frac', ascending=False)

def lig_hbonds(avghb, lig, frac=0.4):
    """Return a DataFrame with the ligand HBonds with avg duration higher than a percentage (default 0.4)"""
    df = pd.read_csv(avghb, delim_whitespace=True)
    hb = df[df.Frac >= frac]
    hb = hb.rename(columns = {'#Acceptor':'Acceptor'})
    hb_lig1 = hb[hb['Acceptor'].str.contains(lig)]
    hb_lig2 = hb[hb['Donor'].str.contains(lig)]
    hb_lig = pd.concat([hb_lig1, hb_lig2], ignore_index=True)
    return hb_lig.sort_values('Frac', ascending=False)

def lig_saltbridges(path1, path2=None):
    first = pd.read_csv(path1,  delim_whitespace=True)
    try:
        second = pd.read_csv(path2,  delim_whitespace=True)
    except:
        second = None
    if second == None:
        df = first
    else:
        df = pd.concat([first, second], ignore_index=True)
    df = df.sort_values('Frac', ascending=False)
    return df.reset_index(drop=True)

def intra_saltbridges(avgsb, frac=0.2):
    """Return a DataFrame with SaltBridges with avg duration higher than a percentage (default 0.4)"""
    #create dataframe
    df = pd.read_csv(avgsb, delim_whitespace=True)
    df = df[df.Frac >= frac]
    df = df.rename(columns = {'#Acceptor':'Acceptor'})
    df = df.loc[~df['Acceptor'].str.endswith('@O')]
    df = df.loc[~df['Donor'].str.endswith('@N')]
    df = df.sort_values('Frames', ascending=False)
    df = df.reset_index(drop=True)
    return df

def solv_bridges(path):
    """Return a DataFrame of the solvent bridges established"""
    #create dataframe
    df_sb = pd.read_csv(path, delimiter=',', header=None, skiprows=1, names=['Residues', 'Frames'])
    df_sb['Frames'] = df_sb['Frames'].str.replace('\sframes.', '')
    df_sb['Residues'] = df_sb['Residues'].str.replace('Bridge Res', '')
    #df_sb = df_sb.loc[df_sb['Residues'].str.contains(lig)]
    df_sb['Frames'] = df_sb['Frames'].astype(float)
    #df_sb['Frac'] = df_sb['Frames']/totFrames
    #df_sb = df_sb[df_sb.Frac >= frac]
    return df_sb.sort_values('Frames', ascending=False)

def contacts(contacts, notnative, frac=0.2, avg=4.5):
    """General parser for contacts.dat with variable Frac and Avg"""
    df = pd.read_csv(contacts, delim_whitespace=True, skiprows=[0,1,notnative])
    df['Frac.'] = pd.to_numeric(df['Frac.'], errors='coerce')
    df['Nframes'] = pd.to_numeric(df['Nframes'], errors='coerce')
    df['Avg'] = pd.to_numeric(df['Avg'], errors='coerce')
    df['Stdev'] = pd.to_numeric(df['Stdev'], errors='coerce')
    df = df.rename(columns = {'Frac.':'Frac'})
    df = df[df['Frac'].notna()]
    df = df[df.Frac >= frac]
    df = df[df.Avg <= avg]
    df = df.drop(['#'], axis=1)
    #Drop row with intraresidue contact
    df = df.reset_index(drop=True)
    todrop = []
    for i, n in df.iterrows():
        contact = re.split((r'@|:'), df.iloc[i]['Contact'])
        if contact[1] == contact[3]:
            todrop.append(i)        
    df.drop(todrop, inplace=True)
    
#     ls_cont = []
#     contatti = list(df['Contact'])
#     for entry in contatti:
#         ls_cont.append(entry.split('@')[1].split(':')[1])
    #return set(ls_cont), df.sort_values('Nframes', ascending=False)
    
    return df.sort_values('Nframes', ascending=False)


def contacts_interhelices(contacts, frac=0.5, cutoff=5):
    """Return a DataFrame if the contacts are longer than frac (def 0.5) and with avg distance lower than cutoff (def 5),
    and if more than 4 residues are between the residues in contact"""
    df = pd.read_csv(contacts, delim_whitespace=True, skiprows=2)
    df['Frac.'] = pd.to_numeric(df['Frac.'], errors='coerce')
    df['Nframes'] = pd.to_numeric(df['Nframes'], errors='coerce')
    df['Avg'] = pd.to_numeric(df['Avg'], errors='coerce')
    df['Stdev'] = pd.to_numeric(df['Stdev'], errors='coerce')
    df = df.rename(columns = {'Frac.':'Frac'})
    df = df[df['Frac'].notna()]
    df = df[df.Frac >= frac]
    df = df[df.Avg <= cutoff]
    df = df.loc[~df['Contact'].str.contains('#')]
    df = df.reset_index(drop=True)
    df = df.drop(['#'], axis=1)
    new = pd.DataFrame(columns=df.columns)

    lista = []
    fails=[]
    for i, n in df.iterrows():
        residues = re.split(':|@', df.iloc[i]['Contact'])
        try:
            if abs(int(residues[1])
                   -int(residues[3])) > 5:
                if 'H' not in residues[2] and 'H' not in residues[4]:
                    lista.append(df.iloc[i])
                    new.loc[i]=df.iloc[i]

        except IndexError:
            fails.append(df.iloc[i])
    return new

def rmsd(rmsd):
    """Return a plot of rmsd.dat"""
    df = pd.read_csv(rmsd, delim_whitespace=True, skiprows=1, header=None, names=['Frame', 'RMSD'])
    return df


# PLOT TIMESERIES

def dfseries_sb(path):
    
    try:
        dfsb1_1 = pd.read_csv(f"{path}/uuhbonds_1.dat", delim_whitespace=True)
    except FileNotFoundError:
        dfsb1_1 = pd.DataFrame()
    try:
        dfsb1_2 = pd.read_csv(f"{path}/uuhbonds_2.dat", delim_whitespace=True)
    except FileNotFoundError:
        dfsb1_2 = pd.DataFrame()
        
    if (dfsb1_1.empty == False) and (dfsb1_2.empty == True):
        dfsb1_1 = dfsb1_1.rename(columns={"#Frame": "Frame"})
        dfsb1 = dfsb1_1
    elif (dfsb1_1.empty == True) and (dfsb1_2.empty == False):
        dfsb1_2 = dfsb1_2.rename(columns={"#Frame": "Frame"})
        dfsb1 = dfsb1_2
    elif (dfsb1_1.empty == False) and (dfsb1_2.empty == False):
        dfsb1_1 = dfsb1_1.rename(columns={"#Frame": "Frame"})
        dfsb1_2 = dfsb1_2.drop(['#Frame'], axis=1)
        dfsb1 = pd.concat([dfsb1_1, dfsb1_2], axis=1)
    else:
        dfsb1=dfsb1_1
        print("Series not present")
    return dfsb1



def dfseries_hp(path):
    dfhp1=True
    dfhp2=True
    try:
        dfhp1_1 =pd.read_csv(f"{path}/series.dat", delim_whitespace=True)
        dfhp1_1 = dfhp1_1.rename(columns={"#Frame": "Frame"})
    except:
        dfhp1 = False
    
    try:
        dfhp1_2 = pd.read_csv(f"{path}/nnseries.dat", delim_whitespace=True)
        dfhp1_2 = dfhp1_2.rename(columns={"#Frame": "Frame"})
    except:
        dfhp2 = False

    if dfhp1 == False:
        dfhp1 = dfhp1_2
    
    elif dfhp1_2 == False:
        dfhp1 = dfhp1_1
    
    else:
        dfhp1_2 = dfhp1_2.drop(['Frame'], axis=1)
        dfhp1 = pd.concat([dfhp1_1, dfhp1_2], axis=1)
   
    return dfhp1



def df_4catplot(df, ligname, interaction, resinfo=None):
    '''Convert normal dataframe(from hbonds, salt_bridge, hydrophobic analyses) to df useful for catplot. Needs type of interction and a df "resinfo" describing resname and numbers (from resinfo command in cpptraj)'''
    
    #Find all the interacting residues
    prot_res= []
    columns = list(df.columns)

    #if series from nativecontact analysis change column names
    if interaction == 'Hydrophobic':
        for col in columns:
            if col != '#Frame' and col != 'Frame':
                lig = col.split('_')[0]
                lignum = lig.split('@')[0]
                latm = lig.split('@')[1]
                res = col.split('_')[1]
                resnum = res.split('@')[0]
                ratm = res.split('@')[1]
                row = resinfo.loc[resinfo['#Res'] == int(resnum[1:])]
                resname = row.iloc[0]['Name']
                row = resinfo.loc[resinfo['#Res'] == int(lignum[1:])]
                lignme = row.iloc[0]['Name']
                new = f'{lignme}_{lignum[1:]}@{latm}-{resname}_{resnum[1:]}@{ratm}'
                df = df.rename(columns={col: new})
                columns = list(df.columns)
    #end    

    
    for col in columns:
        if col != '#Frame' and col != 'Frame':
            res1 = col.split('-')[0]
            res1 = res1.split('@')[0]
            res2 = col.split('-')[1]
            res2 = res2.split('@')[0]
            if res1 != ligname:
                prot_res.append(res1)
#                 new_clmns.append(res1.split('_')[0]+res1.split('_')[1])
            else:
                prot_res.append(res2)
#                 new_clmns.append(res2.split('_')[0]+res2.split('_')[1])





#     new_clmns = list(set(new_clmns))
    prot_res = list(set(prot_res))
    #end

    res_dict = {}

    df_series = pd.DataFrame(columns=['Frame']+list(prot_res))
    if '#Frame' in columns:
        df_series['Frame'] = df['#Frame']
    else:
        df_series['Frame'] = df['Frame']
    df_series = df_series.fillna(0)

    #Iterate over interacting residues and fill a new df with bool for every frame
    for res in prot_res:
        df_res = df.filter(regex=res, axis=1)
        for idx, row in df_res.iterrows():
            if (df_res.loc[idx] == 1).any():
                df_series.loc[idx, res] = True
            else:
                df_series.loc[idx, res] = False
    #end

    
    #Create the final df for the cat plot with seaborn
    df_series = df_series.reset_index()
    df_series = pd.melt(df_series, id_vars=["Frame"], var_name=["Residue"])
    df_series = df_series[df_series["value"] != False]
    df_series.reset_index(inplace=True, drop=True)
    df_series = df_series[df_series.Residue != 'index']
    df_series['interaction'] = interaction
    return df_series

