#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import re

# FUNCTIONS FOR INPUT CREATION

def input_intrahb(parm7, TRAJ_PATH, len_res, OUT_PATH):
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
    
def input_intrasb(parm7, TRAJ_PATH, len_res, OUT_PATH):
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.*
autoimage
hbond donormask :ARG,LYS@N* acceptormask :ASP,GLU@O* dist 4 angle -1 out saltbridges.dat avgout avg_sb.dat series uuseries sb_series.dat
run
quit"""
    return txt
    
def input_intrahp(parm7, TRAJ_PATH, len_res, OUT_PATH):
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.*
autoimage
nativecontacts @/C&!@CA,C writecontacts contacts.dat out numcontacts.dat first savenonnative series seriesout series.dat seriesnnout nnseries.dat
run
quit"""
    return txt
    
def input_lighb(parm7, TRAJ_PATH, len_res, ligname, OUT_PATH):
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
hbond Bridges :1-{len_res},:{ligname} solventdonor :WAT solventacceptor :WAT@O \
solvout solvent_ligavg.dat bridgeout ligbridge.dat
hbond Substrate :1-{len_res},:{ligname} out {ligname}hb.dat avgout avg_{ligname}hb.dat series uuseries {ligname}hb_series.dat nointramol
run
lifetime Substrate[solutehb] out Substrate.lifetime.dat
runanalysis
quit
"""
    return txt
    
def input_ligsb(parm7, TRAJ_PATH, len_res, ligname, OUT_PATH):
    txt = f"""parm {parm7}
trajin {TRAJ_PATH}mdcrd.md*
autoimage
hbond donormask :{ligname}@N* acceptormask :ASP,GLU@O* dist 4 angle -1 out {ligname}_sb1.dat avgout avg_{ligname}sb1.dat series uuseries {ligname}_uuhbonds_1.dat
hbond donormask :ARG,LYS@N* acceptormask :{ligname}@O* dist 4 angle -1 out {ligname}_sb2.dat avgout avg_{ligname}sb2.dat series uuseries {ligname}_uuhbonds_2.dat
run
quit
"""
    return txt
    
def input_lighp(parm7, TRAJ_PATH, len_res, OUT_PATH):
    txt = f"""parm ../../slc25a29_lys1.parm7
trajin {TRAJ_PATH}/mdcrd.md*
autoimage
nativecontacts :{ligname}@/C :1-{len_res}@/C&!@CA,C writecontacts {ligname}_hpcontacts.dat out {ligname}_numhpcontacts.dat first savenonnative series seriesout {ligname}_series.dat seriesnnout {ligname}_nnseries.dat
run
quit"""
    return txt


# FUNCTIONS FOR OUTPUT ANALYSES



def hbonds(avghb, frac=0.4):
    """Return a DataFrame with HBonds with avg duration higher than a percentage (default 0.4)"""
    df = pd.read_csv(avghb, delim_whitespace=True)
    hb = df[df.Frac >= frac]
    hb = hb.rename(columns = {'#Acceptor':'Acceptor'})
    return hb.sort_values('Frac', ascending=False)

def hbonds_lig(avghb, lig, frac=0.4):
    """Return a DataFrame with the ligand HBonds with avg duration higher than a percentage (default 0.4)"""
    df = pd.read_csv(avghb, delim_whitespace=True)
    hb = df[df.Frac >= frac]
    hb = hb.rename(columns = {'#Acceptor':'Acceptor'})
    hb_lig1 = hb[hb['Acceptor'].str.contains(lig)]
    hb_lig2 = hb[hb['Donor'].str.contains(lig)]
    hb_lig = pd.concat([hb_lig1, hb_lig2], ignore_index=True)
    return hb_lig.sort_values('Frac', ascending=False)

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

def contacts(contacts, frac=0.2, avg=4.5):
    """General parser for contacts.dat with variable Frac and Avg"""
    df = pd.read_csv(contacts, delim_whitespace=True, skiprows=2)
    df['Frac.'] = pd.to_numeric(df['Frac.'], errors='coerce')
    df['Nframes'] = pd.to_numeric(df['Nframes'], errors='coerce')
    df['Avg'] = pd.to_numeric(df['Avg'], errors='coerce')
    df['Stdev'] = pd.to_numeric(df['Stdev'], errors='coerce')
    df = df.rename(columns = {'Frac.':'Frac'})
    df = df[df['Frac'].notna()]
    df = df[df.Frac >= frac]
    df = df[df.Avg <= avg]
    df = df.drop(['#'], axis=1)
    ls_cont = []
    contatti = list(df['Contact'])
    for entry in contatti:
        ls_cont.append(entry.split('@')[1].split(':')[1])
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

def salt_bridges(path1, path2):
    first = pd.read_csv(path1,  delim_whitespace=True)
    second = pd.read_csv(path2,  delim_whitespace=True)
    df = pd.concat([first, second], ignore_index=True)
    df = df.sort_values('Frac', ascending=False)
    return df.reset_index(drop=True)


def df_4catplot(df, ligname, interaction, resinfo=None):
    '''Convert normal dataframe(from hbonds, salt_bridge, hydrophobic analyses) to df useful for catplot. Needs type of interction and a df "resinfo" describing resname and numbers (from resinfo command in cpptraj)'''
    
    #Find all the interacting residues
    prot_res= []
#     new_clmns=[]
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
                ligname = row.iloc[0]['Name']
                new = f'{ligname}_{lignum[1:]}@{latm}-{resname}_{resnum[1:]}@{ratm}'
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

