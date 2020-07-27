#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:17:10 2020

@author: gboladekayode
"""
#class vibrational_modes(object):
#
#
#    def __init__(self,outcar):
#        self.outcar=outcar

def read_vibrational_modes(outcar):
    '''Returns DataFrame of Energy Coordinates,vibrational modes,frequencies'''
    import pandas as pd
    from ase.io import read
    atoms=read(outcar)
    elements=atoms.get_chemical_symbols()
    rawData=[]
    f=open(outcar,'r')
    count=0
    lines=f.readlines()
    for line in lines:
        if "THz" in line:
            frequency= line.replace("f/i=","f/i =").split()[3]
            freq_type= line.replace("f/i=","f/i =").split()[1]
            Index=lines.index(line)
            sections=lines[Index+2:Index+len(elements)+2]
            rawData.append([len(elements)])
            rawData.append([count])
            for section in sections:
                section_line=section.split()
                sec_index=sections.index(section)
                rawData.append([elements[sec_index],float(section_line[0])-float(section_line[3]),float(section_line[1])-float(section_line[4]),float(section_line[2])-float(section_line[5]),float(frequency),freq_type])
            count=count+1
            
            
            rawData.append([len(elements)])
            rawData.append([count])
            for section in sections:
                section_line=section.split()
                sec_index=sections.index(section)
                rawData.append([elements[sec_index],float(section_line[0])-float(section_line[3])/2,float(section_line[1])-float(section_line[4])/2,float(section_line[2])-float(section_line[5])/2,float(frequency),freq_type])
            count=count+1
            
            
            rawData.append([len(elements)])
            rawData.append([count])
            for section in sections:
                section_line=section.split()
                sec_index=sections.index(section)
                rawData.append([elements[sec_index],float(section_line[0]),float(section_line[1]),float(section_line[2]),float(frequency),freq_type])
            count=count+1
            
            
            rawData.append([len(elements)])
            rawData.append([count])
            for section in sections:
                section_line=section.split()
                sec_index=sections.index(section)
                rawData.append([elements[sec_index],float(section_line[0])+float(section_line[3])/2,float(section_line[1])+float(section_line[4])/2,float(section_line[2])+float(section_line[5])/2,float(frequency),freq_type])
            count=count+1
            
            
            
            rawData.append([len(elements)])
            rawData.append([count])
            for section in sections:
                section_line=section.split()
                sec_index=sections.index(section)
                rawData.append([elements[sec_index],float(section_line[0])+float(section_line[3]),float(section_line[1])+float(section_line[4]),float(section_line[2])+float(section_line[5]),float(frequency),freq_type])
            count=count+1

    df=pd.DataFrame(rawData,index=None)
    return df


def get_vibrational_frequency(outcar,kind="all"):
    '''Returns Vibrational frequencies in "THz" from an OUTCAR file'''
    import pandas as pd
    import numpy as np
    df=read_vibrational_modes(outcar)
    if kind == "all":
        freq=df[df.columns[4]].unique()
        freq=[x for x in freq if ~np.isnan(x)]
        return freq
    elif kind=="imag":
        freq=df.loc[df[5]=="f/i"]
        freq=freq[freq.columns[4]].unique()
        return freq
    elif kind=="non_imag":
        freq=df.loc[df[5]=="f"]
        freq=freq[freq.columns[4]].unique()
        return freq

def write_xyz(outcar,output_file_name="sample"):
    '''Writes the coordinates and virbrational modes into an xyz file'''
    import pandas as pd
    df=read_vibrational_modes(outcar)
    File=output_file_name+".xyz"
    df.drop([df.columns[4],df.columns[5]], axis=1).to_csv(File, header=None, index=None, sep=' ', mode='a')




def view_vibrational_modes(outcar,frequency_index="all",wrap_center=[0.3,0.3,0.8]):
    '''views the Vibrational modes of an OUTCAR file in the "-", "0" and  "+" states "'''
    from ase.visualize import view
    from ase.io import read
    from ase.geometry import wrap_positions
    from ase import Atoms
    import ase
    import numpy as np
    import os
    if frequency_index == "all":
        output_file_name="sample"
        write_xyz(outcar,output_file_name)
        atoms=read(output_file_name+".xyz@:")
        #print(atoms)
        a=[]
        for atom in atoms:
            atom.wrap(center=wrap_center)
            a.append(atom)
        view(atoms)
        os.remove(output_file_name+".xyz")
    elif frequency_index != "all": # for a list of frequency index
        output_file_name="sample"
        new_df=read_vibrational_modes(outcar)
        Atoms=read(outcar)
        elements=Atoms.get_chemical_symbols()
        for index in frequency_index:
            start=index*3*(len(elements)+2)
            end=3*(len(elements)+2)
            df=new_df[start:start+end]
            File=output_file_name+".xyz"
            df.drop(df.columns[4], axis=1).drop(df.columns[5], axis=1).to_csv(File, header=None, index=None, sep=' ', mode='a')
        atoms=read(output_file_name+".xyz@:")
        a=[]
        for atom in atoms:
            atom.wrap(center=wrap_center)
            a.append(atom)
        view(a)
        os.remove(output_file_name+".xyz")

