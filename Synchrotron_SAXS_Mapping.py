# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:39:14 2023

@author: amart
"""



import numpy as np

import pyFAI

import pandas as pd

import seaborn as sns

import matplotlib.pyplot as plt

import os 

import fabio

from scipy.integrate import simpson
from numpy import trapz
from sklearn.metrics import auc

from scipy.optimize import curve_fit

sns.set_theme()
sns.set_style("ticks")
sns.set_context("paper",font_scale=1.5)
sns.set_palette("tab10")

os.chdir(r"XXXXXXXX") #Input working directory

def recta(x,m,a):
    
    y = 10**a*x**m
    
    return y


z = {}

x = {}

i = 0
j = 0

for file in os.listdir("raw"):
    if file == "mesh_snapdmesh_4029":
        
        for data_file in os.listdir(f"./raw/{file}/data_00/"):
                
            data_dir = f"./raw/{file}/data_00/{data_file}"
            
            
            img = fabio.open(data_dir)
            
            header = img.getheader()

            z_pos = round(float(header["stz"]),2)
            x_pos = round(float(header["sx"]),2)
            
            if str(z_pos) not in z.keys():
                
                z[str(z_pos)] = i
                
                i+=1
            
            if  str(x_pos) not in x.keys():
                
                x[str(x_pos)] = j
                j+=1

    


for mesh in os.listdir("./Processed/"):
    map_data = np.zeros((len(z),len(x)))
    absorption_map = np.zeros((len(z),len(x)))
    if "mesh" in mesh:

        for file in os.listdir(f"./Processed/{mesh}/data/LaVue1D/"):
        
            if "pilatus" in file:
                name = file[:-4]     
                #print(name)
                
                processed_file = f"./Processed/{mesh}/data/LaVue1D/{name}.dat"
                raw_file = f"./raw/{mesh}/data_00/{name}.edf"
                
                raw = fabio.open(raw_file)
                
                header = raw.getheader()
        
                z_pos = round(float(header["stz"]),2)
                x_pos = round(float(header["sx"]),2)
                
                row = z[str(z_pos)]
                col = x[str(x_pos)]
                
                
                Photo = float(header["Photo"])/float(header["Monitor"])
        
                processed = pd.read_table(processed_file,skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
                
                
                
                ROI = processed.loc[(processed.q > 0.43) & (processed.q<0.7),:].reset_index()
                ROI.I = ROI.I/ROI.I.iloc[0]
                
                y1 = np.log10(ROI.iloc[0,2])
                y2 = np.log10(ROI.iloc[-1,2])
                
                x1 = np.log10(ROI.iloc[0,1])
                x2 = np.log10(ROI.iloc[-1,1])
                
                m = (y2-y1)/(x2-x1)
                
                a = y2-m*x2
                
                plt.plot(ROI.q,ROI.I)
                plt.plot(ROI.q,recta(ROI.q,m,a))
                plt.xscale("log")
                plt.yscale("log")
                #plt.show()
                plt.close()
                
                area = auc(ROI.q,ROI.I-recta(ROI.q,m,a))
                map_data[row,col] = area
                absorption_map[row,col] = Photo
                
                
                
        fig = plt.figure(figsize=(6,6))     
        map_data = map_data/np.max(map_data)
        plt.imshow(map_data,aspect=1,cmap = "gnuplot2",vmin=0)  
        plt.colorbar()  
        #plt.axis("Off")
        plt.xlabel("Relative x position to (0,0) ($\mu$m)")
        plt.ylabel("Relative z position to (0,0) ($\mu$m)")
        plt.xticks(np.arange(0,10,1),[int(float(i)*1000+104) for i in x.keys()])
        plt.yticks(np.arange(0,10,1),[int(np.abs(float(i))*1000-27463) for i in z.keys()])
        plt.title("Normalized Concentration Map")
        plt.savefig("Concentration map_log.tif",dpi=300,bbox_inches="tight")
        plt.show()
        plt.close()

        
        fig = plt.figure(figsize=(6,6))     
        absorption_map = absorption_map/np.max(absorption_map)
        plt.imshow(absorption_map,aspect=1,vmin=0.9)  
        plt.colorbar()  
        #plt.axis("Off")
        plt.xlabel("Relative x position to (0,0) ($\mu$m)")
        plt.ylabel("Relative z position to (0,0) ($\mu$m)")
        plt.xticks(np.arange(0,10,1),[int(float(i)*1000+104) for i in x.keys()])
        plt.yticks(np.arange(0,10,1),[int(np.abs(float(i))*1000-27463) for i in z.keys()])
        plt.title("Normalized Transmittance Map")
        plt.savefig("Concentration map_log_absorption.tif",dpi=300,bbox_inches="tight")
        break
    #plt.savefig()

map_data = pd.DataFrame(map_data)


map_data.to_csv("Concentration_map_data_6.0_mg_ml.txt",sep = "\t",index=None)