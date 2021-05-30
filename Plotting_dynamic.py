# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 14:47:10 2021

@author: PradoDomercq
"""
# Import packages
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from celluloid import Camera
import numpy as np
import os


def Multiplots_dynamic(ConcFinal_num_m3, Clist, listRS, stepSize, t_span, compartments, MPformslabels, MPforms, imputPulse, imputFlow,composition):
    
    #define time resolution for the results
    numTstep_hour=(60/stepSize)
    Time_months=t_span[::(int(numTstep_hour*24*30))]
    Time_days=t_span[::(int(numTstep_hour*24))]
    Time_halfMonth=t_span[::(int(numTstep_hour*24*15))]
    
    #Select time resolution for extracting results
    time_extract=Time_months
    
    def extract_SizeBins (t, comp, MPform,): #concentration values for each size bin of the comp-Mpform selected aong the river (RSs)
        Aa=[]
        Ab=[]
        Ac=[]
        Ad=[]
        Ae=[]
        for i in range(len(listRS)):
            Aa.append(ConcFinal_num_m3.values[t, Clist.index("C_"+str(listRS[i])+comp+MPform+"a")])
            Ab.append(ConcFinal_num_m3.values[t, Clist.index("C_"+str(listRS[i])+comp+MPform+"b")])
            Ac.append(ConcFinal_num_m3.values[t, Clist.index("C_"+str(listRS[i])+comp+MPform+"c")])
            Ad.append(ConcFinal_num_m3.values[t, Clist.index("C_"+str(listRS[i])+comp+MPform+"d")]) 
            Ae.append(ConcFinal_num_m3.values[t, Clist.index("C_"+str(listRS[i])+comp+MPform+"e")]) 
        return [Aa, Ab, Ac, Ad, Ae]    


    
    #List of dataframes (len = timepoints). 
    #Each dataframe corresponds to a time point and and organises concentration values of MPs per MP aggregation state and compartment
    #The cell corresponding to a specific agg state and compartment contains a list of concentrations per river section for each size bin (list (len RS) of lists(number of lists = len size bins)) 
    y_multiplot_time= [[] for t in range(len(time_extract))]
    for t in range(len(time_extract)):
        y_multiplot =pd.DataFrame(index=compartments, columns=MPformslabels)
        #Estimate relative abundance of MP types in the whole river and of particles per compartment
        for j in range(len(compartments)):
                for k in range(len(MPforms)):
                    y_multiplot.iloc[j,k] = extract_SizeBins(int(time_extract[t]/stepSize), str(j+1),MPforms[k])  
        y_multiplot_time[t]=y_multiplot               
      
    #####MULTI PLOTS
        
    #create folder to store the figures for the day if it doesnt already exists
    
    os.chdir(results_path)
    Fig_folder= "/Figures"
    os.path.isdir(results_path+Fig_folder)
    
    new_path = os.path.isdir(results_path+Fig_folder)
    if not new_path:
        os.mkdir("Figures")
        print("Created Folder : ", Fig_folder)
    else:
        print(Fig_folder, "folder already exists.")
    
    results_figures_path= results_path+Fig_folder
    
    os.chdir(results_figures_path)

    #### Choose style and colour palette
    palette = plt.get_cmap('Set2')
    plt.style.use('seaborn-white')
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='ivory', alpha=0.5)
    props2 = dict(boxstyle='round', facecolor='white', alpha=0.5)

    #### 
    fig, axs = plt.subplots(len(compartments),len(MPforms), figsize=(30, 15),sharex='col')
    labels = ["Cotton", "Cotton","Polyester", "Polyester", "Rayon"]
    
    #Distance values
    x =[d/1000 for d in RS_cumLength_m]
    compartmentsLabel=["Surface\n Water", "Flowing\n Water", "Stagnant\n Water", "Sediment"]
    
    if imputFlow == 0:
        fig.suptitle( composition +" microplastics fate in generic river over time (pulse= "+ str(imputPulse)+" particles of 1mm)", fontsize=30)
    else:
        fig.suptitle( composition +" microplastics fate in generic river over time (inflow= "+str(imputFlow)+" particles/min of 1mm)", fontsize=30)
    camera = Camera(fig)
    for t in range(len(time_extract)): 
        for j in range(len(compartments)):
            for k in range(len(MPforms)):
                #Plot
                y = y_multiplot_time[t].iloc[j,k]
                axs[j,k].plot(x, y[0], linewidth=2.5,color=palette(0), label='Cotton')
                axs[j,k].plot(x, y[1], linewidth=2.5,color=palette(1), label='Cotton')
                axs[j,k].plot(x, y[2], linewidth=2.5,color=palette(2), label='Polyester')
                axs[j,k].plot(x, y[3], linewidth=2.5,color=palette(3), label='Polyester')
                axs[j,k].plot(x, y[4], linewidth=2.5,color=palette(4), label='Rayon')
                axs[j,k].text(0.90,0.90,str(int(time_extract[t]/60/24)) + ' days', bbox=props, fontsize=12, transform=axs[j,k].transAxes, ha='center', va='center')
                #Set axes labels and title
                if j== 0:
                    axs[j,k].text(0.5,1.1, MPformslabels[k]+ "\n " +str(RelativeAbun_MPtype_t.iloc[t,k+1])+ " %",fontsize=20, bbox=props2, transform= axs[j,k].transAxes, ha='center')        
                else:
                    axs[j,k].set_title(MPformslabels[k]+ " " + composition+ " MPs in "+ compartments[j],fontsize=15)
          
                if k==0:
                    axs[j,k].text(-0.25, 0.5, compartmentsLabel[j]+"\n"+str(RelativeAbun_Comp.iloc[t,j+1])+" %",fontsize=20,rotation=0, va='center',ha='center',bbox=props2, transform=axs[j,k].transAxes)
                    
            
                axs[j,k].set_yscale('log')
                #if j==0:
                    #axs[j,k].set_ylim(0,1.5)
                #elif j==1:
                    #axs[j,k].set_ylim(0,0.004)
                #elif j== 2:
                   # axs[j,k].set_ylim(-0.04,0.04)
                #else:
                    #axs[j,k].set_ylim(-0.04,0.04)
                #axs[j,k].set_xlim(x[0],x[-1])
                axs[j,k].tick_params(axis='x', labelsize=12)
                axs[j,k].tick_params(axis='y', labelsize=12)
                   
        # Axis titles
        plt.text(0.05, 0.5, "Concentration of particles (Num/$m^3$)", fontsize=20, transform=plt.gcf().transFigure, rotation='vertical',ha='center', va='center')
        plt.text(0.5, 0.08, "Distance (km)", fontsize=20, transform=plt.gcf().transFigure, ha='center', va='center')
        #plt.legend(labels,bbox_to_anchor=(0.5, -0.18), loc='center',ncol=5, fontsize=15 )
        camera.snap()
               
           
    animationTitle= composition+ "_MPs_fate_dynamicPlot.gif"
    animation = camera.animate()
    animation.save(animationTitle, writer='matplotlib.animation.PillowWriter', fps=2) 
