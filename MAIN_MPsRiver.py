# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 17:08:03 2020

Last modified on Mon Jan 25 2021

@author: PradoDomercq
"""

#Timer
from datetime import datetime
startTime = datetime.now()

"""IMPORT MODULES AND LIBRARIES NEEDED"""
import os
import pandas as pd
import itertools
import numpy as np
from objectGenerationRiver_func import*
from DynamicModelCore_River_MPs import*
from GlobalConstants import * 
from readImputParam import readProcessparam, microplasticData,readCompartmentData
#from createClistWithBoxes import createClistBox
#import modlue of process rate constant calculations
import RC_GeneratorRiver
from dilutionVol_calculator_func import*
from Multiplots_function import*


#####Set current working directory
cwd = os.getcwd()

"""IMPORT INPUT FILES"""
os.chdir(cwd+"/Inputs")
#Import process parameters
process_df= readProcessparam ("process_paramRiver.txt")
#Import MP parameters (radius, volume, etc) Same way
MP_prop = microplasticData("microplasticsSizeClass.txt")
#Import compartment info
compartments_prop = readCompartmentData("compartmentsGenericRiverSec_prop.txt")
#Add river section depth field
RSdepth = []
for row in range(len(compartments_prop)):
        RSdepth.append(round(sum(compartments_prop.depth_m[0:4]),2))
compartments_prop["depthRS_m"]=RSdepth
os.chdir(cwd)
"""SET UP MODEL"""
#RIVER SECTIONS
numberRS=len (compartments_prop)/4
listRS = [*range(0,int(numberRS),1)]
riverSect = [str(item) for item in listRS]
riverLengths = [str(it) for it in compartments_prop["length_m"]]
riverSectLength= riverLengths[0::4]
RS_cumLength_m =[]
for d in range(len(riverSectLength)):
    if d==0:
        RS_cumLength_m.append(float(riverSectLength[d]))
    else:
        RS_cumLength_m.append(float(riverSectLength[d])+float(RS_cumLength_m[d-1]))


##############
    
#RIVER COMPARTMENTS
riverComp = ["1", "2", "3", "4"]
compartments = ["Surface Water", "Flowing Water", "Stagnant Water", "Sediment"]

#MICROPLASTICS FORMS 
MPforms = ["A", "B", "C", "D"]
MPformslabels = ["Free", "Heteroaggregated", "Biofiolm-covered", "Biofilm-heteroaggregated"]

#SIZE BINS
sizeBin =["a", "b", "c", "d", "e"]
sizeBinLabel = ["Cotton", "Cotton","Polyester", "Polyester", "Rayon"]# Detection limit for MPs via Fourier Transform Infrared Spectroscopy is 20um

#MPS RIVER PROCESSES (FATE AND TRANSPORT) LIST
processList = ["degradation", "fragmentation", "heteroagg", "breakup", "settling","rising", "advection", "mixing", "biofilm", "resusp", "burial","sedTransport"]


"""DEFINE RUN PARAMETERS"""

#define solver --Y steady state versus time resolved
SOLVER = "Dynamic" #options: SteadyState OR Dynamic (SteadyState under developement)
mode = "Standard" #options: Standard OR Monthly (Monthly under developement)
mode2 = "Timelimit" #options: Timelimit OR raw 
#(Timelimit mode sets up a time limit of 30 min on the processes to be accounted for while
# raw leaves the RC as calcualted) The raw version can straing the solver due to time.
record = "True" # if "True" the results and RC and intercations dataframes will be recorded in the Results folder
#define MP composition to study (from microplastics imput file, composition: XX, PP, PS, mixed)
composition = "XX"
#define imput location and MP type: RiverSection, river compartment, MP aggregation state and size bin. 
#Ex: 02Ae (RS=0:1, comp= 2:flowing water, MPtype:A:FreeMP, sizeBin:e:1000um)
imputMP= "02Ae"
imputFlow=100 # (in number of particles per minute)
imputPulse=0

#set starting time of simulation (typically 0)

#Minutes (t_span) or days (t_days) 
t0 = 0 #set starting time of simulation (typically 0)
daysSimulation = 720
tmax = 24*60*daysSimulation #set final simulation time (in min)
min_day = 24*60
t_days = np.linspace(t0, tmax, daysSimulation + 1) #define timespan & time_step, e.g. 1min (min, max, number of elements)
stepSize= 60
timesteps = int(min_day*daysSimulation/stepSize) #number of simulation steps (currently corresponding to a stepsize of 1day = 60*24 min)
t_span = np.linspace(0, tmax, int(timesteps)+1, dtype=int)


#############
from datetime import datetime, timedelta
date_time_str = '2020-01-01 00:00'
DayStart = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M')
LastDay = DayStart + timedelta(minutes=tmax)
date = DayStart

daterun = date.today()
daterun_label = daterun.strftime("%d_%m_%Y")
##############

####Define name for the model run to be stored under
runtitle= "GenericRiver_" + composition +"_runtime"+str(daysSimulation)+"_ts"+str(stepSize)+"min_Imput_"+ imputMP + "_abcde_"+mode2+"_sedTransp"


"""COMBINATIONS LIST"""

#Generate COMBINATIONS
combinations = list(itertools.product(riverSect,riverComp,MPforms,sizeBin))
#Generate raw list of combinations and lists of concentrations (C) and inflows (I)
CombList = []
Ilist = []
Clist =[]
def convertTuple(tup): 
    str =  ''.join(tup) 
    return str
for e in combinations:
    Clist.append("C_" + convertTuple(e))
    Ilist.append("I_" + convertTuple(e))
    CombList.append(convertTuple(e))
    
"""Estimate RATE CONSTANTS"""# of the list of processes of the "combinations" (SPECIES) 
##Estimate rate constants for all "species" = combination of RS*riverComp*MPforms*sizeBins (len (Clist))

#Dataframe to store RC por each SPECIES
RC_df = pd.DataFrame(index=processList , columns= CombList)
#Add a row to store values of volume of the particles and its density
extraRow = pd.DataFrame(index= ["volume_m3","density_kg_m3"] , columns= CombList)
RC_df=pd.concat([RC_df,extraRow])


if numberRS <= 10:
    print("Less than 10 RS")
#if len(Clist) <= 10*len(riverComp)*len(MPforms)*len(sizeBin):
    for spc in Clist:
       
        riverSection = spc[2]
        compartment= spc[3]
        aggState = spc[4]
        sizeBinIdx = spc[5]
    
        #Generate particle OBJECTS MPs and SPM for the corresponding size bin
        
        #generate location index from imput tables
        MP_index = np.where((MP_prop['composition']== composition) & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
        SPM_index = np.where((MP_prop['name']== "SPM1") & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
        comp_index = np.where((compartments_prop['compartment']== int(compartment)) & (compartments_prop['riverSection']== int(riverSection)))[0][0]
        
        
        [MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM]=preProcessElements(MP_prop, MP_index, SPM_index, compartments_prop, comp_index)#add compartments_prop, comp_index to include corresponding SPM concentration
    
        #Create a dictionary with the objects 
        particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
          
        #load compartment objects corresponding to the River Section
        
        [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, comp_index)
        #Create a dictionary with the environmental compartments objects 
        comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
        
        #import data needed for each spcn(t_half, radius...)
        #for processes
        idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
        
    
        RC_df.loc["degradation", spc[2:6]]= RC_GeneratorRiver.degradation(process_df.t_half_d.loc[idx]) 
        RC_df.loc["fragmentation", spc[2:6]]= RC_GeneratorRiver.fragmentation(process_df.k_frag_d.loc[idx], particles_dict[aggState].radius_m,  particles_dict[aggState].volume_m3,  particles_dict[aggState].diameter_um, sizeBinIdx)
        RC_df.loc["heteroagg", spc[2:6]]= RC_GeneratorRiver.heteroagg(process_df.alpha.loc[idx], MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, particles_dict["A"].CSF)
        RC_df.loc["breakup", spc[2:6]]= RC_GeneratorRiver.breakup(process_df, idx, MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
        RC_df.loc["settling", spc[2:6]]= RC_GeneratorRiver.settling(particles_dict["A"].density_kg_m3, particles_dict[aggState].radius_m, surface.depth_m, "CalcFiber", compartment, particles_dict["A"].CSF)
        RC_df.loc["rising", spc[2:6]]= RC_GeneratorRiver.rising(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, flowingWater.depth_m, "CalcFiber", compartment, particles_dict["A"].CSF)
        RC_df.loc["advection", spc[2:6]]= RC_GeneratorRiver.advection(comp_dict, compartment)
        RC_df.loc["mixing", spc[2:6]]= RC_GeneratorRiver.mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
        RC_df.loc["biofilm", spc[2:6]]= RC_GeneratorRiver.biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
        RC_df.loc["resusp", spc[2:6]]= RC_GeneratorRiver.resusp(compartment, comp_dict)
        RC_df.loc["burial", spc[2:6]]= RC_GeneratorRiver.burial(compartment, comp_dict)
        RC_df.loc["sedTransport", spc[2:6]]= RC_GeneratorRiver.sedTransport(compartment,comp_dict)
        RC_df.loc["volume_m3", spc[2:6]]= particles_dict[aggState].volume_m3
        RC_df.loc["density_kg_m3", spc[2:6]]= particles_dict[aggState].density_kg_m3

elif numberRS <= 100:
    print("10 or more RS but less than 100")
    for spc in Clist:
        if spc in Clist[0:10*len(riverComp)*len(MPforms)*len(sizeBin)]:
            riverSection = spc[2]
            compartment= spc[3]
            aggState = spc[4]
            sizeBinIdx = spc[5]
            #spcAll=spc[2:6]
        elif spc in Clist[10*len(riverComp)*len(MPforms)*len(sizeBin):100*len(riverComp)*len(MPforms)*len(sizeBin)]:
        #else:
            riverSection = spc[2:4]
            compartment= spc[4]
            aggState = spc[5]
            sizeBinIdx = spc[6]
            #spcAll=spc[2:7]
        #Generate particle OBJECTS MPs and SPM for the corresponding size bin
        
        #generate location index from imput tables
        MP_index = np.where((MP_prop['composition']== composition) & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
        SPM_index = np.where((MP_prop['name']== "SPM1") & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
        comp_index = np.where((compartments_prop['compartment']== int(compartment)) & (compartments_prop['riverSection']== int(riverSection)))[0][0]
        
        [MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM]=preProcessElements(MP_prop, MP_index, SPM_index, compartments_prop, comp_index)#add compartments_prop, comp_index to include corresponding SPM concentration
    
        #Create a dictionary with the objects 
        particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
          
        #load compartment objects corresponding to the River Section
        
        [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, comp_index)
        #Create a dictionary with the environmental compartments objects 
        comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
        
        #import data needed for each spcn(t_half, radius...)
        #for processes
        idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
        
    
        RC_df.loc["degradation", spc[2:]]= RC_GeneratorRiver.degradation(process_df.t_half_d.loc[idx]) 
        RC_df.loc["fragmentation", spc[2:]]= RC_GeneratorRiver.fragmentation(process_df.k_frag_d.loc[idx], particles_dict[aggState].radius_m,  particles_dict[aggState].volume_m3,  particles_dict[aggState].diameter_um, sizeBinIdx)
        RC_df.loc["heteroagg", spc[2:]]= RC_GeneratorRiver.heteroagg(process_df.alpha.loc[idx], MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment,particles_dict["A"].CSF)
        RC_df.loc["breakup", spc[2:]]= RC_GeneratorRiver.breakup(process_df, idx, MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
        RC_df.loc["settling", spc[2:]]= RC_GeneratorRiver.settling(particles_dict["A"].density_kg_m3, particles_dict[aggState].radius_m, surface.depth_m, "CalcFiber", compartment,particles_dict["A"].CSF)
        RC_df.loc["rising", spc[2:]]= RC_GeneratorRiver.rising(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, flowingWater.depth_m, "CalcFiber", compartment,particles_dict["A"].CSF)
        RC_df.loc["advection", spc[2:]]= RC_GeneratorRiver.advection(comp_dict, compartment)
        RC_df.loc["mixing", spc[2:]]= RC_GeneratorRiver.mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
        RC_df.loc["biofilm", spc[2:]]= RC_GeneratorRiver.biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
        RC_df.loc["resusp", spc[2:]]= RC_GeneratorRiver.resusp(compartment, comp_dict)
        RC_df.loc["burial", spc[2:]]= RC_GeneratorRiver.burial(compartment, comp_dict)
        RC_df.loc["sedTransport", spc[2:]]= RC_GeneratorRiver.sedTransport(compartment,comp_dict)
        RC_df.loc["volume_m3", spc[2:]]= particles_dict[aggState].volume_m3
        RC_df.loc["density_kg_m3", spc[2:]]= particles_dict[aggState].density_kg_m3
        
elif numberRS <= 1000:
    print("100 or more RS but less than 1000")
    for spc in Clist:
        if spc in Clist[0:10*len(riverComp)*len(MPforms)*len(sizeBin)]:
            riverSection = spc[2]
            compartment= spc[3]
            aggState = spc[4]
            sizeBinIdx = spc[5]
        elif spc in Clist[10*len(riverComp)*len(MPforms)*len(sizeBin):100*len(riverComp)*len(MPforms)*len(sizeBin)]:
            riverSection = spc[2:4]
            compartment= spc[4]
            aggState = spc[5]
            sizeBinIdx = spc[6]
        elif spc in Clist[100*len(riverComp)*len(MPforms)*len(sizeBin):1000*len(riverComp)*len(MPforms)*len(sizeBin)]:
            riverSection = spc[2:5]
            compartment= spc[5]
            aggState = spc[6]
            sizeBinIdx = spc[7]
        #Generate particle OBJECTS MPs and SPM for the corresponding size bin
        
        #generate location index from imput tables
        MP_index = np.where((MP_prop['composition']== composition) & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
        SPM_index = np.where((MP_prop['name']== "SPM1") & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
        comp_index = np.where((compartments_prop['compartment']== int(compartment)) & (compartments_prop['riverSection']== int(riverSection)))[0][0]
        
        
        [MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM]=preProcessElements(MP_prop, MP_index, SPM_index, compartments_prop, comp_index)#add compartments_prop, comp_index to include corresponding SPM concentration
    
        #Create a dictionary with the objects 
        particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
          
        #load compartment objects corresponding to the River Section
        
        [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, comp_index)
        #Create a dictionary with the environmental compartments objects 
        comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
        
        #import data needed for each spcn(t_half, radius...)
        #for processes
        idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
        
    
        RC_df.loc["degradation", spc[2:]]= RC_GeneratorRiver.degradation(process_df.t_half_d.loc[idx]) 
        RC_df.loc["fragmentation", spc[2:]]= RC_GeneratorRiver.fragmentation(process_df.k_frag_d.loc[idx], particles_dict[aggState].radius_m,  particles_dict[aggState].volume_m3,  particles_dict[aggState].diameter_um, sizeBinIdx)
        RC_df.loc["heteroagg", spc[2:]]= RC_GeneratorRiver.heteroagg(process_df.alpha.loc[idx], MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment)
        RC_df.loc["breakup", spc[2:]]= RC_GeneratorRiver.breakup(process_df, idx, MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
        RC_df.loc["settling", spc[2:]]= RC_GeneratorRiver.settling(self,particles_dict["A"].density_kg_m3, particles_dict[aggState].radius_m, surface.depth_m, "CalcFiber", compartment)
        RC_df.loc["rising", spc[2:]]= RC_GeneratorRiver.rising(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, flowingWater.depth_m, "Stokes", compartment)
        RC_df.loc["advection", spc[2:]]= RC_GeneratorRiver.advection(comp_dict, compartment,compartments_prop.discharge_m3s[comp_index])
        RC_df.loc["mixing", spc[2:]]= RC_GeneratorRiver.mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
        RC_df.loc["biofilm", spc[2:]]= RC_GeneratorRiver.biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
        RC_df.loc["resusp", spc[2:]]= RC_GeneratorRiver.resusp(compartment, comp_dict)
        RC_df.loc["burial", spc[2:]]= RC_GeneratorRiver.burial(compartment, comp_dict)
        RC_df.loc["sedTransport", spc[2:]]= RC_GeneratorRiver.sedTransport(compartment,comp_dict)
        RC_df.loc["volume_m3", spc[2:]]= particles_dict[aggState].volume_m3
        RC_df.loc["density_kg_m3", spc[2:]]= particles_dict[aggState].density_kg_m3
else:
    print("Number of river sections excedes 1000. Implement new code for RCs")    
#substitute NAN values by 0 in Rate cosntants dataframe 
for x in RC_df:
    RC_df[x] = RC_df[x].fillna(0)
    

    
#OPTIONAL SPEED PROCESS CONTROL
#To speed up the solver eliminate those RC for wich proccesses are faster than 30min (1800sec)--> RC<5.56x10-4
if mode2 == "Timelimit":
    RC_df_templim=pd.DataFrame(index=processList , columns= CombList)
    #Add a row to store values of volume of the particles and its density
    extraRow = pd.DataFrame(index= ["volume_m3","density_kg_m3"] , columns= CombList)
    RC_df_templim=pd.concat([RC_df_templim,extraRow])
    for x in RC_df:
        for y in range(len(RC_df)-2):   
            if type(RC_df[x][y]) is tuple:
                # when there are fragements formed extract the fragmentation array to multiply kfrag by number of fragments formed
                frag = RC_df[x][y]
                fragRC = frag[0]
                if (fragRC !=0) and (fragRC > 0.000556):
                    #print(RC_df.index[y]+" RC for " +x+ " is fastest than 30min")
                    fragRC=0.000556
                    RC_df_templim[x][y]=(fragRC,frag[1])
                else:
                    RC_df_templim[x][y]=RC_df[x][y]
            else:   
                if (RC_df[x][y]!=0) and (RC_df[x][y] > 0.000556):
                    #print(RC_df.index[y]+" RC for " +x+ " is fastest than 30min")
                    RC_df_templim[x][y]=0.000556
                else:
                    #print(RC_df[x][y])
                    RC_df_templim[x][y]=RC_df[x][y]
        for y in range(len(RC_df)-2,len(RC_df)):
            RC_df_templim[x][y]= RC_df[x][y]
    for x in RC_df_templim:
        if (RC_df_templim.loc["breakup",[x]][0] == 0.000556):
            RC_df_templim.loc["breakup",[x]] = 0.000556/10



"""FILL MATRIX OF INTERACTIONS""" 
#Fill in interactions matrix with process rate constants
from fillRCmatrixInteractionsTransport_func import*
from fillRCinteractionMatrices_func import*


if mode2 == "Timelimit":
    interactions_df = fillMatrixInteractionsALL(RC_df_templim, Clist,compartments_prop)
    Interaction_dataframes= fillMatricesOfInteractions(RC_df_templim, Clist,compartments_prop)
    klosses_df= fillMatricesOfIrrevLoss(RC_df_templim, Clist,compartments_prop) 
    
else:
    interactions_df = fillMatrixInteractionsALL(RC_df, Clist,compartments_prop)
    Interaction_dataframes= fillMatricesOfInteractions(RC_df, Clist,compartments_prop)
    klosses_df= fillMatricesOfIrrevLoss(RC_df, Clist,compartments_prop)



"""SOLVE THE MODEL"""      

#Vector of volumes corresponding to the compartments of the river
dilution_vol_m3= volumesVector(Clist,compartments_prop)

#Dataframes for storing starting and final concentrations for all possible species of each size class
  
Concentrations_t0 = pd.DataFrame(index=Clist, columns=['number (#/m3)'])
for p in range(len(Concentrations_t0)):
        Concentrations_t0.iloc[p][0]= 0
Concentrations_t0.loc["C_"+imputMP]=imputPulse
for p in range(len(Concentrations_t0)):
        Concentrations_t0.iloc[p][0]= Concentrations_t0.iloc[p][0]/dilution_vol_m3[p]

#set inflow 
Ilist = []
for C in Clist:
    Ilist.append("I"+ C[1:])
inflow_vector = pd.DataFrame(index=Ilist, columns=["number (#/m3)"])
inflow_vector.loc[:,:] = 0
inflow_vector.loc["I_"+imputMP] = imputFlow# to be read in from emission files
for q in range(len(inflow_vector)):
        inflow_vector.iloc[q][0]= inflow_vector.iloc[q][0]/dilution_vol_m3[q]
 

#########DYNAMIC MODEL CORE#############

#solve model in steady state mode
if SOLVER == "SteadyState":
    q=50 #MPs released to MainWater per second
    SteadyStateResults = SteadyState(MainWater, SurfWater, Sediment, q)#to be implemented (OLD VERSION)
    
    #solve model dynamically
elif SOLVER == "Dynamic":
    
   DynamicResults = solverPrado(Concentrations_t0, interactions_df,inflow_vector, tmax, timesteps, t_span)
   
   ConcFinal_num_m3=DynamicResults[0]
   
   #Substitute values smaller than 10-5 to 0
   ConcFinal_num_m3 = ConcFinal_num_m3.apply(lambda x: [y if y >= 1e-15 else 0 for y in x])
  
   ###Convert particle number to mass concentration
   ConcFinal_kg_m3= pd.DataFrame(index=t_span , columns= Concentrations_t0.index)  
   volume= RC_df.loc["volume_m3"].to_numpy()
   density= RC_df.loc["density_kg_m3"].to_numpy()
   for row in ConcFinal_num_m3.index:
        ConcFinal_kg_m3.loc[row, :]= ConcFinal_num_m3.loc[row, :]*volume*density
     
else:
    print("Error: incorrect solver selected")
    
    
##Estimate relative abundances over the whole river:
    
from relativeAbundance_calculator_func import*

RelativeAbundances=Relative_Abundance_calc (ConcFinal_num_m3, Clist,MPforms,riverComp,sizeBin, MPformslabels,compartments,stepSize,dilution_vol_m3,sizeBinLabel,t_span)
RelativeAbun_MPtype_t=RelativeAbundances[0]
RelativeAbun_Comp=RelativeAbundances[1]
RelativeAbun_Size=RelativeAbundances[2]



#SAVE model RESULTS to the results folder
#change folder
if record == "True":   
    os.chdir(cwd+"/Results")
    
    #create folder for the day if it doesnt already exists
    
    path = cwd+"/Results"
    os.path.isdir(path)
    old_path = (daterun_label)
    new_path = os.path.isdir(old_path)
    if not new_path:
        os.makedirs(old_path)
        print("Created Folder : ", old_path)
    else:
        print(old_path, "folder already exists.")
    
    results_path= cwd+"/Results/"+old_path
    
    os.chdir(results_path)
    #save
    #rate constants 
    dfRC_filename = "RC_df" + runtitle +"_"+ daterun_label + ".csv"
    RC_df_templim.to_csv(dfRC_filename)
    
    #interactions dataframe to results folder
    interactions_filename = "interactionsdf_"+ runtitle +"_"+  daterun_label + ".csv"
    interactions_df.to_csv(interactions_filename)
    
    #Results
    filename= "ConcVStime_particleNum_m3"+ runtitle +"_"+ daterun_label+".csv"
    filename1= "MassVStime_kg"+ runtitle +"_"+ daterun_label+".csv"
    ConcFinal_num_m3.to_csv(filename)
    ConcFinal_kg_m3.to_csv(filename1)  
    RelativeAbun_MPtype_t.to_csv("Relative_abundance_MPtype1_"+composition+daterun_label+".csv")  
    RelativeAbun_Comp.to_csv("Relative_abundance_Compartment1_"+composition+daterun_label+".csv")  
    #Generate and save figures
    Multiplots (ConcFinal_num_m3,t0, tmax, timesteps, results_path, composition, compartments, MPforms, MPformslabels, sizeBin, sizeBinLabel, RS_cumLength_m)   
    #back to working directory
    os.chdir(cwd)
    

print(datetime.now() - startTime)