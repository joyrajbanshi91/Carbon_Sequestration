# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 12:07:54 2022

@author: Joy
"""
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.pylab import rcParams
rcParams['figure.figsize'] = 10,5

df = pd.read_csv(r"D:\GCAM_Output_Analysis\GCAM_Crop_CarbonSeq\Cropland_Soil_Seq\Data\Cropland_Data_NZ_FP_2070.csv")
Soil_C =  df['Soil_C'].iloc[0]
Veg_C = df['Veg_C'].iloc[0]


# Assigning incremental years taking from 1st and last year input values
Year = np.arange(df["Year"].iloc[0],df["Year"].iloc[-1]+1,1)
df2 = pd.DataFrame(Year, columns=['Year'])
# joing the year with the input dataframe and linear interpolating areas based on years
merged_data= df2.merge(df, how='outer',on=["Year"])
merged_data['LU_Area_ha']=merged_data['LU_Area_ha'].interpolate(method='linear')

# Calculating land area difference
merged_data['land_diff'] = merged_data['LU_Area_ha'].diff(periods=-1)
merged_data['land_diff']= merged_data['land_diff'].shift(periods=1, fill_value=0)

# Calculating Soil stock
merged_data['Soil_Carbon_Stock']=merged_data['land_diff']*Soil_C


#####################################  soil Emission #########################

# The first function is for Existing Area Soil Carbon Seq 

def Soil_exponential_Init(start_year,end_year):
    start_year=start_year 
    end_year=end_year
    SoilTimeScale=25  # 25 years for soil carbon recovery following LUC in Indian Region (Source: GCAM)
    halfLife=SoilTimeScale/10
    log2 = math.log( 2.0)
    Constant_K = log2 / halfLife;      
    cumStockDiff_t1 = 0
    yearCounter = 0
    exponential_diff=[]
    for i in range(start_year,end_year+1,1):
        yearCounter+=1
        cumStockDiff_t2 = ( 1.0 - math.exp( -1.0 * Constant_K * yearCounter ))
        exponential_diff.append(cumStockDiff_t2-cumStockDiff_t1)
        cumStockDiff_t1 = cumStockDiff_t2
    exponential_diff=np.array(exponential_diff)
    return exponential_diff

# The second function is for additional area soil carbon seq  

def Soil_exponential(start_year,end_year,carbon_diff):
    start_year=start_year 
    end_year=end_year
    carbon_diff=carbon_diff
    SoilTimeScale=25  # 25 years for soil carbon recovery following LUC in Indian Region (Source: GCAM)
    halfLife=SoilTimeScale/10
    log2 = math.log( 2.0)
    Constant_K = log2 / halfLife;      
    cumStockDiff_t1 = 0
    yearCounter = 0
    exponential_diff=[]
    for i in range(start_year,end_year+1,1):
        yearCounter+=1
        cumStockDiff_t2 = ( 1.0 - math.exp( -1.0 * Constant_K * yearCounter ))
        exponential_diff.append(cumStockDiff_t2-cumStockDiff_t1)
        cumStockDiff_t1 = cumStockDiff_t2
    exponential_diff=np.array(exponential_diff)*carbon_diff
    return exponential_diff


############################

#Existing_Cropland_Soil_C = 24  # Current Soil C stock 24 tC/ha

expo_diff = Soil_exponential_Init(merged_data['Year'][0],21000) # puting a random high years value to get the curve
expo_soil_cropland = pd.DataFrame(expo_diff, columns=['expo_diff'])
Ini_Crop_Area=merged_data['LU_Area_ha'].iloc[0]
Ini_Crop_Soil_C=Ini_Crop_Area*Soil_C
expo_soil_cropland['Soil_C_Seq']=expo_soil_cropland['expo_diff']*Ini_Crop_Soil_C

soil_Age = 20 ## This to define what is the current age of the soil 
                #so that it matches with the existing soil C sequestration
                # Since current carbon stock is 24 tC/ha, this will take some more years to reach the steady-state 
                # conditions of 27 tc/ha.. So we assumed some age for simplicity
Ini_crop_Soil_C_Seq=expo_soil_cropland['Soil_C_Seq'].iloc[soil_Age-1 : soil_Age+(len(merged_data)-1)]
Ini_crop_Soil_C_Seq=np.array(-abs(Ini_crop_Soil_C_Seq)) 


############# Calculating individual arrays for exponential seqestration..... 


vec2=[]
c=0
for i,j in zip (merged_data['land_diff'],merged_data['Soil_Carbon_Stock']):
    sy=merged_data['Year'][c]
    endy=merged_data['Year'][len(merged_data)-1]
    if i==0:
         vec2.append(0)   
    elif i<0:
        b=(Soil_exponential(sy, endy,j))
        vec2.append(b)
    elif i>0:
        b=(Soil_exponential(sy, endy,j))
        vec2.append(b)
    c+=1




###########################  Calculating Soil Emission #################

first_year=merged_data['Year'][0]
empty_list=list(np.zeros(len(merged_data)))
Soil_df=pd.DataFrame({'Year_{}'.format(first_year):empty_list})

c=0
for i in (merged_data['Year'].iloc[0 :]):
    Soil_df.loc[:, 'Year_{}'.format(i)] =pd.Series(vec2[c])
    c+=1

#Soil_df.drop(columns=Soil_df.columns[0], axis=1, inplace=True)

c=0
i = merged_data['Year'][0] #1
for j in Soil_df.columns:
    Soil_df["Year_{}".format(i)] = Soil_df[j].shift(periods=c, fill_value=0)
    i += 1
    c+=1
Soil_df=Soil_df.fillna(0)

Soil_df["Soil_Seq"]=Soil_df.sum(axis=1)
Soil_df["Ini_Crop_Soil_Seq"]=Ini_crop_Soil_C_Seq
Soil_df["Area"]=merged_data['LU_Area_ha']
Soil_df["Area_Change"]=merged_data['land_diff']
Soil_df["Total_Soil_Emission"]=Soil_df["Soil_Seq"]+Soil_df["Ini_Crop_Soil_Seq"]


########## Calculating the carbon density values of individual years...

Initial_Cropland_Soil_C_Stock = merged_data['LU_Area_ha'].iloc[0]*24  # Mean Carbon Stock in India 
Total_Soil_C_Stock=[]
for i in Soil_df["Total_Soil_Emission"]:
    if i <0:
        total_stock=(Initial_Cropland_Soil_C_Stock)+(abs(i))
        Total_Soil_C_Stock.append(total_stock)
        Initial_Cropland_Soil_C_Stock=total_stock
    if i >0:
        total_stock=(Initial_Cropland_Soil_C_Stock)-(i)
        Total_Soil_C_Stock.append(total_stock)
        Initial_Cropland_Soil_C_Stock=total_stock
        
Soil_df["Total_C_Stock"]=np.array(Total_Soil_C_Stock)


############# Preapring Final Dataframe ###########
        
Cropland_C_Seq = pd.DataFrame(Year[0 :],columns=['Year'])
Cropland_C_Seq['Carbon_Seq']=(np.array(Soil_df['Total_Soil_Emission']))/1000000 # Convert to (MtC)  
Cropland_C_Seq['Carbon_Stock']=Soil_df['Total_C_Stock'] / 1000000 # Convert to (MtC)  

Cropland_C_Seq['Carbon_Seq_MtCO2']=(Cropland_C_Seq['Carbon_Seq'])*3.67
Cropland_C_Seq['Carbon_Stock_MtCO2']=(Cropland_C_Seq['Carbon_Stock'])*3.67

Cropland_C_Seq['Carbon_Stock_MtCO2']=(Cropland_C_Seq['Carbon_Stock'])*3.67

Cropland_C_Seq.to_csv("D:\GCAM_Output_Analysis\GCAM_Crop_CarbonSeq\Cropland_Soil_Seq\Output\Cropland_Soil_Emissions_Sequestration_NZ_FP_2070.csv", index=False)

################  Plotting the graphs


fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Line plot for Carbon Sequestration
axs[0].plot(Cropland_C_Seq['Year'], Cropland_C_Seq['Carbon_Seq_MtCO2'], marker='o', color='b')
axs[0].set_title('Carbon Sequestration Over Years')
axs[0].set_xlabel('Year')
axs[0].set_ylabel('Carbon Sequestration (MtCO2)')

# Bar plot for Carbon Stock
axs[1].bar(Cropland_C_Seq['Year'], Cropland_C_Seq['Carbon_Stock_MtCO2'], color='g')
axs[1].set_title('Carbon Stock Over Years')
axs[1].set_xlabel('Year')
axs[1].set_ylabel('Carbon Stock (MtCO2)')

plt.tight_layout()
plt.show()
        

print ("Processing Over!")

