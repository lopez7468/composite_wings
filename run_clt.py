
import numpy as np
from composite_theory_functions import *
from my_clt_functions import *
import pandas as pd


material_lib = pd.read_csv('material_library.csv')
material_lib.set_index('CLT_var', inplace= True)

#schedule and angles have to be in a dictionary of the different possible schedules

all_schedules2 = pd.DataFrame()

all_schedules2['Angles'] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
all_schedules2['schedule'] = [3, 3, 3, 3, 3,3, 3, 3, 3, 3,5,3,3,3,3,3,3, 3, 3, 3, 3,]

all_schedules= {
 1: all_schedules2, 
 }
 


#schedule = [3, 3, 3, 3, 3,3,3,3,3,3]
#angles = [0,0,0,0,0,0,0,0,0,0,0]
surface_area = 576 # in^2 

# get schedule + laminae props

schedul_props = {}

i = 1

for schedule in all_schedules:
    for index, row in all_schedules[schedule].iterrows():

        schedul_props[i] = [float(row['schedule']),
                            float(row['Angles']),
                            0,
                            0,
                            #float(material_lib[str(row['schedule'])]['Alpha_1']),
                            #float(material_lib[str(row['schedule'])]['Alpha_2']),
                            float(material_lib[str(row['schedule'])]['E_11']) *6894.76,
                            float(material_lib[str(row['schedule'])]['E_12'])*6894.76,
                            float(material_lib[str(row['schedule'])]['V_12']),
                            float(material_lib[str(row['schedule'])]['G_12'])*6894.76,
                            float(material_lib[str(row['schedule'])]['thickness'])*25.4,
                            float(material_lib[str(row['schedule'])]['density'] )] 
        #schedul_props[i] = layer_list
        i = i+1


d11 = myCLT(schedul_props)[12]



# max deflection in 

b = 419 # mm or 16.5 in
l = 457.2 # mm or 18 in
p = 88.96 # N or 230 lb

d_m = (p*l**3)/(3*d11*b)

print(d_m/25.4)
