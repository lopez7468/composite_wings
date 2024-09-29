


import numpy as np
from composite_theory_functions import *
import pandas as pd



myE11 = 125744                       # Laminate E1 [MPa]  
myE22 = 10030                        # Laminate E2 [MPa]
myNu12 = 0.271                       # Laminate nu12 [-]   
myG12 = 5555                         # Laminate G12 [MPa]                                                                  
        


myLaminate = [45,-45,0,90,90,0,-45,45]
myPlyNumber = len(myLaminate)
myLaminateThickness = 0.125


A11,A12,A16,A22,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66 = CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber)
        
        
print(A11,A12,A16,A22,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66)