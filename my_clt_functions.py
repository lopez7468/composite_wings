
import numpy as np

def CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber, schedule_dict):

#0   row['schedule'],
#1   row['Angles'],
#2   material_lib[str(row['schedule'])]['Alpha_1'],
#3  material_lib[str(row['schedule'])]['Alpha_2'], 
#4   material_lib[str(row['schedule'])]['E_11'],
#5   material_lib[str(row['schedule'])]['E_12'], 
#6   material_lib[str(row['schedule'])]['V_12'], 
#7   material_lib[str(row['schedule'])]['G_12'], 
#8   material_lib[str(row['schedule'])]['thickness'],
#9   material_lib[str(row['schedule'])]['density'] ]  






#def CLT(schedule_dict):
    #--------------------------------------------------------------------------
    #
    # This function calculates the ABD Composite Stiffness Matrix Components
    # Classical Lamination Theory (CLT)
    #
    #--------------------------------------------------------------------------

        
    
    Z = []
    Z.append(- sum(laminate_thickness)/2)
    
    laminate_thickness = []
    for item in schedule_dict:
        thickness = schedule_dict[item][9]

        laminate_thickness.append(thickness)

    

    for item in laminate_thickness:
        Z.append(Z[-1]+ item)
    
    # Calcuate Reduced Stiffnesses
    
    myQ11 = []
    myQ12 = []
    myQ16 = []
    myQ22 = []
    myQ26 = []
    myQ66 = []
    
    myQ11_S = []
    myQ12_S = []
    myQ16_S = []
    myQ22_S = []
    myQ26_S = []
    myQ66_S = []
    
    for item in schedule_dict:
        e_11 = schedule_dict[item][4]
        e_22 = schedule_dict[item][5]
        v_12 = schedule_dict[item][6]
        g_12 = schedule_dict[item][7]



        myQ11.append((e_11**2)/(e_11-v_12**2*e_22))
        myQ12.append((v_12*e_11*e_22)/(e_11-v_12**2*e_22))
        myQ16.append(0)
        myQ22.append((e_11*e_22)/(e_11-v_12**2*e_22))
        myQ26.append(0)
        myQ66.append(g_12)
    
    i=0    
    for item in schedule_dict:

        angle = schedule_dict[item][1]

        myQ11_S.append(myQ11[i]*np.cos(angle*np.pi/180.0)**4 + 2*(myQ12[i]+2*myQ66[i])*(np.cos(angle*np.pi/180.0))**2*(np.sin(angle*np.pi/180.0))**2 + myQ22[i]*(np.sin(angle*np.pi/180.0))**4)
        myQ12_S.append(myQ12[i]*((np.cos(angle*np.pi/180.0)**4)+(np.sin(angle*np.pi/180.0)**4))+ (myQ11[i]+myQ22[i]-4*myQ66[i])*np.cos(angle*np.pi/180.0)**2*np.sin(angle*np.pi/180.0)**2)
        myQ16_S.append((myQ11[i]-myQ12[i]-2*myQ66[i])*np.cos(angle*np.pi/180.0)**3*np.sin(angle*np.pi/180.0) - (myQ22[i] - myQ12[i] - 2*myQ66[i])*np.cos(angle*np.pi/180.0)*np.sin(angle*np.pi/180.0)**3)
        myQ22_S.append(myQ11[i]*np.sin(angle*np.pi/180.0)**4 + 2*(myQ12[i]+2*myQ66[i])*(np.cos(angle*np.pi/180.0))**2*(np.sin(angle*np.pi/180.0))**2 + myQ22[i]*(np.cos(angle*np.pi/180.0))**4)
        myQ26_S.append((myQ11[i]-myQ12[i]-2*myQ66[i])*np.cos(angle*np.pi/180.0)*np.sin(angle*np.pi/180.0)**3 - (myQ22[i] - myQ12[i] - 2*myQ66[i])*np.cos(angle*np.pi/180.0)**3*np.sin(angle*np.pi/180.0))
        myQ66_S.append((myQ11[i] + myQ22[i] - 2*myQ12[i] - 2*myQ66[i])*np.cos(angle*np.pi/180.0)**2*np.sin(angle*np.pi/180.0)**2 + myQ66[i]*(np.cos(angle*np.pi/180.0)**4+np.sin(angle*np.pi/180.0)**4))
    
    
    
    # Calcualte A Matrix
    #-------------
            
    A11_v = []
    A12_v = []
    A16_v = []
    A22_v = []
    A26_v = []
    A66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        A11_v.append(myQ11_S[j]*(Z[i]-Z[i-1]))
        A12_v.append(myQ12_S[j]*(Z[i]-Z[i-1]))
        A16_v.append(myQ16_S[j]*(Z[i]-Z[i-1]))
        A22_v.append(myQ22_S[j]*(Z[i]-Z[i-1]))
        A26_v.append(myQ26_S[j]*(Z[i]-Z[i-1]))
        A66_v.append(myQ66_S[j]*(Z[i]-Z[i-1]))
        i = i+1
    
    A11 = sum(A11_v)
    A12 = sum(A12_v)
    A16 = sum(A16_v)
    A22 = sum(A22_v)
    A26 = sum(A26_v)
    A66 = sum(A66_v)
    
    # Calcualte B Matrix
    #-------------
            
    B11_v = []
    B12_v = []
    B16_v = []
    B22_v = []
    B26_v = []
    B66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        B11_v.append(0.5*myQ11_S[j]*(Z[i]**2-Z[i-1]**2))
        B12_v.append(0.5*myQ12_S[j]*(Z[i]**2-Z[i-1]**2))
        B16_v.append(0.5*myQ16_S[j]*(Z[i]**2-Z[i-1]**2))
        B22_v.append(0.5*myQ22_S[j]*(Z[i]**2-Z[i-1]**2))
        B26_v.append(0.5*myQ26_S[j]*(Z[i]**2-Z[i-1]**2))
        B66_v.append(0.5*myQ66_S[j]*(Z[i]**2-Z[i-1]**2))
        i = i+1
    
    
    
    
    B11 = sum(B11_v)
    B12 = sum(B12_v)
    B16 = sum(B16_v)
    B22 = sum(B22_v)
    B26 = sum(B26_v)
    B66 = sum(B66_v)
    
    # Calcualte D Matrix
    #-------------
            
    D11_v = []
    D12_v = []
    D16_v = []
    D22_v = []
    D26_v = []
    D66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        D11_v.append((1/3.0)*myQ11_S[j]*(Z[i]**3-Z[i-1]**3))
        D12_v.append((1/3.0)*myQ12_S[j]*(Z[i]**3-Z[i-1]**3))
        D16_v.append((1/3.0)*myQ16_S[j]*(Z[i]**3-Z[i-1]**3))
        D22_v.append((1/3.0)*myQ22_S[j]*(Z[i]**3-Z[i-1]**3))
        D26_v.append((1/3.0)*myQ26_S[j]*(Z[i]**3-Z[i-1]**3))
        D66_v.append((1/3.0)*myQ66_S[j]*(Z[i]**3-Z[i-1]**3))
        i = i+1
    
    
    D11 = sum(D11_v)
    D12 = sum(D12_v)
    D16 = sum(D16_v)
    D22 = sum(D22_v)
    D26 = sum(D26_v)
    D66 = sum(D66_v)
    return A11,A12,A16,A22,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66


#-----------------------------------------------------------------------------------------------

# def CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber):



