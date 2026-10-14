import numpy as np
V_A2 = [0.750,0.700,0.700]
V_B2 = [1.390,1.310,1.310]
V_suma = [np.sum(V_A2), np.sum(V_B2)]
IA2 = 0.04225
IB2 = 0.08500
R_exp = V_suma[0]/IA2 #cambia la corriente
print(R_exp,V_suma[0]) #Muestra el Resultado experimental ademas de al suma del voltaje
R_teo = 63.3
E_R = (np.abs(R_exp-R_teo))/R_teo*100
print(E_R) #muestra el error

V_final = 2.50 #2.50 o 5.00
E_V =(np.abs(V_suma[0]-V_final))/V_final*100
print(E_V)
