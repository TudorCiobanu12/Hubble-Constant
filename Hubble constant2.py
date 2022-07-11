 # -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 15:43:09 2020

@author: Aero
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

'READING FILES SECTION'
parallax, err_par, Period, m, A, err_A =np.loadtxt('MW_Cepheids.txt', unpack=True\
                                                  ,dtype=float,comments="#"\
                                                  ,usecols=(1,2,3,4,5,6))
    
logP1, m1 = np.loadtxt('hst_gal1_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP2, m2 = np.loadtxt('hst_gal2_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP3, m3 = np.loadtxt('hst_gal3_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP4, m4 = np.loadtxt('hst_gal4_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP5, m5 = np.loadtxt('hst_gal5_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP6, m6 = np.loadtxt('hst_gal6_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP7, m7 = np.loadtxt('hst_gal7_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
logP8, m8 = np.loadtxt('hst_gal8_cepheids.txt', unpack=True, dtype=float, comments='#', usecols=(1,2))
Recession, Extinction =np.loadtxt('galaxy.txt',unpack=True, dtype=float, comments='#',usecols=(1,2))

'STEP 1 FUNCTIONS'
def getDistance(parallax): #distance Earth_star
      distance=1000/parallax     
      return distance

def absolute_magnitude(m,A,distance): #absolute magnitude of each Cepheid
    M=m-5*(np.log10(distance))+5-A
    return M


'CURVE FIT FUNCTIONS'
def func(x,a,b):#using curve fit to find out alpha and beta which are the slope and intercept
    return a*(x)+b

def func2(x,b): #using curve fit to find out the hubble constant which is the slope 
    return x*(b)

def func3(x,a):
    return x*(a)

'ERROR PROPOGATION FUNCTIONS'
def err_on_magnitude(err_A,sigma_distance,m,A): #the error on the absolute magnitude involves the error propagation formula because there are independent values in the formula for M   
    sigma_magnitude=np.sqrt(25*(1/(distance*np.log(10))*sigma_distance)**2+err_A**2)
    return sigma_magnitude 

def err_on_distance(err_par,parallax): #error on distance Earth-star for each Cephied 
    err_d_pc=(err_par/parallax**2)*1000 #the formula from the lab manual for error propagation
    return err_d_pc

def error_one(distance1,slopeErr,intErr,alfa,beta): #finding the formula for error on distance for each galaxy using error propagation because the formula contains many independent values
   
    error_mean_m1=np.std(m1)/np.sqrt(len(m1)) #error on the mean for apparent magnitude, where the np.std function computes the standard deviation
    error_mean_logP1=np.std(logP1)/np.sqrt(len(logP1)) #error on the mean for logPs
    
    a=distance1/10**6-alfa*np.mean(logP1)-beta+5-0.0992 #the constant "a" is meant to simplify the the error propagation equation for not making it harder to read
    
    error_d1=np.sqrt((0.465017*distance1/10**6*error_mean_m1)**2+(0.460517*np.exp(-0.460517*a)*slopeErr)**2+(0.460517*np.exp(-0.460517*a)*error_mean_logP1)**2+(0.460517*np.exp(0.460517*a)*intErr)**2)
    return error_d1
    
'RUN THE FUCNTION FOR ALL EIGHT GALAXIES TO FIND OUT THE ERROR ON DISTANCE'
def error_two(distance2,slopeErr,intErr,alfa,beta): 
   
    error_mean_m2=np.std(m2)/np.sqrt(len(m2))
    error_mean_logP2=np.std(logP2)/np.sqrt(len(logP2))      
    b=distance2/10**6-alfa*np.mean(logP2)-beta+5-0.0434     
    error_d2=np.sqrt((0.465017*distance2/10**6*error_mean_m2)**2+(0.460517*np.exp(-0.460517*b)*slopeErr)**2+(0.460517*np.exp(-0.460517*b)*error_mean_logP2)**2+(0.460517*np.exp(0.460517*b)*intErr)**2)
    return error_d2

def error_three(distance3,slopeErr,intErr,alfa,beta): 
   
    error_mean_m3=np.std(m3)/np.sqrt(len(m3))
    error_mean_logP3=np.std(logP3)/np.sqrt(len(logP3))      
    c=distance3/10**6-alfa*np.mean(logP3)-beta+5-0.0775      
    error_d3=np.sqrt((0.465017*distance3/10**6*error_mean_m3)**2+(0.460517*np.exp(-0.460517*c)*slopeErr)**2+(0.460517*np.exp(-0.460517*c)*error_mean_logP3)**2+(0.460517*np.exp(0.460517*c)*intErr)**2)
    return error_d3

def error_four(distance4,slopeErr,intErr,alfa,beta): 
   
    error_mean_m4=np.std(m4)/np.sqrt(len(m4))
    error_mean_logP4=np.std(logP4)/np.sqrt(len(logP4))      
    d=distance4/10**6-alfa*np.mean(logP4)-beta+5-0.0682     
    error_d4=np.sqrt((0.465017*distance4/10**6*error_mean_m4)**2+(0.460517*np.exp(-0.460517*d)*slopeErr)**2+(0.460517*np.exp(-0.460517*d)*error_mean_logP4)**2+(0.460517*np.exp(0.460517*d)*intErr)**2)
    return error_d4

def error_five(distance5,slopeErr,intErr,alfa,beta): 
   
    error_mean_m5=np.std(m5)/np.sqrt(len(m5))
    error_mean_logP5=np.std(logP5)/np.sqrt(len(logP5))      
    e=distance5/10**6-alfa*np.mean(logP5)-beta+5-0.0558   
    error_d5=np.sqrt((0.465017*distance5/10**6*error_mean_m5)**2+(0.460517*np.exp(-0.460517*e)*slopeErr)**2+(0.460517*np.exp(-0.460517*e)*error_mean_logP5)**2+(0.460517*np.exp(0.460517*e)*intErr)**2)
    return error_d5

def error_six(distance6,slopeErr,intErr,alfa,beta): 
   
    error_mean_m6=np.std(m6)/np.sqrt(len(m6))
    error_mean_logP6=np.std(logP6)/np.sqrt(len(logP6))      
    f=distance6/10**6-alfa*np.mean(logP6)-beta+5-0.0806      
    error_d6=np.sqrt((0.465017*distance6/10**6*error_mean_m6)**2+(0.460517*np.exp(-0.460517*f)*slopeErr)**2+(0.460517*np.exp(-0.460517*f)*error_mean_logP6)**2+(0.460517*np.exp(0.460517*f)*intErr)**2)
    return error_d6

def error_seven(distance7,slopeErr,intErr,alfa,beta): 
   
    error_mean_m7=np.std(m7)/np.sqrt(len(m7))
    error_mean_logP7=np.std(logP7)/np.sqrt(len(logP7))      
    g=distance7/10**6-alfa*np.mean(logP7)-beta+5-0.1736     
    error_d7=np.sqrt((0.465017*distance7/10**6*error_mean_m7)**2+(0.460517*np.exp(-0.460517*g)*slopeErr)**2+(0.460517*np.exp(-0.460517*g)*error_mean_logP7)**2+(0.460517*np.exp(0.460517*g)*intErr)**2)
    return error_d7

def error_eight(distance8,slopeErr,intErr,alfa,beta): 
   
    error_mean_m8=np.std(m7)/np.sqrt(len(m8))
    error_mean_logP8=np.std(logP7)/np.sqrt(len(logP8))      
    h=distance8/10**6-alfa*np.mean(logP8)-beta+5-0.0434     
    error_d8=np.sqrt((0.465017*distance8/10**6*error_mean_m8)**2+(0.460517*np.exp(-0.460517*h)*slopeErr)**2+(0.460517*np.exp(-0.460517*h)*error_mean_logP8)**2+(0.460517*np.exp(0.460517*h)*intErr)**2)
    return error_d8

'DISATNCE FUNCTION SECTION'
def distance_one(alfa,logP1,beta,m1):
    distance_modulus1=10**(0.2*(np.mean(m1)-alfa*np.mean(logP1)-beta+5-0.0992)) #the formula for distance, the last number in the equation is the extinction of each galaxy
    return distance_modulus1

def distance_two(alfa,logP2,beta,m2):
    distance_modulus2=10**(0.2*(np.mean(m2)-alfa*np.mean(logP2)-beta+5-0.0434))
    return distance_modulus2 

def distance_three(alfa,logP3,beta,m3):
    distance_modulus3=10**(0.2*(np.mean(m3)-alfa*np.mean(logP3)-beta+5-0.0775))
    return distance_modulus3 
   
def distance_four(alfa,logP4,beta,m4):
    distance_modulus4=10**(0.2*(np.mean(m4)-alfa*np.mean(logP4)-beta+5-0.0682))
    return distance_modulus4 

def distance_five(alfa,logP5,beta,m5):
    distance_modulus5=10**(0.2*(np.mean(m5)-alfa*np.mean(logP5)-beta+5-0.0558))
    return distance_modulus5  

def distance_six(alfa,logP6,beta,m6):
    distance_modulus6=10**(0.2*(np.mean(m6)-alfa*np.mean(logP6)-beta+5-0.0806))
    return distance_modulus6 

def distance_seven(alfa,logP7,beta,m7):
    distance_modulus7=10**(0.2*(np.mean(m7)-alfa*np.mean(logP7)-beta+5-0.1736))
    return distance_modulus7 

def distance_eight(alfa,logP8,beta,m8):
    distance_modulus8=10**(0.2*(np.mean(m8)-alfa*np.mean(logP8)-beta+5-0.0434))
    return distance_modulus8   
                                                                                                                                                                                                                                                                                                                                                                               
print('______________STEP1______________')     
print()
distance=getDistance(parallax)
print("Distance Earth-star: ",distance ,"parsec")
print()
sigma_distance=err_on_distance(err_par, parallax)
print("Error on distance Earth-star: ", sigma_distance,"parsec")
print()
Magnitude=absolute_magnitude(m, A, distance)
print(" Absolute magnitude of each Cepheids: ", Magnitude)
print()
Error_magnitude=err_on_magnitude(err_A, sigma_distance, m, A)
print("Error on absolute magnitude: ", Error_magnitude)
print()

alfa=np.polyfit(np.log10(Period),Magnitude,1)[0] 
beta=np.polyfit(np.log10(Period),Magnitude,1)[1] 
print("alfa: ",alfa)
print("beta",beta)
log_period=np.log10(Period) 

'PLOTTING THE LOG_PERIOD AGAINST MAGNITUDE '
plt.errorbar(log_period,Magnitude,marker='o',linewidth=0,yerr=Error_magnitude,elinewidth=1,markersize=3)
plt.plot(log_period,alfa*log_period+beta, linewidth=2)
plt.xlabel("log_period")
plt.ylabel("Magnitude")
plt.show()
'CHI^2 FOR THE FIRST STEP'
y=alfa*np.log10(Period)+beta #model
y_m=Magnitude #observed value/result
sig_y=Error_magnitude #error on observed value/result
chi2 = np.sum(((y-y_m)**2)/(sig_y**2)) #model fitting check

print()
print("chi^2:",chi2)
print()
'REDUCED CHI^2'
reduced_chi=chi2/10 #chi^2 divided by the number of data points 
print("reduced CHI^2: ",reduced_chi)
print()
popt,pcov = curve_fit(func,log_period,Magnitude)  #finding the slope and intercept using curve_fit
solpe=popt[0] #alpha
bestInt=popt[1] #beta
slopeErr=np.sqrt(pcov[0,0])#error on alpha
intErr=np.sqrt(pcov[1,1]) #error on beta
print("Alpha: ", popt[0], "+/-",slopeErr)
print("Beta: ", popt[1], "+/-",intErr)
print()
print('________________STEP2_________________')
print()
distance1=distance_one(alfa, logP1, beta, m1)
print("distance Earth-NGC3627 IS:", distance1/10**6,"Mpc") #divide each distance by 10^6 to express it in megaparsec
err_d1=error_one(distance1, slopeErr, intErr,alfa,beta)
print("Error on distance1:  +/-", err_d1/10**6 ,"Mpc") #each error is divided by 10^6 because the error of a quantity must be in the same units as the quantity
print()
distance2=distance_two(alfa, logP2, beta, m2)
print("distance Earth-NGC3982 IS:", distance2/10**6 ,"Mpc")
err_d2=error_two(distance2, slopeErr, intErr,alfa,beta)
print("Error on distance2:  +/-", err_d2/10**6 ,"Mpc")
print()
distance3=distance_three(alfa, logP3, beta, m3)
print("distance Earth-NGC4496A IS:", distance3/10**6 ,"Mpc")
err_d3=error_three(distance3, slopeErr, intErr,alfa,beta)
print("Error on distance3:  +/-", err_d3/10**6 ,"Mpc")
print()
distance4=distance_four(alfa, logP4, beta, m4)
print("distance Earth-NGC4527 IS:", distance4/10**6 ,"Mpc")
err_d4=error_four(distance3, slopeErr, intErr,alfa,beta)
print("Error on distance4:  +/-", err_d4/10**6 ,"Mpc")
print()
distance5=distance_five(alfa, logP5, beta, m5)
print("distance Earth-NGC4536 IS:", distance5/10**6 ,"Mpc")
err_d5=error_four(distance5, slopeErr, intErr,alfa,beta)
print("Error on distance5:  +/-", err_d5/10**6 ,"Mpc")
print()
distance6=distance_six(alfa, logP6, beta, m6)
print("distance Earth-NGC4639 IS:", distance6/10**6 ,"Mpc")
err_d6=error_six(distance6, slopeErr, intErr,alfa,beta)
print("Error on distance6:  +/-", err_d6/10**6 ,"Mpc")
print()
distance7=distance_seven(alfa, logP7, beta, m7)
print("distance Earth-NGC5253 IS:", distance7/10**6 ,"Mpc")
err_d7=error_seven(distance7, slopeErr, intErr,alfa,beta)
print("Error on distance7:  +/-", err_d7/10**6 ,"Mpc")
print()
distance8=distance_eight(alfa, logP8, beta, m8)
print("distance Earth-IC4182 IS:", distance8/10**6 ,"Mpc")
err_d8=error_eight(distance8, slopeErr, intErr,alfa,beta)
print("Error on distance8:  +/-", err_d8/10**6 ,"Mpc")
print()
print()

print()
print("____________________STEP3___________________")
print()

Dgal=[] #introducing all distances into an array
Dgal.append(distance1)
Dgal.append(distance2)
Dgal.append(distance3)
Dgal.append(distance4)
Dgal.append(distance5)
Dgal.append(distance6)
Dgal.append(distance7)
Dgal.append(distance8)
Dgal=np.array(Dgal)

Err_gal=[] #introducing all the errors on distance into an array
Err_gal.append(err_d1)
Err_gal.append(err_d2)
Err_gal.append(err_d3)
Err_gal.append(err_d4)
Err_gal.append(err_d5)
Err_gal.append(err_d6)
Err_gal.append(err_d7)
Err_gal.append(err_d8)
Err_gal=np.array(Err_gal)


popt,pcov= curve_fit(func2,Dgal/10**6,Recession)  #fiding the hubble constant using the curve fit for Dgal against Recession velocity
solpe=popt[0]
slopeErr1=np.sqrt(pcov[0,0])
print("Hubble constant: ", popt[0],"+/-", slopeErr1, "km/s/megaparsec")


plt.errorbar(Recession,Dgal,marker='o',linewidth=0,yerr=Err_gal,elinewidth=1,markersize=3)
plt.plot(Recession,Dgal, linewidth=0,markersize=8)
plt.ylabel("Dgal")
plt.xlabel("Recession velocity")
plt.show()

print("____________________STEP4_____________________")
print()
Err_age_universe=slopeErr1/(popt[0])**2# error propagation for the age of the universe using the lab formula
popt1,pcov= curve_fit(func3,Dgal/10**6,Recession)   #finding the Hubble constant
solpe=popt1[0]
print("AGE OF THE UNIVERSE: ", 1/popt1[0],"+/-", Err_age_universe,"sec*megaparsec/km")
print()
print("The unit for the age of the universe is sec*megaparsec/km, therefore a conversion has to be made in order to express it in 10^9 years(billions).")
print()
print("1 sec*megapasec/km= 9.7781310*10^11 years")

age_universe=1/popt[0] #formula for the age of the universe where the popt[0] is the Hubble constant
x=9.7781310*10**11 
Actual_age=age_universe*x #conversion for the age
real_err=Err_age_universe*x #conversion for error on age
print() 
print("AGE OF THE UNIVERSE:", Actual_age/10**9,"+/-",real_err/10**9,"*10^9 years (gigayears)") #divide the age of the universe and the 



