# -*- coding: utf-8 -*-
"""
Created on Thu May 14 19:28:15 2020

@author: Oscar Jackson
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
plt.ioff()



#defining liner function with intercept 0
intercept = 0
def func2(x, slope): 
    line =  slope*x + intercept
    return line

#defining liner function
def func1(x, slope,intercept): 
    line = intercept + slope*x
    return line

print("")
print("")
print("The Cepheid Period-Luminosity Relation:")
print("==============================================================================")
#load the data from MW_Cepheids.dat
parallax, parallax_error, Period, m, A, A_error  = np.loadtxt('MW_Cepheids.dat',\
                                                                               unpack=True,\
                                                                               usecols=(1,2,3,4,5,6),\
                                                                               dtype=float)
#convert period array to log10 period
logPeriod = np.log10(Period)

#convert parallax to parsec
temp =np.log10(1000/parallax)
#calculate absolute magnitudes using the distance, apparent magnitudes and extinctions
M = m - 5*temp + 5 - A
#calculate the error for the absolute magnitudes, by propagting the errors
M_error = ((((5*parallax_error)/(np.log(10)*parallax))**2)+((A_error)**2))**0.5

#Calculate the degree of freedom
dof1 = len(logPeriod)-2
#perform the curvefit
param, cova_mat = opt.curve_fit(f=func1, xdata=logPeriod, ydata=M, sigma=M_error, p0 = (-2,-2), absolute_sigma=True)
#define the best slope and gradient
alpha = param[0]
beta = param[1]
#find vaules corresponding to best slope and intercept
y_best = alpha*logPeriod + beta 


#find the chi2
best_chi2 = np.sum(((M-y_best)**2.0)/(M_error**2.0))
redu_chi2 = best_chi2/dof1

#Check if best fit is consistent with data
redu_chi2_check_low = 1 - ((2/dof1)**0.5)
redu_chi2_check_high = 1 + ((2/dof1)**0.5)


#plot the first plot
plt.subplot(2,1,1)
plt.xlabel('logPeriod')
plt.ylabel('Absolute Magnitude')
plt.title('The Cepheid Period-Luminosity Relation')
plt.plot(logPeriod,M,marker="o", markersize=5, linestyle="None")
plt.errorbar(logPeriod,M,yerr=M_error,fmt=' ')
plt.plot(logPeriod,y_best, '-r')
#invert as magnitudes are negative
plt.gca().invert_yaxis()
plt.show()

#print out vaules
print("Best Slope: "+ str(alpha) + " Error: " + str(cova_mat[0][0]))
print("Best Intercept: "+ str(beta)  + " Error: " + str(cova_mat[1][1]))
print("Best chi2: " + str(best_chi2))
print("Best reduced chi2: " + str(redu_chi2))
print("Is chi2 consistent? From " + str(redu_chi2_check_low) + " to " + str(redu_chi2_check_high))
print("==============================================================================")
print("")
print("")






print("Distances of Nearby Galaxies")
print("==============================================================================")

d = np.zeros(8)
d_error = np.zeros(8)

#get galaxy name
g_name = np.loadtxt('galaxy_data.dat', unpack=True, usecols=(0), dtype=str)
#load galaxies' recession velocity and magnitude
v , A = np.loadtxt('galaxy_data.dat', unpack=True, usecols=(1,2), dtype=float)


for i in range(1,9):
    #loop to load data for each galaxies 
    logP, m  = np.loadtxt('hst_gal'+str(i)+'_cepheids.dat', unpack=True, usecols=(1,2), dtype=float) 
    #calculate absoulte magnitude of the set of cepheids varibles
    M =  alpha*logP + beta
    #calculate distance of the set of cepheids varibles
    Ds = 10**((m-M+5-A[i-1])/5)
    #equation to propagate errors for the set od distances
    temp1 = -0.2*np.log(10)*(np.e**((np.log(10)*M*-1)/5))*10**((m-A[i-1]+5)/5)*(((((cova_mat[0][0]**0.5)*logP)**2)+(((cova_mat[1][1]**0.5)**0.5))))
    #find weighted mean 
    D = np.sum(Ds/temp1)/np.sum(1/temp1)
    #set galaxy distance
    d[i-1] = D
    #set galaxy distance error as the weighted mean of all the errors 
    d_error[i-1] = np.sum(Ds*temp1)/np.sum(Ds)
    #print name, distance and distance error of galaxy
    print(str(g_name[i-1]) + " Distance: " + str(D) + "Error: " + str(d_error[i-1]*-1))


#convert parsecs to mega parsecs
d = d*(10**-6)
d_error = d_error*(10**-6)

print("==============================================================================")
print("")
print("")



print("Calucating Hubble Constant (intercept = 0)")
print("==============================================================================")
#calculate the degrees of freedom.
dof1 = len(v)-2
dof2 = len(v)-1 #for set intercept


#preform curve fit with no forced intercept
param1, cova_mat1 = opt.curve_fit(f=func1, xdata=v, ydata=d,sigma = d_error, p0 = (0.015, 0 ), absolute_sigma=True)
slope1 = param1[0]
intercept1 = param1[1]
#find hubble constant
hubble1 = 1/slope1
#proagate error for the hubble constant
hubble1_error = ((cova_mat1[0,0])**0.5)/(slope1**2)


#find vaules corresponding to best slope and intercept
y_best1 = slope1*v + intercept1
#find the chi2
best_chi2 = np.sum(((d-y_best1)**2.0)/(d_error**2.0))
redu_chi2 = best_chi2/dof1
#Check if best fit is consistent with data
redu_chi2_check_low = 1 - ((2/dof1)**0.5)
redu_chi2_check_high = 1 + ((2/dof1)**0.5)

#print out vaules
print("With no forced intercept: " + str(hubble1) + " Error: " + str(hubble1_error))
print("Best chi2: " + str(best_chi2))
print("Best reduced chi2: " + str(redu_chi2))
print("Is chi2 consistent? From " + str(redu_chi2_check_low) + " to " + str(redu_chi2_check_high))
print("")


#preform curve fit with forced intercept at 0
param2, cova_mat2 = opt.curve_fit(f=func2, xdata=v, ydata=d,sigma = d_error, p0 = (0.015), absolute_sigma=True)
slope2 = param2[0]
#find hubble constant
hubble2 = 1/slope2
#proagate error for the hubble constant
hubble2_error = ((cova_mat2[0,0])**0.5)/(slope2**2)


#find vaules corresponding to best slope
y_best2 = slope2*v 
#find the chi2
best_chi2 = np.sum(((d-y_best2)**2.0)/(d_error**2.0))
redu_chi2 = best_chi2/dof2
#Check if best fit is consistent with data
redu_chi2_check_low = 1 - ((2/dof2)**0.5)
redu_chi2_check_high = 1 + ((2/dof2)**0.5)

#plot line of best fit
y_best3 = hubble2*d 
plt.subplot(2,1,2)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Recession velocity (km/s)')
plt.title('Hubbles Law')
plt.plot(d,v,marker="o", markersize=5, linestyle="None")
plt.errorbar(d,v,xerr=d_error,fmt=' ')
plt.plot(d,y_best3, '-r')
plt.show()

#print out vaules
print("With forced intercept: " + str(hubble2) + " Error: " + str(hubble2_error))
print("Best chi2: " + str(best_chi2))
print("Best reduced chi2: " + str(redu_chi2))
print("Is chi2 consistent? From " + str(redu_chi2_check_low) + " to " + str(redu_chi2_check_high))
print("==============================================================================")
print("")
print("")


print("Age of the Universe")
print("==============================================================================")
#Mpc to km conversion
Mpc_conv = 3.086*10**19
#convert slope (which is 1/h) to seconds
age = slope2 * Mpc_conv
#get error of slope and convert
age_error = Mpc_conv * cova_mat2[0,0]

#print vaules
print(str(age) + " error: " + str(age_error) + " seconds")
#print and convert to billions of years
print(str((age/(365*24*60*60))/1000000000) + " error: " + str((age_error/(365*24*60*60))/1000000000) + " billion years")
print("==============================================================================")
print("")
print("")



