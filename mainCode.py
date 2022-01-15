# Code written by Abhinav Sanjeeva Prasad
# Third Year Electrical and Computer Engineering Student
# Summer Research Student at AP2D Labs
# Finished by: 27th August 2020
# The algorithm below optimises the geometry of a plasmon sensor having a single groove.
# The algorithm varies the depth and the width of the groove and calculates the maximum intensity




# The main code starts from here

import math
import cmath
from cmath import sqrt
dielectricConstantAir = 8.8542*pow(10,-12)

#This function returns the x component of the electric field
# This value is obtained from COMSOL, it will be used to calculate the value of constant D
def get_initial_Ex_value(wavelength):
    wavelength = wavelength/pow(10,-9)
    if wavelength == 532:
        return complex(0.63075, -1.0013)
    elif wavelength == 633:
        return complex(1.5435,-4.3357)
    elif wavelength == 785:
        return complex(0.80339,3.1875)
    elif wavelength == 850:
        return  complex(4.8944,-8.1088)
    elif wavelength == 900:
        return  complex(6.1784,4.5221)
    elif wavelength == 1000:
        return complex(0.48349,-0.024376)


def get_angular_frequency(wavelength):
    speedOfLight = 299792458
    angularFrequency = 2*math.pi*speedOfLight/(wavelength)
    return angularFrequency

def get_dielectric_constant_of_the_metal(metal, wavelength):
    
    wavelength = wavelength/pow(10,-9)
    if metal == "Ag":
        if wavelength == 532:
            return complex(-11.55, 0.37038)
        elif wavelength == 633:
            return complex(-18.295, 0.48085)
        elif wavelength == 785:
            return complex(-29.789, 0.37611)
        elif wavelength == 850:
            return  complex(-35.585, 0.47724)
        elif wavelength == 900:
            return  complex(-40.590, 0.50969)
        elif wavelength == 1000:
          return  complex(-50.629, 0.56924)

    
    elif metal == "Au":
        if wavelength == 532:
            return complex(-4.6810, 2.4266)
        elif wavelength == 633:
            return complex(-11.753, 1.2596)
        elif wavelength == 785:
            return complex(-22.855, 1.4245)
        elif wavelength == 850:
            return  complex(-28.269, 1.7456)
        elif wavelength == 900:
            return  complex(-32.719, 1.9955)
        elif wavelength == 1000:
          return  complex(-41.849, 2.9477)   


def computing_k_zero(wavelength):
    return 2*math.pi/(wavelength)

def computing_neff(grooveWidth, metal, wavelength, grooveDepth):
    angleOfIncidence = 0
    diffractiveModeNum = 1
    m = 1  #m is the mode number and it is one for now
           #m depends on the wavelength and depth of the groove 
           #m is equal to 1 for the wavelength between 532nm to 1000nm and is 1 for the depth between 20nm to 200nm
    theta = 0 #assuming that the angle of incident is 0
    neff = (2*m-1)*wavelength/(4*grooveDepth*math.cos(theta))
    return neff

def computing_beta(grooveWidth, metal, wavelength, grooveDepth):
    
    #This is the first approach to calculate beta
    #neff = computing_neff(grooveWidth, metal, wavelength, grooveDepth)
    #kZero = computing_k_zero(wavelength)
    #return neff*kZero
   
    # This is the second approach
    m = 1 #m is the mode number see "computing_neff" to know more about m.
    return m*cmath.pi/grooveDepth

def computing_k_one(grooveWidth, metal, wavelength, grooveDepth):
    beta = computing_beta(grooveWidth, metal, wavelength, grooveDepth)
    kZero = computing_k_zero(wavelength)
    dielectricConstMetal = get_dielectric_constant_of_the_metal(metal, wavelength)
    return cmath.sqrt(pow(beta,2) - dielectricConstMetal*pow(kZero,2))

def computing_D(grooveWidth, metal, wavelength, grooveDepth):
    initialEx = get_initial_Ex_value(wavelength)
    angularFrequency = get_angular_frequency(wavelength)
    kOne = computing_k_one(grooveWidth, metal, wavelength, grooveDepth)
    dielectricConstMetal = get_dielectric_constant_of_the_metal(metal, wavelength)
    a = grooveWidth/2
    beta = computing_beta(grooveWidth, metal, wavelength,grooveDepth)
    #denominator = beta*cmath.sinh(kOne*a)
    #if denominator == 0:

    #in denominator, it should be sinh instead of cosh, refer to Moein's PHD thesis page 144 on calculation of constant D
    #using sinh, the denominator becomes 0 that causes an error, however, since D is a constant, in terms of optimsing, it should play a significant role
    denominator = (beta*cmath.cosh(kOne*a))

    return initialEx*angularFrequency*kOne*dielectricConstMetal*dielectricConstantAir*a/denominator 

def computing_E2(grooveWidth, metal, wavelength, grooveDepth):
    D = computing_D(grooveWidth, metal, wavelength, grooveDepth)
    angularFrequency = get_angular_frequency(wavelength)
    dielectricConstMetal = get_dielectric_constant_of_the_metal(metal, wavelength)
    kOne = computing_k_one(grooveWidth, metal, wavelength, grooveDepth)
    beta = computing_beta(grooveWidth, metal, wavelength, grooveDepth)
    x = 0 #in the case of a single groove, E2 is calculated at the midpoint of the surface of the groove so x cordinate is 0
    return abs(2*D*cmath.sqrt(pow(kOne*cmath.sinh(kOne*x), 2) + pow(beta*cmath.cosh(kOne*x), 2))/(angularFrequency*dielectricConstantAir*dielectricConstMetal))

def main():
    #retrieving the value of wavelength from the user and checking for errors
    isValidWavelength = False
    wavelength = int(input("Enter one of the following wavelength in nm: 532, 633, 785, 850, 900, 1000 -> "))
    
    while isValidWavelength != True:
        if wavelength != 532 and wavelength != 633 and wavelength != 785 and wavelength != 850 and wavelength != 900 and wavelength != 1000 :
            wavelength = int(input("Invalid wavelength entered! Please try again by entering one of the following wavelength in nm: 532, 633, 785, 850, 900, 1000 -> "))
        else: isValidWavelength = True
    
    wavelength = wavelength*pow(10,-9) #converting wavelength in metres from nm

    #asking the user for which metal to be used and checking for errors
    isValidMetal = False
    metal = input("Choose one of the following metals by entering its short form: Ag or Au -> ")
    while isValidMetal != True:
        if metal != "Ag" and metal != "Au" :
            metal = input("Invalid metal chosen! Please try again by choosing one of the following metals by entering its short form: Ag or Au -> ")
        else: isValidMetal = True
    
    intensity = 0
    maxIntensity = 0
    widthAtMaxInt = 0
    depthAtMaxInt = 0
    #Looping through various geometrical configurations of the sensor by varying the groove width and depth
    for grooveWidth in range(5, 151):
        for grooveDepth in range(20, 201):
            grooveWidth = grooveWidth*pow(10,-9)        #Converting the depth and width in metres
            grooveDepth = grooveDepth*pow(10,-9)
            
            intensity = pow(computing_E2(grooveWidth,metal,wavelength,grooveDepth),2)
            if intensity > maxIntensity :
                maxIntensity = intensity
                widthAtMaxInt = grooveWidth
                depthAtMaxInt = grooveDepth
    
    #denominator = cmath.sinh(computing_k_one(widthAtMaxInt,metal,wavelength,depthAtMaxInt)*widthAtMaxInt/2)
    #print(complex(denominator))
    #beta = computing_beta(widthAtMaxInt, metal, wavelength, grooveDepth)
    #print(float(beta))

    #converting it into nano metres
    widthAtMaxInt = widthAtMaxInt/pow(10,-9)
    depthAtMaxInt = depthAtMaxInt/pow(10, -9)
    wavelength = wavelength/pow(10,-9)
    print("The maximum electricomagnetic field intensity for the metal: " ,metal,  "at the wavelength of: ",wavelength, "nm is: ", maxIntensity, " ")
    print ("The groove width at the maximum intensity is: ", widthAtMaxInt, "nm")
    print("The groove depth at maximum intensity is: ", depthAtMaxInt, "nm")



if __name__ == "__main__":
    main()   

'''
#The code below is used for finding a root for a non linear function

from scipy import optimize

import cmath
import math

def f(x):
    return cmath.tanh(5*pow(10,-8)*cmath.sqrt(pow(x,2) - 4*pow(cmath.pi,2)*3.53*pow(10,12)*11.56*8.85*pow(10,-12))) + cmath.sqrt(pow(x,2) - 4*pow(cmath.pi,2)*3.53*pow(10,12)*8.85*pow(10,-12))*11.56*8.85*pow(10,-12)/(8.85*pow(10,-12)*cmath.sqrt(pow(x,2) - 4*pow(cmath.pi,2)*3.53*pow(10,12)*11.56*8.85*pow(10,-12)))
#x is beta

root = optimize.ridder(f,0, 10000000000)
print(complex(root))
'''




'''
#The code below is another algorithm used for finding a root for a non linear function


#source for concept: http://hplgit.github.io/prog4comp/doc/pub/p4c-sphinx-Python/._pylight007.html
#source for code: https://github.com/hplgit/prog4comp/blob/master/src/py/brute_force_root_finder_flat.py

from numpy import linspace, exp, cos
import numpy
import math

def f(x):
    return math.tanh(5*pow(10,-8)*math.sqrt(pow(x,2) - 4*pow(math.pi,2)*3.53*pow(10,12)*11.56*8.85*pow(10,-12))) + math.sqrt(pow(x,2) - 4*pow(math.pi,2)*3.53*pow(10,12)*8.85*pow(10,-12))*11.56*8.85*pow(10,-12)/(8.85*pow(10,-12)*math.sqrt(pow(x,2) - 4*pow(math.pi,2)*3.53*pow(10,12)*11.56*8.85*pow(10,-12)))

x = linspace(0, 4, 10001)

y = f(x)
x = np.array(x)
y = np.array(y)

root = None  # Initialization
for i in range(len(x)-1):
    if y[i]*y[i+1] < 0:
        root = x[i] - (x[i+1] - x[i])/(y[i+1] - y[i])*y[i]
        break  # Jump out of loop
    elif y[i] == 0:       
        root = x[i]
        break  # Jump out of loop

if root is None:
    print ('Could not find any root in [%g, %g]' % (x[0], x[-1]))
else:
    print('Find (the first) root as x=%g' % root)


'''

#Check this out
# https://stackoverflow.com/questions/15213141/how-to-do-nonlinear-complex-root-finding-in-python
