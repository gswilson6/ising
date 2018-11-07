from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
from math import*


def main():
    N = 15#lattice dimension
    iterations = 10000#monte carlo iterations
    kB = 1#boltzmann's constant
    n = 20#number of temperature and h data points
    J = 1#coupling constant

    #reading to a text file for the phase diagram
    f = open('phasediag.txt', 'w')

    #Generates lists of zeros to append the calculated data to
    En = np.zeros(n)
    Cv = np.zeros(n)
    M = np.zeros(n)
    Sus = np.zeros(n)
    hlist = np.zeros(n)

##    #Generates lists for the equilibrium check
##    steps = np.zeros(iterations)
##    Elist = np.zeros(iterations)

    #h and t lists
    h = np.linspace(-0.5,0.5,n)
    T = np.linspace(0.1,5,n)

    #running over h values
    for q in range(len(h)):

        #running over T values
        for l in range(len(T)):

            #initializing my energy and magnetization to zero
            E1 = 0
            M1 = 0

            config = 2*np.random.randint(2, size=(N,N))-1#creat a NxN lattice of all -1 and +1 values

            #Monte Carlo part
            for m in range(iterations):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                spin =  config[a, b]
                mag = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                e = 2*spin*mag
                if e <= 0:
                    spin *= -1
                elif rand() < np.exp(-e*(1/(kB*T[l]))):
                    spin *= -1
                config[a, b] = spin

                #Energy and Magnetization Calculation  
                energy=0
                for i in range(N):
                    for j in range(N):
                        spin1 = config[i,j]
                        mag1 = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
                        energy = -J*mag1*spin1 - h[q]*spin1 + energy

                mm = np.sum(config)

                E1 = E1 + energy
                M1 = M1 + mm
                En[l] = E1/(N*N)
                M[l] = M1/(N*N)
                
##            #This is writing to the file that will produce my 3d phase diagram 
##            f.write(str(h[q]))
##            f.write(',')
##            f.write(str(T[l]))
##            f.write(',')
##            f.write(str(mm))
##            f.write('\n')

        #Specific Heat and Magnetic Susceptibility Calculations
        for l in range(len(T)):
            if l==0 or l==1 or l==(len(T)-1) or l==(len(T)-2):
                Cv[l]=0
                Sus[l]=0
            else:
                Cv[l] = (En[l+1]-En[l-1])/(2*4.9/n)
                Sus[l] = abs((M[l+1]-M[l-1]))/(2*4.9/n)
            

##    f.close()

    #Plotting everythin and saving the plots to external files
    f1 = plt.figure()
    plt.plot(T, En, color="red", label=' Energy')
    plt.xlabel("Temperature (T)")
    plt.ylabel("Energy")
    f1.savefig("TvsE.pdf")

    f2 = plt.figure()
    plt.plot(T, abs(M), color="blue", label=' Magnetization')
    plt.xlabel("Temperature (T)")
    plt.ylabel("Magnetization")
    f2.savefig("TvsMag.pdf")

    f3 = plt.figure()
    plt.plot(T, Cv, color="green", label=' Specific Heat')
    plt.xlabel("Temperature (T)")
    plt.ylabel("Specific Heat")
    f3.savefig("TvsCv.pdf")

    f4 = plt.figure()
    plt.plot(T, Sus, color="black", label=' Susceptibility')
    plt.xlabel("Temperature (T)")
    plt.ylabel("Susceptibility")
    f4.savefig("TvsSus.pdf")



main()
