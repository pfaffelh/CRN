import numpy as np
import sympy
from sympy import *

def crn_lln(network_name,species,educts,products,scaling_species,scaling_rates):
    reaction_number=shape(educts)[1] #number of reactions                            
    species_number=len(species)  #number of species                                   
    reaction_matrix=(products-educts).T  #reaction matrix                           
    rates=[symbols('k%d' %i) for i in range(reaction_number)] #reaction rates      
    N=symbols('N') #scaling factor
    rates_scaled=[rates[i]*N**(max(0,scaling_rates[i]-1)) for i in range(reaction_number)] #scaled reaction rates

    #Calculates the constant linear combinations of species in the network
    nullspace_reaction_matrix=reaction_matrix.nullspace()
    nullspace_reaction_matrix_as_matrix=Matrix([[nullspace_reaction_matrix[n][i] for i in range(species_number)] for n in range(len(nullspace_reaction_matrix))]).T
    index_fast_species = [i for i in range(len(scaling_species)) if scaling_species[i]==0]           
    index_slow_species = [i for i in range(len(scaling_species)) if scaling_species[i]==1]           
    reaction_matrix_slow=Matrix([[nullspace_reaction_matrix[n][i] for i in range(species_number) if i not in index_fast_species]for n in range(len(nullspace_reaction_matrix))]).T 
    nullspace_reaction_matrix_slow=reaction_matrix_slow.nullspace()
    constant_linear_combinations_fast=[Matrix([(nullspace_reaction_matrix_as_matrix*nullspace_reaction_matrix_slow[n])[i] for i in range(species_number) if i in index_fast_species]) for n in range(len(nullspace_reaction_matrix_slow))]
    M=[symbols('M%d' %i) for i in range(len(nullspace_reaction_matrix))]  #constants for the linear combinations            
    v=[symbols('v%d' %i) for i in range(species_number-len(index_fast_species))] #slow species          
    z=[symbols('z%d' %i) for i in range(len(index_fast_species))] #fast species          
    vz= v + z
    i_fast=0
    i_slow=0
    for j in range(len(vz)):
        if scaling_species[j]==0:
            vz[j]=z[i_fast]
            i_fast=i_fast+1
        else:
            vz[j]=v[i_slow]
            i_slow=i_slow+1
    for i in range(len(constant_linear_combinations_fast)):
        print(f'The constant linear combinations involving only fast species are {((constant_linear_combinations_fast[i].T)*(Matrix([[z[n]] for n in range(len(z))])))[0,0]}')
        if i<len(constant_linear_combinations_fast)-1:
            print('and')
    for i in range(len(nullspace_reaction_matrix)-len(constant_linear_combinations_fast)):
        z1=[i/N for i in z]
        v1z=v + z1
        print(f'In addition, the constant linear combinations that include slow species are {((nullspace_reaction_matrix[i].T)*(Matrix([[v1z[n]] for n in range(len(v1z))])))[0,0]}')
        if i<len(nullspace_reaction_matrix)-len(constant_linear_combinations_fast)-1:
            print('and')
    index_irrelevant_fast_species=[]
    stepvariable=len(index_fast_species)
    for i in range(len(constant_linear_combinations_fast)):
        for j in range(stepvariable-1,-1,-1):
            if ((constant_linear_combinations_fast[i].T)*(Matrix([[z[n]] for n in range(len(z))])))[0,0].coeff(z[j])!=0:
                index_irrelevant_fast_species.append(index_fast_species[j])
                stepvariable=j
                break
    index_relevant_fast_species=[i for i in index_fast_species if i not in index_irrelevant_fast_species]
    reduced_reactíon_matrix_fast=reaction_matrix.col([i for i in index_relevant_fast_species])
    #print(f'The reduced fast network is {reduced_reactíon_matrix_fast}')
    index_irrelevant_slow_species=[]
    #Asks user which slow species can be eliminated from the network.
    for i in range(len(nullspace_reaction_matrix)-len(constant_linear_combinations_fast)):
        while True:
            index_user=input(f'What index u want to eliminate for the {i+1} constant linear combination?')
            if index_user.isdigit():
                index_user=int(index_user)
                if index_user in index_slow_species:
                    if ((nullspace_reaction_matrix[i].T)*(Matrix([[vz[n]] for n in range(len(vz))])))[0,0].coeff(vz[index_user]) != 0:
                        break
                    else:
                        print(f'Please enter an index which occurs in the {i+1} constant linear combinaation')
                else:
                    print("The index must occur in the network") 
            else:
                print("Please enter a number")
        index_irrelevant_slow_species.append(index_user)
    index_relevant_slow_species=[i for i in index_slow_species if i not in index_irrelevant_slow_species]
    reduced_reactíon_matrix_slow=reaction_matrix.col([i for i in index_relevant_slow_species])
    #print(f'The reduced slow network is {reduced_reactíon_matrix_slow}')

    fp=[symbols('fp%d' %i) for i in index_relevant_slow_species] #first derivatives of the slow species
    fpp=[symbols('fpp%d%d' %(i,j)) for i in index_relevant_slow_species for j in index_relevant_slow_species] #second derivatives of the slow species
    u=[symbols('u%d' %i) for i in index_relevant_slow_species] #fluctuations slow species
    #Variables for the ansatz
    a=[symbols('a%d' %i) for i in range(len(fp)*len(index_relevant_fast_species))] 
    b=[symbols('b%d' %i) for i in range(len(fp)*len(index_relevant_fast_species))]
    c=[symbols('c%d' %i) for i in range(len(fpp)*len(index_relevant_fast_species))]
    d=[symbols('d%d' %i) for i in range(len(fpp)*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2))]

    #Generates the vector consisting of the species after eliminating.
    vector_Generator=[None]*species_number
    Mk=[None]*len(nullspace_reaction_matrix)
    for i in range(len(Mk)):
        if i<len(constant_linear_combinations_fast):
            Mk[i]=((constant_linear_combinations_fast[i].T)*(Matrix([[z[n]] for n in range(len(z))])))[0,0]
        else:
            Mk[i]=((nullspace_reaction_matrix[i-len(constant_linear_combinations_fast)].T)*(Matrix([[vz[n]] for n in range(len(vz))])))[0,0]
    stepvariable=len(M)-1
    for i in range(species_number):
        if i in index_relevant_slow_species:
            vector_Generator[i]=vz[i]
        elif i in index_relevant_fast_species:
            vector_Generator[i]=vz[i]
        elif i in index_irrelevant_slow_species:
            abcd=(M[stepvariable]-Mk[stepvariable]+Mk[stepvariable].coeff(vz[i])*vz[i])/(Mk[stepvariable].coeff(vz[i]))
            vector_Generator[i]=lambdify(z,abcd)(*([0]*len(z)))
            stepvariable=stepvariable-1
        elif i in index_irrelevant_fast_species:
            vector_Generator[i]=M[stepvariable]-Mk[stepvariable]+Mk[stepvariable].coeff(vz[i])*vz[i]
            stepvariable=stepvariable-1
    #print(f'The optimised species vector is {vector_Generator}')

    z=[None]*len(index_relevant_fast_species)
    for i in range(len(z)):
        j=index_fast_species.index(index_relevant_fast_species[i])
        z[i]=symbols('z%d' %j)
    #print(z)


    #Calculates the generatorpart G0f for a function f only depending on slow species.
    G0f=0
    for n in range(len(fp)):
        for i in range(reaction_number):
            G0f1=1
            for j in range(species_number):
                if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:
                    for k in range(abs(educts.col(i).T[j])):
                        G0f1=G0f1*(vector_Generator[j])
            G0f+=G0f1*reduced_reactíon_matrix_slow[i,n]*fp[n]*rates[i]
    #print(f'G0f = {G0f}')    


    #Calculates the generatorpart G1g for the ansatzfunction g.
    G1g=0
    for n in range(len(fp)):
        for m in range(len(index_relevant_fast_species)):
            for i in range(reaction_number):
                G1f1=1
                for j in range(species_number):
                    if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:  
                        for k in range(abs(educts.col(i).T[j])):
                            G1f1=G1f1*(vector_Generator[j])
                G1g+=G1f1*reduced_reactíon_matrix_fast[i,m]*a[m+n*len(index_relevant_fast_species)]*fp[n]*rates[i]
    #print(f'G1g = {G1g}')

    #Calculates the limit generator of the slow species (LLN) by solving a linear equation system.
    Gf=expand((G0f+G1g))
    coeff=[Eq((Gf.coeff(fp[j])).coeff(z[i]),0) for j in range(len(fp)) for i in range(len(index_relevant_fast_species))]                     
    sol_LLN=solve(coeff,a)
    a_sol=[sol_LLN[a[i]] for i in range(len(a))]
    mu_LLN=simplify(lambdify(a,Gf)(*a_sol))
    mu_LLN_rates=lambdify(rates,mu_LLN)
    mu_LLN_limit=simplify(limit(mu_LLN_rates(*rates_scaled),N,oo)) #takes the limit if the scaling of the rates is >1  
    print(f'The LLN of the {network_name} is \n {mu_LLN_limit}') 
  

def crn_clt(network_name,species,educts,products,scaling_species,scaling_rates):
    #Same part as in crn_LLN. We need this limit and the solutions for the ansatzfunction g for the CLT.
    reaction_number=shape(educts)[1] #number of reactions                            
    species_number=len(species)  #number of species                                   
    reaction_matrix=(products-educts).T  #reaction matrix                           
    rates=[symbols('k%d' %i) for i in range(reaction_number)] #reaction rates      
    N=symbols('N') #scaling factor
    rates_scaled=[rates[i]*N**(max(0,scaling_rates[i]-1)) for i in range(reaction_number)] #scaled reaction rates

    #Calculates the constant linear combinations of species in the network
    nullspace_reaction_matrix=reaction_matrix.nullspace()
    nullspace_reaction_matrix_as_matrix=Matrix([[nullspace_reaction_matrix[n][i] for i in range(species_number)] for n in range(len(nullspace_reaction_matrix))]).T
    index_fast_species = [i for i in range(len(scaling_species)) if scaling_species[i]==0]           
    index_slow_species = [i for i in range(len(scaling_species)) if scaling_species[i]==1]           
    reaction_matrix_slow=Matrix([[nullspace_reaction_matrix[n][i] for i in range(species_number) if i not in index_fast_species]for n in range(len(nullspace_reaction_matrix))]).T 
    nullspace_reaction_matrix_slow=reaction_matrix_slow.nullspace()
    constant_linear_combinations_fast=[Matrix([(nullspace_reaction_matrix_as_matrix*nullspace_reaction_matrix_slow[n])[i] for i in range(species_number) if i in index_fast_species]) for n in range(len(nullspace_reaction_matrix_slow))]
    M=[symbols('M%d' %i) for i in range(len(nullspace_reaction_matrix))]  #constants for the linear combinations            
    v=[symbols('v%d' %i) for i in range(species_number-len(index_fast_species))] #slow species          
    z=[symbols('z%d' %i) for i in range(len(index_fast_species))] #fast species          
    vz= v + z
    i_fast=0
    i_slow=0
    for j in range(len(vz)):
        if scaling_species[j]==0:
            vz[j]=z[i_fast]
            i_fast=i_fast+1
        else:
            vz[j]=v[i_slow]
            i_slow=i_slow+1
    for i in range(len(constant_linear_combinations_fast)):
        print(f'The constant linear combinations involving only fast species are {((constant_linear_combinations_fast[i].T)*(Matrix([[z[n]] for n in range(len(z))])))[0,0]}')
        if i<len(constant_linear_combinations_fast)-1:
            print('and')
    for i in range(len(nullspace_reaction_matrix)-len(constant_linear_combinations_fast)):
        z1=[i/N for i in z]
        v1z=v + z1
        print(f'In addition, the constant linear combinations that include slow species are {((nullspace_reaction_matrix[i].T)*(Matrix([[v1z[n]] for n in range(len(v1z))])))[0,0]}')
        if i<len(nullspace_reaction_matrix)-len(constant_linear_combinations_fast)-1:
            print('and')
    index_irrelevant_fast_species=[]
    stepvariable=len(index_fast_species)
    for i in range(len(constant_linear_combinations_fast)):
        for j in range(stepvariable-1,-1,-1):
            if ((constant_linear_combinations_fast[i].T)*(Matrix([[z[n]] for n in range(len(z))])))[0,0].coeff(z[j])!=0:
                index_irrelevant_fast_species.append(index_fast_species[j])
                stepvariable=j
                break
    index_relevant_fast_species=[i for i in index_fast_species if i not in index_irrelevant_fast_species]
    reduced_reaction_matrix_fast=reaction_matrix.col([i for i in index_relevant_fast_species])
    #print(f'The reduced fast network is {reduced_reactíon_matrix_fast}')
    index_irrelevant_slow_species=[]
    #Asks user which slow species can be eliminated from the network.
    for i in range(len(nullspace_reaction_matrix)-len(constant_linear_combinations_fast)):
        while True:
            index_user=input(f'What index u want to eliminate for the {i+1} constant linear combination?')
            if index_user.isdigit():
                index_user=int(index_user)
                if index_user in index_slow_species:
                    if ((nullspace_reaction_matrix[i].T)*(Matrix([[vz[n]] for n in range(len(vz))])))[0,0].coeff(vz[index_user]) != 0:
                        break
                    else:
                        print(f'Please enter an index which occurs in the {i+1} constant linear combinaation')
                else:
                    print("The index must occur in the network") 
            else:
                print("Please enter a number")
        index_irrelevant_slow_species.append(index_user)
    index_relevant_slow_species=[i for i in index_slow_species if i not in index_irrelevant_slow_species]
    reduced_reaction_matrix_slow=reaction_matrix.col([i for i in index_relevant_slow_species])
    #print(f'The reduced slow network is {reduced_reactíon_matrix_slow}')

    fp=[symbols('fp%d' %i) for i in index_relevant_slow_species] #first derivatives of the slow species
    fpp=[symbols('fpp%d%d' %(i,j)) for i in index_relevant_slow_species for j in index_relevant_slow_species] #second derivatives of the slow species
    u=[symbols('u%d' %i) for i in index_relevant_slow_species] #fluctuations slow species
    #Variables for the ansatz
    a=[symbols('a%d' %i) for i in range(len(fp)*len(index_relevant_fast_species))] 
    b=[symbols('b%d' %i) for i in range(len(fp)*len(index_relevant_fast_species))]
    c=[symbols('c%d' %i) for i in range(len(fpp)*len(index_relevant_fast_species))]
    d=[symbols('d%d' %i) for i in range(len(fpp)*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2))]

    #Generates the vector consisting of the species after eliminating.
    vector_Generator=[None]*species_number
    Mk=[None]*len(nullspace_reaction_matrix)
    for i in range(len(Mk)):
        if i<len(constant_linear_combinations_fast):
            Mk[i]=((constant_linear_combinations_fast[i].T)*(Matrix([[z[n]] for n in range(len(z))])))[0,0]
        else:
            Mk[i]=((nullspace_reaction_matrix[i-len(constant_linear_combinations_fast)].T)*(Matrix([[vz[n]] for n in range(len(vz))])))[0,0]
    stepvariable=len(M)-1
    for i in range(species_number):
        if i in index_relevant_slow_species:
            vector_Generator[i]=vz[i]
        elif i in index_relevant_fast_species:
            vector_Generator[i]=vz[i]
        elif i in index_irrelevant_slow_species:
            abcd=(M[stepvariable]-Mk[stepvariable]+Mk[stepvariable].coeff(vz[i])*vz[i])/(Mk[stepvariable].coeff(vz[i]))
            vector_Generator[i]=lambdify(z,abcd)(*([0]*len(z)))
            stepvariable=stepvariable-1
        elif i in index_irrelevant_fast_species:
            vector_Generator[i]=M[stepvariable]-Mk[stepvariable]+Mk[stepvariable].coeff(vz[i])*vz[i]
            stepvariable=stepvariable-1
    #print(f'The optimised species vector is {vector_Generator}')

    z=[None]*len(index_relevant_fast_species)
    for i in range(len(z)):
        j=index_fast_species.index(index_relevant_fast_species[i])
        z[i]=symbols('z%d' %j)
    #print(z)


    #Calculates the generatorpart G0f for a function f only depending on slow species.
    G0f=0
    for n in range(len(fp)):
        for i in range(reaction_number):
            G0f1=1
            for j in range(species_number):
                if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:
                    for k in range(abs(educts.col(i).T[j])):
                        G0f1=G0f1*(vector_Generator[j])
            G0f+=G0f1*reduced_reaction_matrix_slow[i,n]*fp[n]*rates[i]
    #print(f'G0f = {G0f}')    


    #Calculates the generatorpart G1g for the ansatzfunction g.
    G1g=0
    for n in range(len(fp)):
        for m in range(len(index_relevant_fast_species)):
            for i in range(reaction_number):
                G1f1=1
                for j in range(species_number):
                    if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:  
                        for k in range(abs(educts.col(i).T[j])):
                            G1f1=G1f1*(vector_Generator[j])
                G1g+=G1f1*reduced_reaction_matrix_fast[i,m]*a[m+n*len(index_relevant_fast_species)]*fp[n]*rates[i]
    #print(f'G1g = {G1g}')

    #Calculates the limit generator of the slow species (LLN) by solving a linear equation system.
    Gf=expand((G0f+G1g))
    coeff=[Eq((Gf.coeff(fp[j])).coeff(z[i]),0) for j in range(len(fp)) for i in range(len(index_relevant_fast_species))]                     
    sol_LLN=solve(coeff,a)
    a_sol=[sol_LLN[a[i]] for i in range(len(a))]
    mu_LLN=simplify(lambdify(a,Gf)(*a_sol))
    mu_LLN_rates=lambdify(rates,mu_LLN)
    mu_LLN_limit=simplify(limit(mu_LLN_rates(*rates_scaled),N,oo)) #takes the limit if the scaling of the rates is >1  
    #print(f'The LLN of the {network_name} is \n {mu_LLN_limit}')  
    
    #Calculates the generatorpart L0f for a function f only depending on fluctuations of the slow species
    L0f=0
    for n in range(len(fpp)):
        for i in range(reaction_number):
            L0f1=1
            for j in range(species_number):
                if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:
                    for k in range(abs(educts.col(i).T[j])):
                        L0f1=L0f1*(vector_Generator[j])
            if n%len(index_relevant_slow_species)==n//len(index_relevant_slow_species):
                L0f+=0.5*L0f1*reduced_reaction_matrix_slow[i,n//len(index_relevant_slow_species)]**2*fpp[n]*rates[i]
            else:
                L0f+=0.5*L0f1*reduced_reaction_matrix_slow[i,n//len(index_relevant_slow_species)]*reduced_reaction_matrix_slow[i,n%len(index_relevant_slow_species)]*fpp[n]*rates[i]
        if n%len(index_relevant_slow_species)==n//len(index_relevant_slow_species):
            for l in range(1,max(reduced_reaction_matrix_slow)+1,1):
                L0f+=G0f.coeff(v[n//len(index_relevant_slow_species)]**l)*u[n//len(index_relevant_slow_species)]*l*v[n//len(index_relevant_slow_species)]**(l-1)
    #print(f'L0f = {L0f}')  


    #Calculates the generatorpart L1g for the ansatzfunction g. 
    L1g=0
    for n in range(len(fpp)):
        for m in range(len(index_relevant_fast_species)):
            for i in range(reaction_number):
                L1g1=1
                for j in range(species_number):
                    if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:  
                        for k in range(abs(educts.col(i).T[j])):
                            L1g1=L1g1*(vector_Generator[j])
                if reduced_reaction_matrix_slow[i,n//len(index_relevant_slow_species)]!=0:
                    L1g+=L1g1*reduced_reaction_matrix_slow[i,n//len(index_relevant_slow_species)]*(z[m]+reduced_reaction_matrix_fast[i,m])*a[m+n//len(index_relevant_slow_species)*len(index_relevant_fast_species)]*fpp[n]*rates[i]
            L1g-=mu_LLN_limit.coeff(fp[n//len(index_relevant_slow_species)])*fpp[n]*a[m+n//len(index_relevant_slow_species)*len(index_relevant_fast_species)]*z[m]
        if n%len(index_relevant_slow_species)==n//len(index_relevant_slow_species):
            for l in range(1,max(reduced_reaction_matrix_slow)+1,1):
                L1g+=G1g.coeff(v[n//len(index_relevant_slow_species)]**l)*u[n//len(index_relevant_slow_species)]*l*v[n//len(index_relevant_slow_species)]**(l-1)
    #print(f'L1g = {L1g}')


    #Calculates the generatorpart L2h for the ansatzfunction h.
    L2h=0
    for n in range(len(fpp)):
        for m in range(len(index_relevant_fast_species)):
            for i in range(reaction_number):
                L2h1=1
                for j in range(species_number):
                    if (educts.col(i).T*diag(vector_Generator, unpack=True))[j]!=0:  
                        for k in range(abs(educts.col(i).T[j])):
                            L2h1=L2h1*(vector_Generator[j])
                if n%len(index_relevant_slow_species)==n//len(index_relevant_slow_species):
                    L2h+=L2h1*reduced_reaction_matrix_fast[i,m]*b[m+n//len(index_relevant_slow_species)*len(index_relevant_fast_species)]*fp[n//len(index_relevant_slow_species)]*rates[i]
                L2h+=L2h1*reduced_reaction_matrix_fast[i,m]*c[m+n*len(index_relevant_fast_species)]*fpp[n]*rates[i]
                for k in range(m,len(index_relevant_fast_species),1):
                    if(reduced_reaction_matrix_fast[i,m]==0):
                        if(reduced_reaction_matrix_fast[i,k]==1):
                            L2h+=L2h1*reduced_reaction_matrix_fast[i,k]*d[len(index_relevant_fast_species)-1+k+n*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2)]*(z[m]+z[k])*fpp[n]*rates[i]
                            break
                        if(reduced_reaction_matrix_fast[i,k]==-1):
                            L2h+=L2h1*reduced_reaction_matrix_fast[i,k]*d[len(index_relevant_fast_species)-1+k+n*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2)]*(z[m]+z[k]-1)*fpp[n]*rates[i]
                            break
                    if(reduced_reaction_matrix_fast[i,m]==1):
                        if(reduced_reaction_matrix_fast[i,k]==1):
                            L2h+=L2h1*reduced_reaction_matrix_fast[i,m]*d[m+n*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2)]*(z[m])*fpp[n]*rates[i]
                        if(reduced_reaction_matrix_fast[i,k]==0):
                            L2h+=L2h1*reduced_reaction_matrix_fast[i,m]*d[len(index_relevant_fast_species)-1+k+n*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2)]*(z[m]+z[k])*fpp[n]*rates[i]
                    if(reduced_reaction_matrix_fast[i,m]==-1):
                        if(reduced_reaction_matrix_fast[i,k]==-1):
                            L2h+=L2h1*reduced_reaction_matrix_fast[i,m]*d[m+n*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2)]*(z[m]-1)*fpp[n]*rates[i]
                        if(reduced_reaction_matrix_fast[i,k]==0):
                            L2h+=L2h1*reduced_reaction_matrix_fast[i,m]*d[len(index_relevant_fast_species)-1+k+n*(len(index_relevant_fast_species)*(len(index_relevant_fast_species)+1)//2)]*(z[m]+z[k]-1)*fpp[n]*rates[i]
    #print(f'L2h = {L2h}')


    #Calculates the limit generator of the fluctuations of the slow species (CLT) by solving a linear equation system.
    Lf=expand(lambdify(a,L0f+L1g+L2h)(*a_sol))
    d=[d[n] for n in range(len(d)) if Lf.coeff(d[n])!=0]
    coeff_fp=[Eq(simplify((Lf.coeff(fp[j])).coeff(z[i])),0) for j in range(len(fp)) for i in range(len(index_relevant_fast_species))]                       
    coeff_fpp=[Eq(simplify((Lf.coeff(fpp[j])).coeff(z[i])),0) for j in range(len(fpp)) for i in range(len(index_relevant_fast_species))]
    coeff_fpp_2=[Eq(simplify((Lf.coeff(fpp[j])).coeff(z[i]*z[k])),0) for j in range(len(fpp)) for i in range(len(index_relevant_fast_species)) for k in range(len(index_relevant_fast_species)) if k>=i]
    sol_CLT=solve(coeff_fp+coeff_fpp+coeff_fpp_2,b+c+d)
    b_sol=[sol_CLT[b[i]] for i in range(len(b))]
    c_sol=[sol_CLT[c[i]] for i in range(len(c))]
    d_sol=[sol_CLT[d[i]] for i in range(len(d))]
    mu_CLT=simplify(lambdify(b+c+d,Lf)(*(b_sol+c_sol+d_sol)))
    mu_CLT_rates=lambdify(rates,mu_CLT)
    mu_CLT_limit=simplify(limit(mu_CLT_rates(*rates_scaled),N,oo))  #takes the limit if the scaling of the rates is >1 
    mu_CLT_limit=mu_CLT_limit.expand()
    #returns the drift part of the limit generator for the slow fluctuation of every species
    for n in range(len(fp)):
        print(f'The drift-proportion of {fp[n]} is {simplify(mu_CLT_limit.coeff(fp[n]))}.')
    #returns the diffusion part of the limit generator for the slow fluctuation of every species 
    for n in range(len(fpp)):
        if n%len(index_relevant_slow_species)==n//len(index_relevant_slow_species):
            print(f'The sigma-proportion of {fpp[n]} is {simplify(mu_CLT_limit.coeff(fpp[n]))}.')
        else:
            if n//len(index_relevant_slow_species)<n%len(index_relevant_slow_species):
                print(f'The sigma-proportion of {fpp[n]} is {simplify(mu_CLT_limit.coeff(fpp[n])+mu_CLT_limit.coeff(fpp[n*n%len(index_relevant_slow_species)+n//len(index_relevant_slow_species)]))}.')