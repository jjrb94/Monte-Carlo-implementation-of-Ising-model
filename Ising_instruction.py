# This is an instruction to the whole project, which shows the basic structure of code and guides the following work.
# for the functions below, we need to write arguments, return variables and a brief description for each of them.
# unfinished yet

def initialize_lattice(n): # the lattice will be represnted by a n x n matrix 
''' 
The initialization of lattice should be the first step of simulation. In this function, 
we give the lattice a random initial state by built-in function: np.random.random(), the

'''
    return np.random.choice([1, -1], size=(n, n))


def spin_flip(): # I dont think we need to do this in a seperate function 
'''
Flip a spin, namely change the state of one lattice point from +1 to -1 or from -1 to +1, 
once a time to generate a new state of lattice.
'''
    return


def Monte_Carlo(): # This can also be done in the Ising_simulation function 
'''
If the old state of lattice and a new state of lattice are given, use Metropolis Algorithm to determine whether
we should accept it or not.
'''
    return

def Energy(sigma_k,sum_sigma_i,J):  # Double check if this is the correct equation to find the energy for each site.               
'''
This function will calculate the energy for each random spin site.
The energy equation is written down, where sum_sigma_i represent the sum of all the neighboring spins of spin k as stated on the readMe.
'''
    return -J*sigma_k*sum_sigma_i

def Magnetization(lat):
'''
This fucntion calculate the average magnetization. lat = lattice created in the ising fucntion. 
'''
    return np.sum(lat)/(len(lat))


def Ising_simulation(n, steps, J, T):
    '''
    '''
    
    k_b = 1 # Set the actuall bolztman constant if needed
    
    lattice = lattice_gerator(n)
    energies = []
    
    for step in range(steps):
        # We will use this random generator to obtain the random indexes for the random spin site on the lattice
        i = np.random.randint(n)
        j = np.random.randint(n)
        
        s_k = lattice[i][j] # This is our random chosen spin site 
        s_i_sum = lattice[i+1][j] + lattice[i-1][j] + lattice[i][j+1] + lattice[i][j-1] # This is the sum of the neighborin spins for the specific site
        
        E = Energy(s_k,s_i_sum,J)
        
        delta_E = E-(-E) # The energy is given by the defference between the energy of the spin original configuration 
                         # and the energy if the spin was flip i.e changed in sign. 
        
        if delta_E < 0 or np.random.random() < np.exp(-delta_E/(k_b*T)): # If any of this two conditions is met, then the spin is flipped.
            lattice[i][j] = -lattice[i][j]
            
        energies.append(E) # This line will add the energy values for each spin site to a list which will then use to find the avarge energy
        
    # Advcice if we should use separete function to do the calculation of the evarage_energy and the evarage_energy^2.
    average_energy = np.sum(np.exp(-energies/(k_b*T))*energies)/np.sum(np.exp(-energies/(k_b*T)))
    
    average_energy_2 = np.sum(np.exp(-(energies**2)/(k_b*T))*(energies**2)/np.sum(np.exp(-(energies**2)/(k_b*T))
    
    specific_heat = (average_energy_2 - average_energy**2)/(T**2)
                                                                                 
     M = Magnetization(lattice)    
     
     #We need to add the correlation calculations                                                                             
                              
    return specific_heat, M 

