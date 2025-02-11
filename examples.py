import numpy as np
from sympy as sp
import json

with open("MM.json", "r") as f:
    data = json.load(f) 
print(data)

productM = Matrix([[r["products"][i] for i in r["species"] if i in r["products"].keys()]  for r in MM["reactions"]])

network_name_MM='Michaelis-Menten-Network'
species_MM=['S_0','S_1','E_0','E_1']                               
educts_MM=Matrix([[1,0,1,0],[0,0,0,1],[0,0,0,1]]).T                
products_MM=Matrix([[0,0,0,1],[1,0,1,0],[0,1,1,0]]).T              
scaling_species_MM=[1,1,0,0]                                       
scaling_rates_MM=[0,1,1]   

network_name_Hill='Hill_network'                        
species_Hill=['S_0','S_1','E_0','E_1','E_2']                                                                            
educts_Hill=Matrix([[1,0,1,0,0],[0,0,0,1,0],[1,0,0,1,0],[0,0,0,0,1],[0,0,0,0,1]]).T               
products_Hill=Matrix([[0,0,0,1,0],[1,0,1,0,0],[0,0,0,0,1],[1,0,0,1,0],[0,1,0,1,0]]).T               
scaling_species_Hill=[1,1,0,0,0]                                      
scaling_rates_Hill=[1,2,2,1,1] 

network_name_allos='Allosteric Activation'                         
species_allos=['S','A','P','E','EA','EAS']                               
educts_allos=Matrix([[0,1,0,1,0,0],[0,0,0,0,1,0],[1,0,0,0,1,0],[0,0,0,0,0,1],[0,0,0,0,0,1]]).T                  
products_allos=Matrix([[0,0,0,0,1,0],[0,1,0,1,0,0],[0,0,0,0,0,1],[1,0,0,0,1,0],[0,0,1,0,1,0]]).T            
scaling_species_allos=[1,1,1,0,0,0]                                       
scaling_rates_allos=[0,1,0,1,1] 

network_name_gen='gene expression'                         
species_gen=['off','on','R','P']                               
educts_gen=Matrix([[1,0,0,0],[0,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,1,0],[0,0,0,1]]).T          #neutral model        
products_gen=Matrix([[0,1,0,0],[1,0,0,0],[0,1,1,0],[0,0,0,0],[0,0,1,1],[0,0,0,0]]).T        #neutral model     
#educts=Matrix([[1,0,0,0],[0,1,0,1],[0,1,0,0],[0,0,1,0],[0,0,1,0],[0,0,0,1]]).T         #negative model        
#products=Matrix([[0,1,0,0],[1,0,0,1],[0,1,1,0],[0,0,0,0],[0,0,1,1],[0,0,0,0]]).T       #negative model   
#educts=Matrix([[1,0,0,1],[0,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,1,0],[0,0,0,1]]).T         #positive model        
#products=Matrix([[0,1,0,1],[1,0,0,0],[0,1,1,0],[0,0,0,0],[0,0,1,1],[0,0,0,0]]).T       #positive model    
scaling_species_gen=[0,0,0,1]                                       
scaling_rates_gen=[1,1,1,1,1,0] 