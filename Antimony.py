# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:48:59 2019

@author: Paolo
"""

# plotting utility from Tellurium documentation
def my_plot(r, result, sizeX=6, sizeY=4):

    import pylab as p

    if result is None:
        raise Exception("no simulation result")

    # assume result is a standard numpy array

    selections = r.timeCourseSelections

    if len(result.shape) != 2 or result.shape[1] != len(selections):
        raise Exception("simulation result columns not equal to number of selections,"
                        "likely a simulation has not been run")

    times = result[:,0]

    p.figure(figsize=(sizeX, sizeY))
    
    for i in range(1, len(selections)):
        series = result[:,i]
        name = selections[i]
        p.plot(times, series, label=str(name))
        
    p.legend()
 
    p.show()


import tellurium as te

MM = '''
model Michaelis_Menten
  
  // Species initializations:
  R5P = 5e-7
  RibA = 2e-7
  R5P_RibA = 0
  D2B4P = 5e-7
  RibH = 2e-7
  D2B4P_RibH = 0
  D8RL = 5e-7
  RibE = 2e-7
  D8RL_RibE = 0
  Riboflavin = 0
  RibC = 2e-7
  Riboflavin_RibC =0
  FMN = 0
  FMN_RibC = 0
  FAD =0
  
  GTP = 5e-7
  GTP_RibA = 0
  D6P4 = 0
  RibD = 2e-7
  RibD_D6P4 = 0
  A6U = 0
  
  // Reactions:
      
  Synthesis1: -> R5P;    (k_R5P);
  Synthesis2: -> GTP;     (k_GTP);
  Biomass1:  Riboflavin ->;   (k_Riboflavin * Riboflavin);
  Biomass2:  FAD ->;   (k_FAD * FAD);
  
  BindingR1: R5P  + RibA  -> R5P_RibA;      (k1 * RibA * R5P);                      
  DissociationR1: R5P_RibA     -> R5P + RibA;    (k1r * R5P_RibA );
  ConversionR1: R5P_RibA     -> RibA + D2B4P;  (k2 * R5P_RibA);
  
  BindingR2: D2B4P +RibH  -> D2B4P_RibH;    (k1 * RibH * D2B4P);
  DissociationR2: D2B4P_RibH   -> RibH + D2B4P;  (k1r * D2B4P_RibH);
  ConversionR2: D2B4P_RibH   -> RibH + D8RL;   (k2 * D2B4P_RibH);
  
  BindingR3: D8RL+RibE    -> D8RL_RibE;    (k1 * RibE * D8RL);
  DissociationR3: D8RL_RibE    -> RibE + D8RL;  (k1r * D8RL_RibE);
  ConversionR3: D8RL_RibE    -> RibE + Riboflavin; (k2 * D8RL_RibE);
  
  BindingR4: RibC + Riboflavin -> Riboflavin_RibC;   (k1 * RibC * Riboflavin);
  DissociationR4: Riboflavin_RibC   -> RibC + Riboflavin; (k1r * Riboflavin_RibC);
  ConversionR4:  Riboflavin_RibC   -> RibC + FMN;        (k2 * Riboflavin_RibC);
  
  BindingR5: RibC + FMN -> FMN_RibC;                 (k1 * RibC * FMN);
  DissociationR5: FMN_RibC   -> RibC + FMN;               (k1r * FMN_RibC);
  ConversionR5: FMN_RibC   -> RibC + FAD;               (k2 * FMN_RibC);
  
  BindingR6: GTP + RibA   -> GTP_RibA;      (k1 * RibA * GTP);
  DissociationR6: GTP_RibA     -> GTP + RibA;    (k1r * GTP_RibA);
  ConversionR6: GTP_RibA     -> RibA + D6P4;   (k2 * GTP_RibA);
  
  BindingR7: D6P4 + RibD  -> RibD_D6P4;     (k1 * RibD * D6P4);
  DissociationR7: RibD_D6P4    -> D6P4 + RibD;   (k1r * RibD_D6P4);
  ConversionR7: RibD_D6P4    -> RibD + A6U;   (k2 * RibD_D6P4);





   






  // Variable initialization:
  k1  = 1e6; 
  k1r = 1e-4; 
  k2  = 0.1;
  k_Riboflavin = 0;
  k_R5P = 0;
  k_GTP = 0;
  k_FAD = 0;
end'''

r = te.loada(MM)
result = r.simulate(0, 100, 1000)

# plot the simulation output 
#r.plot()
my_plot(r, result, sizeX=12, sizeY=8)

# export to SBML
r.exportToSBML('FluxBalanceAnalysis.xml', current = False)



