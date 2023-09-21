"""
@author: Karl Seltzer
contributors: Karl Seltzer, Havala Pye: USEPA
"""

import rdkit
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from datetime import date,datetime

startTime = datetime.now()

def dfappend_cb6r3_ae7_tracer(dfin):
  '''
    Compound to mechanism mapping
  '''

  mech = 'CB6R3_AE7_TRACER'
  dfin['CB6R3_AE7_TRACER']=''  

  # mech4import prep
  column_names = ['mechanism','SPECIES_ID','mech_species','moles_ratio']
  dfmech4import = pd.DataFrame(columns=column_names)
  
  # Perform mapping, append info to data frame, create mech4import df
  for idx, row in dfin.iterrows():
    Pvap        = row['VP_Pascal_OPERA']
    molwght     = row['SPEC_MW']
    smiles      = row['Smiles Notation']
    koh         = row['ATMOSPHERIC_HYDROXYLATION_RATE_(AOH)_CM3/MOLECULE*SEC_OPERA']
    spec_id     = row['SPECIES_ID']
    log10cstar  = np.log10(molwght * Pvap / 8.31451 / 298.15 * 1000000)
    if ( pd.isnull(smiles) or smiles=='-' or pd.isnull(koh) ):
                             mechspecies = 'NONBAF' # Unreactive
    elif ( spec_id==279 ): mechspecies = 'ALD2_PRIMARY'
    elif ( spec_id==465 ): mechspecies = 'FORM_PRIMARY'
    else: mechspecies = get_cb6r3_ae7_tracer_roc(smiles,log10cstar)
    
    mech4import = pd.Series(data={'mechanism':mech,'SPECIES_ID':row['SPECIES_ID'],
                                  'mech_species':mechspecies,'moles_ratio':1})
    dfmech4import = dfmech4import.append(mech4import,ignore_index=True)

  # write mech4import df to file
  today = date.today()
  dfmech4import.to_csv('./mechanism_forImport_'+mech+'_speciate5_2_'+str(today)+'.csv',index=False,header=False)
  
  return 

def get_cb6r3_ae7_tracer_roc(smiles,log10cstar):
  '''
  Function maps input reactive organic carbon (ROC) species to CB6R3_AE7_TRACER species.
  Uses functional group and molecule info from RDKit http://www.rdkit.org/
  Function inputs, for ONE compound (make loop outside this function):
      smiles string (should be canonical for explicit species, some alt added) 
      log10(Cstar in micrograms/m3) 
  '''
  
  # Prep inputs
  m       = Chem.MolFromSmiles(smiles)
  smiles  = smiles.upper()
  ndbl    = smiles.count('=')
  ntrip   = smiles.count('#')
  nC      = smiles.count('C')-smiles.count('CL')
  nO      = smiles.count('O')
  nN      = smiles.count('N')
  nSi     = smiles.count('SI')
  nS      = smiles.count('S')
  nH      = 0
  for atom in m.GetAtoms():
      nH += atom.GetTotalNumHs()

  # Count functional groups (http://rdkit.org/docs/source/rdkit.Chem.Fragments.html)
  nacid     = rdkit.Chem.Fragments.fr_COO(m,countUnique=True)     # carboxylic acid
  nketone   = rdkit.Chem.Fragments.fr_ketone(m,countUnique=True)
  naldehyde = rdkit.Chem.Fragments.fr_aldehyde(m,countUnique=True)
  ncarbonyl = nketone + naldehyde
  nbenzene  = rdkit.Chem.Fragments.fr_benzene(m,countUnique=True)
  nalcohol  = rdkit.Chem.Fragments.fr_Al_OH(m,countUnique=True) + \
              rdkit.Chem.Fragments.fr_Ar_OH(m,countUnique=True)      # aliphatic and aromatic
  nfuran    = rdkit.Chem.Fragments.fr_furan(m,countUnique=True) # number of furan rings
  for atom in range(len(m.GetAtoms())):
      if m.GetAtomWithIdx(atom).IsInRing(): 
          gotring = 1 # 0 = no ring, 1 = ring
          break
      else: gotring = 0
  
  if   ( nC <= 6 ):                 mechspecies = 'NONBAF'
  elif ( nO>0 or nN>0 or nSi>0 or nS>0 or ndbl>0 or ntrip>0): 
                                    mechspecies = 'NONBAF'
  elif ( nbenzene>1 ):              mechspecies = 'NONBAF'
  elif ( nC>6 and gotring==1):      mechspecies = 'SOAALK'
  elif ( nC>8):                     mechspecies = 'SOAALK'
  else: mechspecies = 'NONBAF'
    
  return mechspecies
  # end of function

############################################################################################
#### Get input dataset.
df = pd.read_excel('./export_species_properties.xlsx')
############################################################################################

#df = df.iloc[0:50,:]
dfappend_cb6r3_ae7_tracer(df)
############################################################################################

print("Time to generate mechanism for import file: ",datetime.now() - startTime)