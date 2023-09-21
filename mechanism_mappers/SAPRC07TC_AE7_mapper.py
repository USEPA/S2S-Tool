"""
@author: Karl Seltzer
contributors: Karl Seltzer, Havala Pye: USEPA
"""

import rdkit
import itertools
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from datetime import date,datetime

startTime = datetime.now()

def dfappend_saprc07tc_ae7(dfin):
  '''
    Compound to mechanism mapping
  '''

  mech = 'SAPRC07TC_AE7'
  dfin['SAPRC07TC_AE7']=''  

  # mech4import prep
  column_names = ['mechanism','SPECIES_ID','mech_species','moles_ratio']
  dfmech4import = pd.DataFrame(columns=column_names)
  
  # Perform mapping, append info to data frame, create mech4import df
  for idx, row in dfin.iterrows():
    Pvap        = row['VP_Pascal_OPERA']
    molwght     = row['SPEC_MW']
    smiles      = row['Smiles Notation']
    koh         = row['ATMOSPHERIC_HYDROXYLATION_RATE_(AOH)_CM3/MOLECULE*SEC_OPERA']
    log10cstar  = np.log10(molwght * Pvap / 8.31451 / 298.15 * 1000000)
    if ( pd.isnull(smiles) or smiles=='-' or pd.isnull(koh) ):
                                    mechspecies = 'UNR' # Unreactive
    else: mechspecies = get_saprc07tc_ae7_roc(smiles,log10cstar,koh)

    mech4import = pd.Series(data={'mechanism':mech,'SPECIES_ID':row['SPECIES_ID'],
                                  'mech_species':mechspecies,'moles_ratio':1})
    dfmech4import = dfmech4import.append(mech4import,ignore_index=True)

  # write mech4import df to file
  today = date.today()
  dfmech4import.to_csv('./mechanism_forImport_'+mech+'_speciate5_2_'+str(today)+'.csv',index=False,header=False)
  
  return 

def get_saprc07tc_ae7_roc(smiles,log10cstar,koh):
  '''
  Function maps input reactive organic carbon (ROC) species to SAPRC07TC_AE7 species.
  Uses functional group and molecule info from RDKit http://www.rdkit.org/
  Function inputs, for ONE compound (make loop outside this function):
      smiles string (should be canonical for explicit species, some alt added) 
      log10(Cstar in micrograms/m3) 
  C* may not be used if the compound can be mapped without it. 
  '''

  # Prep inputs
  m       = Chem.MolFromSmiles(smiles)
  smiles  = smiles.upper()

  # Count C=C and atoms
  nCdblC  = smiles.count('=C')-smiles.count('O=C')
  nCtripC = smiles.count('#C')
  nC      = smiles.count('C')-smiles.count('CL')
  nO      = smiles.count('O')
  nN      = smiles.count('N')
  nSi     = smiles.count('SI')
  nHalo   = smiles.count('CL') + smiles.count('BR') + smiles.count('F') - smiles.count('FE') + \
            smiles.count('I') - smiles.count('SI') - smiles.count('NI') - smiles.count('TI') - \
            smiles.count('LI')
  nBranch = smiles.count('(C') + smiles.count('CC1') + smiles.count('C1C') + \
            smiles.count('(=C') + smiles.count('(/C')
  ntermalke = smiles[:2].count('C=') + smiles[:2].count('C#') + smiles[-2:].count('=C') + \
              smiles[-2:].count('#C') + smiles.count('(=C)') + smiles[:2].count('C\\')
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
  npyrrole  = rdkit.Chem.Fragments.fr_Nhpyrrole(m,countUnique=True) # number of pyrrole rings
  naroN     = rdkit.Chem.Fragments.fr_Ar_N(m,countUnique=True) # number of aromatic nitrogens
  for atom in range(len(m.GetAtoms())):
      if m.GetAtomWithIdx(atom).IsInRing(): 
          gotring = 1 # 0 = no ring, 1 = ring
          break
      else: gotring = 0
  nnitrate = smiles.count('O[N+](=O)[O-]') + smiles.count('O[N+]([O-])=O')
  nperoxide = smiles.count('COO') +  smiles.count('OOC') - smiles.count('COOC')
  tfmonoterpene = (nC == 10 and nH == 18 and nO == 1) or (nC == 10 and nH == 16) 

  # Mapper is for ROC only and not elemental carbon
  if   ( nC <= 0 ):                 mechspecies = 'NROG'  # Nonreactive
  elif ( smiles == '[C]' ):         mechspecies = 'NROG'  # Nonreactive

  # Explicit species
  elif ( smiles == 'CC(O)=O' ):     mechspecies = 'AACD'  # Acetic Acid
  elif ( smiles == 'CC(=O)C' or smiles == 'CC(C)=O' ): 
                                    mechspecies = 'ACET'  # Acetone
  elif ( smiles == 'C=CC=O' ):      mechspecies = 'ACRO'  # Acrolein
  elif ( smiles == 'C#C' ):         mechspecies = 'ACYE'  # Acetylene; Ethyne
  elif ( smiles == 'CC1=CC(C)=C(C)C=C1' ):
                                    mechspecies = 'B124'  # 1,2,4-Trimethylbenzene
  elif ( smiles == 'CC(=O)C(C)=O' ): 
                                    mechspecies = 'BACL'  # Biacetyl
  elif ( smiles == 'C=CC=C' ):      mechspecies = 'BDE13' # 1,3-Butadiene
  elif ( nC==6 and nH==6 and nO==0 and nbenzene==1 ):   
                                    mechspecies = 'BENZ'  # Benzene   
  elif ( smiles == 'CC=O' ):        mechspecies = 'CCHO'  # Acetaldehyde
  elif ( smiles == "C"  ):          mechspecies = 'CH4'   # Methane
  elif ( smiles == 'C=C' ):         mechspecies = 'ETHE'  # Ethene
  elif ( smiles == 'CCO' ):         mechspecies = 'ETOH'  # Ethanol
  elif ( smiles == 'C(=O)O' or smiles == 'OC=O' ):       
                                    mechspecies = 'FACD'  # Formic Acid
  elif ( smiles == 'C=O' ):         mechspecies = 'HCHO'  # Formaldehyde  
  elif ( smiles == 'O=CC=O' ):      mechspecies = 'GLY'   # Glyoxal
  elif ( smiles == 'CC(=C)C=C' ):   mechspecies = 'ISOP'  # Isoprene
  elif ( smiles == 'CC(=C)C=O' ):   mechspecies = 'MACR'  # Methacrolein
  elif ( smiles == 'CO' ):          mechspecies = 'MEOH'  # Methanol 
  elif ( smiles == 'CC(=O)C=O' ):   mechspecies = 'MGLY'  # Methylglyoxal
  elif ( smiles == 'CC(=O)C=C' ):   mechspecies = 'MVK'   # Methyl Vinyl Ketone
  elif ( nC==10 and nH==8 and nO==0 and nbenzene==2 ):   
                                    mechspecies = 'NAPH'  # Naphthalene   
  elif ( smiles == 'CCC' ):         mechspecies = 'PRPE'  # Propane 
  elif ( tfmonoterpene and nCdblC==1 ): 
                                    mechspecies = 'APIN'  # a-pinene
  elif ( smiles == 'CC1=CC=CC=C1' ):
                                    mechspecies = 'TOLU'  # Toluene
  elif ( smiles == 'CC1=CC(C)=CC=C1' ):
                                    mechspecies = 'MXYL'  # m-Xylene
  elif ( smiles == 'CC1=C(C)C=CC=C1' ):
                                    mechspecies = 'OXYL'  # o-Xylene
  elif ( smiles == 'CC1=CC=C(C)C=C1' ):
                                    mechspecies = 'PXYL'  # p-Xylene
  # Lumped species
  # Volatility defined
  elif ( log10cstar < 2.5 ):        mechspecies = 'NVOL'  # Nonvolatile
  elif ( log10cstar < 6.5 ):        mechspecies = 'IVOC'  # IVOCs
  elif ( tfmonoterpene ):           mechspecies = 'TERP'  # Monoterpenes
  elif ( nC == 15 and nH == 24 ) and ( nCdblC>=1 ):
                                    mechspecies = 'SESQ'  # Sesquiterpenes
  elif ( nC==4 and nH==6 and nO==1 and nketone==1 ): 
                                    mechspecies = 'IPRD'  # Isoprene products
  elif ( nC>=2 and nacid>=1 ):
                                    mechspecies = 'PACD'  # Higher acids
  # Single-ring aromatics
  elif ( nbenzene==1 ):
       if ( nC>=6 and nC<=7 and nalcohol==1 ):    
                                    mechspecies = 'CRES'  # Phenols and Cresols
       elif ( naldehyde>=1 ):       mechspecies = 'BALD'  # Aromatic Aldehydes
       else:                                              # ARO Series
           air_density = 1013. / 298. / 1.380658E-19
           molec_percm3_2_ppm = (1. / 1E6) * air_density
           koh_temp = koh * molec_percm3_2_ppm * 60       # Units: ppm-1 min-1
           if ( koh_temp<=2.0e4 ):  mechspecies = 'ARO1'
           else:                    mechspecies = 'ARO2'
  # Ketones and other non-aldehyde oxygenated products
  elif ( nketone>=1 and naldehyde==0 ):
       if ( koh>=5.0e-12 ):         mechspecies = 'PRD2'
       elif ( koh<=5.0e-12 and koh>=5.0e-13 ):
                                    mechspecies = 'MEK'
       else:                        mechspecies = 'NROG'
  elif ( nC>=2 and naldehyde>=1 ):  mechspecies = 'RCHO'  # Propionaldehyde and higher aldehydes
  elif ( nnitrate>=1 ):             mechspecies = 'RNO3'  # Lumped Organic Nitrates
  # Olefins
  elif ( nCdblC>=1 or nCtripC>=1 ):
       air_density = 1013. / 298. / 1.380658E-19
       molec_percm3_2_ppm = (1. / 1E6) * air_density
       koh_temp = koh * molec_percm3_2_ppm * 60       # Units: ppm-1 min-1
       if ( koh_temp<=7.0e4 ):      mechspecies = 'OLE1'
       else:                        mechspecies = 'OLE2'
  # ALK Series
  else:
       air_density = 1013. / 298. / 1.380658E-19
       molec_percm3_2_ppm = (1. / 1E6) * air_density
       koh_temp = koh * molec_percm3_2_ppm * 60       # Units: ppm-1 min-1
       if ( koh_temp<=2.0e2 ):      mechspecies = 'NROG'                              
       elif ( koh_temp<=5.0e2 and koh_temp>2.0e2 ):
                                    mechspecies = 'ALK1'
       elif ( koh_temp<=2.5e3 and koh_temp>5.0e2 ):
                                    mechspecies = 'ALK2'
       elif ( koh_temp<=5.0e3 and koh_temp>2.5e3 ):
                                    mechspecies = 'ALK3'
       elif ( koh_temp<=1.0e4 and koh_temp>5.0e3 ):
                                    mechspecies = 'ALK4'
       elif ( koh_temp>1.0e4 ):     mechspecies = 'ALK5'
       else:                        mechspecies = 'NROG' 

  return mechspecies
  # end of function

############################################################################################
#### Get input dataset.
df = pd.read_excel('./export_species_properties.xlsx')
############################################################################################

dfappend_saprc07tc_ae7(df)
############################################################################################

print("Time to generate mechanism for import file: ",datetime.now() - startTime)