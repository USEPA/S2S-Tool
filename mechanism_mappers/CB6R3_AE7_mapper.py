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

def dfappend_cb6r3_ae7(dfin):
  '''
    Compound to mechanism mapping
  '''

  mech = 'CB6R3_AE7'
  dfin['CB6R3_AE7']=''  

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
                                    mechspecies, mole_ratio = 'UNR', 1 # Unreactive
    else: mechspecies, mole_ratio = get_cb6r3_ae7_roc(smiles,log10cstar,koh)
    if np.count_nonzero(mechspecies) > 1:
        for i in range(len(mechspecies)):
            mech4import = pd.Series(data={'mechanism':mech,'SPECIES_ID':row['SPECIES_ID'],
                                          'mech_species':mechspecies[i],'moles_ratio':mole_ratio[i]})
            dfmech4import = dfmech4import.append(mech4import,ignore_index=True)
    else:
        mech4import = pd.Series(data={'mechanism':mech,'SPECIES_ID':row['SPECIES_ID'],
                                      'mech_species':mechspecies,'moles_ratio':mole_ratio})
        dfmech4import = dfmech4import.append(mech4import,ignore_index=True)

  # write mech4import df to file
  today = date.today()
  dfmech4import.to_csv('./mechanism_forImport_'+mech+'_speciate5_2_'+str(today)+'.csv',index=False,header=False)
  
  return 

def get_cb6r3_ae7_roc(smiles,log10cstar,koh):
  '''
  Function maps input reactive organic carbon (ROC) species to CB6R3_AE7 species.
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
  if   ( nC <= 0 ):                 mechspecies, mole_ratio = 'UNR', 1   # Unreactive
  elif ( smiles == '[C]' ):         mechspecies, mole_ratio = 'UNR', 1   # Unreactive

  # Explicit species
  elif ( smiles == 'CC(O)=O' ):     mechspecies, mole_ratio = 'AACD', 1  # Acetic Acid
  elif ( smiles == 'CC(=O)C' or smiles == 'CC(C)=O' ): 
                                    mechspecies, mole_ratio = 'ACET', 1  # Acetone
  elif ( smiles == 'CC=O' ):        mechspecies, mole_ratio = 'ALD2', 1  # Acetaldehyde
  elif ( nC==6 and nH==6 and nO==0 and nbenzene==1 ):   
                                    mechspecies, mole_ratio = 'BENZ', 1  # Benzene   
  elif ( smiles == "C"  ):          mechspecies, mole_ratio = 'CH4',  1  # Methane
  elif ( smiles == 'C=C' ):         mechspecies, mole_ratio = 'ETH',  1  # Ethene
  elif ( smiles == "CC" ):          mechspecies, mole_ratio = 'ETHA', 1  # Ethane  
  elif ( smiles == 'C#C' ):         mechspecies, mole_ratio = 'ETHY', 1  # Ethyne
  elif ( smiles == 'CCO' ):         mechspecies, mole_ratio = 'ETOH', 1  # Ethanol
  elif ( smiles == 'C(=O)O' or smiles == 'OC=O' ):       
                                    mechspecies, mole_ratio = 'FACD', 1  # Formic Acid
  elif ( smiles == 'C=O' ):         mechspecies, mole_ratio = 'FORM', 1  # Formaldehyde  
  elif ( smiles == 'O=CC=O' ):      mechspecies, mole_ratio = 'GLY',  1  # Glyoxal
  elif ( smiles == 'OCC=O' ):       mechspecies, mole_ratio = 'GLYD', 1  # Glycolaldehyde
  elif ( smiles == 'CC(=C)C=C' ):   mechspecies, mole_ratio = 'ISOP', 1  # Isoprene
  elif ( smiles == 'CO' ):          mechspecies, mole_ratio = 'MEOH', 1  # Methanol 
  elif ( smiles == 'CC(=O)C=O' ):   mechspecies, mole_ratio = 'MGLY', 1  # Methylglyoxal
  elif ( nC==10 and nH==8 and nO==0 and nbenzene==2 ):   
                                    mechspecies, mole_ratio = 'NAPH', 1  # Naphthalene   
  elif ( smiles == 'CCC' ):         mechspecies, mole_ratio = 'PRPA', 1  # Propane 
  elif ( tfmonoterpene and nCdblC==1 ): 
                                    mechspecies, mole_ratio = 'APIN', 1  # a-pinene

  # Lumped species
  # Volatility defined
  elif ( log10cstar < 2.5 ):        mechspecies, mole_ratio = 'NVOL', 1   # Nonvolatile
  elif ( log10cstar < 6.5 ):        mechspecies, mole_ratio = 'IVOC', 1   # IVOCs
  # Low reactivity --> Unreactive
  elif ( koh<=1.1e-12 ):            mechspecies, mole_ratio = 'UNR', 1
  # Monoterpenes
  elif ( tfmonoterpene ): mechspecies, mole_ratio = 'TERP', 1
  # Isoprene products; methacrolein, methyl vinyl ketone
  elif ( nC==4 and nH==6 and nO==1 and nketone==1 ): 
                                    mechspecies, mole_ratio = 'ISPD', 1
  # Cyclodienes
  elif ( gotring==1 and nCdblC==2 and nbenzene==0):
      mechspecies, mole_ratio = 'IOLE', 1
      carbon_count = nC - 4
      if carbon_count > 0:
          if ( nCdblC>2 and ntermalke>0 ):
              nCdblC = nCdblC - 1
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'OLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
              else:
                  mechspecies = mechspecies,'OLE',
                  mole_ratio  = mole_ratio,ntermalke
              carbon_count = carbon_count - 2 * ntermalke
          if ( nCdblC>2 and carbon_count>=4 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],1
              else:
                  mechspecies = mechspecies,'IOLE',
                  mole_ratio  = mole_ratio,1
              carbon_count = carbon_count - 4
          if ( nketone>0 and carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'KET',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
              else:
                  mechspecies = mechspecies,'KET',
                  mole_ratio  = mole_ratio,nketone
              carbon_count = carbon_count - nketone
          if ( carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'PAR',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
              else:
                  mechspecies = mechspecies,'PAR',
                  mole_ratio  = mole_ratio,carbon_count
  # furans or pyrroles
  elif ( nfuran>=1 or npyrrole>=1 ):
      mechspecies, mole_ratio = 'OLE', 2
      carbon_count = nC - 4
      if carbon_count > 0:
          if ( nCdblC>2 and ntermalke>0 ):
              nCdblC = nCdblC - 1
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'OLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
              else:
                  mechspecies = mechspecies,'OLE',
                  mole_ratio  = mole_ratio,ntermalke
              carbon_count = carbon_count - 2 * ntermalke
          if ( nCdblC>2 and carbon_count>=4 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],1
              else:
                  mechspecies = mechspecies,'IOLE',
                  mole_ratio  = mole_ratio,1
              carbon_count = carbon_count - 4
          if ( nketone>0 and carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'KET',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
              else:
                  mechspecies = mechspecies,'KET',
                  mole_ratio  = mole_ratio,nketone
              carbon_count = carbon_count - nketone
          if ( carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'PAR',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
              else:
                  mechspecies = mechspecies,'PAR',
                  mole_ratio  = mole_ratio,carbon_count
  # Heterocyclic aromatic compounds w/ 2 non-carbon atoms
  elif ( (naroN>=2) or (naroN>=1 and nO>=1) and nbenzene==0):
      mechspecies, mole_ratio = 'OLE', 1
      carbon_count = nC - 2
      if carbon_count > 0:
          if ( nCdblC>2 and ntermalke>0 ):
              nCdblC = nCdblC - 1
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'OLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
              else:
                  mechspecies = mechspecies,'OLE',
                  mole_ratio  = mole_ratio,ntermalke
              carbon_count = carbon_count - 2 * ntermalke
          if ( nCdblC>2 and carbon_count>=4 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],1
              else:
                  mechspecies = mechspecies,'IOLE',
                  mole_ratio  = mole_ratio,1
              carbon_count = carbon_count - 4
          if ( nketone>0 and carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'KET',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
              else:
                  mechspecies = mechspecies,'KET',
                  mole_ratio  = mole_ratio,nketone
              carbon_count = carbon_count - nketone
          if ( carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'PAR',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
              else:
                  mechspecies = mechspecies,'PAR',
                  mole_ratio  = mole_ratio,carbon_count
  # Triple bonds; contains no other functional groups
  elif ( nCtripC>=1 and nBranch==0 and nCdblC==0 and nacid==0 and ncarbonyl==0 and nalcohol==0):
      mechspecies, mole_ratio = 'OLE', 1
      carbon_count = nC - 2
      if carbon_count > 0:
          mechspecies = mechspecies,'PAR',
          mole_ratio  = mole_ratio,carbon_count
  # Triple bonds; contains other functional groups
  elif ( nCtripC>=1 and (nBranch>=1 or nCdblC>=1 or nacid>=1 or ncarbonyl>=1 or nalcohol>=1)):
      mechspecies, mole_ratio = 'PAR', nCtripC
      carbon_count = nC - nCtripC
      if carbon_count > 0:
          if ( (nCdblC>1 or ntermalke>0)  and carbon_count>=2):
              nCdblC = nCdblC - 1
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'OLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
              else:
                  mechspecies = mechspecies,'OLE',
                  mole_ratio  = mole_ratio,ntermalke
              carbon_count = carbon_count - 2 * ntermalke
          if ( nCdblC>2 and carbon_count>=4 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],1
              else:
                  mechspecies = mechspecies,'IOLE',
                  mole_ratio  = mole_ratio,1
              carbon_count = carbon_count - 4
          if ( nketone>0 and carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'KET',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
              else:
                  mechspecies = mechspecies,'KET',
                  mole_ratio  = mole_ratio,nketone
              carbon_count = carbon_count - nketone
          if ( carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'PAR',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
              else:
                  mechspecies = mechspecies,'PAR',
                  mole_ratio  = mole_ratio,carbon_count
  # Single-ring aromatics
  elif ( nbenzene==1 ): # Single-ring aromatics
       if ( nC>=7 and nalcohol>=2 ):      mechspecies, mole_ratio = 'CAT1', 1  # Methyl-catechols
       elif ( nC>=7 and nalcohol>=1 and nN>=1 ):
                                          mechspecies, mole_ratio = 'CRON', 1  # Nitro-cresols
       elif ( nC==7 and nalcohol==1 ):    mechspecies, mole_ratio = 'CRES', 1  # Cresols
       elif ( nHalo >= 4):                mechspecies, mole_ratio = 'UNR', 6   # Halogenated Benzenes
       elif ( nHalo > 0):                                                      # Halogenated Benzenes
           mechspecies = 'PAR','UNR',
           mole_ratio  = 1, 5
       elif ( nC>=7 and nBranch==1 ): # Toluene and other monoalkyl aromatics
           mechspecies, mole_ratio = 'TOL', 1
           carbon_count = nC - 7
           if carbon_count > 0:
               if len(mechspecies)==2:
                   mechspecies = mechspecies[0],mechspecies[1],'PAR',
                   mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
               else:
                   mechspecies = mechspecies,'PAR',
                   mole_ratio  = mole_ratio,carbon_count
       elif ( nC>=8 and nBranch>1 ): # Xylene and other polyalkyl aromatics
           mechspecies, mole_ratio = 'XYLMN', 1
           carbon_count = nC - 8
           if carbon_count > 0:
               if len(mechspecies)==2:
                   mechspecies = mechspecies[0],mechspecies[1],'PAR',
                   mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
               else:
                   mechspecies = mechspecies,'PAR',
                   mole_ratio  = mole_ratio,carbon_count
       else:
           carbon_count = nC
           if carbon_count > 0:
               if ( nCdblC==4 and ntermalke>0 ):
                   nCdblC = nCdblC - 1
                   try: mechspecies
                   except NameError:
                       mechspecies = 'OLE'
                       mole_ratio  = ntermalke
                   else:
                       if len(mechspecies)==2:
                           mechspecies = mechspecies[0],mechspecies[1],'OLE',
                           mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
                       else:
                           mechspecies = mechspecies,'OLE',
                           mole_ratio  = mole_ratio,ntermalke
                   carbon_count = carbon_count - 2 * ntermalke
               if ( nCdblC==4 and carbon_count>=4 ):
                   try: mechspecies
                   except NameError:
                       mechspecies = 'IOLE'
                       mole_ratio  = 1
                   else:
                       if len(mechspecies)==2:
                           mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                           mole_ratio  = mole_ratio[0],mole_ratio[1],1
                       else:
                           mechspecies = mechspecies,'IOLE',
                           mole_ratio  = mole_ratio,1
                   carbon_count = carbon_count - 4
               if ( nketone>0 and carbon_count>0 ):
                   try: mechspecies
                   except NameError:
                       mechspecies = 'KET'
                       mole_ratio  = nketone
                   else:
                       if len(mechspecies)==2:
                           mechspecies = mechspecies[0],mechspecies[1],'KET',
                           mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
                       else:
                           mechspecies = mechspecies,'KET',
                           mole_ratio  = mole_ratio,nketone
                   carbon_count = carbon_count - nketone
               if ( carbon_count>0 ):
                   try: mechspecies
                   except NameError:
                       mechspecies = 'PAR'
                       mole_ratio  = carbon_count
                   else:
                       if len(mechspecies)==2:
                           mechspecies = mechspecies[0],mechspecies[1],'PAR',
                           mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
                       else:
                           mechspecies = mechspecies,'PAR',
                           mole_ratio  = mole_ratio,carbon_count
  # Multi-ring aromatics
  elif ( nbenzene > 1 and nO/nC == 0 ):
      if ( nCdblC>5 and ntermalke>0 ):
          nCdblC = nCdblC - 1
          try: mechspecies
          except NameError:
              mechspecies = 'OLE'
              mole_ratio  = ntermalke
          else:
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'OLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
              else:
                  mechspecies = mechspecies,'OLE',
                  mole_ratio  = mole_ratio,ntermalke
          carbon_count = carbon_count - 2 * ntermalke
      if ( nCdblC>5 and carbon_count>=4 ):
          try: mechspecies
          except NameError:
              mechspecies = 'IOLE'
              mole_ratio  = 1
          else:
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],1
              else:
                  mechspecies = mechspecies,'IOLE',
                  mole_ratio  = mole_ratio,1
          carbon_count = carbon_count - 4
      if ( nketone>0 and carbon_count>0 ):
          try: mechspecies
          except NameError:
              mechspecies = 'KET'
              mole_ratio  = nketone
          else:
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'KET',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
              else:
                  mechspecies = mechspecies,'KET',
                  mole_ratio  = mole_ratio,nketone
          carbon_count = carbon_count - nketone
      if ( carbon_count>0 ):
          try: mechspecies
          except NameError:
              mechspecies = 'PAR'
              mole_ratio  = carbon_count
          else:
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'PAR',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
              else:
                  mechspecies = mechspecies,'PAR',
                  mole_ratio  = mole_ratio,carbon_count  
  # Higher aldehydes
  elif ( nC>=2 and naldehyde>=1 ):  # Propionaldehyde and higher aldehydes
      mechspecies, mole_ratio = 'ALDX', 1
      carbon_count = nC - 2
      if carbon_count > 0:
          if ( nCdblC>0 and ntermalke>0 ):
              nCdblC = nCdblC - 1
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'OLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
              else:
                  mechspecies = mechspecies,'OLE',
                  mole_ratio  = mole_ratio,ntermalke
              carbon_count = carbon_count - 2 * ntermalke
          if ( nCdblC>0 and carbon_count>=4 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],1
              else:
                  mechspecies = mechspecies,'IOLE',
                  mole_ratio  = mole_ratio,1
              carbon_count = carbon_count - 4
          if ( nketone>0 and carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'KET',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
              else:
                  mechspecies = mechspecies,'KET',
                  mole_ratio  = mole_ratio,nketone
              carbon_count = carbon_count - nketone
          if ( carbon_count>0 ):
              if len(mechspecies)==2:
                  mechspecies = mechspecies[0],mechspecies[1],'PAR',
                  mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
              else:
                  mechspecies = mechspecies,'PAR',
                  mole_ratio  = mole_ratio,carbon_count
  # Higher acids
  elif ( nC>=2 and nacid>=1 ): # Peroxyacetic and higher peroxycarboxylic acids
      mechspecies, mole_ratio = 'PACD', 1
      carbon_count = nC - 2
      if carbon_count > 0:
          if len(mechspecies)==2:
              mechspecies = mechspecies[0],mechspecies[1],'PAR',
              mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
          else:
              mechspecies = mechspecies,'PAR',
              mole_ratio  = mole_ratio,carbon_count
 # Other
  else:
      carbon_count = nC
      if carbon_count > 0:
          if ( nCdblC>0 and ntermalke>0 ):
              nCdblC = nCdblC - ntermalke
              try: mechspecies
              except NameError:
                  mechspecies = 'OLE'
                  mole_ratio  = ntermalke
              else:
                  if len(mechspecies)==2:
                      mechspecies = mechspecies[0],mechspecies[1],'OLE',
                      mole_ratio  = mole_ratio[0],mole_ratio[1],ntermalke
                  else:
                      mechspecies = mechspecies,'OLE',
                      mole_ratio  = mole_ratio,ntermalke
              carbon_count = carbon_count - 2 * ntermalke
          if ( nCdblC>0 and carbon_count>=4 ):
              try: mechspecies
              except NameError:
                  mechspecies = 'IOLE'
                  mole_ratio  = 1
              else:
                  if len(mechspecies)==2:
                      mechspecies = mechspecies[0],mechspecies[1],'IOLE',
                      mole_ratio  = mole_ratio[0],mole_ratio[1],1
                  else:
                      mechspecies = mechspecies,'IOLE',
                      mole_ratio  = mole_ratio,1
              carbon_count = carbon_count - 4
          if ( nketone>0 and carbon_count>0 ):
              try: mechspecies
              except NameError:
                  mechspecies = 'KET'
                  mole_ratio  = nketone
              else:
                  if len(mechspecies)==2:
                      mechspecies = mechspecies[0],mechspecies[1],'KET',
                      mole_ratio  = mole_ratio[0],mole_ratio[1],nketone
                  else:
                      mechspecies = mechspecies,'KET',
                      mole_ratio  = mole_ratio,nketone
              carbon_count = carbon_count - nketone
          if ( carbon_count>0 ):
              try: mechspecies
              except NameError:
                  mechspecies = 'PAR'
                  mole_ratio  = carbon_count
              else:
                  if len(mechspecies)==2:
                      mechspecies = mechspecies[0],mechspecies[1],'PAR',
                      mole_ratio  = mole_ratio[0],mole_ratio[1],carbon_count
                  else:
                      mechspecies = mechspecies,'PAR',
                      mole_ratio  = mole_ratio,carbon_count   
  # Limit functional groups to 1 with this priority: TOL > XYL > IOLE > OLE > ALDX > KET; per Greg Yarwood
  countfg   = 0
  countole  = 0
  if 'IOLE' in mechspecies: countfg += 1
  for ele in mechspecies:
      if (ele == 'OLE'):
          countole += 1
          if (ele == 'IOLE'): countole -= 1
  if countole > 0:  countfg += 1
  if 'ALDX' in mechspecies:  countfg += 1
  if 'KET' in mechspecies:  countfg += 1
  
  if countfg > 1:
      carbon_count = nC
      if 'IOLE' in mechspecies:
          ind = 0
          for ele in mechspecies:
              if (ele == 'IOLE'):
                  temp_mechspecies = 'IOLE'
                  temp_mole_ratio  = mole_ratio[ind]
                  carbon_count -= 4
                  if ( carbon_count>0 ):
                      temp_mechspecies = temp_mechspecies,'PAR',
                      temp_mole_ratio  = temp_mole_ratio,carbon_count
              else: pass
              ind += 1
      elif 'OLE' in mechspecies:
          ind = 0
          for ele in mechspecies:
              if (ele == 'OLE'):
                  temp_mechspecies = 'OLE'
                  temp_mole_ratio  = mole_ratio[ind]
                  carbon_count -= 2
                  if ( carbon_count>0 ):
                      temp_mechspecies = temp_mechspecies,'PAR',
                      temp_mole_ratio  = temp_mole_ratio,carbon_count
              else: pass
              ind += 1
      elif 'KET' in mechspecies:
          ind = 0
          for ele in mechspecies:
              if (ele == 'KET'):
                  temp_mechspecies = 'KET'
                  temp_mole_ratio  = mole_ratio[ind]
                  carbon_count -= 1
                  if ( carbon_count>0 ):
                      temp_mechspecies = temp_mechspecies,'PAR',
                      temp_mole_ratio  = temp_mole_ratio,carbon_count
              else: pass
              ind += 1
      elif 'ALDX' in mechspecies:
          ind = 0
          for ele in mechspecies:
              if (ele == 'ALDX'):
                  temp_mechspecies = 'ALDX'
                  temp_mole_ratio  = mole_ratio[ind]
                  carbon_count -= 2
                  if ( carbon_count>0 ):
                      temp_mechspecies = temp_mechspecies,'PAR',
                      temp_mole_ratio  = temp_mole_ratio,carbon_count
              else: pass
              ind += 1
      mechspecies = temp_mechspecies
      mole_ratio  = temp_mole_ratio
  # If OLE/PAR > 1, OLE = (PAR + 2 OLE)/3 and then PAR = OLE; per Greg Yarwood
  countpar  = 0
  countole  = 0
  if 'PAR' in mechspecies and 'IOLE' in mechspecies:
      pass
  elif 'PAR' in mechspecies and 'OLE' in mechspecies:
      ind = 0
      for ele in mechspecies:
          if (ele == 'PAR'):
              countpar = mole_ratio[ind]
          ind += 1
      ind = 0
      for ele in mechspecies:
          if (ele == 'OLE'):
              countole = mole_ratio[ind]
          ind += 1
  elif 'IOLE' in mechspecies:
      pass
  elif 'OLE' in mechspecies:
      countole = mole_ratio
  if countpar == 0: countpar = 0.1
  if countole / countpar > 1:
      if countpar == 0.1: countpar = 0
      temp_mechspecies = 'OLE'
      temp_mole_ratio  = (countpar + 2 * countole) / 3
      temp_mechspecies = temp_mechspecies,'PAR',
      temp_mole_ratio  = temp_mole_ratio,((countpar + 2 * countole) / 3)
      mechspecies = temp_mechspecies
      mole_ratio  = temp_mole_ratio
  # If KET/PAR > 0.333, KET = (PAR + KET)/4 and then PAR = 3*KET; per Greg Yarwood
  countpar  = 0
  countket  = 0
  if 'PAR' in mechspecies and 'KET' in mechspecies:
      ind = 0
      for ele in mechspecies:
          if (ele == 'PAR'):
              countpar = mole_ratio[ind]
          ind += 1
      ind = 0
      for ele in mechspecies:
          if (ele == 'KET'):
              countket = mole_ratio[ind]
          ind += 1
  elif 'KET' in mechspecies:
      countket = mole_ratio
  if countpar == 0: countpar = 0.1
  if countket / countpar > 0.3334:
      if countpar == 0.1: countpar = 0
      temp_mechspecies = 'KET'
      temp_mole_ratio  = (countpar + countket) / 4
      temp_mechspecies = temp_mechspecies,'PAR',
      temp_mole_ratio  = temp_mole_ratio,((countpar + countket) / 4 * 3)
      mechspecies = temp_mechspecies
      mole_ratio  = temp_mole_ratio
  return mechspecies, mole_ratio
  # end of function

############################################################################################
#### Get input dataset.
df = pd.read_excel('./export_species_properties.xlsx')
############################################################################################

dfappend_cb6r3_ae7(df)
############################################################################################

print("Time to generate mechanism for import file: ",datetime.now() - startTime)