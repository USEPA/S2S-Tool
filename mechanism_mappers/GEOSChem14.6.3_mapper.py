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

def dfappend_geoschem1463(dfin):
  '''
    Compound to mechanism mapping
  '''

  mech = 'GEOSChem14.6.3'
  dfin['GEOSChem14.6.3']=''  

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
                                    mechspecies = 'UNK' # Unknown
    else: mechspecies = get_geoschem_roc(smiles,koh,log10cstar)

    mech4import = pd.Series(data={'mechanism':mech,'SPECIES_ID':row['SPECIES_ID'],
                                  'mech_species':mechspecies,'moles_ratio':1})
    dfmech4import = pd.concat([dfmech4import,pd.DataFrame([mech4import])],ignore_index=True)
    #dfmech4import = dfmech4import.concat(mech4import,ignore_index=True)

  # write mech4import df to file
  today = date.today()
  dfmech4import.to_csv('./mechanism_forImport_'+mech+'_speciate5_2_'+str(today)+'.csv',index=False,header=False)
  
  return 

def get_geoschem_roc(smiles,koh,log10cstar,phase=None):
    '''
    Function maps input reactive organic carbon (ROC) species to GEOS-Chem v14.6.3 ROC species.
    Uses functional group and molecule info from RDKit http://www.rdkit.org/
    Function inputs, for ONE compound (make loop outside this function):
        smiles string (should be canonical for explicit species, some alt added) 
        kOH (in cm3/molec-s)
        log10(Cstar in micrograms/m3)
    kOH and C* may not be used if the compound can be mapped without it. 
    '''
    
    # Prep inputs
    m            = Chem.MolFromSmiles(smiles)
    smiles_upper = smiles.upper()
    
    # Count C=C and atoms
    nCdblC  = smiles_upper.count('=C')-smiles_upper.count('O=C')
    if nCdblC < 0:
        nCdblC = 0
    nC      = smiles_upper.count('C')-smiles_upper.count('CL')
    nO      = smiles_upper.count('O')
    nN      = smiles_upper.count('N')
    nSi     = smiles_upper.count('SI')
    nH      = 0
    for atom in m.GetAtoms():
        nH += atom.GetTotalNumHs()
    # O:C ratio
    if nC > 0:
        OtoC = nO/nC
    
    # Count functional groups (http://rdkit.org/docs/source/rdkit.Chem.Fragments.html)
    nacid     = rdkit.Chem.Fragments.fr_COO(m,countUnique=True)     # carboxylic acid
    nketone   = rdkit.Chem.Fragments.fr_ketone(m,countUnique=True)
    naldehyde = rdkit.Chem.Fragments.fr_aldehyde(m,countUnique=True)
    ncarbonyl = nketone + naldehyde
    nbenzene  = rdkit.Chem.Fragments.fr_benzene(m,countUnique=True)
    nalcohol  = rdkit.Chem.Fragments.fr_Al_OH(m,countUnique=True) + \
                  rdkit.Chem.Fragments.fr_Ar_OH(m,countUnique=True)      # aliphatic and aromatic
    nfuran    = rdkit.Chem.Fragments.fr_furan(m,countUnique=True) # number of furan rings
    # gotring variable is never used so this could be removed
    #     it may be useful to keep it here commented out in case it is needed in the future
    #for atom in range(len(m.GetAtoms())):
        #if m.GetAtomWithIdx(atom).IsInRing():
        #    gotring = 1 # 0 = no ring, 1 = ring
        #    break
        #else: gotring = 0
    nnitrate = smiles_upper.count('O[N+](=O)[O-]') + smiles_upper.count('O[N+]([O-])=O')
    nperoxide = smiles_upper.count('COO') +  smiles_upper.count('OOC') - smiles_upper.count('COOC')
    #namine = smiles_upper.count('CN') + smiles_upper.count('NC') + smiles_upper.count('C(N') + smiles_upper.count('N(C')
    tfmonoterpene = (nC == 10 and nH == 18 and nO == 1) or (nC == 10 and nH == 16) 

    # Mapper is for ROC only and not elemental carbon
    if   ( nC <= 0 ):                   mechspecies = 'UNK'     # Unknown
    elif ( smiles == '[C]' ):           mechspecies = 'UNK'     # Map elemental carbon to UNK
    elif ( smiles == 'C#[O+]' ):        mechspecies = 'UNK'     # Map CO to UNK

    # Explicit species
    elif ( smiles == 'CC(C)=O' ):       mechspecies = 'ACET'    # acetone
    elif ( smiles == 'C=CC=O' ):        mechspecies = 'ACR'     # acrolein
    elif ( smiles == 'CC(O)=O' ):       mechspecies = 'ACTA'    # acetic acid
    elif ( smiles == 'CC=O' ):          mechspecies = 'ALD2'    # acetaldehyde
    elif ( nC==6 and nH==6 and nO==0 and nbenzene==1 ):
                                        mechspecies = 'BENZ'    # benzene    
    elif ( smiles == 'C#C' ):           mechspecies = 'C2H2'    # acetylene aka ethyne
    elif ( smiles == 'C=C' ):           mechspecies = 'C2H4'    # ethylene aka ethene
    elif ( smiles == 'CC' ):            mechspecies = 'C2H6'    # ethane
    elif ( smiles == 'CCC' ):           mechspecies = 'C3H8'    # propane
    elif ( smiles == 'C=CC=C' ):        mechspecies = 'C4H6'    # 1,3 butadiene   
    elif ( smiles == 'ClC(Cl)(Cl)Cl' ): mechspecies = 'CCl4'    # carbon tetrachloride
    elif ( smiles == 'BrCBr' ):         mechspecies = 'CH2Br2'  # dibromomethane aka methylene bromide
    elif ( smiles == 'ClCCl' ):         mechspecies = 'CH2Cl2'  # dichloromethane aka methylene chloride
    elif ( smiles == 'C=O' ):           mechspecies = 'CH2O'    # formaldehyde
    elif ( smiles == 'CBr' ):           mechspecies = 'CH3Br'   # methyl bromide aka bromomethane
    elif ( smiles == 'CCl' ):           mechspecies = 'CH3Cl'   # methyl chloride aka chloromethane
    elif ( smiles == 'CI' ):            mechspecies = 'CH3I'    # methyl iodide
    elif ( smiles == 'C' ):             mechspecies = 'CH4'     # methane
    elif ( smiles == 'ClC(Cl)Cl' ):     mechspecies = 'CHCl3'   # chloroform aka trichlormethane
    elif ( smiles == 'CCC1=CC=CC=C1' ): mechspecies = 'EBZ'     # ethylbenzene
    elif ( smiles == 'CCO' ):           mechspecies = 'EOH'     # ethanol
    elif ( smiles == 'O=CC=O' ):        mechspecies = 'GLYX'    # glyoxal
    elif ( smiles == 'OCC(O)=O' ):      mechspecies = 'HACTA'   # glycolic acid aka hydroxyacetic acid
    elif ( smiles == 'OC=O' ):          mechspecies = 'HCOOH'   # formic acid
    elif ( smiles == 'CC(=C)C=C' ):     mechspecies = 'ISOP'    # isoprene
    elif ( smiles == 'CC(=C)C1CCC(C)=CC1' ):
                                        mechspecies = 'LIMO'    # limonene
    elif ( smiles == 'CC(=C)C=O' ):     mechspecies = 'MACR'    # methacrolein aka 2-methyl-2-propenal
    elif ( smiles == 'CCC(C)=O' ):      mechspecies = 'MEK'     # methyl ethyl ketone
    elif ( smiles == 'CC(=O)C=O' ):     mechspecies = 'MGLY'    # methylglyoxal
    elif ( smiles == 'CO'):             mechspecies = 'MOH'     # methanol
    elif ( smiles == 'CC(=O)C=C' ):     mechspecies = 'MVK'     # methyl vinyl ketone
    elif ( smiles == 'C1=CC2=CC=CC=C2C=C1' ):
                                        mechspecies = 'NAP'     # naphthalene
    elif ( smiles == 'O=C=S' ):         mechspecies = 'OCS'     # glyoxal
    elif ( smiles == 'OC1=CC=CC=C1' ):  mechspecies = 'PHEN'    # phenol aka carbolic acid
    elif ( smiles == 'C=CC1=CC=CC=C1'): mechspecies = 'STYR'    # styrene
    elif ( smiles == 'CC1=CC=CC=C1' ):  mechspecies = 'TOLU'    # toluene

    # Send SVOCs and lower volatility ROC to UNR
    elif ( log10cstar < 2.5 ):          mechspecies = 'UNR'     # SVOCs
    
    # Low reactivity species
    elif ( koh < 3.5e-13 ):             mechspecies = 'UNR'     # unreactive; low reactivity gas

    # IVOC species with nSi > 0
    elif ( nSi > 0 ):                   mechspecies = 'IVOC'    # Any silanes/siloxanes
                                        
    # Multi-ring aromatics (PAH and NAPH can be collapsed together if necessary)  
    elif ( nbenzene > 1 ):              mechspecies = 'NAP'     # Naphthalene-like, PAH with 2 rings

    # Lumped terpene species 
    elif ( smiles == 'CC1=CCC2CC1C2(C)C' ): 
                                        mechspecies = 'MTPA'    # alpha-pinene
    elif ( smiles == 'CC1(C)C2CC1C(=C)CC2' ): 
                                        mechspecies = 'MTPA'    # beta-pinene
    elif ( smiles == 'CC(C)C12CCC(=C)C1C2' ): 
                                        mechspecies = 'MTPA'    # sabinene
    elif ( smiles == 'CC1=CCC2C(C1)C2(C)C' ): 
                                        mechspecies = 'MTPA'    # carene
    elif ( tfmonoterpene ):             mechspecies = 'MTPO'    # other monoterpenes

    # Furans and dienes other than 1,3 BDE
    elif ( nfuran > 0 or nCdblC==2):    mechspecies = 'FURA'    # furans and other dienes
    
    # IVOC species with nC > 10 binned
    elif ( log10cstar < 6.5 and nC>10):  
                                        mechspecies = 'IVOC'    # IVOCs with nC > 10

    # Single-ring aromatics (excluding explicit species)
    elif ( nbenzene > 0 ): # Single-ring aromatics
        if ( naldehyde > 0 ):           mechspecies = 'BALD'    # benzaldehyde, tolualdehyde, and arom. aldehydes
        elif ( nC>=7 and nalcohol>=1 ): mechspecies = 'CSL'     # cresol
        elif ( nC==9 ):                 mechspecies = 'TMB'     # trimethylbenzenes
        # any single-ring aromatics that have not been mapped by rules above
        else:                           mechspecies = 'XYLE'    # xylenes and other aromatics

    # Species with double bonds, not aromatic
    elif ( nCdblC==1 and nC>=3 ):       mechspecies = 'PRPE'    # >= C3 alkenes
    
    # Oxygenated species without double bonds
    elif ( ncarbonyl>=1 and nC>=3 ):    mechspecies = 'RCHO'    # >= C3 aldehydes and ketones
    elif ( nalcohol>=1 and nC>2 ):      mechspecies = 'ROH'     # > C2 alcohols 
    elif ( nacid>=1 and nC>2 ):         mechspecies = 'RCOOH'   # > C2 organic acids

    # "Alkane" Series
    elif ( nC>3 and nC<=5 and OtoC==0 ):   
                                        mechspecies = 'ALK4'
    elif ( nC>=6 and nC<=10 and OtoC==0 ):  
                                        mechspecies = 'ALK6'
    # Catch remaining IVOCs
    elif ( log10cstar < 6.5 ):          mechspecies = 'IVOC'    # IVOCs

    else:                               mechspecies = 'UNK'     # Species is unknown
    
    return mechspecies
    # end of function

############################################################################################
#### Get input dataset.
df = pd.read_excel('./export_species_properties.xlsx')
############################################################################################

dfappend_geoschem1463(df)
############################################################################################

print("Time to generate mechanism for import file: ",datetime.now() - startTime)