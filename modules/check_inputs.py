import sys
import pandas as pd

####################################################################################################
### These functions perform QA checks on input files. 
####################################################################################################

####################################################################################################
def check_basis_output(MECH_BASIS,OUTPUT):
    if OUTPUT=='VOC':
        if MECH_BASIS=='CB6R3_AE7' or MECH_BASIS=='CB6R3_AE8' or MECH_BASIS=='CB6R5_AE7' or \
           MECH_BASIS=='CB6R5_AE8' or MECH_BASIS=='CRACMMv1.0' or MECH_BASIS=='SAPRC07TC_AE7' or \
           MECH_BASIS=='SAPRC07TC_AE8'  or MECH_BASIS=='CB6R4_CF2' or MECH_BASIS=='SAPRC07_CF2' or \
           MECH_BASIS=='CB6R3_AE7_TRACER': pass
        else:
            print('MECH_BASIS for OUTPUT==VOC is not allowed.')
            sys.exit('Only CB6R3_AE7, CB6R3_AE8, CB6R5_AE7, CB6R5_AE8, CRACMMv0.3 are accepted.')
    elif OUTPUT=='PM':
        if MECH_BASIS=='PM-AE6' or MECH_BASIS=='PM-AE8' or MECH_BASIS=='PM-CR1':
            pass
        else:
            print('MECH_BASIS for OUTPUT==PM is not allowed.')
            sys.exit('Only PM-AE6, PM-AE8, and PM-CR1 are accepted.')
    else: sys.exit('OUTPUT entered is not recognized. Only VOC and PM are allowed.')
####################################################################################################

####################################################################################################
def check_inputs(carbons,mech4import,mechPM,tbl_tox,OUTPUT):
    if carbons.empty and OUTPUT=='VOC':
        sys.exit('The MECH_BASIS entered is not in the carbons file.')
    else:
        pass
    if mech4import.empty and OUTPUT=='VOC':
        sys.exit('The MECH_BASIS entered is not in the mechanism_forImport file.')
    else:
        pass
    if mechPM.empty and OUTPUT=='PM':
        sys.exit('The MECH_BASIS and AQM entered is not in the mech_pm file.')
    else:
        pass
    if tbl_tox.empty:
        sys.exit('The AQM entered is not in the tbl_tox file.')
    else:
        pass
####################################################################################################

####################################################################################################
def check_species_profiles(profiles,species):
    prof_sum  = species.groupby('PROFILE_CODE',as_index=False)['WEIGHT_PERCENT'].sum()  # Sum WEIGHT_PERCENT for each profile
    temp_prof = pd.merge(profiles,prof_sum[['PROFILE_CODE','WEIGHT_PERCENT']],on='PROFILE_CODE',how ='left') # append WEIGHT_PERCENT
    zero_sums = temp_prof[temp_prof['WEIGHT_PERCENT'].isna()]
    zero_sums = zero_sums.reset_index(drop=True) # reset index
    if zero_sums.empty:
        pass
    else:
        print('Profile(s) are in export_profiles but not in export_species. These include:')
        print(zero_sums.loc[:,'PROFILE_CODE'])
        sys.exit()
####################################################################################################

####################################################################################################
def check_species_properties(species_props,species):
    temp_spec = pd.merge(species,species_props[['SPECIES_ID','SPEC_MW']],on='SPECIES_ID',how ='left').fillna(0) # append MolWght
    zero_MWs  = temp_spec.loc[temp_spec['SPEC_MW'] == 0.0]
    zero_MWs  = zero_MWs.reset_index(drop=True) # reset index
    if zero_MWs.empty:
        pass
    else:
        print('Species are in export_species but not in export_species_properties. These include:')
        print(zero_MWs.loc[:,'SPECIES_ID'])
        sys.exit()
####################################################################################################

####################################################################################################
def check_tox_properties(species_props,tbl_tox):
    temp_tox  = pd.merge(tbl_tox,species_props[['SPECIES_ID','SPEC_MW']],on='SPECIES_ID',how ='left').fillna(0) # append MolWght
    zero_MWs  = temp_tox.loc[temp_tox['SPEC_MW'] == 0.0]
    zero_MWs  = zero_MWs.reset_index(drop=True) # reset index
    if zero_MWs.empty:
        pass
    else:
        print('Species are in tbl_tox but not in export_species_properties. These include:')
        print(zero_MWs.loc[:,'SPECIES_ID'])
        sys.exit()
####################################################################################################

####################################################################################################
def check_species_mech4import(mech4import,species,MECH_BASIS):
    spec_sum  = mech4import.groupby('SPECIES_ID',as_index=False)['Moles'].sum()  # Sum moles for each SPECIES_ID
    temp_spec = pd.merge(species,spec_sum[['SPECIES_ID','Moles']],on='SPECIES_ID',how ='left').fillna(0) # append moles
    zero_mol  = temp_spec.loc[temp_spec['Moles'] == 0.0]
    zero_mol  = zero_mol.reset_index(drop=True) # reset index
    if zero_mol.empty:
        pass
    else:
        print('Species are in export_species but not in the mechanism_forImport file. These include:')
        print(zero_mol.loc[:,'SPECIES_ID'])
        print('For MECH_BASIS: '+MECH_BASIS)
        sys.exit()
####################################################################################################

####################################################################################################
def check_profiles_volatility(profiles,poa_volatility):
    poa_copy  = poa_volatility.copy()
    poa_copy['Total'] = poa_copy['-2'] + poa_copy['-1'] + poa_copy['0'] + poa_copy['1'] + poa_copy['2']
    temp_prof = pd.merge(profiles,poa_copy,on=['CATEGORY_LEVEL_1_Generation_Mechanism','CATEGORY_LEVEL_2_Sector_Equipment'], \
                         how ='left').fillna(0)
    no_vol    = temp_prof.loc[temp_prof['Total'] == 0.0]
    no_vol    = no_vol.reset_index(drop=True) # reset index
    if no_vol.empty:
        pass
    else:
        print('PM profiles are in export_profiles whose OM volatility is not specified in poa_volatility. These include:')
        print(no_vol.loc[:,'PROFILE_CODE'])
        sys.exit()
####################################################################################################

####################################################################################################
def check_volatility(poa_volatility):
    poa_copy  = poa_volatility.copy()
    poa_copy['Total'] = poa_copy['-2'] + poa_copy['-1'] + poa_copy['0'] + poa_copy['1'] + poa_copy['2']
    wrong_vol = poa_copy.loc[poa_copy['Total'] != 1.0]
    wrong_vol = wrong_vol.reset_index(drop=True) # reset index
    if wrong_vol.empty:
        pass
    else:
        print('The POA volatility specified for some categories does not equal 1. These include:')
        print(wrong_vol.loc[:,'CATEGORY_LEVEL_1_Generation_Mechanism']+', '+wrong_vol.loc[:,'CATEGORY_LEVEL_2_Sector_Equipment'])
        sys.exit()
####################################################################################################