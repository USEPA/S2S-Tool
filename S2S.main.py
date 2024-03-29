import sys
from datetime import date,datetime
import pandas as pd

startTime = datetime.now()
today     = date.today()

####################################################################################################
### S2S-Tool: U.S. EPA's SPECIATE-to-SMOKE Tool generates GSPRO and GSCNV files that are used by 
### SMOKE to generate gridded emissions for photochemical modeling. All profiles are in the U.S. 
### EPA's SPECIATE database, and the type of GSPRO/GSCNV files generated by the S2S-Tool are 
### dependent on the speciation methods employed for the desired modeling platform. SPECIATE is 
### a database of organic gas, particle, and mercury speciation profiles. These profiles provide 
### greater specificity than what is needed for a chemical mechanism within a photochemical model.
### The S2S-Tool bridges this gap and translates SPECIATE data into a format that is chemical 
### mechanism specific.
####################################################################################################

####################################################################################################
### User Input
### Photochemical mechanism; options include CB6R3_AE7, CB6R5_AE7, CB7_AE7, CRACMMv1.0, SAPRC07TC_AE7, CB6R3_AE7_TRACER, CB6R4_CF2, CB7_CF2, SAPRC07_CF2, PM-AE6, PM-CR1
MECH_BASIS = 'CB6R3_AE7'
### Output type; options include VOC, PM
OUTPUT     = 'VOC'
### Select run type; options include CRITERIA, INTEGRATE, NOINTEGRATE
RUN_TYPE   = 'CRITERIA'
### Select air quality model; options include CMAQ, CAMX
AQM        = 'CMAQ'
### Assign acceptable deviation from 100% allowable. Applies only to gas profiles.
TOLERANCE  = 0.05 # 0.05 = 5%
### If applicable, input pollutant for toxics INTEGRATE scenario
TOX_IN     = 'NONHAPTOG' # NONHAPTOG ; TOM ; RESID_PM
### CAMx FCRS file name & path:
FCRS_FILE  = './input/camx_fcrs.profile.csv'
### MW file name & path:
MW_FILE    = './input/mechanism_mw.csv'
### mechanism_forImport file name & path:
M4I_FILE   = './input/mechanism_forImport_SPECIATEv5_3.csv'
### pm_mech file name & path:
PMM_FILE   = './input/mech_pm_ae5_ae6_cr1.csv'
### tbl_tox file name & path:
TOX_FILE   = './input/tbl_tox_NBAFM.csv' # 'tbl_tox_NBAFM.csv' ; 'tbl_tox_TOM.csv' ; 'tbl_tox_RESID_PM.csv' ; 'tbl_tox_MOVES_HAPS.csv'
####################################################################################################

####################################################################################################
### Fixed input data; can be updated by user as needed.
### GSPRO/GSCNV output file names
if RUN_TYPE == 'INTEGRATE':
    CNV_OUT    = './output/gscnv.'+MECH_BASIS+'_'+RUN_TYPE+'_'+TOX_IN+'_'+OUTPUT+'.'+AQM+'.'+str(today)+'.txt'
    PRO_OUT = './output/gspro.'+MECH_BASIS+'_'+RUN_TYPE+'_'+TOX_IN+'_'+OUTPUT+'.'+AQM+'.'+str(today)+'.txt'
else:
    CNV_OUT    = './output/gscnv.'+MECH_BASIS+'_'+RUN_TYPE+'_'+OUTPUT+'.'+AQM+'.'+str(today)+'.txt'
    PRO_OUT = './output/gspro.'+MECH_BASIS+'_'+RUN_TYPE+'_'+OUTPUT+'.'+AQM+'.'+str(today)+'.txt'
### Location of modules:
sys.path.append('./modules/')
### Import oxygen/metals ratios:
oxygen_metals    = pd.read_csv('./input/oxygen_metal_Ratios.csv')
### Import SPECIATE profiles table to be processed:
profiles         = pd.read_csv('./input/export_profiles.csv',converters={'PROFILE_CODE': str})
### Import SPECIATE species table to be processed:
species          = pd.read_csv('./input/export_species.csv',converters={'PROFILE_CODE': str})
### Import SPECIATE species properties table to facilitate processing:
species_props    = pd.read_csv('./input/export_species_properties.csv')
### Import list of SPECIATE profiles where FPRM --> FCRS:
camx_fcrs        = pd.read_csv(FCRS_FILE)
### Import list of profiles to append with variable IN/OUT pollutants:
gscnv_append     = pd.read_csv('./input/gscnv_append.csv')
### Import MW file:
molwght          = pd.read_csv(MW_FILE)
### Import mechanism_forImport file:
mech4import      = pd.read_csv(M4I_FILE)
### Import pm_mech file:
mechPM           = pd.read_csv(PMM_FILE)
### Import toxics table file:
tbl_tox          = pd.read_csv(TOX_FILE)
### Import POA volatility bin assignments:
poa_volatility   = pd.read_csv('./input/POA_VolatilityBins.csv')
### Import mechanism-specific mapping of OM/OC/NCOM:
poa_mapping      = pd.read_csv('./input/POA_mapping.csv')
####################################################################################################

####################################################################################################
### Filter imported data for the target MECH_BASIS, AQM, etc.
molwght          = molwght.loc[molwght['Mechanism'] == MECH_BASIS] # filter molwght for relevant MECH_BASIS
molwght          = molwght.reset_index(drop=True) # reset index
mech4import      = mech4import.loc[mech4import['Mechanism'] == MECH_BASIS] # filter mech4import for relevant MECH_BASIS
mech4import      = mech4import.reset_index(drop=True) # reset index
mechPM           = mechPM.loc[mechPM['Mechanism'] == MECH_BASIS] # filter mechPM for relevant MECH_BASIS
mechPM           = mechPM.reset_index(drop=True) # reset index
mechPM           = mechPM.loc[mechPM['AQM'] == AQM] # filter mechPM for relevant AQM
mechPM           = mechPM.reset_index(drop=True) # reset index
if OUTPUT == 'VOC':
    profiles     = profiles.loc[profiles['PROFILE_TYPE'] == 'GAS'] # filter profiles for relevant PROFILE_TYPEs
    profiles     = profiles.reset_index(drop=True) # reset index
elif OUTPUT == 'PM':
    PMLIST       = ['PM','PM-AE6','PM-AE8','PM-CR1'] # removes PM-Simplified
    profiles     = profiles.loc[profiles['PROFILE_TYPE'].isin(PMLIST)] # filter profiles for relevant PROFILE_TYPEs
    profiles     = profiles.reset_index(drop=True) # reset index
else: sys.exit('There is an issue with your OUTPUT entry. Only VOC and PM are accepted.')
tbl_tox          = tbl_tox.loc[tbl_tox['AQM'] == AQM] # filter tbl_tox for relevant AQM
tbl_tox          = tbl_tox.reset_index(drop=True) # reset index
temp_prof        = profiles.loc[:,'PROFILE_CODE'] # pull PROFILE_CODE list
species          = species.loc[species['PROFILE_CODE'].isin(temp_prof)] # filter species array for PROFILE_CODE
####################################################################################################

####################################################################################################
### This module contains several functions that perform QA checks on the input files.
import check_inputs
### This module generates the gscnv file for the target MECH_BASIS.
import gscnv
### This module generates the gspro file for the target MECH_BASIS.
import gspro
####################################################################################################

####################################################################################################
print('NOTICE: Type of Output is '+OUTPUT)
print('NOTICE: AQM is '+AQM)
print('NOTICE: Type of run is '+RUN_TYPE)
print('NOTICE: Mechanism basis is '+MECH_BASIS)
print('NOTICE: Profile Tolerance is ',TOLERANCE)
####################################################################################################

### QA check: ensure MECH_BASIS and OUTPUT make sense
check_inputs.check_basis_output(MECH_BASIS,OUTPUT)
### QA check: ensure mech4import, mechPM, etc are not empty following filters above
check_inputs.check_inputs(molwght,mech4import,mechPM,tbl_tox,OUTPUT)
### QA check: all profiles in export_profiles are in export_species
check_inputs.check_species_profiles(profiles,species)
### QA check: all species in export_species are in export_species_properties
check_inputs.check_species_properties(species_props,species)
### QA check: all species in tbl_tox are in export_species_properties
check_inputs.check_tox_properties(species_props,tbl_tox)
### QA check: all species in molwght are in mech4import
check_inputs.check_molwght_mech4import(molwght,mech4import,OUTPUT)
### QA check: all species in export_species are in mech4import
if OUTPUT=='VOC':
    check_inputs.check_species_mech4import(mech4import,species,MECH_BASIS)
elif OUTPUT == 'PM':
    pass
else: sys.exit('There is an issue with your OUTPUT entry. Only VOC and PM are accepted.')
### QA check: all PM profile types are in the SV-POA list
if OUTPUT=='VOC':
    pass
elif OUTPUT == 'PM':
    check_inputs.check_profiles_volatility(profiles,poa_volatility)
else: sys.exit('There is an issue with your OUTPUT entry. Only VOC and PM are accepted.')
### QA check: all SV-POA entries sum to 1
if OUTPUT=='VOC':
    pass
elif OUTPUT == 'PM':
    check_inputs.check_volatility(poa_volatility)
else: sys.exit('There is an issue with your OUTPUT entry. Only VOC and PM are accepted.')

### Generate GSCNV file:
if OUTPUT=='VOC':
    if RUN_TYPE=='CRITERIA' or RUN_TYPE=='INTEGRATE' or RUN_TYPE=='NOINTEGRATE':
        print('NOTICE: Beginning generation of GSCNV file.')
        ### Create gscnv file for the target MECH_BASIS
        gscnv.gen_gscnv(profiles,species,species_props,tbl_tox,gscnv_append,MECH_BASIS,RUN_TYPE,TOLERANCE,CNV_OUT)
        ### Formats and adds header to gscnv file
        gscnv.format_and_header(MECH_BASIS,RUN_TYPE,AQM,MW_FILE,TOX_FILE,CNV_OUT)
        print('NOTICE: Finished generating GSCNV file.')
    else: sys.exit('There is an issue with your RUN_TYPE entry. Only CRITERIA, INTEGRATE, and NOINTEGRATE are accepted.')
elif OUTPUT == 'PM': 
    print('NOTICE: Type of Output is '+OUTPUT+', no GSCNV generated.')
    pass
else: sys.exit('There is an issue with your OUTPUT entry. Only VOC and PM are accepted.')

### Generate GSPRO file:
if OUTPUT=='VOC':
    if RUN_TYPE=='CRITERIA' or RUN_TYPE=='INTEGRATE' or RUN_TYPE=='NOINTEGRATE':
        print('NOTICE: Beginning generation of GSPRO file.')
        ### Create gspro file for the target MECH_BASIS
        gspro.gen_gspro_voc(profiles,species,species_props,molwght,mech4import,tbl_tox,MECH_BASIS,RUN_TYPE,TOLERANCE,TOX_IN,PRO_OUT)
        ### Formats and adds header to gspro file
        gspro.format_and_header(tbl_tox,TOX_IN,MECH_BASIS,RUN_TYPE,AQM,MW_FILE,FCRS_FILE,TOX_FILE,PRO_OUT)
        print('NOTICE: Finished generating GSPRO file.')
    else: sys.exit('There is an issue with your RUN_TYPE entry. Only CRITERIA, INTEGRATE, and NOINTEGRATE are accepted.')
elif OUTPUT == 'PM':
    if AQM=='CMAQ' or (AQM=='CAMX' and MECH_BASIS=='PM-AE6'):
        if RUN_TYPE=='CRITERIA' or RUN_TYPE=='INTEGRATE':
            print('NOTICE: Beginning generation of GSPRO file.')
            ### Create gspro file for the target MECH_BASIS
            gspro.gen_gspro_pm(profiles,species,species_props,mechPM,tbl_tox,poa_volatility,poa_mapping,camx_fcrs,oxygen_metals,MECH_BASIS,RUN_TYPE,AQM,TOX_IN,PRO_OUT)
            ### Formats and adds header to gspro file
            gspro.format_and_header(tbl_tox,TOX_IN,MECH_BASIS,RUN_TYPE,AQM,MW_FILE,FCRS_FILE,TOX_FILE,PRO_OUT)
            print('NOTICE: Finished generating GSPRO file.')
        elif RUN_TYPE=='NOINTEGRATE': sys.exit('NOINTEGRATE is not an acceptable RUN_TYPE when OUTPUT is PM.')
        else: sys.exit('There is an issue with your RUN_TYPE entry. Only CRITERIA, INTEGRATE, and NOINTEGRATE are accepted.')
    else: sys.exit('There is an issue with your AQM and/or MECH_BASIS entry. Note that MECH_BASIS = PM-AE6 is only available when AQM = CAMX.')
else: sys.exit('There is an issue with your OUTPUT entry. Only VOC and PM are accepted.')

print("Time to run the S2S Tool: ",datetime.now() - startTime)