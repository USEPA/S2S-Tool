import sys
import numpy as np
import pandas as pd
from datetime import date

today     = date.today()

####################################################################################################
### This function generates a gspro file for the target MECH_BASIS.
####################################################################################################

####################################################################################################
def gen_gspro_voc(profiles,species,species_props,carbons,mech4import,tbl_tox,MECH_BASIS,RUN_TYPE,TOLERANCE,TOX_IN,PRO_OUT):
    
    ### gscnv file columns
    column_names  = ['PROFILE','INPUT.POLL','MODEL.SPECIES','MASS.FRACTION','MOLECULAR.WGHT','MASS.FRACTION1']
    dfgspro       = pd.DataFrame(columns=column_names)

    for i in range(len(profiles)):
        prof      = profiles.loc[i,'PROFILE_CODE'] # target profile for this iteration
        temp_spec = species.loc[species['PROFILE_CODE'] == prof] # pull target speciation profile
        temp_spec = temp_spec.reset_index(drop=True) # reset index
        
        if temp_spec['WEIGHT_PERCENT'].sum() / 100 < (1 - TOLERANCE) or temp_spec['WEIGHT_PERCENT'].sum() / 100 > (1 + TOLERANCE):
            print('Profile '+prof+' is outside the specified TOLERANCE and was not processed.')
            continue
        
        if RUN_TYPE=='CRITERIA':
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / temp_spec['WEIGHT_PERCENT'].sum()  # Renormalize Wght%
            i_poll    = 'TOG'
            METHANE   = pd.Series(data={'SPECIES_ID':529})
            nmog_perc = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(METHANE)] # pull target speciation profile minus METHANE
            nmog_perc = nmog_perc['WEIGHT_PERCENT'].sum()
        elif RUN_TYPE=='NOINTEGRATE':
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / temp_spec['WEIGHT_PERCENT'].sum()  # Renormalize Wght%
            TOXLIST   = tbl_tox.loc[:,'SPECIES_ID'] # pull integrated SPECIES_ID list
            TOXLIST   = TOXLIST.drop_duplicates() # drop duplicates
            TOXLIST   = TOXLIST.reset_index(drop=True) # reset index
            temp_spec = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(TOXLIST)] # pull target speciation profile minus TOXLIST
            i_poll    = 'TOG'
            METHANE   = pd.Series(data={'SPECIES_ID':529})
            nmog_perc = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(METHANE)] # pull target speciation profile minus METHANE
            nmog_perc = nmog_perc['WEIGHT_PERCENT'].sum()
        else: # if RUN_TYPE=='INTEGRATE':
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / temp_spec['WEIGHT_PERCENT'].sum()  # Renormalize Wght%
            TOXLIST   = tbl_tox.loc[:,'SPECIES_ID'] # pull integrated SPECIES_ID list
            TOXLIST   = TOXLIST.drop_duplicates() # drop duplicates
            TOXLIST   = TOXLIST.reset_index(drop=True) # reset index
            temp_spec = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(TOXLIST)] # pull target speciation profile minus TOXLIST
            i_poll    = TOX_IN
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / temp_spec['WEIGHT_PERCENT'].sum()  # Renormalize Wght%
            METHANE   = pd.Series(data={'SPECIES_ID':529})
            nmog_perc = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(METHANE)] # pull target speciation profile minus METHANE
            nmog_perc = nmog_perc['WEIGHT_PERCENT'].sum()

        temp_spec = temp_spec.loc[temp_spec['WEIGHT_PERCENT'] > 0.0 ] # Remove species where WEIGHT_PERCENT == 0.0

        if temp_spec.empty:
            print('Profile '+prof+' is empty following pollutant integration and therefore not processed.')
            continue

        spec_list = temp_spec.loc[:,'SPECIES_ID'] # pull unique species in profile
        temp_m4i  = mech4import.loc[mech4import['SPECIES_ID'].isin(spec_list)] # filter mech4import for spec_list
        temp_m4i  = temp_m4i.reset_index(drop=True) # reset index
        
        temp_m4i  = pd.merge(temp_m4i,species_props[['SPECIES_ID','SPEC_MW']],on='SPECIES_ID',how ='left') # append MolWght
        temp_m4i  = pd.merge(temp_m4i,temp_spec[['SPECIES_ID','WEIGHT_PERCENT']],on='SPECIES_ID',how ='left') # append Wght %
        temp_m4i  = pd.merge(temp_m4i,carbons[['Species','nC']],on='Species',how ='left') # append nC per 1 model species

        temp_m4i.loc[:,'totC'] = temp_m4i.loc[:,'Moles'] * temp_m4i.loc[:,'nC'] # calculate total nC per model species
        temp_m4i.loc[:,'moleSpec_massVOC'] = temp_m4i.loc[:,'WEIGHT_PERCENT'] * temp_m4i.loc[:,'Moles'] / \
                                             temp_m4i.loc[:,'SPEC_MW'] # calculate mole model species / mass VOC

        nC        = temp_m4i.groupby('SPECIES_ID',as_index=False)['totC'].sum() # calculate total nC for each compound
        molesplit = temp_m4i.groupby('Species',as_index=False)['moleSpec_massVOC'].sum() # calculate total mole model species / mass VOC
        
        temp_m4i  = pd.merge(temp_m4i,nC[['SPECIES_ID','totC']],on='SPECIES_ID',how ='left',suffixes=('_1','_2')) # append total nC for each compound
        
        temp_m4i.loc[:,'modSpec_MW'] = temp_m4i.loc[:,'SPEC_MW'] * temp_m4i.loc[:,'nC'] / temp_m4i.loc[:,'totC_2'] # calculate model species grams per mole
        temp_m4i['moleSpec_massVOC'] /= temp_m4i.groupby('Species')['moleSpec_massVOC'].transform(sum) # normalize mole model species / mass VOC
        temp_m4i.loc[:,'frac_modSpec_MW'] = temp_m4i.loc[:,'modSpec_MW'] * temp_m4i.loc[:,'moleSpec_massVOC'] # calculate fraction of model species grams per mole

        final_MW  = temp_m4i.groupby('Species',as_index=False)['frac_modSpec_MW'].sum() # calculate final MW for each model species
        molesplit = pd.merge(molesplit,final_MW[['Species','frac_modSpec_MW']],on='Species',how ='left') # append final MW for each model species

        molesplit.loc[:,'MASS.FRACTION']  = molesplit.loc[:,'moleSpec_massVOC'] * molesplit.loc[:,'frac_modSpec_MW'] # calculate final split factors
        molesplit.loc[:,'MASS.FRACTION1'] = molesplit.loc[:,'moleSpec_massVOC'] * molesplit.loc[:,'frac_modSpec_MW'] # calculate final split factors
        molesplit.loc[:,'PROFILE']        = prof # add profile to dataframe
        molesplit.loc[:,'INPUT.POLL']     = i_poll # add input pollutant to dataframe
        molesplit = molesplit.rename(columns={'Species': 'MODEL.SPECIES'}) # rename Species column
        molesplit = molesplit.rename(columns={'frac_modSpec_MW': 'MOLECULAR.WGHT'}) # rename frac_modSpec_MW column
        molesplit = molesplit.drop(['moleSpec_massVOC'],axis=1) # drop moleSpec_massVOC column
        molesplit = molesplit[column_names] # reorder columns
        
        if nmog_perc == 0.0:
            pass
        else:
            nmog = pd.Series(data={'PROFILE':prof,'INPUT.POLL':i_poll,'MODEL.SPECIES':'NMOG','MASS.FRACTION':nmog_perc,\
                                   'MOLECULAR.WGHT':1,'MASS.FRACTION1':nmog_perc})
            molesplit = molesplit.append(nmog,ignore_index=True)

        dfgspro   = dfgspro.append(molesplit) # append profile gspro to final gspro

    dfgspro['MASS.FRACTION']  = dfgspro['MASS.FRACTION'].astype(float).apply(lambda x: '%.6E' % x)
    dfgspro['MOLECULAR.WGHT'] = dfgspro['MOLECULAR.WGHT'].astype(float).apply(lambda x: '%.6E' % x)
    dfgspro['MASS.FRACTION1'] = dfgspro['MASS.FRACTION1'].astype(float).apply(lambda x: '%.6E' % x)

    ### Output gspro df to file
    dfgspro.to_csv(PRO_OUT,index=False,header=False)
####################################################################################################

####################################################################################################
def gen_gspro_pm(profiles,species,species_props,mechPM,tbl_tox,poa_volatility,poa_mapping,camx_fcrs,oxygen_metals,MECH_BASIS,RUN_TYPE,AQM,TOX_IN,PRO_OUT):
    
    ### gscnv file columns
    column_names  = ['PROFILE','INPUT.POLL','MODEL.SPECIES','MASS.FRACTION','MOLECULAR.WGHT','MASS.FRACTION1']
    dfgspro       = pd.DataFrame(columns=column_names)
    
    for i in range(len(profiles)):
        prof      = profiles.loc[i,'PROFILE_CODE'] # target profile for this iteration
        temp_spec = species.loc[species['PROFILE_CODE'] == prof] # pull target speciation profile
        temp_spec = temp_spec.reset_index(drop=True) # reset index
        
    ### if PROFILE_TYPE == PM, translate data into a "PM-ready" profile.
        if profiles.loc[i,'PROFILE_TYPE'] == 'PM': # if PROFILE_TYPE == PM, translate data into a "PM-ready" profile.
            #continue
            temp_spec_pm  = pd.DataFrame(columns=['PROFILE_CODE','SPECIES_ID']) # setup temp_spec_pm array
            temp_spec_pm['SPECIES_ID']   = [626,797,699,700,613,784,2669,488,292,694,715,2303,329,2772,525,2302,669,526,785,696,337,795,\
                                            2668,2671,666,767,347,379,612,380,778,468,298,693,689,697,779,586,649,695,328,487,\
                                            714,296,300,519,1861,528,520,2670] # hardcode "PM-ready" species
            temp_spec_pm['PROFILE_CODE'] = prof # assign PROFILE_CODE
            temp_spec_pm  = pd.merge(temp_spec_pm,temp_spec[['SPECIES_ID','WEIGHT_PERCENT']],on='SPECIES_ID',how ='left') # append WEIGHT_PERCENT
            temp_spec_pm['WEIGHT_PERCENT'] = temp_spec_pm['WEIGHT_PERCENT'].fillna(0) # replace all NaN with zeros
            if temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2303, 'WEIGHT_PERCENT'].iloc[0] > 0:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'] = 0.0
            if temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2772, 'WEIGHT_PERCENT'].iloc[0] > 0:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'] = 0.0
            if temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2302, 'WEIGHT_PERCENT'].iloc[0] > 0:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'] = 0.0
            if temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 785, 'WEIGHT_PERCENT'].iloc[0] > 0:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'] = 0.0
            if temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 337, 'WEIGHT_PERCENT'].iloc[0] > 0:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 795, 'WEIGHT_PERCENT'] = 0.0
            # Add particulate water, PH2O
            if profiles.loc[i,'CATEGORY_LEVEL_1_Generation_Mechanism'] == 'Combustion' or profiles.loc[i,'CATEGORY_LEVEL_1_Generation_Mechanism'] == 'Miscellaneous':
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2668, 'WEIGHT_PERCENT'] = (temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 699, 'WEIGHT_PERCENT'].iloc[0] + \
                                                                                          temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 784, 'WEIGHT_PERCENT'].iloc[0]) * 0.24
            else:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2668, 'WEIGHT_PERCENT'] = 0.0
            # Add PSO4 from sulfur, if necessary
            if temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 699, 'WEIGHT_PERCENT'].iloc[0] > 0:
                pass
            elif temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 700, 'WEIGHT_PERCENT'].iloc[0] == 0:
                pass
            else:
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 699, 'WEIGHT_PERCENT'] = temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 700, 'WEIGHT_PERCENT'].iloc[0] * (96/32)
            temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 700, 'WEIGHT_PERCENT'] = 0
            # Add MO, metal bound oxygen, if necessary
            oxygen_metals_temp  = pd.merge(oxygen_metals,temp_spec_pm[['SPECIES_ID','WEIGHT_PERCENT']],on='SPECIES_ID',how ='left') # append WEIGHT_PERCENT
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'].iloc[0] == 0 or oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2303, 'WEIGHT_PERCENT'].iloc[0] == 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'] = 0.0
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2303, 'WEIGHT_PERCENT'] = 0.0
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'].iloc[0] == 0 or oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2302, 'WEIGHT_PERCENT'].iloc[0] == 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'] = 0.0
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2302, 'WEIGHT_PERCENT'] = 0.0
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'].iloc[0] == 0 or oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2772, 'WEIGHT_PERCENT'].iloc[0] == 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'] = 0.0
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2772, 'WEIGHT_PERCENT'] = 0.0
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0] == 0 or oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 785, 'WEIGHT_PERCENT'].iloc[0] == 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'] = 0.0
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 785, 'WEIGHT_PERCENT'] = 0.0
            oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'] = oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'].iloc[0] - \
                                                                                      oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2303, 'WEIGHT_PERCENT'].iloc[0]
            oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'] = oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'].iloc[0] - \
                                                                                      oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2302, 'WEIGHT_PERCENT'].iloc[0]
            oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'] = oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'].iloc[0] - \
                                                                                      oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 2772, 'WEIGHT_PERCENT'].iloc[0]
            oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'] = oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0] - \
                                                                                      oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 785, 'WEIGHT_PERCENT'].iloc[0]
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'].iloc[0] <= 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'] = 0.0
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'].iloc[0] <= 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'] = 0.0
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'].iloc[0] <= 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'] = 0.0
            if oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0] <= 0:
                oxygen_metals_temp.loc[oxygen_metals_temp['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'] = 0.0
            oxygen_metals_temp['MOunadjusted'] = oxygen_metals_temp['oxy/metal_ratio'] * oxygen_metals_temp['WEIGHT_PERCENT'] # Calculate MOunadjusted                
            MOunadjusted = oxygen_metals_temp['MOunadjusted'].sum()
            NeutralizedSO4 = (0.5 * 96 / 18) * temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 784, 'WEIGHT_PERCENT'].iloc[0]
            NonNeutralizedSO4 = temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 699, 'WEIGHT_PERCENT'].iloc[0] - NeutralizedSO4
            if NonNeutralizedSO4 <= 0:
                MOadjusted = MOunadjusted
            else:
                MOadjusted = MOunadjusted - NonNeutralizedSO4 * (16 / 96)
            if MOadjusted <= 0:
                MOadjusted = 0.0
            else: pass
            temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2670, 'WEIGHT_PERCENT'] = MOadjusted
            # Add PNCOM, if necessary
            if 'mobile' in profiles.loc[i,'CATEGORY_LEVEL_2_Sector_Equipment'].lower():
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2669, 'WEIGHT_PERCENT'] = temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 626, 'WEIGHT_PERCENT'].iloc[0] * 0.25
            elif 'biomass burning' in profiles.loc[i,'CATEGORY_LEVEL_2_Sector_Equipment'].lower():
                temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 2669, 'WEIGHT_PERCENT'] = temp_spec_pm.loc[temp_spec_pm['SPECIES_ID'] == 626, 'WEIGHT_PERCENT'].iloc[0] * 0.7
            else:
                pass
            temp_spec_pm  = temp_spec_pm.loc[temp_spec_pm['WEIGHT_PERCENT'] > 0.0 ] # Remove species where WEIGHT_PERCENT == 0.0
            # Renormalize, if necessary
            if temp_spec_pm['WEIGHT_PERCENT'].sum() > 101.0:
                temp_spec_pm['WEIGHT_PERCENT'] = temp_spec_pm['WEIGHT_PERCENT'] / temp_spec_pm['WEIGHT_PERCENT'].sum() * 100 # Renormalize Wght%
            else:
                pass
            temp_spec = temp_spec_pm

        if RUN_TYPE=='CRITERIA':
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / 100  # Sum of Wght% = 1
            i_poll    = 'PM2_5'
        else: # if RUN_TYPE=='INTEGRATE':
            TOXLIST   = tbl_tox.loc[:,'SPECIES_ID'] # pull integrated SPECIES_ID list
            TOXLIST   = TOXLIST.drop_duplicates() # drop duplicates
            TOXLIST   = TOXLIST.reset_index(drop=True) # reset index
            temp_spec = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(TOXLIST)] # pull target speciation profile minus TOXLIST
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / temp_spec['WEIGHT_PERCENT'].sum()  # Renormalize Wght%
            i_poll    = TOX_IN
        
        temp_mech = pd.merge(mechPM,temp_spec[['SPECIES_ID','WEIGHT_PERCENT']],on='SPECIES_ID',how ='left') # append WEIGHT_PERCENT
        temp_mech['WEIGHT_PERCENT'] = temp_mech['WEIGHT_PERCENT'].fillna(0) # replace all NaN with zeros

        poa_mech  = pd.merge(poa_mapping,temp_spec[['SPECIES_ID','WEIGHT_PERCENT']],on='SPECIES_ID',how ='left') # append WEIGHT_PERCENT
        poa_mech['WEIGHT_PERCENT'] = poa_mech['WEIGHT_PERCENT'].fillna(0) # replace all NaN with zeros
        poa_mech  = poa_mech.loc[poa_mech['WEIGHT_PERCENT'] > 0.0 ] # Remove species where WEIGHT_PERCENT == 0.0
        poa_mech  = poa_mech.reset_index(drop=True) # reset index
        
        ### Extract SV-POA profile assigned to source type
        temp_poa  = poa_volatility.loc[(poa_volatility['CATEGORY_LEVEL_1_Generation_Mechanism']==profiles.loc[i,'CATEGORY_LEVEL_1_Generation_Mechanism']) & \
                                       (poa_volatility['CATEGORY_LEVEL_2_Sector_Equipment']==profiles.loc[i,'CATEGORY_LEVEL_2_Sector_Equipment'])]
        temp_poa  = temp_poa.reset_index(drop=True) # reset index

        ### Process OM/OC/NCOM from temp_spec for appropriate MECH_BASIS
        if MECH_BASIS == 'PM-AE6':
            if profiles.loc[i,'PROFILE_TYPE'] == 'PM-AE6' or profiles.loc[i,'PROFILE_TYPE'] == 'PM':
                pass # OC/PNCOM already appropriately mapped in base profile.
            
            elif profiles.loc[i,'PROFILE_TYPE'] == 'PM-AE8':
                temp_mech.loc[temp_mech['Species'] == 'POC', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC', 'WEIGHT_PERCENT'].sum() # add POC
                temp_mech.loc[temp_mech['Species'] == 'PNCOM', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum() # add PNCOM
            
            elif profiles.loc[i,'PROFILE_TYPE'] == 'PM-CR1':
                temp_mech.loc[temp_mech['Species'] == 'POC', 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * 1 / \
                                                                                 profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'] # add POC
                temp_mech.loc[temp_mech['Species'] == 'PNCOM', 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * (1 - 1 / \
                                                                                   profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO']) # add PNCOM
            else: sys.exit('PROFILE_TYPE = '+profiles.loc[i,'PROFILE_TYPE']+' for profile '+prof+' is not recognized.')

        elif MECH_BASIS == 'PM-AE8':
            if profiles.loc[i,'PROFILE_TYPE'] == 'PM-AE6' or profiles.loc[i,'PROFILE_TYPE'] == 'PM':
                temp_mech.loc[temp_mech['Species'] == 'POCN2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'-2'][0] # add POCN2
                temp_mech.loc[temp_mech['Species'] == 'POCN1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'-1'][0] # add POCN1
                temp_mech.loc[temp_mech['Species'] == 'POCP0', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'0'][0] # add POCP0
                temp_mech.loc[temp_mech['Species'] == 'POCP1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'1'][0] # add POCP1
                temp_mech.loc[temp_mech['Species'] == 'POCP2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'2'][0] # add POCP2
                temp_mech.loc[temp_mech['Species'] == 'PNCOMN2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'-2'][0] # add PNCOMN2
                temp_mech.loc[temp_mech['Species'] == 'PNCOMN1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'-1'][0] # add PNCOMN1
                temp_mech.loc[temp_mech['Species'] == 'PNCOMP0', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'0'][0] # add PNCOMP0
                temp_mech.loc[temp_mech['Species'] == 'PNCOMP1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'1'][0] # add PNCOMP1
                temp_mech.loc[temp_mech['Species'] == 'PNCOMP2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum() * \
                                                                                   temp_poa.loc[:,'2'][0] # add PNCOMP2

            elif profiles.loc[i,'PROFILE_TYPE'] == 'PM-AE8':
                pass # SV-POA already appropriately mapped in base profile.

            elif profiles.loc[i,'PROFILE_TYPE'] == 'PM-CR1':
                temp_mech.loc[temp_mech['Species'] == 'POCN2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == -2, 'WEIGHT_PERCENT'].sum() * \
                                                                                   1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'] # add POCN2
                temp_mech.loc[temp_mech['Species'] == 'PNCOMN2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == -2, 'WEIGHT_PERCENT'].sum() * \
                                                                                     (1 - 1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO']) # add PNCOMN2
                temp_mech.loc[temp_mech['Species'] == 'POCN1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == -1, 'WEIGHT_PERCENT'].sum() * \
                                                                                   1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'] # add POCN1
                temp_mech.loc[temp_mech['Species'] == 'PNCOMN1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == -1, 'WEIGHT_PERCENT'].sum() * \
                                                                                     (1 - 1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO']) # add PNCOMN1
                temp_mech.loc[temp_mech['Species'] == 'POCP0', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 0, 'WEIGHT_PERCENT'].sum() * \
                                                                                   1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'] # add POCP0
                temp_mech.loc[temp_mech['Species'] == 'PNCOMP0', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 0, 'WEIGHT_PERCENT'].sum() * \
                                                                                     (1 - 1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO']) # add PNCOMP0
                temp_mech.loc[temp_mech['Species'] == 'POCP1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 1, 'WEIGHT_PERCENT'].sum() * \
                                                                                   1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'] # add POCP1
                temp_mech.loc[temp_mech['Species'] == 'PNCOMP1', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 1, 'WEIGHT_PERCENT'].sum() * \
                                                                                     (1 - 1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO']) # add PNCOMP1
                temp_mech.loc[temp_mech['Species'] == 'POCP2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 2, 'WEIGHT_PERCENT'].sum() * \
                                                                                   1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'] # add POCP2
                temp_mech.loc[temp_mech['Species'] == 'PNCOMP2', 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 2, 'WEIGHT_PERCENT'].sum() * \
                                                                                     (1 - 1 / profiles.loc[i,'ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO']) # add PNCOMP2

            else: sys.exit('PROFILE_TYPE = '+profiles.loc[i,'PROFILE_TYPE']+' for profile '+prof+' is not recognized.')

        elif MECH_BASIS == 'PM-CR1':
            if profiles.loc[i,'PROFILE_TYPE'] == 'PM-AE6' or profiles.loc[i,'PROFILE_TYPE'] == 'PM':
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3394., 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * temp_poa.loc[:,'-2'][0] # add POCN2
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3395., 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * temp_poa.loc[:,'-1'][0] # add POCN1
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3396., 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * temp_poa.loc[:,'0'][0] # add POCP0
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3397., 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * temp_poa.loc[:,'1'][0] # add POCP1
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3398., 'WEIGHT_PERCENT'] = poa_mech.loc[:,'WEIGHT_PERCENT'].sum() * temp_poa.loc[:,'2'][0] # add POCP2

            elif profiles.loc[i,'PROFILE_TYPE'] == 'PM-AE8':
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3394., 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == -2, 'WEIGHT_PERCENT'].sum() # add POCN2
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3395., 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == -1, 'WEIGHT_PERCENT'].sum() # add POCN1
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3396., 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 0, 'WEIGHT_PERCENT'].sum() # add POCP0
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3397., 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 1, 'WEIGHT_PERCENT'].sum() # add POCP1
                temp_mech.loc[temp_mech['SPECIES_ID'] == 3398., 'WEIGHT_PERCENT'] = poa_mech.loc[poa_mech['log10Cstar'] == 2, 'WEIGHT_PERCENT'].sum() # add POCP2
                
            elif profiles.loc[i,'PROFILE_TYPE'] == 'PM-CR1':
                pass # SV-POA already appropriately mapped in base profile.
            
            else: sys.exit('PROFILE_TYPE = '+profiles.loc[i,'PROFILE_TYPE']+' for profile '+prof+' is not recognized.')
        
        else: sys.exit('MECH_BASIS is not recognized.')
        
        ### Chlorine: If 337 is absent, but 795 is present, use 795
        if temp_mech.loc[temp_mech['SPECIES_ID'] == 337, 'WEIGHT_PERCENT'].iloc[0] == 0.0 and 795 in temp_spec['SPECIES_ID'].values:
            temp_mech.loc[temp_mech['Species'] == 'PCL', 'WEIGHT_PERCENT'] = temp_spec.loc[temp_spec['SPECIES_ID'] == 795, 'WEIGHT_PERCENT'].iloc[0]
        else: pass

        ### Calcium: If 2303 is absent, but 329 is present, use 329  
        if AQM=='CMAQ':
            if temp_mech.loc[temp_mech['SPECIES_ID'] == 2303, 'WEIGHT_PERCENT'].iloc[0] == 0.0 and 329 in temp_spec['SPECIES_ID'].values:
                temp_mech.loc[temp_mech['Species'] == 'PCA', 'WEIGHT_PERCENT'] = temp_spec.loc[temp_spec['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'].iloc[0]
            else: pass
        else: pass
		
        ### Magnesium: If 2772 is absent, but 525 is present, use 525
        if AQM=='CMAQ':
            if temp_mech.loc[temp_mech['SPECIES_ID'] == 2772, 'WEIGHT_PERCENT'].iloc[0] == 0.0 and 525 in temp_spec['SPECIES_ID'].values:
                temp_mech.loc[temp_mech['Species'] == 'PMG', 'WEIGHT_PERCENT'] = temp_spec.loc[temp_spec['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'].iloc[0]
            else: pass
        else: pass
		
        ### Potassium: If 2302 is absent, but 669 is present, use 669
        if AQM=='CMAQ':
            if temp_mech.loc[temp_mech['SPECIES_ID'] == 2302, 'WEIGHT_PERCENT'].iloc[0] == 0.0 and 669 in temp_spec['SPECIES_ID'].values:
                temp_mech.loc[temp_mech['Species'] == 'PK', 'WEIGHT_PERCENT'] = temp_spec.loc[temp_spec['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'].iloc[0]
            else: pass
        else: pass
		
        ### Sodium: If 785 is absent, but 696 is present, use 696
        if AQM=='CMAQ':
            if temp_mech.loc[temp_mech['SPECIES_ID'] == 785, 'WEIGHT_PERCENT'].iloc[0] == 0.0 and 696 in temp_spec['SPECIES_ID'].values:
                temp_mech.loc[temp_mech['Species'] == 'PNA', 'WEIGHT_PERCENT'] = temp_spec.loc[temp_spec['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0]
            else: pass
        else: # if AQM=='CAMX'
            temp_mech.loc[temp_mech['SPECIES_ID'] == 785., 'Species'] = 'NA'
            if temp_mech.loc[temp_mech['SPECIES_ID'] == 785, 'WEIGHT_PERCENT'].iloc[0] == 0.0 and 696 in temp_spec['SPECIES_ID'].values:
                temp_mech.loc[temp_mech['Species'] == 'NA', 'WEIGHT_PERCENT'] = temp_spec.loc[temp_spec['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0]
            else: pass
        
        ### Calculate PMOTHR or FPRM/FCRS for CMAQ and CAMX, respectively
        if temp_mech['WEIGHT_PERCENT'].sum() < 1.0:
            if AQM=='CMAQ':
                if TOX_IN=='TOM':
                    pass
                else:
                    temp_mech.loc[temp_mech['Species'] == 'PMOTHR', 'WEIGHT_PERCENT'] = 1.0 - temp_mech['WEIGHT_PERCENT'].sum()
            else: # if AQM=='CAMX'
                if TOX_IN=='TOM':
                    pass
                else:
                    temp_fcrs = camx_fcrs.loc[:,'PROFILE_CODE'] # pull integrated SPECIES_ID list
                    temp_fcrs = camx_fcrs.loc[camx_fcrs['PROFILE_CODE'].astype(str) == prof] # pull prof from camx_fcrs, if available
                    if temp_fcrs.empty:
                        temp_mech.loc[temp_mech['Species'] == 'FPRM', 'WEIGHT_PERCENT'] = 1.0 - temp_mech['WEIGHT_PERCENT'].sum()
                    else:
                        temp_mech.loc[temp_mech['Species'] == 'FCRS', 'WEIGHT_PERCENT'] = 1.0 - temp_mech['WEIGHT_PERCENT'].sum()
        else: pass

        # Remove added species, if necessary, and renormalize if RUN_TYPE=='INTEGRATE'
        if RUN_TYPE=='INTEGRATE':
            TOXLIST   = tbl_tox.loc[:,'SPECIES_ID'] # pull integrated SPECIES_ID list
            temp_mech = temp_mech.loc[~temp_mech['SPECIES_ID'].isin(TOXLIST)] # remove TOXLIST species that might have been added
            temp_mech.loc[:,'WEIGHT_PERCENT'] = temp_mech.loc[:,'WEIGHT_PERCENT'] / temp_mech.loc[:,'WEIGHT_PERCENT'].sum()  # Renormalize Wght%
        else: # if RUN_TYPE=='CRITERIA':
            pass
        
        ### Calculate POA for CAMX and drop PNCOM
        temp_mech = temp_mech.reset_index(drop=True) # reset index
        if AQM=='CAMX':
            poa_row   = {'AQM': AQM, 'Mechanism': temp_mech.loc[0,'Mechanism'], 'SPECIES_ID': 9999, 'Species': 'POA', \
                         'WEIGHT_PERCENT':temp_mech.loc[temp_mech['Species'] == 'POC', 'WEIGHT_PERCENT'].sum() + \
                         temp_mech.loc[temp_mech['Species'] == 'PNCOM', 'WEIGHT_PERCENT'].sum()}
            temp_mech = temp_mech.append(poa_row, ignore_index = True)
            temp_mech.drop(temp_mech.index[temp_mech['Species'] == 'PNCOM'], inplace=True)
        else: pass

        temp_mech = temp_mech.groupby('Species').sum().reset_index() # Sum model species with multiple WEIGHT_PERCENT
        temp_mech = temp_mech.loc[temp_mech['WEIGHT_PERCENT'] > 0.0 ] # Remove species where WEIGHT_PERCENT == 0.0

        if temp_mech.empty:
            print('Profile '+prof+' is empty following pollutant integration and therefore not processed.')
            continue
        
        temp_mech = temp_mech.rename(columns={'Species': 'MODEL.SPECIES'}) # rename Species column
        temp_mech = temp_mech.rename(columns={'WEIGHT_PERCENT': 'MASS.FRACTION'}) # rename Species column
        temp_mech = temp_mech.drop(['SPECIES_ID'],axis=1) # drop SPECIES_ID column
        temp_mech.loc[:,'MASS.FRACTION1'] = temp_mech.loc[:,'MASS.FRACTION'] # calculate final split factors
        temp_mech.loc[:,'MOLECULAR.WGHT'] = 1.0
        temp_mech.loc[:,'PROFILE'] = prof
        temp_mech.loc[:,'INPUT.POLL'] = i_poll
        temp_mech = temp_mech[column_names] # reorder columns
        temp_mech = temp_mech.reset_index(drop=True) # reset index

        dfgspro   = dfgspro.append(temp_mech) # append profile gspro to final gspro
        
    dfgspro['MASS.FRACTION']  = dfgspro['MASS.FRACTION'].astype(float).apply(lambda x: '%.6E' % x)
    dfgspro['MOLECULAR.WGHT'] = dfgspro['MOLECULAR.WGHT'].astype(float).apply(lambda x: '%.6E' % x)
    dfgspro['MASS.FRACTION1'] = dfgspro['MASS.FRACTION1'].astype(float).apply(lambda x: '%.6E' % x)

    ### Output gspro df to file
    dfgspro.to_csv(PRO_OUT,index=False,header=False)
####################################################################################################

####################################################################################################
def format_and_header(tbl_tox,TOX_IN,MECH_BASIS,RUN_TYPE,AQM,CAR_FILE,FCRS_FILE,TOX_FILE,PRO_OUT):

    ################################################################################################
    ### Import gscnv csv file.
    f1_gspro = np.genfromtxt(PRO_OUT,delimiter=',',dtype='str')
    ################################################################################################
    
    ################################################################################################
    headerline1   = '#S2S_AQM             '+AQM
    headerline2   = '#S2S_CAMX_FCRS       '+FCRS_FILE
    headerline3   = '#S2S_CARBONS         '+CAR_FILE
    headerline4   = '#S2S_MECH_BASIS      '+MECH_BASIS
    headerline5   = '#S2S_RUN_TYPE        '+RUN_TYPE
    headerline6   = '#S2S_RUN_DATE        '+str(today)
    if RUN_TYPE == 'CRITERIA':
        headerline7   = '#S2S_TBL_TOX         Not Applicable'
    else: # RUN_TYPE == 'INTEGRATE' or RUN_TYPE == 'NOINTEGRATE':
        headerline7   = '#S2S_TBL_TOX         '+TOX_FILE
    headerline    = '\n'.join([headerline1,headerline2,headerline3,headerline4,headerline5,headerline6,headerline7])
    if RUN_TYPE=='INTEGRATE':
        temp_tox  = tbl_tox.copy()
        temp_tox.loc[:,'AQM'] = '#NHAP'
        temp_tox.loc[:,'SPECIES_ID'] = TOX_IN
        temp_tox = temp_tox.to_string(header=None,index=False,justify='left')
        headerline = '\n'.join([headerline,temp_tox])
    else: pass # if RUN_TYPE=='CRITERIA' or RUN_TYPE=='NOINTEGRATE'

    np.savetxt(PRO_OUT,f1_gspro[:],fmt='%-20s %-20s %-10s %-13s %-13s %-13s',header=headerline,comments='')
    ################################################################################################

####################################################################################################