import numpy as np
import pandas as pd
from datetime import date

today     = date.today()

####################################################################################################
### This function generates a gscnv file for the target MECH_BASIS.
####################################################################################################

####################################################################################################
def gen_gscnv(profiles,species,species_props,tbl_tox,gscnv_append,MECH_BASIS,RUN_TYPE,TOLERANCE,CNV_OUT):
    
    ### gscnv file columns
    column_names  = ['INPUT.POLL','OUTPUT.POLL','PROFILE','OUTPUT.MASS/INPUT.MASS']
    dfgscnv       = pd.DataFrame(columns=column_names)

    if RUN_TYPE=='CRITERIA' or RUN_TYPE=='NOINTEGRATE':
        for i in range(len(profiles)):
            i_poll    = 'VOC'
            o_poll    = 'TOG'
            prof      = profiles.iloc[i,0]
            temp_spec = species.loc[species['PROFILE_CODE'] == prof] # pull target speciation profile

            # Check if sum of WEIGHT_PERCENT is outside specified TOLERANCE. If YES, continue. Printed comment provided in gspro module.            
            if temp_spec.loc[:,'WEIGHT_PERCENT'].sum() / 100 < (1 - TOLERANCE) or temp_spec.loc[:,'WEIGHT_PERCENT'].sum() / 100 > (1 + TOLERANCE):
                continue

            temp_spec = pd.merge(temp_spec,species_props[['SPECIES_ID','NonVOCTOG']],on='SPECIES_ID',how ='left')

            if temp_spec.loc[temp_spec['NonVOCTOG'] == 0, 'WEIGHT_PERCENT'].sum() == 0.:
                ratio = 0.0
            else:
                ratio = 100/temp_spec.loc[temp_spec['NonVOCTOG'] == 0, 'WEIGHT_PERCENT'].sum()
            gscnv_row = pd.Series(data={'INPUT.POLL':i_poll,'OUTPUT.POLL':o_poll,
                                        'PROFILE':prof,'OUTPUT.MASS/INPUT.MASS':ratio})
            dfgscnv   = dfgscnv.append(gscnv_row,ignore_index=True)

        if RUN_TYPE=='CRITERIA':
            gscnv_append = pd.merge(gscnv_append,dfgscnv[['PROFILE','OUTPUT.MASS/INPUT.MASS']],on='PROFILE',how ='left') # append TOG/VOC ratio
            dfgscnv   = dfgscnv.append(gscnv_append,ignore_index=True)
        else: pass

    else: # if RUN_TYPE=='INTEGRATE'
        for i in range(len(profiles)):
            i_poll    = 'NONHAPVOC'
            o_poll    = 'NONHAPTOG'
            prof      = profiles.iloc[i,0]
            temp_spec = species.loc[species['PROFILE_CODE'] == prof] # pull target speciation profile
            
            # Check if sum of WEIGHT_PERCENT is outside specified TOLERANCE. If YES, continue. Printed comment provided in gspro module.
            if temp_spec.loc[:,'WEIGHT_PERCENT'].sum() / 100 < (1 - TOLERANCE) or temp_spec.loc[:,'WEIGHT_PERCENT'].sum() / 100 > (1 + TOLERANCE):
                continue

            TOXLIST   = tbl_tox.loc[:,'SPECIES_ID'] # pull integrated SPECIES_ID list
            TOXLIST   = TOXLIST.drop_duplicates() # drop duplicates
            TOXLIST   = TOXLIST.reset_index(drop=True) # reset index
            temp_spec = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(TOXLIST)] # pull target speciation profile minus TOXLIST
            temp_spec = pd.merge(temp_spec,species_props[['SPECIES_ID','NonVOCTOG']],on='SPECIES_ID',how ='left')
            if temp_spec.loc[temp_spec['NonVOCTOG'] == 0, 'WEIGHT_PERCENT'].sum() == 0.:
                ratio = 0.0
            else:
                ratio = temp_spec['WEIGHT_PERCENT'].sum()/temp_spec.loc[temp_spec['NonVOCTOG'] == 0, 'WEIGHT_PERCENT'].sum()
            gscnv_row = pd.Series(data={'INPUT.POLL':i_poll,'OUTPUT.POLL':o_poll,
                                        'PROFILE':prof,'OUTPUT.MASS/INPUT.MASS':ratio})
            dfgscnv   = dfgscnv.append(gscnv_row,ignore_index=True)
    
    dfgscnv['OUTPUT.MASS/INPUT.MASS'] = dfgscnv['OUTPUT.MASS/INPUT.MASS'].astype(float).apply(lambda x: '%.8F' % x)
    
    ### Output gscnv df to file
    dfgscnv.to_csv(CNV_OUT,index=False,header=False)

####################################################################################################
def format_and_header(MECH_BASIS,RUN_TYPE,AQM,MW_FILE,TOX_FILE,CNV_OUT):

    ################################################################################################
    ### Import gscnv csv file.
    f1_gscnv = np.genfromtxt(CNV_OUT,delimiter=',',dtype='str')
    ################################################################################################
    
    ################################################################################################
    headerline1   = '#S2S_AQM             '+AQM
    headerline2   = '#S2S_CAMX_FCRS       Not Applicable'
    headerline3   = '#S2S_MW              '+MW_FILE
    headerline4   = '#S2S_MECH_BASIS      '+MECH_BASIS
    headerline5   = '#S2S_RUN_TYPE        '+RUN_TYPE
    headerline6   = '#S2S_RUN_DATE        '+str(today)
    if RUN_TYPE == 'CRITERIA':
        headerline7   = '#S2S_TBL_TOX         Not Applicable'
    else: # RUN_TYPE == 'INTEGRATE' or RUN_TYPE == 'NOINTEGRATE':
        headerline7   = '#S2S_TBL_TOX         '+TOX_FILE
    headerline8   = '#BY PROFILE'
    headerline    = '\n'.join([headerline1,headerline2,headerline3,headerline4,headerline5,headerline6,headerline7,headerline8])
    np.savetxt(CNV_OUT,f1_gscnv[:],fmt='%-20s',header=headerline,comments='')
    ################################################################################################

####################################################################################################