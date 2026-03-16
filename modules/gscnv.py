import numpy as np
import pandas as pd
from datetime import date

today     = date.today()

####################################################################################################
### This function generates a gscnv file for the target MECH_BASIS.
####################################################################################################

####################################################################################################
def gen_gscnv(profiles, species, species_props, tbl_tox, gscnv_append, MECH_BASIS, RUN_TYPE, TOLERANCE, CNV_OUT):
    """
    Generate gscnv file for the target MECH_BASIS (VOC output only).
    """
    # Ensure explicit columns are used (don’t rely on positions)
    if 'PROFILE_CODE' not in profiles.columns:
        raise ValueError("profiles must contain 'PROFILE_CODE'")

    # Pre-merge species_props once to get NonVOCTOG
    sp = species.merge(species_props[['SPECIES_ID', 'NonVOCTOG']], on='SPECIES_ID', how='left')

    # Compute total weight per profile for tolerance check
    total_wt = sp.groupby('PROFILE_CODE', as_index=True)['WEIGHT_PERCENT'].sum()

    # Tolerance check: sum should be within [100*(1-TOL), 100*(1+TOL)]
    low, high = 100.0 * (1 - TOLERANCE), 100.0 * (1 + TOLERANCE)
    valid_profiles = total_wt.index[(total_wt >= low) & (total_wt <= high)]

    # Pulls unique profiles that meet tolerance check
    prof_codes = pd.Index(profiles['PROFILE_CODE'].unique(), dtype=valid_profiles.dtype)
    prof_codes = prof_codes.intersection(valid_profiles)

    if RUN_TYPE in {'CRITERIA', 'NOINTEGRATE'}:
        # Denominator: 
        denom = sp.loc[sp['NonVOCTOG'] == 0].groupby('PROFILE_CODE')['WEIGHT_PERCENT'] \
                 .sum().reindex(prof_codes, fill_value=0.0)

        # Ratio: 100 / denom (0 if denom == 0; aka NonVOCTOG.sum == 0)
        ratio = (pd.Series(100.0, index=denom.index).div(denom.replace(0.0, np.nan))
                .fillna(0.0).to_numpy())

        dfgscnv = pd.DataFrame({'INPUT.POLL': 'VOC', 'OUTPUT.POLL': 'TOG',
                                'PROFILE': prof_codes.values, 'OUTPUT.MASS/INPUT.MASS': ratio})

        if RUN_TYPE == 'CRITERIA':
            # Append additional rows by merging ratios into gscnv_append on PROFILE
            dfgscnv = dfgscnv.copy()
            app = gscnv_append.merge(dfgscnv[['PROFILE', 'OUTPUT.MASS/INPUT.MASS']],
                on='PROFILE', how='left')
            dfgscnv = pd.concat([dfgscnv, app], ignore_index=True)

    elif RUN_TYPE == 'INTEGRATE':
        # Remove toxic species
        tox_ids = pd.Index(tbl_tox['SPECIES_ID'].drop_duplicates())
        sp_nontox = sp.loc[~sp['SPECIES_ID'].isin(tox_ids)]

        # Numerator: total remaining weight per profile
        num = sp_nontox.groupby('PROFILE_CODE')['WEIGHT_PERCENT'] \
                       .sum().reindex(prof_codes, fill_value=0.0)
        # Denominator: sum of WEIGHT_PERCENT for species with NonVOCTOG == 0
        denom = sp_nontox.loc[sp_nontox['NonVOCTOG'] == 0] \
                          .groupby('PROFILE_CODE')['WEIGHT_PERCENT'] \
                          .sum().reindex(prof_codes, fill_value=0.0)

        # Ratio: num / denom (0 if denom == 0; aka NonVOCTOG.sum == 0)        
        ratio = (num.astype(float).div(denom.replace(0.0, np.nan))
                .fillna(0.0).to_numpy())

        dfgscnv = pd.DataFrame({'INPUT.POLL': 'NONHAPVOC', 'OUTPUT.POLL': 'NONHAPTOG',
                                'PROFILE': prof_codes.values, 'OUTPUT.MASS/INPUT.MASS': ratio})
    else:
        raise ValueError("RUN_TYPE must be one of: CRITERIA, INTEGRATE, or NOINTEGRATE")

    # Format numeric column with 8 decimal places
    dfgscnv['OUTPUT.MASS/INPUT.MASS'] = dfgscnv['OUTPUT.MASS/INPUT.MASS'].astype(float).map(lambda x: f"{x:.8f}")
    ### Output gscnv df to file
    dfgscnv.to_csv(CNV_OUT, index=False, header=False)
####################################################################################################
        
####################################################################################################
def format_and_header(MECH_BASIS,RUN_TYPE,AQM,MW_FILE,TOX_FILE,CNV_OUT):
    """
    Prepend a header block to the existing GSCNV at CNV_OUT.
    """
    ### Import gscnv csv file.
    f1_gscnv = np.genfromtxt(CNV_OUT,delimiter=',',dtype='str')
    ### Add header lines.
    header_lines = [
        f"#S2S_AQM             {AQM}",
        "#S2S_CAMX_FCRS       Not Applicable",
        f"#S2S_MW              {MW_FILE}",
        f"#S2S_MECH_BASIS      {MECH_BASIS}",
        f"#S2S_RUN_TYPE        {RUN_TYPE}",
        f"#S2S_RUN_DATE        {today}",
        f"#S2S_TBL_TOX         {'Not Applicable' if RUN_TYPE == 'CRITERIA' else TOX_FILE}",
        "#BY PROFILE"]

    total_header = "\n".join(header_lines)
    np.savetxt(CNV_OUT,f1_gscnv[:],fmt='%-20s',header=total_header,comments='')
####################################################################################################