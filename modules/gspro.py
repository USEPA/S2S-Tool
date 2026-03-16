import sys
import numpy as np
import pandas as pd
from datetime import date

today     = date.today()

####################################################################################################
### This function generates a gspro file for the target MECH_BASIS.
####################################################################################################

####################################################################################################
def gen_gspro_voc(profiles: pd.DataFrame, species: pd.DataFrame, molwght: pd.DataFrame, 
                   mech4import: pd.DataFrame, tbl_tox: pd.DataFrame, MECH_BASIS: str, 
                   RUN_TYPE: str, TOLERANCE: float, TOX_IN: str, PRO_OUT: str):

    # Compute total weight per profile for tolerance check
    total_wt = species.groupby('PROFILE_CODE', as_index=True)['WEIGHT_PERCENT'].sum()

    # Tolerance check: sum should be within [100*(1-TOL), 100*(1+TOL)]
    low, high = 100.0 * (1 - TOLERANCE), 100.0 * (1 + TOLERANCE)
    valid_profiles = total_wt.index[(total_wt >= low) & (total_wt <= high)]

    # Unique profiles list that pass tolerance
    prof_codes = profiles['PROFILE_CODE'].unique()
    prof_codes = pd.Index(prof_codes).intersection(valid_profiles)

    # Filter species DataFrame for valid profiles and drop zero weights
    sp = species.loc[species['PROFILE_CODE'].isin(prof_codes), ['PROFILE_CODE', 'SPECIES_ID', 'WEIGHT_PERCENT']].copy()
    sp = sp.loc[sp['WEIGHT_PERCENT'] > 0.0]

    # Normalize to fraction within each profile
    sp['WEIGHT_PERCENT'] = sp['WEIGHT_PERCENT'] / sp.groupby('PROFILE_CODE')['WEIGHT_PERCENT'].transform('sum')

    # Remove toxics for NOINTEGRATE / INTEGRATE
    if RUN_TYPE in {'NOINTEGRATE', 'INTEGRATE'}:
        tox_ids = set(tbl_tox['SPECIES_ID'].drop_duplicates().tolist())
        sp = sp.loc[~sp['SPECIES_ID'].isin(tox_ids)]

    # Normalize to fraction within each profile for CRITERIA and INTEGRATE
    if RUN_TYPE in {'CRITERIA', 'INTEGRATE'}:
        sp['WEIGHT_PERCENT'] = sp['WEIGHT_PERCENT'] / sp.groupby('PROFILE_CODE')['WEIGHT_PERCENT'].transform('sum')

    # Compute sum per profile, drop and report profiles that are empty following pollutant integration
    sum_by_profile = (sp.groupby('PROFILE_CODE')['WEIGHT_PERCENT'].sum()
                    .reindex(prof_codes, fill_value=0.0))
    empty_profiles = sum_by_profile.index[sum_by_profile.eq(0.0)]
    keep_profiles  = sum_by_profile.index[sum_by_profile.gt(0.0)]

    for pc in empty_profiles:
         print(f'Profile "{pc}" is empty following pollutant integration and therefore not processed.')

    sp = sp.loc[sp['PROFILE_CODE'].isin(keep_profiles)].copy()

    # NMOG fraction per profile (exclude SPECIES_ID == 529; aka methane)
    nmog_perc = sp.loc[sp['SPECIES_ID'] != 529].groupby('PROFILE_CODE')['WEIGHT_PERCENT'].sum()
    nmog_perc = nmog_perc.reindex(prof_codes, fill_value=0.0)

    # Set input pollutant name
    i_poll = 'TOG' if RUN_TYPE in {'CRITERIA', 'NOINTEGRATE'} else TOX_IN

    # mech_map: mech4import with model species MW attached (MW_model), per row
    mech_map = mech4import.merge(molwght[['Species', 'SPEC_MW']], on='Species', how='left'
                                 ).rename(columns={'SPEC_MW': 'MW_model'})

    # Sum over mapping rows per compound SPECIES_ID
    eff_mw_by_specid = (mech_map['Moles'] * mech_map['MW_model']).groupby(mech_map['SPECIES_ID']).sum()

    # Inner-join is safe and avoids NaNs
    sp_mech = sp.merge(mech_map[['SPECIES_ID', 'Species', 'Moles', 'MW_model']],
                       on='SPECIES_ID',how='inner')
    
    # Attach effective MW per compound SPECIES_ID
    sp_mech['eff_MW'] = sp_mech['SPECIES_ID'].map(eff_mw_by_specid)

    # Calculate mole model species / mass VOC
    sp_mech['moleSpec_massVOC'] = sp_mech['WEIGHT_PERCENT'] * sp_mech['Moles'] / sp_mech['eff_MW']

    grouped = sp_mech.groupby(['PROFILE_CODE', 'Species'], as_index=False).agg(
        moleSpec_massVOC=('moleSpec_massVOC', 'sum'), MW_model=('MW_model', 'first'))
    
    # Final mass fractions for model species
    grouped['MASS.FRACTION'] = grouped['moleSpec_massVOC'] * grouped['MW_model']
    
    dfgspro = grouped.rename(columns={'PROFILE_CODE': 'PROFILE', 'Species': 'MODEL.SPECIES',
                                      'MW_model': 'MOLECULAR.WGHT'})
    [['PROFILE', 'MODEL.SPECIES', 'MASS.FRACTION', 'MOLECULAR.WGHT']]

    # Add INPUT.POLL and MASS.FRACTION1
    dfgspro.insert(1, 'INPUT.POLL', i_poll)
    dfgspro['MASS.FRACTION1'] = dfgspro['MASS.FRACTION']

    # Add NMOG
    if MECH_BASIS != 'CB6R3_AE7_TRACER':
        add = (nmog_perc.rename('MASS.FRACTION').rename_axis('PROFILE').reset_index())
        add['INPUT.POLL'] = i_poll
        add['MODEL.SPECIES'] = 'NMOG'
        add['MOLECULAR.WGHT'] = 1.0
        add['MASS.FRACTION1'] = add['MASS.FRACTION']
        add = add[['PROFILE', 'INPUT.POLL', 'MODEL.SPECIES', 'MASS.FRACTION', 'MOLECULAR.WGHT', 'MASS.FRACTION1']]
        dfgspro = pd.concat([dfgspro, add], ignore_index=True)
    else: # Remove NONBAF from CB6R3_AE7_TRACER GSPROs
        dfgspro = dfgspro.loc[dfgspro['MODEL.SPECIES'] != 'NONBAF']

    # Order columns and format numbers
    dfgspro = dfgspro[['PROFILE', 'INPUT.POLL', 'MODEL.SPECIES', 'MASS.FRACTION', 'MOLECULAR.WGHT', 'MASS.FRACTION1']]
    dfgspro['MASS.FRACTION']  = dfgspro['MASS.FRACTION'].astype(float).map(lambda x: f"{x:.6E}")
    dfgspro['MOLECULAR.WGHT'] = dfgspro['MOLECULAR.WGHT'].astype(float).map(lambda x: f"{x:.6E}")
    dfgspro['MASS.FRACTION1'] = dfgspro['MASS.FRACTION1'].astype(float).map(lambda x: f"{x:.6E}")

    # Write
    dfgspro = dfgspro.sort_values(['PROFILE'], kind='mergesort').reset_index(drop=True) # Ordered by Profile
    dfgspro.to_csv(PRO_OUT, index=False, header=False)
####################################################################################################

####################################################################################################
def build_pm_ready_profile(p: pd.Series, temp_spec: pd.DataFrame, oxygen_metals: pd.DataFrame) -> pd.DataFrame:
    """
    Translate a base SPECIATE PM profile (ptype == 'PM') into a 'PM-ready' profile.

    Parameters
    ----------
    p : pd.Series
        Row from the profiles DataFrame for the current profile. Must include
        PROFILE_CODE, CATEGORY_LEVEL_1_Generation_Mechanism, CATEGORY_LEVEL_2_Sector_Equipment.
    temp_spec : pd.DataFrame
        Species rows for the current PROFILE_CODE with columns ['PROFILE_CODE','SPECIES_ID','WEIGHT_PERCENT'].
        WEIGHT_PERCENT is expected to be percent values (sum near 100).
    oxygen_metals : pd.DataFrame
        DataFrame with metal oxygen ratios. Must include columns ['SPECIES_ID', 'oxy/metal_ratio'].

    Returns
    -------
    pd.DataFrame
        A PM-ready species DataFrame with columns ['PROFILE_CODE','SPECIES_ID','WEIGHT_PERCENT'].
        WEIGHT_PERCENT is in percent (not fraction), potentially renormalized.
    """
    prof = p['PROFILE_CODE']  # Pull profile code
    l1   = p['CATEGORY_LEVEL_1_Generation_Mechanism']  # Pull L1 category
    l2   = p['CATEGORY_LEVEL_2_Sector_Equipment']  # Pull L2 category

    # Hardcoded list of "PM-ready" species IDs
    pm_ready_ids = [
        626, 797, 699, 700, 613, 784, 2669, 488, 292, 694, 715, 2303, 329, 2772, 525, 2302, 669, 526, 785,
        696, 337, 795, 2668, 2671, 666, 767, 347, 379, 612, 380, 778, 468, 298, 693, 689, 697, 779, 586,
        649, 695, 328, 487, 714, 296, 300, 519, 1861, 528, 520, 2670]
    
    temp_pm = pd.DataFrame({'PROFILE_CODE': prof, 'SPECIES_ID': pm_ready_ids})
    temp_pm = temp_pm.merge(temp_spec[['SPECIES_ID', 'WEIGHT_PERCENT']], on='SPECIES_ID', how='left')
    temp_pm['WEIGHT_PERCENT'] = temp_pm['WEIGHT_PERCENT'].fillna(0.0)

    # Helpers
    def get_w(spec_id: int) -> float:
        return float(temp_pm.loc[temp_pm['SPECIES_ID'] == spec_id, 'WEIGHT_PERCENT'].iloc[0])

    def set_w(spec_id: int, value: float) -> None:
        temp_pm.loc[temp_pm['SPECIES_ID'] == spec_id, 'WEIGHT_PERCENT'] = float(value)

    # If ionic species present, zero out the elemental counterpart
    if get_w(2303) > 0.0:  # Ca++ present
        set_w(329, 0.0)    # zero Ca
    if get_w(2772) > 0.0:  # Mg++ present
        set_w(525, 0.0)    # zero Mg
    if get_w(2302) > 0.0:  # K+ present
        set_w(669, 0.0)    # zero K
    if get_w(785) > 0.0:   # Na+ present
        set_w(696, 0.0)    # zero Na
    if get_w(337) > 0.0:   # Cl- present
        set_w(795, 0.0)    # zero Cl

    # Particulate water (PH2O, SPECIES_ID 2668)
    if l1 in ('Combustion', 'Miscellaneous'):
        ph2o = (get_w(699) + get_w(784)) * 0.24
        set_w(2668, ph2o)
    else:
        set_w(2668, 0.0)

    # PSO4 from S if needed (699 = PSO4, 700 = S)
    if get_w(699) > 0.0:
        pass
    elif get_w(700) == 0.0:
        pass
    else:
        # Convert S to PSO4 (molecular weight ratio: 96/32)
        set_w(699, get_w(700) * (96.0 / 32.0))
    set_w(700, 0.0)

    # Metal-bound oxygen (MO, SPECIES_ID 2670)
    om_temp = oxygen_metals.merge(temp_pm[['SPECIES_ID', 'WEIGHT_PERCENT']],
                                  on='SPECIES_ID',how='left')
    om_temp['WEIGHT_PERCENT'] = om_temp['WEIGHT_PERCENT'].fillna(0.0)

    # For pairs (elemental vs ionic), if either is absent (0), set both to 0 to avoid negative oxygen calc
    # 329 (Ca) vs 2303 (Ca++)
    if get_w(329) == 0.0 or get_w(2303) == 0.0:
        om_temp.loc[om_temp['SPECIES_ID'].isin([329, 2303]), 'WEIGHT_PERCENT'] = 0.0
    # 669 (K) vs 2302 (K+)
    if get_w(669) == 0.0 or get_w(2302) == 0.0:
        om_temp.loc[om_temp['SPECIES_ID'].isin([669, 2302]), 'WEIGHT_PERCENT'] = 0.0
    # 525 (Mg) vs 2772 (Mg++)
    if get_w(525) == 0.0 or get_w(2772) == 0.0:
        om_temp.loc[om_temp['SPECIES_ID'].isin([525, 2772]), 'WEIGHT_PERCENT'] = 0.0
    # 696 (Na) vs 785 (Na+)
    if get_w(696) == 0.0 or get_w(785) == 0.0:
        om_temp.loc[om_temp['SPECIES_ID'].isin([696, 785]), 'WEIGHT_PERCENT'] = 0.0

    # Subtract ionic from elemental for oxygen calc
    def subtract_pair(elem: int, ion: int):
        e = float(om_temp.loc[om_temp['SPECIES_ID'] == elem, 'WEIGHT_PERCENT'].iloc[0])
        i = float(om_temp.loc[om_temp['SPECIES_ID'] == ion, 'WEIGHT_PERCENT'].iloc[0])
        v = max(e - i, 0.0)
        om_temp.loc[om_temp['SPECIES_ID'] == elem, 'WEIGHT_PERCENT'] = v

    subtract_pair(329, 2303)  # Ca - Ca++
    subtract_pair(669, 2302)  # K  - K+
    subtract_pair(525, 2772)  # Mg - Mg++
    subtract_pair(696, 785)   # Na - Na+

    # Compute unadjusted MO: sum(oxy/metal_ratio * metal_weight)
    if 'oxy/metal_ratio' not in om_temp.columns:
        raise ValueError("oxygen_metals must include column 'oxy/metal_ratio'")
    MOunadjusted = float((om_temp['oxy/metal_ratio'] * om_temp['WEIGHT_PERCENT']).sum())

    # Adjust MO for non-neutralized sulfate (based on NH4 and PSO4)
    NeutralizedSO4    = (0.5 * 96.0 / 18.0) * get_w(784)  # based on NH4 weight
    NonNeutralizedSO4 = get_w(699) - NeutralizedSO4
    if NonNeutralizedSO4 > 0.0:
        MOadjusted = MOunadjusted - NonNeutralizedSO4 * (16.0 / 96.0)
    else:
        MOadjusted = MOunadjusted
    MOadjusted = max(MOadjusted, 0.0)
    set_w(2670, MOadjusted)  # MO

    # PNCOM (SPECIES_ID 2669) additions by sector
    l2_lower = str(l2).lower()
    if 'mobile' in l2_lower:
        set_w(2669, get_w(626) * 0.25)
    elif 'biomass burning' in l2_lower:
        set_w(2669, get_w(626) * 0.7)
    else:
        set_w(2669, get_w(626) * 0.4)
        
    # Keep only positive weights
    temp_pm = temp_pm.loc[temp_pm['WEIGHT_PERCENT'] > 0.0].copy()

    # Renormalize if sum exceeds 101%
    s = float(temp_pm['WEIGHT_PERCENT'].sum())
    if s > 101.0:
        temp_pm['WEIGHT_PERCENT'] = temp_pm['WEIGHT_PERCENT'] / s * 100.0

    temp_pm.reset_index(drop=True, inplace=True)

    return temp_pm[['PROFILE_CODE', 'SPECIES_ID', 'WEIGHT_PERCENT']]
####################################################################################################

####################################################################################################
def gen_gspro_pm(profiles,species,mechPM,tbl_tox,poa_volatility,poa_mapping,camx_fcrs,oxygen_metals,MECH_BASIS,RUN_TYPE,AQM,TOX_IN,PRO_OUT):
    
    ### gspro file columns
    column_names  = ['PROFILE','INPUT.POLL','MODEL.SPECIES','MASS.FRACTION','MOLECULAR.WGHT','MASS.FRACTION1']
    dfgspro       = pd.DataFrame(columns=column_names)

    TOXLIST = set(tbl_tox['SPECIES_ID'].drop_duplicates().tolist()) # pull integrated SPECIES_ID list
    is_camx = (AQM == 'CAMX')
    is_cmaq = (AQM == 'CMAQ')

    # Set a weight for a species by name or SPECIES_ID
    def set_weight(temp_mech: pd.DataFrame, key, value, by: str = 'Species'):
        temp_mech.loc[temp_mech[by] == key, 'WEIGHT_PERCENT'] = float(value)
        
    for i in range(len(profiles)):
        p     = profiles.iloc[i]
        prof  = p['PROFILE_CODE']
        ptype = p['PROFILE_TYPE']
        l1    = p['CATEGORY_LEVEL_1_Generation_Mechanism']
        l2    = p['CATEGORY_LEVEL_2_Sector_Equipment']

        # Pull target speciation profile
        temp_spec = species.loc[species['PROFILE_CODE'] == prof].copy().reset_index(drop=True)

        # Translate data into a "PM-ready" profile
        if ptype == 'PM':
            temp_spec = build_pm_ready_profile(p, temp_spec, oxygen_metals)

        # Integrate vs criteria normalization
        if RUN_TYPE == 'CRITERIA':
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / 100.0
            i_poll = 'PM2_5'
        else: # if RUN_TYPE=='INTEGRATE'
            # # pull target speciation profile minus TOXLIST
            temp_spec = temp_spec.loc[~temp_spec['SPECIES_ID'].isin(TOXLIST)].copy()
            # Renormalize Wght%
            temp_spec['WEIGHT_PERCENT'] = temp_spec['WEIGHT_PERCENT'] / temp_spec['WEIGHT_PERCENT'].sum()
            i_poll = TOX_IN

        # Append WEIGHT_PERCENT and then replace all NaN with zeros
        temp_mech = mechPM.merge(temp_spec[['SPECIES_ID', 'WEIGHT_PERCENT']].drop_duplicates(),
                                 on='SPECIES_ID', how='left')
        temp_mech['WEIGHT_PERCENT'] = temp_mech['WEIGHT_PERCENT'].fillna(0.0)

        # Append WEIGHT_PERCENT, replace all NaN with zeros, 
        # remove species where WEIGHT_PERCENT == 0.0, and calculate total POA weight percent
        poa_mech = poa_mapping.merge(temp_spec[['SPECIES_ID', 'WEIGHT_PERCENT']].drop_duplicates(),
                                     on='SPECIES_ID', how='left').fillna({'WEIGHT_PERCENT': 0.0})
        poa_mech = poa_mech.loc[poa_mech['WEIGHT_PERCENT'] > 0.0].reset_index(drop=True)
        poa_total = float(poa_mech['WEIGHT_PERCENT'].sum())

        # Extract SV-POA profile assigned to source type
        temp_poa = poa_volatility.loc[(poa_volatility['CATEGORY_LEVEL_1_Generation_Mechanism'] == l1) &
                                      (poa_volatility['CATEGORY_LEVEL_2_Sector_Equipment'] == l2)
                                      ].reset_index(drop=True)
        
        # Mechanism-specific mapping of OM/OC/NCOM from temp_spec
        if MECH_BASIS == 'PM-AE6':
            if ptype in ('PM-AE6', 'PM'):
                pass  # OC/PNCOM already appropriately mapped in base profile
            elif ptype == 'PM-AE8':
                set_weight(temp_mech, 'POC', poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC',   'WEIGHT_PERCENT'].sum(), by='Species') # add POC
                set_weight(temp_mech, 'PNCOM', poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum(), by='Species') # add PNCOM
            elif ptype in ('PM-CR1', 'PM-CR2'):
                if poa_total > 0.0:
                    ratio = float(p['ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'])
                    set_weight(temp_mech, 'POC',   poa_total * (1.0 / ratio), by='Species')
                    set_weight(temp_mech, 'PNCOM', poa_total * (1.0 - 1.0 / ratio), by='Species')
            else:
                sys.exit(f'PROFILE_TYPE = {ptype} for profile {prof} is not recognized.')

        elif MECH_BASIS == 'PM-CR1':
            if ptype in ('PM-AE6', 'PM'):
                A_N2 = float(temp_poa['N2ALK'].iat[0] + temp_poa['N2OXY8'].iat[0] + temp_poa['N2OXY4'].iat[0] + temp_poa['N2OXY2'].iat[0])
                A_N1 = float(temp_poa['N1ALK'].iat[0] + temp_poa['N1OXY6'].iat[0] + temp_poa['N1OXY3'].iat[0] + temp_poa['N1OXY1'].iat[0])
                A_P0 = float(temp_poa['P0ALK'].iat[0] + temp_poa['P0OXY4'].iat[0] + temp_poa['P0OXY2'].iat[0])
                A_P1 = float(temp_poa['P1ALK'].iat[0] + temp_poa['P1OXY3'].iat[0] + temp_poa['P1OXY1'].iat[0])
                A_P2 = float(temp_poa['P2ALK'].iat[0] + temp_poa['P2OXY2'].iat[0])
                # Allocate by SPECIES_ID
                set_weight(temp_mech, 3394., poa_total * A_N2, by='SPECIES_ID')  # AROCN2ALK
                set_weight(temp_mech, 3395., poa_total * A_N1, by='SPECIES_ID')  # AROCN1ALK
                set_weight(temp_mech, 3396., poa_total * A_P0, by='SPECIES_ID')  # AROCP0ALK
                set_weight(temp_mech, 3397., poa_total * A_P1, by='SPECIES_ID')  # AROCP1ALK
                set_weight(temp_mech, 3398., poa_total * A_P2, by='SPECIES_ID')  # AROCP2ALK

            elif ptype == 'PM-AE8':
                # AE8 bins mapped to SPECIES_IDs
                set_weight(temp_mech, 3394., poa_mech.loc[poa_mech['log10Cstar'] == -2, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3395., poa_mech.loc[poa_mech['log10Cstar'] == -1, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3396., poa_mech.loc[poa_mech['log10Cstar'] ==  0, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3397., poa_mech.loc[poa_mech['log10Cstar'] ==  1, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3398., poa_mech.loc[poa_mech['log10Cstar'] ==  2, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')

            elif ptype == 'PM-CR1':
                pass  # Already mapped

            elif ptype == 'PM-CR2':
                # Use AE8-style bins for CR2 profile
                set_weight(temp_mech, 3394., poa_mech.loc[poa_mech['log10Cstar'] == -2, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3395., poa_mech.loc[poa_mech['log10Cstar'] == -1, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3396., poa_mech.loc[poa_mech['log10Cstar'] ==  0, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3397., poa_mech.loc[poa_mech['log10Cstar'] ==  1, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
                set_weight(temp_mech, 3398., poa_mech.loc[poa_mech['log10Cstar'] ==  2, 'WEIGHT_PERCENT'].sum(), by='SPECIES_ID')
            else:
                sys.exit(f'PROFILE_TYPE = {ptype} for profile {prof} is not recognized.')

        elif MECH_BASIS == 'PM-CR2':
            if ptype in ('PM-AE6', 'PM'):
                # Direct allocation by volatility bins
                set_weight(temp_mech, 3394., poa_total * float(temp_poa['N2ALK'].iat[0]),  by='SPECIES_ID') # add AROCN2ALK
                set_weight(temp_mech, 3395., poa_total * float(temp_poa['N1ALK'].iat[0]),  by='SPECIES_ID') # add AROCN1ALK
                set_weight(temp_mech, 3396., poa_total * float(temp_poa['P0ALK'].iat[0]),  by='SPECIES_ID') # add AROCP0ALK
                set_weight(temp_mech, 3397., poa_total * float(temp_poa['P1ALK'].iat[0]),  by='SPECIES_ID') # add AROCP1ALK
                set_weight(temp_mech, 3398., poa_total * float(temp_poa['P2ALK'].iat[0]),  by='SPECIES_ID') # add AROCP2ALK

                set_weight(temp_mech, 3524., poa_total * float(temp_poa['N2OXY8'].iat[0]), by='SPECIES_ID') # add AROCN2OXY8
                set_weight(temp_mech, 3523., poa_total * float(temp_poa['N2OXY4'].iat[0]), by='SPECIES_ID') # add AROCN2OXY4
                set_weight(temp_mech, 3522., poa_total * float(temp_poa['N2OXY2'].iat[0]), by='SPECIES_ID') # add AROCN2OXY2
                set_weight(temp_mech, 3527., poa_total * float(temp_poa['N1OXY6'].iat[0]), by='SPECIES_ID') # add AROCN1OXY6
                set_weight(temp_mech, 3526., poa_total * float(temp_poa['N1OXY3'].iat[0]), by='SPECIES_ID') # add AROCN1OXY3
                set_weight(temp_mech, 3525., poa_total * float(temp_poa['N1OXY1'].iat[0]), by='SPECIES_ID') # add AROCN1OXY1

                set_weight(temp_mech, 3529., poa_total * float(temp_poa['P0OXY4'].iat[0]), by='SPECIES_ID') # add AROCP0OXY4
                set_weight(temp_mech, 3528., poa_total * float(temp_poa['P0OXY2'].iat[0]), by='SPECIES_ID') # add AROCP0OXY2
                set_weight(temp_mech, 3531., poa_total * float(temp_poa['P1OXY3'].iat[0]), by='SPECIES_ID') # add AROCP1OXY3
                set_weight(temp_mech, 3530., poa_total * float(temp_poa['P1OXY1'].iat[0]), by='SPECIES_ID') # add AROCP1OXY1
                set_weight(temp_mech, 3532., poa_total * float(temp_poa['P2OXY2'].iat[0]), by='SPECIES_ID') # add AROCP2OXY2

            elif ptype == 'PM-AE8':
                def split_to(id_map, family_cols):
                    fam_total = float(sum(float(temp_poa[c].iat[0]) for c in family_cols))
                    weight = poa_mech.loc[poa_mech['log10Cstar'] == id_map['cstar'], 'WEIGHT_PERCENT'].sum()
                    if fam_total > 0.0 and weight > 0.0:
                        for spec_id, col in id_map['targets']:
                            frac = float(temp_poa[col].iat[0]) / fam_total
                            set_weight(temp_mech, spec_id, weight * frac, by='SPECIES_ID')

                # N2 family
                split_to({'cstar': -2, 'targets': [(3394., 'N2ALK'), (3524., 'N2OXY8'), (3523., 'N2OXY4'), (3522., 'N2OXY2')]},
                         ['N2ALK', 'N2OXY8', 'N2OXY4', 'N2OXY2'])
                # N1 family
                split_to({'cstar': -1, 'targets': [(3395., 'N1ALK'), (3527., 'N1OXY6'), (3526., 'N1OXY3'), (3525., 'N1OXY1')]},
                         ['N1ALK', 'N1OXY6', 'N1OXY3', 'N1OXY1'])
                # P0 family
                split_to({'cstar': 0, 'targets': [(3396., 'P0ALK'), (3529., 'P0OXY4'), (3528., 'P0OXY2')]},
                         ['P0ALK', 'P0OXY4', 'P0OXY2'])
                # P1 family
                split_to({'cstar': 1, 'targets': [(3397., 'P1ALK'), (3531., 'P1OXY3'), (3530., 'P1OXY1')]},
                         ['P1ALK', 'P1OXY3', 'P1OXY1'])
                # P2 family
                split_to({'cstar': 2, 'targets': [(3398., 'P2ALK'), (3532., 'P2OXY2')]},
                         ['P2ALK', 'P2OXY2'])

            elif ptype == 'PM-CR1':
                pass  # already mapped in base profile
                
            elif ptype == 'PM-CR2':
                pass  # already mapped in base profile
                
            else:
                sys.exit(f'PROFILE_TYPE = {ptype} for profile {prof} is not recognized.')

        elif MECH_BASIS == 'PM-GC':
            # OC/PNCOM mapping uses OM/OC/NCOM in poa_mapping
            oc_w   = float(poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'OC',   'WEIGHT_PERCENT'].sum())
            pncom_w = float(poa_mech.loc[poa_mech['OM/OC/NCOM'] == 'NCOM', 'WEIGHT_PERCENT'].sum())
            if ptype in ('PM-AE6', 'PM', 'PM-AE8'):
                set_weight(temp_mech, 'OC', oc_w, by='Species')
                set_weight(temp_mech, 'PNCOM', pncom_w, by='Species')
            elif ptype in ('PM-CR1', 'PM-CR2'):
                if poa_total > 0.0:
                    ratio = float(p['ORGANIC_MATTER_to_ORGANIC_CARBON_RATIO'])
                    set_weight(temp_mech, 'OC',    poa_total * (1.0 / ratio), by='Species')
                    set_weight(temp_mech, 'PNCOM', poa_total * (1.0 - 1.0 / ratio), by='Species')
            else:
                sys.exit(f'PROFILE_TYPE = {ptype} for profile {prof} is not recognized.')

        else:
            sys.exit('MECH_BASIS is not recognized.')

        # Ion substitutions
        # Chlorine: If 337 absent but 795 present, use 795 -> PCL
        if float(temp_mech.loc[temp_mech['SPECIES_ID'] == 337, 'WEIGHT_PERCENT'].iloc[0]) == 0.0 and 795 in temp_spec['SPECIES_ID'].values:
            set_weight(temp_mech, 'PCL', float(temp_spec.loc[temp_spec['SPECIES_ID'] == 795, 'WEIGHT_PERCENT'].iloc[0]), by='Species')

        if is_cmaq:
            # Calcium: 2303 absent, but 329 present -> PCA
            if float(temp_mech.loc[temp_mech['SPECIES_ID'] == 2303, 'WEIGHT_PERCENT'].iloc[0]) == 0.0 and 329 in temp_spec['SPECIES_ID'].values:
                set_weight(temp_mech, 'PCA', float(temp_spec.loc[temp_spec['SPECIES_ID'] == 329, 'WEIGHT_PERCENT'].iloc[0]), by='Species')
            # Magnesium: 2772 absent, 525 present -> PMG
            if float(temp_mech.loc[temp_mech['SPECIES_ID'] == 2772, 'WEIGHT_PERCENT'].iloc[0]) == 0.0 and 525 in temp_spec['SPECIES_ID'].values:
                set_weight(temp_mech, 'PMG', float(temp_spec.loc[temp_spec['SPECIES_ID'] == 525, 'WEIGHT_PERCENT'].iloc[0]), by='Species')
            # Potassium: 2302 absent, 669 present -> PK
            if float(temp_mech.loc[temp_mech['SPECIES_ID'] == 2302, 'WEIGHT_PERCENT'].iloc[0]) == 0.0 and 669 in temp_spec['SPECIES_ID'].values:
                set_weight(temp_mech, 'PK', float(temp_spec.loc[temp_spec['SPECIES_ID'] == 669, 'WEIGHT_PERCENT'].iloc[0]), by='Species')
            # Sodium: 785 absent, 696 present -> PNA
            if float(temp_mech.loc[temp_mech['SPECIES_ID'] == 785,  'WEIGHT_PERCENT'].iloc[0]) == 0.0 and 696 in temp_spec['SPECIES_ID'].values:
                set_weight(temp_mech, 'PNA', float(temp_spec.loc[temp_spec['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0]), by='Species')
        else: # if AQM=='CAMX'
            temp_mech.loc[temp_mech['SPECIES_ID'] == 785., 'Species'] = 'NA'
            if float(temp_mech.loc[temp_mech['SPECIES_ID'] == 785,  'WEIGHT_PERCENT'].iloc[0]) == 0.0 and 696 in temp_spec['SPECIES_ID'].values:
                set_weight(temp_mech, 'NA', float(temp_spec.loc[temp_spec['SPECIES_ID'] == 696, 'WEIGHT_PERCENT'].iloc[0]), by='Species')

        # Calculate PMOTHR or FPRM/FCRS for CMAQ and CAMX, respectively
        total_w = float(temp_mech['WEIGHT_PERCENT'].sum())
        if total_w < 1.0:
            if is_cmaq:
                if TOX_IN != 'TOM':
                    set_weight(temp_mech, 'PMOTHR', 1.0 - total_w, by='Species')
            else: # if AQM=='CAMX'
                if TOX_IN != 'TOM':
                    # CAMx: FCRS vs FPRM
                    temp_fcrs = camx_fcrs.loc[camx_fcrs['PROFILE_CODE'].astype(str) == str(prof)]
                    species_name = 'FCRS' if not temp_fcrs.empty else 'FPRM'
                    set_weight(temp_mech, species_name, 1.0 - total_w, by='Species')

        # INTEGRATE: remove toxics that might have been added and renormalize
        if RUN_TYPE == 'INTEGRATE':
            # remove TOXLIST species that might have been added
            temp_mech = temp_mech.loc[~temp_mech['SPECIES_ID'].isin(TOXLIST)]
            sum_w = float(temp_mech['WEIGHT_PERCENT'].sum())
            if sum_w > 0.0: # Renormalize Wght%
                temp_mech['WEIGHT_PERCENT'] = temp_mech['WEIGHT_PERCENT'] / sum_w
        
        # Calculate POA for CAMX and drop PNCOM
        temp_mech = temp_mech.reset_index(drop=True)
        if is_camx:
            poc = float(temp_mech.loc[temp_mech['Species'] == 'POC', 'WEIGHT_PERCENT'].sum())
            pncom = float(temp_mech.loc[temp_mech['Species'] == 'PNCOM', 'WEIGHT_PERCENT'].sum())
            poa_row = pd.Series({'AQM': AQM, 'Mechanism': temp_mech.loc[0, 'Mechanism'],
                                 'SPECIES_ID': 9999, 'Species': 'POA', 'WEIGHT_PERCENT': poc + pncom})
            temp_mech = pd.concat([temp_mech, pd.DataFrame([poa_row])], ignore_index=True)
            temp_mech = temp_mech.loc[temp_mech['Species'] != 'PNCOM']

        # Sum model species with multiple WEIGHT_PERCENT
        temp_mech = temp_mech.groupby('Species', as_index=False)['WEIGHT_PERCENT'].sum()
        # Remove species where WEIGHT_PERCENT == 0.0
        temp_mech = temp_mech.loc[temp_mech['WEIGHT_PERCENT'] > 0.0]

        if temp_mech.empty:
            print(f'Profile "{prof}" is empty following pollutant integration and therefore not processed.')
            continue

        # Final formatting to output schema
        temp_mech = temp_mech.rename(columns={'Species': 'MODEL.SPECIES', 'WEIGHT_PERCENT': 'MASS.FRACTION'})
        temp_mech['MASS.FRACTION1'] = temp_mech['MASS.FRACTION']
        temp_mech['MOLECULAR.WGHT'] = 1.0
        temp_mech['PROFILE'] = prof
        temp_mech['INPUT.POLL'] = i_poll
        temp_mech = temp_mech[column_names].reset_index(drop=True)
        
        # append profile gspro to final gspro
        dfgspro = pd.concat([dfgspro, temp_mech], ignore_index=True)

    # Format numeric columns
    if not dfgspro.empty:
        dfgspro['MASS.FRACTION']  = dfgspro['MASS.FRACTION'].astype(float).map(lambda x: f"{x:.6E}")
        dfgspro['MOLECULAR.WGHT'] = dfgspro['MOLECULAR.WGHT'].astype(float).map(lambda x: f"{x:.6E}")
        dfgspro['MASS.FRACTION1'] = dfgspro['MASS.FRACTION1'].astype(float).map(lambda x: f"{x:.6E}")

    dfgspro.to_csv(PRO_OUT, index=False, header=False)
####################################################################################################

####################################################################################################
def format_and_header(tbl_tox,TOX_IN,MECH_BASIS,RUN_TYPE,AQM,MW_FILE,FCRS_FILE,TOX_FILE,PRO_OUT):
    """
    Prepend a header block to the existing GSPRO at PRO_OUT.
    """
    ### Import gscnv csv file.
    f1_gspro = np.genfromtxt(PRO_OUT,delimiter=',',dtype='str')
    ### Add header lines.    
    header_lines = [
        f"#S2S_AQM             {AQM}",
        f"#S2S_CAMX_FCRS       {FCRS_FILE}",
        f"#S2S_MW              {MW_FILE}",
        f"#S2S_MECH_BASIS      {MECH_BASIS}",
        f"#S2S_RUN_TYPE        {RUN_TYPE}",
        f"#S2S_RUN_DATE        {today}",
        f"#S2S_TBL_TOX         {'Not Applicable' if RUN_TYPE == 'CRITERIA' else TOX_FILE}"]

    # Add NHAP tox lines only for INTEGRATE
    if RUN_TYPE == 'INTEGRATE':
        if 'Inv.Species' not in tbl_tox.columns:
            raise ValueError("tbl_tox must contain the column 'Inv.Species' to build NHAP header lines.")
        # Append each tox species as its own line (keeps the original order)
        # Format: " #NHAP  {TOX_IN:<12} {Inv.Species}"
        for _, row in tbl_tox.iterrows():
            species_name = str(row['Inv.Species']).strip()
            header_lines.append(f" #NHAP  {TOX_IN:<12s} {species_name}")

    # Final header string
    total_header = "\n".join(header_lines)
    
    # Ensure f1_gspro is 2D for savetxt (handles single-row or empty files)
    if f1_gspro.ndim == 1:
        if f1_gspro.size == 0:
            # If file is empty, create an empty array with 6 columns
            f1_gspro = np.empty((0, 6), dtype=str)
        else:
            f1_gspro = f1_gspro.reshape(1, -1)
            
    fmt = '%-20s %-20s %-10s %-13s %-13s %-13s'
    np.savetxt(PRO_OUT, f1_gspro[:], fmt=fmt, header=total_header, comments='')
####################################################################################################