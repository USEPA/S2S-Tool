from typing import Sequence
import pandas as pd
import numpy as np

####################################################################################################
### QA checks for S2S-Tool inputs
####################################################################################################

####################################################################################################
def _require_columns(df: pd.DataFrame, required: Sequence[str], df_name: str) -> None:
    """
    Ensure that a DataFrame contains all required columns.
    Raises ValueError with a helpful message if any are missing.
    """
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Required columns missing in {df_name}: {missing}. "
            f"Columns present: {list(df.columns)}"
        )
####################################################################################################

####################################################################################################
def check_basis_output(MECH_BASIS: str, OUTPUT: str) -> None:
    """
    Validate that the MECH_BASIS is allowed for the selected OUTPUT.
    """
    allowed_voc = {
        'CB6R3_AE7', 'CB6R5_AE7', 'CB7_AE7', 'CB6R4_CF2', 'CB7_CF2',
        'CB7VCP_CF2', 'CB7VCP_CF2', 'CB6R3_AE7_TRACER', 'CRACMMv1.0',
        'CRACMMv2.0','SAPRC07TC_AE7', 'SAPRC07_CF2', 'GEOSChem14.6.3'}
    allowed_pm = {'PM-AE6', 'PM-CR1', 'PM-CR2', 'PM-GC'}

    if OUTPUT == 'VOC':
        if MECH_BASIS not in allowed_voc:
            raise ValueError(
                "MECH_BASIS for OUTPUT==VOC is not allowed. "
                f"Allowed: {sorted(allowed_voc)}. Got: {MECH_BASIS}"
            )
    elif OUTPUT == 'PM':
        if MECH_BASIS not in allowed_pm:
            raise ValueError(
                "MECH_BASIS for OUTPUT==PM is not allowed. "
                f"Allowed: {sorted(allowed_pm)}. Got: {MECH_BASIS}"
            )
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")
####################################################################################################

####################################################################################################
def check_inputs(molwght: pd.DataFrame, mech4import: pd.DataFrame, mechPM: pd.DataFrame,
                 profiles: pd.DataFrame, tbl_tox: pd.DataFrame, OUTPUT: str) -> None:
    """
    Validate that key filtered inputs are present (non-empty) given OUTPUT.
    """
    if OUTPUT == 'VOC':
        if molwght.empty:
            raise ValueError("The MECH_BASIS entered is not in the molwght file (after filtering).")
        if mech4import.empty:
            raise ValueError("The MECH_BASIS entered is not in the mechanism_forImport file (after filtering).")
        if profiles.empty:
            raise ValueError("profiles file empty (after filtering).")
    elif OUTPUT == 'PM':
        if mechPM.empty:
            raise ValueError("The MECH_BASIS and AQM entered is not in the mech_pm file (after filtering).")
        if profiles.empty:
            raise ValueError("profiles file empty (after filtering).")
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")

    if tbl_tox.empty:
        raise ValueError("The AQM entered is not in the tbl_tox file (after filtering).")
####################################################################################################

####################################################################################################
def check_species_profiles(profiles: pd.DataFrame, species: pd.DataFrame) -> None:
    """
    Ensure all profiles in export_profiles are represented in export_species.
    """
    _require_columns(profiles, ['PROFILE_CODE'], 'profiles')
    _require_columns(species, ['PROFILE_CODE'], 'species')

    profile_codes = set(profiles['PROFILE_CODE'])
    species_codes = set(species['PROFILE_CODE'])

    missing = sorted(profile_codes - species_codes)
    if missing:
        raise ValueError(
            "Profiles are in export_profiles but not in export_species. "
            f"Missing PROFILE_CODEs: {missing}"
        )
####################################################################################################

####################################################################################################
def check_species_properties(species_props: pd.DataFrame, species: pd.DataFrame) -> None:
    """
    Ensure all species in export_species have properties (e.g., molecular weight).
    """
    _require_columns(species_props, ['SPECIES_ID'], 'species_props')
    _require_columns(species, ['SPECIES_ID'], 'species')

    species_ids = set(species['SPECIES_ID'])
    props_ids = set(species_props['SPECIES_ID'])

    missing = sorted(species_ids - props_ids)
    if missing:
        raise ValueError(
            "Species found in export_species but missing from export_species_properties. "
            f"Missing SPECIES_IDs: {missing}"
        )
####################################################################################################

####################################################################################################
def check_tox_properties(species_props: pd.DataFrame, tbl_tox: pd.DataFrame) -> None:
    """
    Ensure all species in tbl_tox have entries in species_props.
    """
    _require_columns(species_props, ['SPECIES_ID'], 'species_props')
    _require_columns(tbl_tox, ['SPECIES_ID'], 'tbl_tox')

    tox_ids = set(tbl_tox['SPECIES_ID'])
    props_ids = set(species_props['SPECIES_ID'])

    missing = sorted(tox_ids - props_ids)
    if missing:
        raise ValueError(
            "Species found in tbl_tox but missing from export_species_properties. "
            f"Missing SPECIES_IDs: {missing}"
        )
####################################################################################################

####################################################################################################
def check_molwght_mech4import(molwght: pd.DataFrame, mech4import: pd.DataFrame, OUTPUT: str) -> None:
    """
    When OUTPUT=='VOC', ensure all mechanism species in mech4import exist in molwght.
    """
    if OUTPUT == 'VOC':
        _require_columns(molwght, ['Species'], 'molwght')
        _require_columns(mech4import, ['Species'], 'mech4import')

        mw_species = set(molwght['Species'])
        mech_species = set(mech4import['Species'])

        missing = sorted(mech_species - mw_species)
        if missing:
            raise ValueError(
                "Species present in mech4import are missing from molwght. "
                f"Missing species: {missing}"
            )
    elif OUTPUT == 'PM':
        return
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")
####################################################################################################

####################################################################################################
def check_species_mech4import(mech4import: pd.DataFrame, species: pd.DataFrame, MECH_BASIS: str, OUTPUT: str) -> None:
    """
    Ensure VOC species present in export_species are present in mech4import.
    """
    if OUTPUT == 'VOC':
        _require_columns(mech4import, ['SPECIES_ID'], 'mech4import')
        _require_columns(species, ['SPECIES_ID'], 'species')
        
        species_ids = set(species['SPECIES_ID'])
        mech_species = set(mech4import['SPECIES_ID'])

        missing = sorted(species_ids - mech_species)
        if missing:
            raise ValueError(
                "Species are in export_species but not in mechanism_forImport . "
                f"SPECIES_IDs: {missing}. For MECH_BASIS: {MECH_BASIS}"
            )
    elif OUTPUT == 'PM':
        return
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")
####################################################################################################

####################################################################################################
def check_profiles_volatility(profiles: pd.DataFrame, poa_volatility: pd.DataFrame, OUTPUT: str) -> None:
    """
    Ensure PM profiles have default OM volatility profiles based on L1 and L2 categories.
    """
    if OUTPUT == 'PM':
        _require_columns(profiles, ['PROFILE_CODE',
                                    'CATEGORY_LEVEL_1_Generation_Mechanism',
                                    'CATEGORY_LEVEL_2_Sector_Equipment'], 'profiles')
        _require_columns(poa_volatility, ['CATEGORY_LEVEL_1_Generation_Mechanism',
                                          'CATEGORY_LEVEL_2_Sector_Equipment'], 'poa_volatility')
        
        l1, l2 = 'CATEGORY_LEVEL_1_Generation_Mechanism', 'CATEGORY_LEVEL_2_Sector_Equipment'

        # Build sets of unique (L1, L2) pairs from each DataFrame
        profiles_pairs = set(
            profiles[[l1, l2]].drop_duplicates().itertuples(index=False, name=None)
        )
        volatility_pairs = set(
            poa_volatility[[l1, l2]].drop_duplicates().itertuples(index=False, name=None)
        )

        missing_pairs = sorted(profiles_pairs - volatility_pairs)
        if missing_pairs:
            # Find which PROFILE_CODEs are impacted by the missing category pairs
            missing_df = pd.DataFrame(missing_pairs, columns=[l1, l2])
            affected = profiles.merge(missing_df, on=[l1, l2], how='inner')
            affected_codes = sorted(affected['PROFILE_CODE'].unique().tolist())

            raise ValueError(
                "PM profiles have unspecified OM volatility in poa_volatility. "
                f"Missing category pairs: {missing_pairs}. "
                f"Affected PROFILE_CODEs: {affected_codes}. "
            )
    elif OUTPUT == 'VOC':
        return
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")
####################################################################################################

####################################################################################################
def check_volatility(poa_volatility: pd.DataFrame, OUTPUT: str, tol: float = 1e-6) -> None:
    """
    Ensure each category’s POA volatility bins sum to 1 within tolerance.
    """
    if OUTPUT == 'PM':
        # Volatility bin columns used across checks
        VOL_BIN_COLS = ['N2ALK', 'N1ALK', 'P0ALK', 'P1ALK', 'P2ALK','N2OXY8', 'N2OXY4', 'N2OXY2', 
                        'N1OXY6', 'N1OXY3', 'N1OXY1','P0OXY4', 'P0OXY2', 'P1OXY3', 'P1OXY1', 'P2OXY2']
    
        _require_columns(poa_volatility, ['CATEGORY_LEVEL_1_Generation_Mechanism',
                                          'CATEGORY_LEVEL_2_Sector_Equipment'], 'poa_volatility')
        vol_bins = poa_volatility.reindex(columns=VOL_BIN_COLS, fill_value=0)
        total = vol_bins.sum(axis=1)

        ok = np.isclose(total.values, 1.0, atol=tol)
        if not np.all(ok):
            wrong = poa_volatility.loc[~ok, ['CATEGORY_LEVEL_1_Generation_Mechanism',
                                             'CATEGORY_LEVEL_2_Sector_Equipment']].copy()
            # Build readable labels
            labels = (wrong['CATEGORY_LEVEL_1_Generation_Mechanism'].astype(str) + ', ' +
                      wrong['CATEGORY_LEVEL_2_Sector_Equipment'].astype(str)).tolist()
            raise ValueError(
                "The POA volatility specified for some categories does not sum to 1 "
                f"(tolerance={tol}). Categories: {labels}"
            )
    elif OUTPUT == 'VOC':
        return
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")
####################################################################################################

####################################################################################################
def check_mech_in_mapping(poa_mapping: pd.DataFrame, MECH_BASIS: str, OUTPUT: str, 
                          case_sensitive: bool = True, strip: bool = True) -> None:
    """
    Ensure MECH_BASIS features values in poa_mapping file.
    """
    if OUTPUT == 'PM':
        cols = list(poa_mapping.columns)
        cols = cols[3:]

        # Normalize string-like column labels (ignore non-strings)
        def _norm(x):
            if not isinstance(x, str):
                return x
            s = x.strip() if strip else x
            return s if case_sensitive else s.lower()

        normalized_cols = [_norm(c) for c in cols]
        target = _norm(MECH_BASIS)

        if target not in normalized_cols:
            # Show available columns as strings for a clear message
            available = [str(c) for c in cols]
            raise ValueError(
                f"Mechanism '{MECH_BASIS}' is not present in poa_mapping file. "
                f"Available columns: {available}"
            )    
    elif OUTPUT == 'VOC':
        return
    else:
        raise ValueError("OUTPUT entered is not recognized. Only 'VOC' and 'PM' are allowed.")
####################################################################################################