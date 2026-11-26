"""Helpers for retrieving nuclear data from the NNDC ENSDF API.

API documentation: https://www.nndc.bnl.gov/ensdf-api/

The endpoint paths and field names are kept in module-level constants so they
can be adjusted if the service layout differs.
"""
import re

import pandas as pd
import requests

BASE_URL = "https://www.nndc.bnl.gov/ensdf-api"
LEVELS_ENDPOINT = "levels"
GAMMAS_ENDPOINT = "gammas"
TIMEOUT = 30

_NUCLIDE_RE = re.compile(r"^\s*(\d+)?\s*([A-Za-z]{1,2})[-\s]?(\d+)?\s*(m\d*)?\s*$")


def parse_nuclide(name):
    """Parse e.g. 'Co60', 'Co-60', 'Co60m', '60Co', 'Ag110m2'.

    Returns (symbol, mass_number, isomer_index) where isomer_index is 0 for the
    ground state, 1 for 'm'/'m1', 2 for 'm2', etc.
    """
    match = _NUCLIDE_RE.match(name)
    if not match:
        raise ValueError(f"Cannot parse nuclide name: {name!r}")
    mass_prefix, symbol, mass_suffix, meta = match.groups()
    mass = mass_prefix or mass_suffix
    if not mass:
        raise ValueError(f"No mass number in nuclide name: {name!r}")
    isomer = 0
    if meta:
        isomer = int(meta[1:]) if len(meta) > 1 else 1
    return symbol.capitalize(), int(mass), isomer


def _fetch(endpoint, nuclide, session=None):
    """Request records for a nuclide and return them as a DataFrame."""
    symbol, mass, _ = parse_nuclide(nuclide)
    getter = session.get if session is not None else requests.get
    response = getter(
        f"{BASE_URL}/{endpoint}",
        params={"nuc": f"{mass}{symbol}"},
        timeout=TIMEOUT,
    )
    response.raise_for_status()
    data = response.json()
    if isinstance(data, dict):
        for key in ("data", "results", "levels", "gammas"):
            if isinstance(data.get(key), list):
                data = data[key]
                break
        else:
            data = [data]
    return pd.DataFrame(data)


def get_levels(nuclide, session=None):
    """DataFrame of adopted levels for a nuclide."""
    return _fetch(LEVELS_ENDPOINT, nuclide, session)


def get_gammas(nuclide, session=None):
    """DataFrame of gamma records for a nuclide."""
    return _fetch(GAMMAS_ENDPOINT, nuclide, session)


def select_isomer(levels, isomer_index):
    """Select the level row for the ground state (0) or the n-th isomer.

    Isomers are the long-lived excited levels (those with a half-life), ordered
    by energy. Expects columns 'energy' and 'half_life_seconds'.
    """
    levels = levels.sort_values("energy").reset_index(drop=True)
    if isomer_index == 0:
        return levels.iloc[0]
    isomers = levels[(levels["energy"] > 0) & levels["half_life_seconds"].notna()]
    if len(isomers) < isomer_index:
        raise LookupError(f"Isomer m{isomer_index} not found")
    return isomers.iloc[isomer_index - 1]


def get_half_life(nuclide, session=None):
    """Half-life in seconds of a nuclide, e.g. 'Co60' or 'Co60m'."""
    _, _, isomer = parse_nuclide(nuclide)
    levels = get_levels(nuclide, session)
    return float(select_isomer(levels, isomer)["half_life_seconds"])


def get_gamma_energies_intensities(nuclide, session=None):
    """DataFrame with 'energy' (keV) and 'intensity' columns, by energy."""
    gammas = get_gammas(nuclide, session)
    cols = [c for c in ("energy", "intensity") if c in gammas.columns]
    return gammas[cols].sort_values("energy").reset_index(drop=True)
