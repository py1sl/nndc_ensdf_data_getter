"""Helpers for retrieving data from the NNDC ENSDF API."""

import json
import re
from urllib.parse import urljoin

import pandas as pd
import requests

BASE_URL = "https://www.nndc.bnl.gov/ensdf-api/"
TIMEOUT = 60
GAMMA_EXPORT_PATH = "/gammas/export"
GAMMA_EXPORT_PARAMS = {"document-type": "adopted"}
LEVEL_EXPORT_PATH = "/levels/export"
LEVEL_EXPORT_PARAMS = {
    "document-type": "adopted",
    "energy-min": 0.001,
    "has-halflife": True,
}

_NUCLIDE_RE = re.compile(r"^\s*(\d+)?\s*([A-Za-z]{1,2})[-\s]?(\d+)?\s*(m\d*)?\s*$")

session = requests.Session()
session.headers.update({"Accept": "application/json"})


def parse_nuclide(name):
    """Parse Co60, Co-60, 60Co, Co60m, or Ag110m2."""
    match = _NUCLIDE_RE.match(name)
    if not match:
        raise ValueError(f"Cannot parse nuclide name: {name!r}")
    mass_prefix, symbol, mass_suffix, meta = match.groups()
    mass = mass_prefix or mass_suffix
    if not mass:
        raise ValueError(f"No mass number in nuclide name: {name!r}")
    isomer = int(meta[1:]) if meta and len(meta) > 1 else int(bool(meta))
    return symbol.capitalize(), int(mass), isomer


def api_url(path):
    """Build an API URL from an absolute- or relative-style API path."""
    return urljoin(BASE_URL, path.lstrip("/"))


def request_json(method, path, *, params=None, json_body=None, timeout=TIMEOUT,
                 http_session=None):
    """Make an API request and return its decoded JSON response."""
    client = http_session or session
    response = client.request(
        method,
        api_url(path),
        params=params,
        json=json_body,
        timeout=timeout,
    )
    if not response.ok:
        try:
            detail = response.json()
        except ValueError:
            detail = response.text
        raise RuntimeError(
            f"{response.status_code} {response.url}\n"
            f"{json.dumps(detail, indent=2)}"
        )
    return response.json()


def get(path, params=None, *, http_session=None):
    """Perform a GET request against the API."""
    return request_json("GET", path, params=params, http_session=http_session)


def results(page):
    """Return records from an API page, or an empty list if absent."""
    return page.get("results", [])


def frame(page_or_rows, columns=None):
    """Convert API results or row records into a pandas DataFrame."""
    rows = (
        results(page_or_rows)
        if isinstance(page_or_rows, dict) and "results" in page_or_rows
        else page_or_rows
    )
    table = rows.copy() if isinstance(rows, pd.DataFrame) else pd.json_normalize(rows)
    if columns:
        table = table.reindex(columns=columns)
    return table


def get_decays(nuclide, params=None, *, http_session=None):
    """Get decay records for a nuclide, returned as a DataFrame."""
    symbol, mass, _ = parse_nuclide(nuclide)
    path = f"/nuclides/{mass}{symbol}/decays"
    return frame(get(path, params=params, http_session=http_session))


def stream_ndjson(path, *, params, max_rows=None, timeout=300,
                  http_session=None):
    """Stream newline-delimited JSON records from the API."""
    if max_rows is not None and max_rows <= 0:
        return []
    client = http_session or session
    rows = []
    with client.get(
        api_url(path),
        params=params,
        stream=True,
        timeout=timeout,
        headers={"Accept": "application/x-ndjson"},
    ) as response:
        response.raise_for_status()
        for line in response.iter_lines():
            if not line:
                continue
            rows.append(json.loads(line))
            if max_rows is not None and len(rows) >= max_rows:
                break
    return rows


def get_gammas(max_rows=None, *, params=None, timeout=300, http_session=None):
    """Get adopted gamma records from the API's NDJSON export."""
    export_params = dict(GAMMA_EXPORT_PARAMS)
    if params:
        export_params.update(params)
    return frame(
        stream_ndjson(
            GAMMA_EXPORT_PATH,
            params=export_params,
            max_rows=max_rows,
            timeout=timeout,
            http_session=http_session,
        )
    )


def get_levels(max_rows=None, *, params=None, timeout=300, http_session=None):
    """Get adopted levels with half-lives from the levels NDJSON export.

    The export is not nuclide-filtered; it can contain many records.
    """
    export_params = dict(LEVEL_EXPORT_PARAMS)
    if params:
        export_params.update(params)
    return frame(
        stream_ndjson(
            LEVEL_EXPORT_PATH,
            params=export_params,
            max_rows=max_rows,
            timeout=timeout,
            http_session=http_session,
        )
    )


def get_excited_state_half_lives(levels=None, max_rows=None, *,
                                 params=None, http_session=None):
    """Return adopted excited levels with positive half-lives as a DataFrame.

    Columns include ``levelId``, ``datasetId``, ``energy_keV``, ``halflife_s``,
    and ``spinParity``. If ``levels`` is omitted, the API's global levels
    export is streamed; it does not select a particular nuclide.
    """
    records = (
        get_levels(max_rows, params=params, http_session=http_session)
        if levels is None
        else frame(levels)
    )
    rows = []
    for level in records.to_dict(orient="records"):
        energy = level.get("energy")
        energy_value = level.get("energy.value")
        energy_unit = level.get("energy.unit")
        if isinstance(energy, dict):
            energy_value = energy.get("value")
            energy_unit = energy.get("unit")
        dataset = level.get("datasetSummary")
        dataset_id = level.get("datasetId", level.get("datasetSummary.id"))
        if isinstance(dataset, dict):
            dataset_id = dataset.get("id")
        spin_parity = level.get(
            "spinParity", level.get("spinParity.stringRepresentation")
        )
        if isinstance(spin_parity, dict):
            spin_parity = spin_parity.get("stringRepresentation")
        rows.append({
            "levelId": level.get("id", level.get("levelId")),
            "datasetId": dataset_id,
            "energy_keV": energy_value,
            "energy_unit": energy_unit,
            "halflife_s": level.get("halflifeSeconds", level.get("halflife_s")),
            "spinParity": spin_parity,
        })
    table = pd.DataFrame(
        rows,
        columns=[
            "levelId", "datasetId", "energy_keV", "energy_unit",
            "halflife_s", "spinParity",
        ],
    )
    table["energy_keV"] = pd.to_numeric(table["energy_keV"], errors="coerce")
    table["halflife_s"] = pd.to_numeric(table["halflife_s"], errors="coerce")
    return table.loc[
        table["energy_unit"].eq("keV")
        & table["energy_keV"].notna()
        & table["energy_keV"].gt(0)
        & table["halflife_s"].gt(0)
    ].reset_index(drop=True)


def get_gamma_energies_intensities(gammas=None, max_rows=None,
                                   *, http_session=None):
    """Return gamma energy (keV) and available intensity data as a DataFrame.

    If ``gammas`` is omitted, adopted gamma records are streamed from the API.
    The export shown in the API example is not nuclide-filtered.
    """
    records = get_gammas(max_rows, http_session=http_session) if gammas is None else frame(gammas)
    rows = []
    for record in records.to_dict(orient="records"):
        energy = record.get("energy")
        energy_value = record.get("energy.value")
        energy_unit = record.get("energy.unit")
        if isinstance(energy, dict):
            energy_value = energy.get("value")
            energy_unit = energy.get("unit")
        if energy_value is None:
            continue
        if energy_unit != "keV":
            continue
        intensity = record.get("intensity")
        if intensity is None:
            intensity = record.get("intensity.value")
        if not isinstance(intensity, dict):
            intensity = intensity if intensity is not None else record.get("relativeIntensity")
        if isinstance(intensity, dict):
            intensity = intensity.get("value")
        rows.append({"energy": energy_value, "intensity": intensity})
    return pd.DataFrame(rows, columns=["energy", "intensity"])


def get_half_life(nuclide, *, params=None, http_session=None):
    """Return an exposed half-life field from the nuclide's decay records.

    The API example provided documents the decays route but not a half-life
    field or metastable path syntax, so this helper does not infer either.
    """
    _, _, isomer = parse_nuclide(nuclide)
    if isomer:
        raise NotImplementedError(
            "The API example does not document how to address metastable nuclides"
        )
    decays = get_decays(nuclide, params=params, http_session=http_session)
    if decays.empty:
        raise LookupError(f"No decay records found for {nuclide}")
    half_life_fields = (
        "halfLife.value", "half-life.value", "half_life.value",
        "halfLife", "half-life", "half_life",
    )
    for field in half_life_fields:
        if field in decays:
            values = decays[field].dropna()
            if not values.empty:
                return values.iloc[0]
    raise LookupError("The decay records do not include a half-life field")
