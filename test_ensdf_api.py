import pytest
import pandas as pd

import ensdf_api as api


class FakeResponse:
    def __init__(self, payload=None, lines=None, ok=True, text=""):
        self.payload = payload
        self.lines = lines or []
        self.ok = ok
        self.text = text
        self.status_code = 400
        self.url = "https://example.test/"

    def json(self):
        if self.payload is None:
            raise ValueError("not JSON")
        return self.payload

    def raise_for_status(self):
        if not self.ok:
            raise RuntimeError("HTTP error")

    def iter_lines(self):
        return self.lines

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False


class FakeSession:
    def __init__(self, payload=None, lines=None):
        self.payload = payload
        self.lines = lines or []
        self.calls = []

    def request(self, method, url, **kwargs):
        self.calls.append((method, url, kwargs))
        return FakeResponse(self.payload)

    def get(self, url, **kwargs):
        self.calls.append(("GET", url, kwargs))
        return FakeResponse(lines=self.lines)


@pytest.mark.parametrize("name,expected", [
    ("Co60", ("Co", 60, 0)),
    ("Co-60", ("Co", 60, 0)),
    ("60Co", ("Co", 60, 0)),
    ("Co60m", ("Co", 60, 1)),
    ("Ag110m2", ("Ag", 110, 2)),
])
def test_parse(name, expected):
    assert api.parse_nuclide(name) == expected


def test_parse_invalid():
    with pytest.raises(ValueError):
        api.parse_nuclide("Co")


def test_get_decays_uses_nuclide_route_and_results():
    client = FakeSession({"results": [{"id": "decay-1", "energy": {"value": 1.2}}]})
    table = api.get_decays("Co60", {"limit": 200}, http_session=client)
    assert list(table["id"]) == ["decay-1"]
    method, url, kwargs = client.calls[0]
    assert method == "GET"
    assert url.endswith("/nuclides/60Co/decays")
    assert kwargs["params"] == {"limit": 200}


def test_request_error_includes_api_response_details():
    client = FakeSession({"detail": "invalid"}, lines=[])
    client.request = lambda *_args, **_kwargs: FakeResponse(
        payload={"detail": "invalid"}, ok=False
    )
    with pytest.raises(RuntimeError, match="invalid"):
        api.get("/nuclides/60Co", http_session=client)


def test_stream_ndjson_uses_export_accept_header_and_row_limit():
    client = FakeSession(lines=[b'{"energy":{"value":1}}', b'{"energy":{"value":2}}'])
    records = api.stream_ndjson(
        "/gammas/export", params={"document-type": "adopted"},
        max_rows=1, http_session=client,
    )
    assert records == [{"energy": {"value": 1}}]
    _, url, kwargs = client.calls[0]
    assert url.endswith("/gammas/export")
    assert kwargs["headers"] == {"Accept": "application/x-ndjson"}
    assert kwargs["params"] == {"document-type": "adopted"}


def test_gamma_energies_extracts_keV_and_nested_intensity():
    records = [
        {"energy": {"value": 1332.5, "unit": "keV"},
         "intensity": {"value": 99.98}},
        {"energy": {"value": 1.0, "unit": "MeV"},
         "intensity": {"value": 2.0}},
        {"energy": {"value": 1173.2, "unit": "keV"}},
    ]
    table = api.get_gamma_energies_intensities(records)
    assert list(table["energy"]) == [1332.5, 1173.2]
    assert list(table["intensity"][:1]) == [99.98]
    assert pd.isna(table["intensity"].iloc[1])


def test_gamma_energy_helper_accepts_normalized_api_frame():
    table = api.get_gamma_energies_intensities([
        {"energy.value": 1332.5, "energy.unit": "keV", "intensity.value": 99.98}
    ])
    assert table.to_dict(orient="records") == [
        {"energy": 1332.5, "intensity": 99.98}
    ]
