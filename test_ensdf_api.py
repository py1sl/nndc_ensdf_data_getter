import pytest

import ensdf_api as api


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self.payload


class FakeSession:
    def __init__(self, payload):
        self.payload = payload
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append((url, params))
        return FakeResponse(self.payload)


LEVELS = [
    {"energy": 58.6, "half_life_seconds": 628.0},
    {"energy": 0.0, "half_life_seconds": 1.663e8},
]


@pytest.mark.parametrize("name,expected", [
    ("Co60", ("Co", 60, 0)), ("Co-60", ("Co", 60, 0)), ("60Co", ("Co", 60, 0)),
    ("Co60m", ("Co", 60, 1)), ("Ag110m2", ("Ag", 110, 2)),
])
def test_parse(name, expected):
    assert api.parse_nuclide(name) == expected


def test_parse_invalid():
    with pytest.raises(ValueError):
        api.parse_nuclide("Co")


def test_half_life_ground_and_meta():
    s = FakeSession(LEVELS)
    assert api.get_half_life("Co60", s) == pytest.approx(1.663e8)
    assert api.get_half_life("Co60m", s) == pytest.approx(628.0)
    assert s.calls[0][1] == {"nuc": "60Co"}


def test_missing_isomer():
    with pytest.raises(LookupError):
        api.get_half_life("Co60m2", FakeSession(LEVELS))


def test_gammas():
    s = FakeSession({"data": [{"energy": 1332.5, "intensity": 99.98},
                              {"energy": 1173.2, "intensity": 99.85}]})
    df = api.get_gamma_energies_intensities("Co60", s)
    assert list(df["energy"]) == [1173.2, 1332.5]
