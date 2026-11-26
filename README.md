# nndc_ensdf_data_getter

Python helpers for retrieving data from the NNDC ENSDF API. The module
`ensdf_api.py` uses pandas DataFrames for tabular results.

`get_levels()` streams adopted level records from `/levels/export`, requesting
levels with half-lives and energies above 0.001 keV. The export is global and
is not filtered to a particular nuclide. To extract positive half-lives for
excited levels:

```python
import ensdf_api

excited_levels = ensdf_api.get_excited_state_half_lives(max_rows=100_000)
```

The resulting DataFrame contains `levelId`, `datasetId`, `energy_keV`,
`energy_unit`, `halflife_s`, and `spinParity`. Provide a list of level records
to `get_excited_state_half_lives(levels=...)` to filter already retrieved data
without making another API request.
