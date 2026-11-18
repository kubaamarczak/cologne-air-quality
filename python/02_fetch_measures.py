import time
from pathlib import Path

import pandas as pd

from uba_api import get_json

# -----------------------------------------------------------------------
# fetch air quality data of selected stations, scope & components #
# -----------------------------------------------------------------

# station ids as constants
CW, RK, CR, TS = 1108, 1114, 1254, 1255 # chorweiler, rodenkirchen, clevischer ring, turiner straße
STATIONS = [CW, RK, CR, TS]

# scope 
SCOPE_HOURLY = 2 # 1SMW

# components
PM10, NO2 = 1, 5

# function to get data per station, component & year
def fetch_one(station, component, year):
    # restrict with station, component & year
    params = {
        "date_from": f"{year}-01-01", "time_from": 1,
        "date_to": f"{year}-12-31", "time_to": 24,
        "station": station, "component": component,
        "scope": SCOPE_HOURLY,
    }
    
    # get measures as json
    r = get_json(endpoint="measures", params=params)

    # parse json
    cols = ["component", "scope", "value", "date_end", "index"]
    rows = [
        {"station": sta, "date_start": timestmp, **dict(zip(cols, vals))}
        for sta, series in r["data"].items()
        for timestmp, vals in series.items()
    ]
    
    df_part = pd.DataFrame(rows)
    df_part["value"] = pd.to_numeric(df_part["value"])
    df_part["date_start"] = pd.to_datetime(df_part["date_start"])
    df_part["date_end"] = df_part["date_start"] + pd.Timedelta(hours=1) # avoids wrong time format
    
    time.sleep(0.5) # less strain on api 
    print(f"{station} {component} {year}")

    return df_part

# function to iterate over stations, components & years and create a data frame
def fetch_measures(stations, components, years):
    results = [fetch_one(s, c, y) for s in stations for c in components for y in years]
    return pd.concat(results, ignore_index=True)

# fetch measures
df_measures = fetch_measures(stations=STATIONS, components=[PM10, NO2], years=range(2016, 2025))

# peek at the final data extraction 
print(df_measures)

# plausibility check
assert len(df_measures) == 631_296, f"unexpected number of rows: {len(df_measures)}"
assert not df_measures.duplicated(["station", "component", "date_start"]).any()

# save data frame as parquet file for further analysis
OUT_DIR = Path(__file__).parent / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)
df_measures.to_parquet(OUT_DIR / "cologne_measures_2016_2024.parquet")