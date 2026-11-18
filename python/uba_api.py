import pandas as pd
import requests

# -----------------------------------------------------------------------
# api base url #
# --------------

BASE = "https://luftdaten.umweltbundesamt.de/api-proxy"


# -----------------------------------------------------------------------
# function to get json outputs #
# ------------------------------

def get_json(endpoint, params):
    response = requests.get(f"{BASE}/{endpoint}/json", params=params, timeout=30)
    response.raise_for_status()
    return response.json()


# -----------------------------------------------------------------------
# function to parse a table #
# ---------------------------

def parse_table(data, cols, numeric=(), dates=(), exp_len=None):
    # get rows
    rows = [dict(zip(cols, v)) for v in data.values()
            if exp_len is None or len(v) == exp_len]
    df = pd.DataFrame(rows)

    # fix types
    for c in numeric:
        df[c] = pd.to_numeric(df[c])
    for c in dates:
        df[c] = pd.to_datetime(df[c])

    return df


