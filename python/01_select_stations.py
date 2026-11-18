import pandas as pd

from uba_api import get_json, parse_table

# -----------------------------------------------------------------------
# get all air quality stations as json #
# --------------------------------------

stations = get_json(endpoint="stations", params={"lang": "en"})


# -----------------------------------------------------------------------
# first overview #
# ----------------

# print(stations.keys())
# print(stations["indices"])
# print(stations["data"])


# -----------------------------------------------------------------------
# parse json to data frame #
# --------------------------

df = parse_table(data=stations["data"], 
                 cols=stations["indices"], 
                 numeric=["station id", "station longitude", "station latitude", 
                          "network id", "station setting id", "station type id"],
                 dates=["station active from", "station active to"])
# print(df)


# -----------------------------------------------------------------------
# filter for cologne stations which are still active #
# ----------------------------------------------------

df_nrw_col = df[(df["station code"].str.startswith("DENW")) &
            (df["station active from"] <= "2016-01-01") & (df["station active to"].isna()) &
            (df["station city"] == "Köln")]
# print(df_nrw_col)
# print(df_nrw_col["station type name"].value_counts())


# -----------------------------------------------------------------------------------------------------
# get measures as json #
# ----------------------

meta = get_json(endpoint="meta", 
                params={"use": "measure", "date_from": "2016-01-01", "date_to": "2024-12-31", "lang": "en"})


# -----------------------------------------------------------------------
# link scopes, components & stations #
# ------------------------------------

# parse json from measure/limits
df_lim = parse_table(data=meta["limits"], 
                     cols=["scope id", "component id", "station id", "date started", "date ended"], 
                     numeric=["scope id", "component id", "station id"],
                     dates=["date started", "date ended"],
                     exp_len=5)
# print(df_lim)

# filter measures for cologne stations
cologne_ids = list(df_nrw_col["station id"])
df_cologne = df_lim[df_lim["station id"].isin(cologne_ids)]
# print(df_cologne)


# -----------------------------------------------------------------------
# build component & scope tables #
# --------------------------------

# components
df_components = parse_table(data=meta["components"],
                            cols=["component id", "component code", "symbol", "unit", "translated name"],
                            numeric=["component id"])
# print(df_components)

# scopes
df_scopes = parse_table(data=meta["scopes"],
                            cols=["scope id", "scope code", "time base", "time scope", "flag"],
                            numeric=["scope id"])
# print(df_scopes)

# join tables together
df_join_comp = pd.merge(df_cologne, df_components, on="component id", how="left")
df_joined = pd.merge(df_join_comp, df_scopes, on="scope id", how="left")
df_filtered = df_joined[["scope code", "symbol", "station id", "date started", "date ended"]]
# print(df_filtered)

# get more info about remaining stations
df_more = pd.merge(df_filtered, df_nrw_col, on="station id", how="left")
df_extended = df_more[["station id", "station name", "scope code", "symbol", "station type name", "date started", "date ended"]]
# print(df_extended)

# show selected stations, scope & components for analysis (scope: hourly average; components: PM₁₀, NO₂)
mask = (
      (df_extended["scope code"] == "1SMW")
      & df_extended["symbol"].isin(["PM₁₀", "NO₂"])
  )
print(df_extended[mask])