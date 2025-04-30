# little script to extract the data I want from a json
# and put it in a csv to plot in R
# simply because I already know how to process a json in python

import pandas as pd
import json

with open("visualization.vl.json", "r") as f:
    data = json.load(f)
if isinstance(data, dict):
    print("Top-level keys:", data.keys())

print(data)

df = pd.DataFrame(data["datasets"]['data-a059397d6b190f2bbcfca262e538e064'])
df.to_csv("qiime_da.csv", index=False)