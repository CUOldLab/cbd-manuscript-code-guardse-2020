import pandas as pd
import numpy as np
import os



def combine_replicates(df, sensor):
    newdat = []; headers = []

    for i in range(1,df.shape[1]):
        cur_header = str(df.columns[i])
        if cur_header.endswith(".1"):
            cur_header = cur_header.replace(".1","_R2")
        else:
            cur_header = cur_header+"_R1"
        headers.append(sensor+"_"+str(cur_header))

    df["T"] = np.linspace(-2,72,df.shape[0])
    df.columns = ["Time"]+headers

    return(df)
        



root = "./sheets/"
fcount = 0
final_df = ""


for f in os.listdir(root):
    print(f)
    floc = root+f
    sensor = "_".join(f.replace(" ","_").split("_")[1:]).split(".")[0]
    if os.path.isfile(floc):
        df = pd.read_csv(floc)
        df = combine_replicates(df, sensor)
        if fcount==0:
            final_df = df
        else:
            final_df = final_df.merge(df, left_on="Time", right_on="Time", how="inner")

        fcount+=1

final_df.to_csv("results/combined_sensors.csv", index=False)      


