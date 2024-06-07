import json
import seaborn as sns
import pandas as pd
import glob
import numpy as np


def signif(x, p):
    return np.format_float_positional(
        x, precision=p, unique=False, fractional=False, trim="k"
    )


def extract_results(prefix, file):

    with open(prefix + file, "r") as f:
        d = json.load(f)
        dcg_df = pd.DataFrame(d["dcg_times"])
        timings = d.pop("dcg_times")
        if len(d["cut_times"]) < len(d["num_active_cuts"]):
            d["cut_times"] += [0]
        df = pd.DataFrame(d)
        df["file"] = file
        df = pd.concat([df, dcg_df], axis=1)
    return df


def extract_results_from_file_prefix(prefix):

    files = glob.glob(prefix + "*")
    files = [
        file.split("\\")[-1]
        for file in files
        if "precompile" not in file and "logs" not in file
    ]
    dfs = []
    for file in files:
        dfs.append(extract_results(prefix, file))

    return pd.concat(dfs).reset_index(drop=True)

def cuts_at_root_node(prefix='data\\experiment_outputs\\cuts_at_root_node\\'):
    cut_progress = extract_results_from_file_prefix(prefix)

    ## GET TIMINGS
    time_cols = ["master_problem", "fire_subproblems", "crew_subproblems", "cut_times"]
    aggs = {i : "sum" for i in time_cols}
    aggs["times"] = "max"
    time_df = cut_progress.groupby("file", as_index=False).agg(aggs)
    time_df["num_fires"] = time_df["file"].apply(lambda x: x.split("_")[-1][:-5]).astype(int)
    time_df.sort_values(["num_fires", "file"])


    ## GET BOUND PROGRESS
    min_lb = cut_progress.groupby("file")["objectives"].transform("min")
    cut_progress["lower_bound_increase"] = cut_progress["objectives"] / min_lb - 1
    cut_progress["strategy"] = cut_progress["file"].apply(lambda x: "_".join(x.split("_")[:-3]))
    cut_progress["num_fires"] = cut_progress["file"].apply(lambda x: x.split("_")[-1][:-5]).astype(int)
    cut_progress.rename(columns={"times" : "time (secs)"}, inplace=True)

    ## MAKE PLOT
    fig = sns.relplot(
        data=cut_progress, x="time (secs)", 
        y="lower_bound_increase",
        col="num_fires",
        hue="strategy",
        kind="line",
        facet_kws={'sharey': False, 'sharex': False}
    )

    fig.figure.savefig("experiments\\figures\\cuts_at_root_node.png") 

cuts_at_root_node()