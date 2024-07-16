import json
import seaborn as sns
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("ggplot")

sns.set_context(
    rc={
        "font.size": 13,
        "axes.titlesize": 20,
        "axes.labelsize": 17,
    }
)


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


def extract_bb_results(prefix):

    files = glob.glob(prefix + "*")
    files = [
        file.split("\\")[-1]
        for file in files
        if "precompile" not in file and "logs" not in file
    ]
    dfs = []
    for file in files:
        with open(prefix + file, "r") as f:
            d = json.load(f)
            df = pd.DataFrame(d)
            df["file"] = file
            dfs.append(df)
    return pd.concat(dfs)


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


def dcg_root_node_stats(
    price_and_cut_path="data\\experiment_outputs\\cuts_at_root_node\\",
    heuristic_benchmark_path="data\\experiment_outputs\\triage_then_route\\",
    exact_benchmark_path="data\\experiment_outputs\\network_flow_direct\\",
    cut_strategy="cglp_adaptive",
):

    df = extract_results_from_file_prefix(price_and_cut_path)
    df = df[df["file"].apply(lambda x: cut_strategy in x)]

    ## GET TIMINGS
    time_cols = [
        "num_plans",
        "num_routes",
        "master_problem",
        "fire_subproblems",
        "crew_subproblems",
    ]
    aggs = {i: "first" for i in time_cols}
    aggs["times"] = "first"
    time_df = df.groupby("file", as_index=False).agg(aggs)
    time_df["num_fires"] = (
        time_df["file"].apply(lambda x: x.split("_")[-1][:-5]).astype(int)
    )
    time_df["num_crews"] = (time_df["num_fires"] * 10 / 3).astype(int)
    time_df.rename(
        columns={
            "times": "DCG",
            "fire_subproblems": "FSP",
            "crew_subproblems": "CSP",
            "master_problem": "RMP",
        },
        inplace=True,
    )

    time_df.sort_values(["num_fires", "file"])
    time_df["RMP"] = np.round(time_df["RMP"] + 1e-12, 4)
    time_df["FSP"] = np.round(time_df["FSP"] + 1e-12, 4)
    time_df["CSP"] = np.round(time_df["CSP"] + 1e-12, 4)
    time_df["DCG"] = np.round(time_df["DCG"] + 1e-12, 4)

    # with open(heuristic_benchmark_path + "output.json", "r") as f:
    #     d = json.load(f)
    #     heuristic_benchmark = pd.DataFrame(d).T.reset_index()
    #     heuristic_benchmark.columns = ["num_crews", "TTR_time", "TTR_cost"]
    #     heuristic_benchmark["num_crews"] = heuristic_benchmark["num_crews"].astype(int)
    # time_df = pd.merge(time_df, heuristic_benchmark, how="left", on="num_crews")

    with open(exact_benchmark_path + "output.json", "r") as f:
        d = json.load(f)
        exact_benchmark = pd.DataFrame(d["linear"]).T.reset_index()
        exact_benchmark.rename(
            columns={"index": "num_crews", "time": "NAT"}, inplace=True
        )
        exact_benchmark["num_crews"] = exact_benchmark["num_crews"].astype(int)

    time_df = pd.merge(time_df, exact_benchmark, how="left", on="num_crews")

    time_df = time_df.sort_values("num_fires")
    cols = [
        "num_fires",
        "num_crews",
        "RMP",
        "FSP",
        "CSP",
        "DCG",
        "NAT",
        "num_plans",
        "num_routes",
    ]
    time_df = time_df[cols]
    latex_table = time_df.to_latex(
        index=False, escape=True, float_format="%.3f", column_format="|cc|ccc|cc|cc|"
    )
    latex_table = (
        latex_table.replace("\\bottomrule", "\\hline")
        .replace("\\midrule", "\\hline")
        .replace("\\toprule", "\\hline")
    )

    with open("experiments\\figures\\gurobi_linear_relaxation_benchmark.txt", "w") as f:
        print(latex_table, file=f)


def cuts_at_root_node(prefix="data\\experiment_outputs\\cuts_at_root_node\\"):
    cut_progress = extract_results_from_file_prefix(prefix)

    ## GET BOUND PROGRESS
    cut_progress["cut_times"] = cut_progress.groupby("file")["cut_times"].transform(
        "cumsum"
    )
    min_lb = cut_progress.groupby("file")["objectives"].transform("min")
    cut_progress["Lower bound increase (prop.)"] = (
        cut_progress["objectives"] / min_lb - 1
    )
    cut_progress["Cut method"] = cut_progress["file"].apply(
        lambda x: "_".join(x.split("_")[:-3])
    )
    cut_progress["num_fires"] = (
        cut_progress["file"].apply(lambda x: x.split("_")[-1][:-5]).astype(int)
    )
    cut_progress.rename(columns={"cut_times": "Time (sec.) generating cuts"}, inplace=True)
    cut_progress = cut_progress[
        cut_progress["Cut method"].isin(
            ["gub_plus_strengthen", "cglp_adaptive", "gub_only"]
        )
    ]
    cut_progress["Cut method"].replace("cglp_adaptive", "CGLP", inplace=True)
    cut_progress["Cut method"].replace("gub_only", "GUB", inplace=True)
    cut_progress["Cut method"].replace("gub_plus_strengthen", "GUB+S", inplace=True)
    cut_progress = cut_progress[cut_progress["num_fires"].isin([9, 12, 15])]

    ## MAKE PLOT
    fig = sns.relplot(
        data=cut_progress,
        x="Time (sec.) generating cuts",
        y="Lower bound increase (prop.)",
        col="num_fires",
        hue="Cut method",
        palette=["red", "darkgrey", "dimgrey"],
        kind="line",
        facet_kws={"sharey": False, "sharex": False, "legend_out" : True},
        linewidth=3,
    )

    number_of_fires = sorted(list(cut_progress["num_fires"].unique()))
    for i, fire in enumerate(number_of_fires):
        fig.axes[0, i].set_title(f"{fire} Fires, {int(10 * fire / 3)} Crews")

        # sns.move_legend(
        #     fig,
        #     loc="upper center",
        #     bbox_to_anchor=(0.5, -0.05),
        #     fancybox=True,
        #     shadow=True,
        #     ncol=3,
        # )

    fig.figure.savefig("experiments\\figures\\cuts_at_root_node.png")


def value_of_cuts_and_branching(
    prefix="data\\experiment_outputs\\branch_price_and_cut\\",
):
    df = extract_bb_results(prefix)
    df["upper_bounds"] = df["upper_bounds"].fillna(np.inf)
    file_extract = df["file"].apply(lambda x: x.split("_"))
    param_list = file_extract.str[0]
    param_list = param_list.apply(lambda x: x.split("+"))
    df["Crews"] = param_list.str[0].astype(int)
    df["Line per crew"] = param_list.str[1].astype(int)
    df["Travel speed"] = param_list.str[2].astype(float)
    df["Cuts"] = file_extract.str[2] == "true"
    df["Branching strategy"] = file_extract.apply(
        lambda x: "_".join(x[x.index("branch") + 2 : x.index("heuristic")])
    )
    df["Heuristic"] = file_extract.str[-1] == "true.json"
    df["Pct. gap"] = (
        (df["upper_bounds"] - df["lower_bounds"]) / df["lower_bounds"] * 100
    )
    df["Pct. gap"] = np.round(df["Pct. gap"] + 1e-12, 4)
    df["Fires"] = (0.3 * df["Crews"]).astype(int)
    df["heuristic_times"] = df.groupby("file")["heuristic_times"].transform("cumsum")
    cols = [
        "Crews",
        "Fires",
        "Line per crew",
        "Travel speed",
        "Cuts",
        "Branching strategy",
        "Heuristic",
        "times",
        "heuristic_times",
        "explored_nodes",
        "Pct. gap",
        "upper_bounds",
        "lower_bounds",
    ]
    tbl = df.drop_duplicates("file", keep="last")[cols].reset_index(drop=True)
    still_gap = tbl["Pct. gap"] > 0
    tbl.loc[still_gap, "times"] = 1200.0
    tbl["times"] = tbl["times"].apply(lambda x: signif(x, 5))
    tbl["Pct. gap"] = tbl["Pct. gap"].apply(lambda x: signif(x, 3))
    tbl.rename(
        inplace=True,
        columns={
            "upper_bounds": "UB",
            "lower_bounds": "LB",
            "times": "Time",
            "explored_nodes": "Nodes",
        },
    )
    latex_table = tbl.to_latex(
        index=False, escape=True, float_format="%.4f", column_format="|ccccc|ccccc|"
    )
    latex_table = (
        latex_table.replace("\\bottomrule", "\\hline")
        .replace("\\midrule", "\\hline")
        .replace("\\toprule", "\\hline")
    )

    with open("experiments\\figures\\value_of_cuts_and_branching.txt", "w") as f:
        print(latex_table, file=f)


def cglp_methods(prefix="data\\experiment_outputs\\cuts_at_root_node\\"):
    cut_progress = extract_results_from_file_prefix(prefix)

    ## GET BOUND PROGRESS
    cut_progress["cut_times"] = cut_progress.groupby("file")["cut_times"].transform(
        "cumsum"
    )
    min_lb = cut_progress.groupby("file")["objectives"].transform("min")
    cut_progress["Lower bound increase (prop.)"] = (
        cut_progress["objectives"] / min_lb - 1
    )
    cut_progress["Cut method"] = cut_progress["file"].apply(
        lambda x: "_".join(x.split("_")[:-3])
    )
    cut_progress["num_fires"] = (
        cut_progress["file"].apply(lambda x: x.split("_")[-1][:-5]).astype(int)
    )
    cut_progress.rename(columns={"cut_times": "Time (sec.) generating cuts"}, inplace=True)
    cut_progress = cut_progress[
        cut_progress["Cut method"].isin(
            ["cglp_adaptive", "cglp_cutting_plane", "cglp_enumerate"]
        )
    ]
    cut_progress["Cut method"].replace("cglp_adaptive", "CGLP-A", inplace=True)
    cut_progress["Cut method"].replace("cglp_cutting_plane", "CGLP-CP", inplace=True)
    cut_progress["Cut method"].replace("cglp_enumerate", "CGLP-E", inplace=True)
    cut_progress = cut_progress[cut_progress["num_fires"].isin([9, 12, 15])]

    ## MAKE PLOT
    fig = sns.relplot(
        data=cut_progress,
        x="Time (sec.) generating cuts",
        y="Lower bound increase (prop.)",
        col="num_fires",
        hue="Cut method",
        kind="line",
        facet_kws={"sharey": False, "sharex": False, "legend_out" : True},
        linewidth=3,
    )

    number_of_fires = sorted(list(cut_progress["num_fires"].unique()))
    for i, fire in enumerate(number_of_fires):
        fig.axes[0, i].set_title(f"{fire} Fires, {int(10 * fire / 3)} Crews")

        # sns.move_legend(
        #     fig,
        #     loc="upper center",
        #     bbox_to_anchor=(0.5, -0.05),
        #     fancybox=True,
        #     shadow=True,
        #     ncol=3,
        # )

    fig.figure.savefig("experiments\\figures\\cglp_methods.png")


dcg_root_node_stats()
cglp_methods()
value_of_cuts_and_branching()
cuts_at_root_node()
