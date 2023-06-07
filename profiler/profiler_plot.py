import argparse
import re
from collections import OrderedDict
from typing import Dict, List, Union

import matplotlib.pyplot as plt
import pandas as pd

from data import FileBasedProfileDB, ProfileData, improve_profileDB_key


def collect_tiledb_stats(data: ProfileData) -> Dict[str, Union[int, float]]:
    """Extract all TileDB stats as dictionary"""
    result = {}
    tiledb_stats = data.tiledb_stats
    lines = tiledb_stats.split("\n")
    for idx, line in enumerate(lines):
        perf_match = re.match(r"\s+\"(.+)\": (\d+\.\d+)\s*,", line)
        if perf_match:
            value = float(perf_match.groups()[1])
            metric = perf_match.groups()[0]
            result[metric] = value
        perf_match = re.match(r"\s+\"(.*)\": (\d+)\s*,", line)
        if perf_match:
            value = int(perf_match.groups()[1])
            metric = perf_match.groups()[0]
            result[metric] = value
    return result


def extract_tiledb_data(data: ProfileData, metric: str) -> Union[int, float, None]:
    """Read the tile DB stats from a profile stored data
    Parse it and extract the expected metric
    """
    tiledb_stats = getattr(data, str("tiledb_stats"))
    lines = tiledb_stats.split("\n")
    for idx, line in enumerate(lines):
        perf_match = re.match(rf"\s+\"{metric}\": (\d+\.\d+)\s*,", line)
        if perf_match:
            value = float(perf_match.groups()[0])
            print(f"value for tiledb_stats {metric} = {value}")
            return value
        perf_match = re.match(rf"\s+\"{metric}\": (\d+)\s*,", line)
        if perf_match:
            value = int(perf_match.groups()[0])
            return value
    return None


def extract_context_data(data: ProfileData, metric: str) -> Union[int, float, None]:
    """Read the context data from a profile stored data and return the expected metric"""
    context = data.context
    if metric in context.keys():
        return context[metric]
    else:
        raise Exception(f"context does not have the following metric {metric}")
    return None


def create_pandas_df(data: List[ProfileData]) -> pd.DataFrame:
    """Create panda datat frame for all the runs of a given process
    Columns are metric names and rows are the runs
    """
    table = []
    # collecting metrics
    tmp = data[0]
    metrics = [m for m in dir(tmp) if not m.startswith("__")]
    tdb_metrics = collect_tiledb_stats(tmp).keys()
    context_metrics = tmp.context.keys()
    all_metrics = metrics + list(tdb_metrics) + list(context_metrics)

    # collecting data
    for run_data in data:
        row = (
            [getattr(run_data, metric) for metric in metrics]
            + [collect_tiledb_stats(run_data)[key] for key in tdb_metrics]
            + [run_data.context[key] for key in context_metrics]
        )
        table.append(row)
    row_labels = [runData.now for runData in data]

    # Creating dataframe
    df = pd.DataFrame(table, columns=all_metrics, index=row_labels)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "process",
        nargs="+",
        help="The main process name (and optional arguments) to be plotted",
    )
    parser.add_argument(
        "-pd",
        "--pandas",
        required=False,
        help="create pandas dataframe and prints it",
        action="store_true",
    )
    parser.add_argument(
        "-m",
        "--metric",
        required=False,
        help="the metric to be plotted",
    )
    parser.add_argument(
        "-tdm",
        "--tiledb_metric",
        required=False,
        help="the tiledb metric to be plotted",
    )
    parser.add_argument(
        "-cm",
        "--context_metric",
        required=False,
        help="the context metric to be plotted",
    )

    args = parser.parse_args()
    print(f"Process to be plotted: {args.process}")

    # extract process run data
    db = FileBasedProfileDB()
    data: ProfileData = db.find(" ".join(args.process))

    # prepare the extracted data for plotting
    pdata = {}
    for profile_data in data:
        if args.pandas:
            df: pd.DataFrame = create_pandas_df(data)
            file_name: str = (
                improve_profileDB_key(" ".join(args.process)) + "_df.pickle"
            )
            df.to_pickle(file_name)
            print(f"The panda DataFrame is stored in {file_name}")
            print(df)
            return
        elif args.metric:
            pdata[profile_data.now] = getattr(profile_data, str(args.metric))
        elif args.tiledb_metric:
            pdata[profile_data.now] = extract_tiledb_data(
                profile_data, str(args.tiledb_metric)
            )
            if not pd[profile_data.now]:
                raise Exception(f"TileDB stat {args.tiledb_metric} not found!")
        elif args.context_metric:
            pdata[profile_data.now] = extract_context_data(
                profile_data, str(args.context_metric)
            )
            if not pdata[profile_data.now]:
                raise Exception(f"TileDB stat {args.context_metric} not found!")
        else:
            raise Exception("No metric or TileDB or context metric specified!")

    plot_data = OrderedDict(pd)
    plt.xticks(rotation=15, ha="right")
    plt.plot(plot_data.keys(), plot_data.values())
    plt.title(args.process)
    plt.show()


if __name__ == "__main__":
    main()
