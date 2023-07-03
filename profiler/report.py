import argparse
import re
from collections import OrderedDict
from sys import stderr
from typing import Dict, List, Union

import attr
import matplotlib.pyplot as plt
import pandas as pd

from .data import FileBasedProfileDB, ProfileData


def collect_tiledb_stats(data: ProfileData) -> Dict[str, Union[int, float]]:
    """Extract all TileDB stats as dictionary"""
    result = {}
    tiledb_stats = data.tiledb_stats

    if tiledb_stats is None:
        return result

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


def create_pandas_df(profile_datas: List[ProfileData]) -> pd.DataFrame:
    """Create pandas dataframe for all the runs of a given command
    Columns are metric names and rows are the runs
    """
    return pd.DataFrame.from_records(data=[attr.asdict(d) for d in profile_datas])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        nargs="+",
        help="The command to be profiled (and optional arguments) to be plotted",
    )
    parser.add_argument(
        "-j",
        "--json",
        required=False,
        help="Displays results as JSON (Pandas DataFrame \"columns\" format)",
        action="store_true",
    )
    parser.add_argument(
        "-m",
        "--metric",
        required=False,
        help="The metric to be plotted",
    )
    parser.add_argument(
        "-tdm",
        "--tiledb_metric",
        required=False,
        help="The tiledb metric to be plotted",
    )
    parser.add_argument(
        "-cm",
        "--context_metric",
        required=False,
        help="the context metric to be plotted",
    )

    args = parser.parse_args()
    print(f"Profiling command to be plotted: {args.command}", file=stderr)

    # extract profiling command run data
    db = FileBasedProfileDB()
    profile_datas: List[ProfileData] = db.find(" ".join(args.command))

    if args.json:
        output_as_json(profile_datas)
        return

    # prepare the extracted data for plotting
    plot_data = {}
    for profile_data in profile_datas:
        if args.metric:
            plot_data[profile_data.datetime] = getattr(profile_data, args.metric)
        elif args.tiledb_metric:
            plot_data[profile_data.datetime] = extract_tiledb_data(
                profile_data, str(args.tiledb_metric)
            )
            if not pd[profile_data.datetime]:
                raise RuntimeError(f"TileDB stat {args.tiledb_metric} not found!")
        elif args.context_metric:
            plot_data[profile_data.datetime] = extract_context_data(
                profile_data, str(args.context_metric)
            )
            if not plot_data[profile_data.datetime]:
                raise RuntimeError(f"TileDB stat {args.context_metric} not found!")
        else:
            raise RuntimeError("No metric or TileDB or context metric specified!")

    plot_data = OrderedDict(plot_data)
    plt.xticks(rotation=15, ha="right")
    plt.plot(plot_data.keys(), plot_data.values())
    plt.title(args.command)
    plt.show()


def output_as_json(profile_datas: List[ProfileData]) -> None:
    print(create_pandas_df(profile_datas).to_json(orient="columns"))


if __name__ == "__main__":
    main()
