import argparse
import re
from collections import OrderedDict

import matplotlib.pyplot as plt

from data import FileBasedProfileDB, ProfileData


def extract_tiledb_data(data: ProfileData, metric: str) -> str:
    """Read the tile DB stats from a profile stored data
    Parse it and extract the expected metric"""
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "process",
        nargs="+",
        help="The main process name (and optional arguments) to be plotted",
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

    args = parser.parse_args()
    print(f"Process to be plotted: {args.process}")

    # extract process run data
    db = FileBasedProfileDB()
    data: ProfileData = db.find(" ".join(args.process))

    # prepare the extracted data for plotting
    pd = {}
    for profile_data in data:
        if args.metric:
            pd[profile_data.now] = getattr(profile_data, str(args.metric))
        elif args.tiledb_metric:
            pd[profile_data.now] = extract_tiledb_data(
                profile_data, str(args.tiledb_metric)
            )
            if not pd[profile_data.now]:
                raise Exception(f"TileDB stat {args.tiledb_metric} not found!")
        else:
            raise Exception("No metric or TileDB metric specified!")

    plot_data = OrderedDict(pd)
    plt.xticks(rotation=15, ha="right")
    plt.plot(plot_data.keys(), plot_data.values())
    plt.title(args.process)
    plt.show()


if __name__ == "__main__":
    main()
