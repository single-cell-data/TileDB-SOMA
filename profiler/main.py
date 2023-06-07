import argparse
import json
import os
import re
import subprocess
from datetime import date, datetime
from subprocess import PIPE
from typing import Dict, Optional

import somacore

import tiledbsoma
from data import FileBasedProfileDB, ProfileData, ProfileDB

TILEDB_STATS_FILE_PATH = "./tiledb_stats.txt"


def read_context() -> Dict[str, str]:
    """Run the context generator process and collects the JSON printed output"""
    context_process = subprocess.Popen(
        ["python", "context_generator.py"], stdout=PIPE, stderr=PIPE
    )
    stdout, stderr = context_process.communicate()
    stats = json.loads(stdout.decode("utf-8"))
    return stats


def build_profile_data(
    process: list[str], output: str, prof1: Optional[str], prof2: Optional[str]
) -> ProfileData:
    """Parse the time utility output to extract performance and memory metrics"""
    lines = output.split("\n")
    lines = lines[-20:]
    for idx, line in enumerate(lines):
        perf_match = re.match(
            "\s+(\d+\.\d+) real\s+(\d+\.\d+) user\s+(\d+\.\d+) sys\s*", line
        )
        if perf_match:
            real_time = float(perf_match.groups()[0])
            usr_time = float(perf_match.groups()[1])
            sys_time = float(perf_match.groups()[2])
            print(
                f"real_time = {real_time} user_time = {usr_time} sys_time = {sys_time}"
            )
        max_resident_set_size_match = re.match(
            "\s+(\d+)\s+maximum resident set size\s*", line
        )
        if max_resident_set_size_match:
            max_resident_set_size = int(max_resident_set_size_match.groups()[0])
            print(f"max_resident_set_size = {max_resident_set_size}")
        page_reclaims_match = re.match("\s+(\d+)\s+page reclaims\s*", line)
        if page_reclaims_match:
            page_reclaims = int(page_reclaims_match.groups()[0])
            print(f"page_reclaims = {page_reclaims}")
        page_faults_match = re.match("\s+(\d+)\s+page faults\s*", line)
        if page_faults_match:
            page_faults = int(page_faults_match.groups()[0])
            print(f"page_faults = {page_faults}")
        cycles_elapsed_match = re.match("\s+(\d+)\s+cycles elapsed\s*", line)
        if cycles_elapsed_match:
            cycles_elapsed = int(cycles_elapsed_match.groups()[0])
            print(f"cycles_elapsed = {cycles_elapsed}")
        peak_memory_footprint_match = re.match(
            "\s+(\d+)\s+peak memory footprint\s*", line
        )
        if peak_memory_footprint_match:
            peak_memory_footprint = int(peak_memory_footprint_match.groups()[0])
            print(f"peak_memory_footprint = {peak_memory_footprint}")
    tiledb_stats = None
    if os.path.isfile(TILEDB_STATS_FILE_PATH):
        with open(TILEDB_STATS_FILE_PATH, "r") as f:
            print("TileDB stats found")
            tiledb_stats = f.read()
    custom_out = [prof1, prof2]
    context: Dict[str, str] = read_context()
    data: ProfileData = ProfileData(
        process=" ".join(process),
        custom_out=custom_out,
        rt=real_time,
        ut=usr_time,
        st=sys_time,
        max_set_size=max_resident_set_size,
        page_reclaims=page_reclaims,
        page_faults=page_faults,
        cycles_elapsed=cycles_elapsed,
        peak_memory=peak_memory_footprint,
        date=str(date.today()),
        now=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),
        tiledb_stats=tiledb_stats,
        somacore_version=somacore.__version__,
        tiledbsoma_version=tiledbsoma.__version__,
        context=context,
    )
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "process", nargs="+", help="The main process and its arguments to run"
    )
    parser.add_argument(
        "-p1",
        "--prof1",
        required=False,
        help="Second profiler and its args to use (like python and R) to produce flamegraph of calls",
    )
    parser.add_argument(
        "-p1out",
        "--prof1_output",
        required=False,
        help="The flamegraph output produced by prof1",
    )
    parser.add_argument(
        "-p2",
        "--prof2",
        required=False,
        help="Third profiler and its args to use (like python and R) to produce flamegraph of calls",
    )
    parser.add_argument(
        "-p2out",
        "--prof2_output",
        required=False,
        help="The flamegraph output produced by prof2",
    )
    args = parser.parse_args()

    print(f"Process to be run: {args.process}")
    # Running the main process using time -v to get detailed memory and time"""
    p = subprocess.Popen(
        ["/usr/bin/time", "-al"] + args.process, stdout=PIPE, stderr=PIPE
    )

    print(f"Running main process, PID = {p.pid}")
    # Running additional profilers to extract flame graphs for the run
    p1 = None
    p2 = None
    if args.prof1 is not None:
        if args.prof1_output is not None:
            command = args.prof1.replace("PID", str(p.pid))
            p1 = subprocess.Popen(command.split())
            # potential example of a custom profiler here is:
            # p1 = subprocess.Popen(["py-spy", "record", "-o", "profile.svg", "--", "python", "tests/objects.py"])
            print(f"Running second profiler {command.split()}")
        else:
            print(
                f"Second profiler {args.prof1} missing output flamegraph file location"
            )

    if args.prof2 is not None:
        if args.prof2_output is not None:
            command = args.prof2.replace("PID", str(p.pid))
            p2 = subprocess.Popen(command.split())
            print(f"Running third profiler {command}")
        else:
            print(
                f"Third profiler {args.prof2} missing output flamegraph file location"
            )

    stdout, stderr = p.communicate()
    if p1 is not None:
        p1.wait()
    if p2 is not None:
        p2.wait()

    o: str = stdout.decode("utf-8")
    print(f"The benchmarked process output:\n {o}")
    # Parse the generated output from the time utility
    data: ProfileData = build_profile_data(
        args.process, stderr.decode("utf-8"), args.prof1_output, args.prof2_output
    )
    # Add the run data to DB
    db: ProfileDB = FileBasedProfileDB()
    db.add(data)
    print("Printing DB:\n")
    print(db)
    db.close()


if __name__ == "__main__":
    main()
