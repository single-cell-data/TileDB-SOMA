import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from subprocess import PIPE
from sys import stderr
from typing import Any, Dict, Optional

import somacore

import tiledbsoma

from .context_generator import host_context

# import context_generator
from .data import ProfileData, ProfileDB, S3ProfileDB

GNU_TIME_FORMAT = (
    'Command being timed: "%C"\n'
    "User time (seconds): %U\n"
    "System time (seconds): %S\n"
    "Percent of CPU this job got: %P\n"
    "Elapsed (wall clock) time (seconds): %e\n"
    "Average shared text size (kbytes): %X\n"
    "Average unshared data size (kbytes): %D\n"
    "Average stack size (kbytes): %p\n"
    "Average total size (kbytes): %K\n"
    "Maximum resident set size (kbytes): %M\n"
    "Average resident set size (kbytes): %t\n"
    "Major (requiring I/O) page faults: %F\n"
    "Minor (reclaiming a frame) page faults: %R\n"
    "Voluntary context switches: %w\n"
    "Involuntary context switches: %c\n"
    "Swaps: %W\n"
    "File system inputs: %I\n"
    "File system outputs: %O\n"
    "Socket messages sent: %s\n"
    "Socket messages received: %r\n"
    "Signals delivered: %k\n"
    "Page size (bytes): %Z\n"
    "Exit status: %x"
)

GNU_TIME_OUTPUT_REGEXP = re.compile(
    r".*Command being timed: \"(?P<command>.+)\"\n"
    r"User time \(seconds\): (?P<user_time_sec>.+)\n"
    r"System time \(seconds\): (?P<system_time_sec>.+)\n"
    r"Percent of CPU this job got: (?P<pct_of_cpu>.+)%\n"
    r"Elapsed \(wall clock\) time \(seconds\): (?P<elapsed_time_sec>.+)\n"
    r"Average shared text size \(kbytes\): (?P<avg_shared_text_sz_kb>.+)\n"
    r"Average unshared data size \(kbytes\): (?P<avg_unshared_text_sz_kb>.+)\n"
    r"Average stack size \(kbytes\): (?P<avg_stack_sz_kb>.+)\n"
    r"Average total size \(kbytes\): (?P<avg_total_sz_kb>.+)\n"
    r"Maximum resident set size \(kbytes\): (?P<max_res_set_sz_kb>.+)\n"
    r"Average resident set size \(kbytes\): (?P<avg_res_set_sz_kb>.+)\n"
    r"Major \(requiring I/O\) page faults: (?P<major_page_faults>.+)\n"
    r"Minor \(reclaiming a frame\) page faults: (?P<minor_page_faults>.+)\n"
    r"Voluntary context switches: (?P<voluntary_context_switches>.+)\n"
    r"Involuntary context switches: (?P<involuntary_context_switches>.+)\n"
    r"Swaps: (?P<swaps>.+)\n"
    r"File system inputs: (?P<file_system_inputs>.+)\n"
    r"File system outputs: (?P<file_system_outputs>.+)\n"
    r"Socket messages sent: (?P<socket_messages_sent>.+)\n"
    r"Socket messages received: (?P<socket_messages_received>.+)\n"
    r"Signals delivered: (?P<signals_delivered>.+)\n"
    r"Page size \(bytes\): (?P<page_size_bytes>.+)\n"
    r"Exit status: (?P<exit_status>.+)\n"
    r".*"
)

# parameterize as cmd line arg
TILEDB_STATS_FILE_PATH = "./tiledb_stats.json"


def build_profile_data(
    stderr_: str, stdout_: str, prof1: Optional[str], prof2: Optional[str]
) -> ProfileData:
    """Parse the time utility output to extract performance and memory metrics"""
    gnu_time_output_values = GNU_TIME_OUTPUT_REGEXP.search(stderr_)
    assert gnu_time_output_values

    gnu_time_output_values = gnu_time_output_values.groupdict()
    # cast all dict int values
    gnu_time_output_values.update(
        {k: int(v) for k, v in gnu_time_output_values.items() if v.isdigit()}
    )

    data = ProfileData(
        stdout=stdout_,
        stderr=stderr_,
        timestamp=datetime.utcnow().timestamp(),
        tiledb_stats=read_tiledb_stats_output(),
        somacore_version=somacore.__version__,
        tiledbsoma_version=tiledbsoma.__version__,
        host_context=host_context(),
        custom_out=[prof1, prof2],
        **gnu_time_output_values,
    )

    return data


def read_tiledb_stats_output() -> Dict[str, Any]:
    if not os.path.isfile(TILEDB_STATS_FILE_PATH):
        return {}

    with open(TILEDB_STATS_FILE_PATH, "r") as f:
        print("TileDB stats found", file=stderr)
        return json.load(f)


def main():
    data_columns = ", ".join([a for a in dir(ProfileData) if a[0] != "_"])
    parser = argparse.ArgumentParser(
        epilog=f"The list of collected metrics by the generic profiler: {data_columns}"
    )
    parser.add_argument(
        "command",
        help="The command and its arguments to be profiled (as quoted, single-argument)",
    )
    parser.add_argument(
        "db_path",
        help="FileDB Path",
    )
    parser.add_argument(
        "-t",
        "--gtime-cmd",
        required=False,
        default="/usr/bin/time",
        help="Path to the gtime command",
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

    args = parser.parse_args(sys.argv[1:])

    print(f"Command to be run: {args.command}", file=stderr)
    # Run the command, using `time -v` to get detailed memory and time"""
    p = subprocess.Popen(
        [args.gtime_cmd, "--format", GNU_TIME_FORMAT] + args.command.split(" "),
        stdout=PIPE,
        stderr=PIPE,
    )

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
    p_stdout, p_stderr = p.communicate()
    if p1 is not None:
        p1.wait()
    if p2 is not None:
        p2.wait()
    print("Profiler run done", file=stderr)
    print("Profiler run done", file=stderr)

    p_stdout = p_stdout.decode("utf-8")
    print(f"The benchmarked process output:\n {p_stdout}", file=stderr)
    # Parse the generated output from the time utility
    data: ProfileData = build_profile_data(
        p_stderr.decode("utf-8"), p_stdout, args.prof1_output, args.prof2_output
    )
    # Add the run data to DB
    db: ProfileDB = S3ProfileDB(args.db_path)
    db_record_file = db.add(data)
    db.close()

    print(
        f"{data.command_key=}, {data.command=}, {data.exit_status=}, {db_record_file=}",
        file=stderr,
    )


if __name__ == "__main__":
    main()
