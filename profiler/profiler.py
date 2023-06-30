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


def read_host_context() -> Dict[str, str]:
    """Run the context generator process and collects the JSON printed output"""
    context_process = subprocess.Popen(
        ["python", "context_generator.py"], stdout=PIPE, stderr=PIPE
    )
    stdout, stderr = context_process.communicate()
    if context_process.returncode != 0:
        raise RuntimeError(f"error reading host context information: {stderr.decode('utf-8')}")

    stats = json.loads(stdout.decode("utf-8"))
    return stats


def build_profile_data(output: str) -> ProfileData:
    """Parse the time utility output to extract performance and memory metrics"""
    print(f"OUTPUT=\n{output}\n")
    perf_match = re.search(r""".*Command being timed: \"(?P<command>.+)\"\n\s+User time \(seconds\): (?P<user_time_sec>.+)\n\s+System time \(seconds\): (?P<system_time_sec>.+)\n\s+Percent of CPU this job got: (?P<pct_of_cpu>.+)%\n\s+Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (?P<elapsed_time>.+)\n\s+Average shared text size \(kbytes\): (?P<avg_shared_text_sz_kb>.+)\n\s+Average unshared data size \(kbytes\): (?P<avg_unshared_text_sz_kb>.+)\n\s+Average stack size \(kbytes\): (?P<avg_stack_sz_kb>.+)\n\s+Average total size \(kbytes\): (?P<avg_total_sz_kb>.+)\n\s+Maximum resident set size \(kbytes\): (?P<max_res_set_sz_kb>.+)\n\s+Average resident set size \(kbytes\): (?P<avg_res_set_sz_kb>.+)\n\s+Major \(requiring I/O\) page faults: (?P<major_page_faults>.+)\n\s+Minor \(reclaiming a frame\) page faults: (?P<minor_page_faults>.+)\n\s+Voluntary context switches: (?P<voluntary_context_switches>.+)\n\s+Involuntary context switches: (?P<involuntary_context_switches>.+)\n\s+Swaps: (?P<swaps>.+)\n\s+File system inputs: (?P<file_system_inputs>.+)\n\s+File system outputs: (?P<file_system_outputs>.+)\n\s+Socket messages sent: (?P<socket_messages_sent>.+)\n\s+Socket messages received: (?P<socket_messages_received>.+)\n\s+Signals delivered: (?P<signals_delivered>.+)\n\s+Page size \(bytes\): (?P<page_size_bytes>.+)\n\s+Exit status: (?P<exit_status>.+)\n.*""", output)
    assert perf_match

    profile_data = perf_match.groupdict()
    # cast all dict int values
    profile_data.update({k: int(v) for k, v in profile_data.items() if v.isdigit()})

    # cast all dict elapsed time values to int (seconds)
    def elapsed_to_seconds(elapsed: str) -> float:
        elapsed_parts = re.match(r"(((\d+):)?((\d+):))?(\d+)(\.\d+)", elapsed).groups()
        return int(elapsed_parts[2] or "0") * 3600 + int(elapsed_parts[4] or "0") * 60 + int(elapsed_parts[5]) + float(elapsed_parts[6])

    profile_data["user_time_sec"] = elapsed_to_seconds(profile_data["user_time_sec"])
    profile_data["system_time_sec"] = elapsed_to_seconds(profile_data["system_time_sec"])
    profile_data["elapsed_time"] = elapsed_to_seconds(profile_data["elapsed_time"])

    tiledb_stats = None
    if os.path.isfile(TILEDB_STATS_FILE_PATH):
        with open(TILEDB_STATS_FILE_PATH, "r") as f:
            print("TileDB stats found")
            tiledb_stats = f.read()
    # custom_out = [prof1, prof2]
    context: Dict[str, str] = read_host_context()
    data: ProfileData = dict(
        date=str(date.today()),
        now=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),
        tiledb_stats=tiledb_stats,
        somacore_version=somacore.__version__,
        tiledbsoma_version=tiledbsoma.__version__,
        context=context,
        **perf_match.groupdict()
    )
    return data


def main():
    parser = argparse.ArgumentParser(
        epilog="""The list of collected metrics by the generic profiler:
         
         process: The process and its parameters to be profile\n
         custom_out: list of custom profilers to be stored\n
         date\n
         time\n
           
         rt: Real time\n
         ut: User time\n
         st: System time\n
                        
         max_set_size\n
         page_reclaims\n
         page_faults\n
         cycles_elapsed\n
         peak_memory\n
         tiledb_stats\n
         somacore_version\n
         tiledbsoma_version\n
                        
         uname: uname -a\n
         total_virtual_mem\n
         total_physical_mem\n
         swap_mem\n 
         cpu_count\n
         python_version\n
             """
    )
    parser.add_argument(
        "process", nargs="+", help="The main process and its arguments to run"
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
    args = parser.parse_args()

    print(f"Process to be run: {args.process}")
    # Running the main process using time -v to get detailed memory and time"""
    p = subprocess.Popen(
        [args.gtime_cmd, "-v"] + args.process, stdout=PIPE, stderr=PIPE
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
    data: ProfileData = build_profile_data(stderr.decode("utf-8"))
    # Add the run data to DB
    db: ProfileDB = FileBasedProfileDB()
    db.add(data)
    print("Printing DB:\n")
    print(db)
    db.close()


if __name__ == "__main__":
    main()
