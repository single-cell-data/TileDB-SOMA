import errno
import json
import os
import re
from abc import ABC, abstractmethod
from typing import Dict, List, Optional

import attr


@attr.define
class ProfileData:
    """This class represents the data stored per run"""
    command: str
    datetime: str
    tiledb_stats: Optional[str]
    somacore_version: str
    tiledbsoma_version: str
    host_context: Dict[str, str]
    user_time_sec: float
    system_time_sec: float
    pct_of_cpu: float
    elapsed_time_sec: float
    avg_shared_text_sz_kb: int
    avg_unshared_text_sz_kb: int
    avg_stack_sz_kb: int
    avg_total_sz_kb: int
    max_res_set_sz_kb: int
    avg_res_set_sz_kb: int
    major_page_faults: int
    minor_page_faults: int
    voluntary_context_switches: int
    involuntary_context_switches: int
    swaps: int
    file_system_inputs: int
    file_system_outputs: int
    socket_messages_sent: int
    socket_messages_received: int
    signals_delivered: int
    page_size_bytes: int
    exit_status: int
    custom_out: List[Optional[str]]


def extract_key_from_filename(filename: str) -> str:
    """Extracts DB key for stats from the corresponding filename"""
    match = re.match(r"(.+)\.run", filename)
    assert match
    return match.groups()[0]


PROFILE_PATH = "./profiling_runs"


def command_key(command: str) -> str:
    """Remove space characters from profileDB command keys."""
    return command.replace(" ", "_").replace("/", "_").replace(".", "_")


class ProfileDB(ABC):
    """Base class for profiler runs database"""

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def find(self, command):
        pass

    @abstractmethod
    def add(self, data: ProfileData):
        pass

    @abstractmethod
    def close(self):
        pass


class FileBasedProfileDB(ProfileDB):
    """Represents a file-based implementation of a ProfileDB
    runs database. Each run is stored as a separate file under a subdirectory structured as `<command>/<timestamp>`.
    """

    def __init__(self, path: str = PROFILE_PATH):
        self.path = path
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __str__(self):
        result = ""
        if os.path.exists(self.path):
            for command in os.listdir(self.path):
                result += (
                    f"{command}: "
                    + str(len(os.listdir(f"{self.path}/{command}")))
                    + "\n"
                )
            return result
        return ""

    def find(self, command) -> List[ProfileData]:
        key = command_key(command)
        if not os.path.exists(f"{self.path}/{key}"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), key)
        dir_list = os.listdir(f"{self.path}/{key}")
        result = []
        for filename in dir_list:
            with open(f"{self.path}/{key}/{filename}", "r") as file:
                result.append(ProfileData(**json.load(file)))
        return result

    def add(self, data: ProfileData):
        key = command_key(data.command)
        os.makedirs(f"{self.path}/{key}", exist_ok=True)
        key2 = data.datetime
        filename = f"{self.path}/{key}/{key2}.run"
        with open(filename, "w") as f:
            json.dump(attr.asdict(data), f)

    def close(self):
        pass
