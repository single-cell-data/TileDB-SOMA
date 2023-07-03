import errno
import glob
import hashlib
import json
import os
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any

import attr


@attr.define
class ProfileData:
    """This class represents the data stored per run"""
    command: str
    timestamp: float
    tiledb_stats: Dict[str, Any]
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


DEFAULT_PROFILE_DB_PATH = "./profiling_db"


def _command_key(command: str) -> str:
    """Remove space characters from profileDB command keys."""
    return hashlib.md5(command.encode('utf-8')).hexdigest()


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

    def __init__(self, path: str = DEFAULT_PROFILE_DB_PATH):
        self.path = path
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __str__(self):
        result = ""
        if os.path.exists(self.path):
            for command_hash in glob.glob(self.path + "/*"):
                with open(os.path.join(command_hash, "command.txt"), "r") as f:
                    command = f.read()
                n_runs = len(glob.glob(os.path.join(command_hash, "*.json")))
                result += f"[{command_hash.split('/')[-1]}] \"{command}\": {n_runs} runs\n"
            return result
        return ""

    def find(self, command) -> List[ProfileData]:
        key = _command_key(command)
        if not os.path.exists(f"{self.path}/{key}"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), key)
        result = []
        for filename in glob.glob(f"{self.path}/{key}/*.json"):
            with open(filename, "r") as file:
                result.append(ProfileData(**json.load(file)))
        return result

    def add(self, data: ProfileData):
        key = _command_key(data.command)
        os.makedirs(f"{self.path}/{key}", exist_ok=True)
        with open(f"{self.path}/{key}/command.txt", "w") as f:
            f.write(data.command.strip())

        key2 = data.timestamp

        filename = f"{self.path}/{key}/{key2}.json"
        with open(filename, "w") as f:
            json.dump(attr.asdict(data), f)

    def close(self):
        pass
