import errno
import json
import os
import re
from abc import ABC, abstractmethod
from typing import Dict, List, Any

"""This class represents the data stored per run"""


ProfileData = Dict[str, Any]


def extract_key_from_filename(filename: str) -> str:
    """Extracts DB key for stats from the corresponding filename"""
    match = re.match("(.+)\.run", filename)
    assert match
    return match.groups()[0]


PROFILE_PATH = "./profiling_runs"


def improve_profileDB_key(process: str) -> str:
    """Remove space characters from profileDB keys."""
    name: str = process.replace(" ", "_").replace("/", "_").replace(".", "_")
    print(f"Profiler key = {name}")
    return name


class ProfileDB(ABC):
    """Base class for profiler runs database"""

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def find(self, process):
        pass

    @abstractmethod
    def add(self, data: ProfileData):
        pass

    @abstractmethod
    def close(self):
        pass


class FileBasedProfileDB(ProfileDB):
    """Represents a file-based implementation of a ProfileDB
    runs database. Each run is stored as a separate file under a subdirectory structured as `<process_name>/<timestamp>`.
    """

    def __init__(self, path: str = PROFILE_PATH):
        self.path = path
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __str__(self):
        result = ""
        if os.path.exists(self.path):
            for process in os.listdir(self.path):
                result += (
                    f"{process}: "
                    + str(len(os.listdir(f"{self.path}/{process}")))
                    + "\n"
                )
            return result
        return ""

    def find(self, process) -> List[ProfileData]:
        key = improve_profileDB_key(process)
        if not os.path.exists(f"{self.path}/{key}"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), key)
        dir_list = os.listdir(f"{self.path}/{key}")
        result = []
        for filename in dir_list:
            with open(f"{self.path}/{key}/{filename}", "r") as file:
                data: ProfileData = json.load(file)
                result.append(data)
        return result

    def add(self, data: ProfileData):
        key = improve_profileDB_key(data["command"])
        os.makedirs(f"{self.path}/{key}", exist_ok=True)
        key2 = data["now"]
        filename = f"{self.path}/{key}/{key2}.run"
        with open(filename, "w") as f:
            json.dump(data, f)

    def close(self):
        pass
