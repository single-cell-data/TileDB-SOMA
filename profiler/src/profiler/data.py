from __future__ import annotations

import errno
import glob
import hashlib
import json
import os
import pathlib
import uuid
from abc import ABC, abstractmethod
from typing import Any

import attr
import boto3
from botocore.client import ClientError


@attr.define
class ProfileData:
    """This class represents the data stored per run"""

    command: str
    timestamp: float
    stdout: str
    stderr: str
    tiledb_stats: dict[str, Any]
    tiledbsoma_version: str
    host_context: dict[str, str]
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
    custom_out: list[str | None]

    command_key: str = attr.field()

    @command_key.default
    def _command_key_factory(self):
        return _command_key(self.command)


DEFAULT_PROFILE_DB_PATH = "./profiling_db"


def _command_key(command: str) -> str:
    """Remove space characters from profileDB command keys."""
    return hashlib.md5(command.encode("utf-8")).hexdigest()


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
        if not pathlib.Path(self.path).exists():
            os.mkdir(self.path)

    def __str__(self):
        result = ""
        if pathlib.Path(self.path).exists():
            for command_hash in glob.glob(self.path + "/*"):
                with open(os.path.join(command_hash, "command.txt")) as f:
                    command = f.read()
                n_runs = len(glob.glob(os.path.join(command_hash, "*.json")))
                result += f'[{command_hash.split("/")[-1]}] "{command}": {n_runs} runs\n'
            return result
        return ""

    def find(self, command) -> list[ProfileData]:
        key = _command_key(command)
        if not pathlib.Path(f"{self.path}/{key}").exists():
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), key)
        result = []
        for filename in glob.glob(f"{self.path}/{key}/*.json"):
            with open(filename) as file:
                result.append(ProfileData(**json.load(file)))
        return result

    def add(self, data: ProfileData) -> str:
        key = _command_key(data.command)
        os.makedirs(f"{self.path}/{key}", exist_ok=True)
        filename = f"{self.path}/{key}/command.txt"

        if not pathlib.Path(filename).exists():
            with open(filename, "w") as f:
                f.write(data.command.strip())
        key2 = data.timestamp

        filename = f"{self.path}/{key}/{key2}.json"
        with open(filename, "w") as f:
            json.dump(attr.asdict(data), f)

        return filename

    def close(self):
        pass


class S3ProfileDB(ProfileDB):
    """Represents a S3-based implementation of a ProfileDB database.
    Each run is stored as a separate S3 object under a key with the structure `<bucket>/<command>/<timestamp>`.
    """

    def read_object_keys(self, prefix: str, suffix: str) -> list[str]:
        # return all the objects kets starting with prefix and ending with suffix
        result = self.s3.list_objects(Bucket=self.bucket_name, Prefix=prefix)
        keys: list[str] = []
        for o in result.get("Contents"):
            object_key = o.get("Key")
            if object_key.endswith(suffix):
                keys.append(object_key)
        return keys

    def read_s3_text(self, key: str) -> str:  # noqa: RET503
        # Assume the key is associated with one object. Otherwise, return the first object
        result = self.s3.list_objects(Bucket=self.bucket_name, Prefix=key)
        for o in result.get("Contents"):
            data = self.s3.get_object(Bucket=self.bucket_name, Key=o.get("Key"))
            return data["Body"].read().decode("utf-8")

    def bucket_exist_and_accessible(self):
        try:
            self.s3.head_bucket(Bucket=self.bucket.name)
        except ClientError:
            return False
        return True

    def __init__(self, bucket_name: str):
        # Initialize bucket's info
        self.s3 = boto3.client("s3")
        self.bucket_name = bucket_name
        self.bucket = boto3.resource("s3").Bucket(self.bucket_name)
        # Check if the bucket exists
        if not self.bucket_exist_and_accessible():
            raise (Exception(f"Bucket {self.bucket_name} does not exist or access is not granted."))

    def __str__(self):
        result = ""
        if self.bucket_exist_and_accessible():
            for command_file_key in self.read_object_keys("", "/command.txt"):
                # Extract runs prefix
                runs_prefix = command_file_key.replace("/command.txt", "")
                for _data_file_key in self.read_object_keys(runs_prefix, ".json"):
                    # Read command.
                    command = self.read_s3_text(command_file_key)
                    # Get the number of runs.
                    n_runs = len(self.read_object_keys(runs_prefix, ".json"))
                    result = +f'[{command_file_key}] "{command}": {n_runs} runs\n'
            return result
        return Exception(f"Bucket {self.bucket_name} does not exist or access is not granted.")

    def find(self, command) -> list[ProfileData]:
        key = _command_key(command)
        result = []
        # Extract all data files associated with this command
        for file_key in self.read_object_keys(key, ".json"):
            text = self.read_s3_text(file_key)
            if len(text) > 0:
                result.append(ProfileData(**json.loads(text)))
        return result

    def add(self, data: ProfileData) -> str:
        key = _command_key(data.command)
        s3_command_key = f"{key}/command.txt"

        filename = str(uuid.uuid4())

        with open(filename, "w") as fp:
            fp.write(data.command.strip())
        self.s3.upload_file(pathlib.Path(str(fp.name)).resolve(), self.bucket_name, s3_command_key)
        pathlib.Path(filename).unlink()

        key2 = data.timestamp
        data_key = f"{key}/{key2}.json"
        filename = str(uuid.uuid4())

        with open(filename, "w") as fp:
            json.dump(attr.asdict(data), fp)
        self.s3.upload_file(pathlib.Path(str(fp.name)).resolve(), self.bucket_name, data_key)
        pathlib.Path(filename).unlink()

        return data_key

    def close(self):
        pass
