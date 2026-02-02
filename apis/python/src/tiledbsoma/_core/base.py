# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Definitions of the most fundamental types used by the SOMA project.

SOMA users should ordinarily not need to import this module directly; relevant
members will be exported to the ``somacore`` namespace.
"""

from __future__ import annotations

import abc
from collections.abc import MutableMapping
from typing import Any, ClassVar

from typing_extensions import LiteralString, Self

from . import options, types


class SOMAObject(metaclass=abc.ABCMeta):
    """The base type for all SOMA objects, containing common behaviors."""

    __slots__ = ("__weakref__",)

    @classmethod
    @abc.abstractmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode = "r",
        *,
        context: Any | None = None,  # noqa: ANN401
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Opens the SOMA object of this type at the given URI.

        Args:
            uri: The URI of the object to open.
            mode: The mode to open this in, either `r` or `w`.
            context: The Context value to use when opening the object.
            platform_config: Platform configuration options specific to
                this open operation.
        Returns: The SOMA object, opened for reading.
        Lifecycle: maturing
        """
        raise NotImplementedError

    @classmethod
    @abc.abstractmethod
    def exists(cls, uri: str, *, context: Any | None = None) -> bool:  # noqa: ANN401
        """Checks whether a SOMA object of this type is stored at the URI.

        Args:
            uri: The URI to check.
            context: The Context value to use when checking existence.

        Returns:
            True if the object exists and is of the correct type.
            False if the object does not exist, or is of a different type.
        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def uri(self) -> str:
        """The URI of this SOMA object.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    def context(self) -> types.ContextBase | None:
        """A value storing implementation-specific configuration information.

        This contains long-lived (i.e., not call-specific) information that is
        used by the SOMA implementation to access storage. This may include
        things like credentials, endpoint locations, or database connections.

        End users should treat this as an opaque value. While it may be passed
        from an existing SOMA object to be used in the creation of a new SOMA
        object, it should not be inspected.

        Lifecycle: maturing
        """
        return None

    @property
    @abc.abstractmethod
    def metadata(self) -> MutableMapping[str, Any]:
        """The metadata of this SOMA object.

        The returned value directly references the stored metadata; reads from
        and writes to it (provided the object is opened) are reflected in
        storage.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def mode(self) -> options.OpenMode:
        """Returns the mode this object was opened in, either ``r`` or ``w``.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def closed(self) -> bool:
        """True if this object has been closed; False if still open.

        Lifecycle: maturing
        """
        raise NotImplementedError

    soma_type: ClassVar[LiteralString]
    """A string describing the SOMA type of this object. This is constant.

    This uses ClassVar since you can't do abstract class properties.
    This is the equivalent, just without abc-based automatic verification.
    Overrides are marked Final with an ignore[misc] because mypy by default
    wants this to be mutable, and doesn't like overriding the mutable member
    with a Final member.
    """

    # Context management

    def close(self) -> None:
        """Releases any external resources held by this object.

        For objects opened for write, this also finalizes the write operation
        and ensures that all writes are completed before returning.

        This is also called automatically by the Python interpreter via
        ``__del__`` when this object is garbage collected, so the implementation
        must be idempotent.

        Lifecycle: maturing
        """
        # Default implementation does nothing.
        return

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: Any) -> None:  # noqa: ANN401
        self.close()

    def __del__(self) -> None:
        self.close()
        super_del = getattr(super(), "__del__", lambda: None)
        super_del()

    # Explicitly use Python's identity-based equality/hash checks.
    # These will show up in the `__mro__` before any other classes
    # provided a SOMAObject base is put first:
    #
    #    class SubType(SomeSOMAObject, MutableMapping):
    #        ...
    #
    #    # sub_type_inst.__eq__ uses object.__eq__ rather than
    #    # MutableMapping.__eq__.

    __eq__ = object.__eq__
    __hash__ = object.__hash__
