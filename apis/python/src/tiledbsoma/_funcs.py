"""Utilities and decorators around modifying functions."""

import inspect
from typing import (
    Any,
    Callable,
    Collection,
    List,
    TypeVar,
)

from typing_extensions import ParamSpec

_Params = ParamSpec("_Params")
_T = TypeVar("_T")
_CT = TypeVar("_CT", bound=Callable[..., Any])


# Define a typeguard_ignore function so that we can use the `@typeguard_ignore`
# decorator without having to depend upon typeguard at runtime.
def typeguard_ignore(f: Callable[_Params, _T]) -> Callable[_Params, _T]:
    """No-op. Returns the argument unchanged."""
    return f


try:
    import typeguard

    typeguard_ignore = typeguard.typeguard_ignore  # noqa: F811
except ImportError:
    pass


def forwards_kwargs_to(
    dst: Callable[..., Any], *, exclude: Collection[str] = ()
) -> Callable[[_CT], _CT]:
    """Decorator function to update the signature with ``dst``'s kwargs.

    Example::

        def _internal(__it, a, b, c=3, *d, e=6, **f) -> None:
            ...

        @forwards_kwargs_to(_internal, exclude=("b",))
        def external(a, param1, param2, **kwargs) -> None:
            _internal(a="x", b=1, **kwargs)
            # Do other stuff

        # inspect.signature(external) is now:
        #     (a, param1, param2, *, c=3, e=6, **f) -> None

    This gives online help (like ``help()`` and autocompletion in ipython)
    information about the entire expected signature of the function without
    having to repeat the entire thing. However, it cannot (currently) provide
    information useful for static typing.
    """
    dst_params = tuple(inspect.Signature.from_callable(dst).parameters.values())

    def wrap(me: _CT) -> _CT:
        my_sig = inspect.Signature.from_callable(me)
        merged: List[inspect.Parameter] = []
        claimed_names = set(exclude)
        for param in my_sig.parameters.values():
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                # This is our **kwargs parameter, so substitute in
                # the destination function's arguments.
                for dst_param in dst_params:
                    if dst_param.kind == inspect.Parameter.VAR_KEYWORD:
                        # This is dst's **kwargs parameter. Pass it directly.
                        merged.append(dst_param)
                    elif _can_be_kwarg(dst_param):
                        # In this case:
                        #     def internal(a, b): ...
                        #     @forwards_kwargs_to(internal)
                        #     def public(b, **kwargs): ...
                        # `b` cannot be forwarded to `internal`, so we skip it.
                        if dst_param.name not in claimed_names:
                            # ...however, `a` can.
                            merged.append(
                                dst_param.replace(kind=inspect.Parameter.KEYWORD_ONLY)
                            )
                            claimed_names.add(dst_param.name)
                    # else: this param is positional-only.
            else:
                # This is one of our other params; add it to the signature.
                if _can_be_kwarg(param):
                    # This is one of our kwargs, so claim its name.
                    claimed_names.add(param.name)
                merged.append(param)

        # Lastly, make sure that we don't have any duplicate names.
        # The name of positional-only, *args, and **kwargs arguments
        # don't matter, so we can change them to make the signature unique.
        for idx, param in enumerate(merged):
            if _can_be_kwarg(param):
                continue  # Is already a named arg; can't be changed.

            # Ensure this isn't a dupe by adding `_`s until it's not.
            while param.name in claimed_names:
                param = param.replace(name=param.name + "_")
            claimed_names.add(param.name)
            merged[idx] = param

        me.__signature__ = my_sig.replace(  # type: ignore[attr-defined]
            parameters=merged
        )
        return me

    return wrap


def _can_be_kwarg(param: inspect.Parameter) -> bool:
    return not param.name.startswith("__") and param.kind in (
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY,
    )
