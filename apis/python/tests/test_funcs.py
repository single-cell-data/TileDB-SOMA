import inspect
import textwrap

import pytest

from tiledbsoma import _funcs


@pytest.mark.parametrize(
    ("dst_sig", "outer_sig", "want"),
    [
        ("(a: int) -> float", "(no_kwargs: int) -> int", "(no_kwargs: int) -> int"),
        (
            "(*, two: int = 2, three: float, **more)",
            "(one = 1, **kwargs) -> str",
            "(one=1, *, two: int = 2, three: float, **more) -> str",
        ),
        (
            "(__pos_only_dst, *args, hello)",
            "(__pos_only_outer, __args, *rest, kw1, **kwargs)",
            "(__pos_only_outer, __args, *rest, kw1, hello)",
        ),
        (
            "(dont_shadow: int, *ignore, do_rename_: str, **do_rename: frozenset)",
            "(__pos_only, *dont_shadow, do_rename: bytes = b'', **kwargs) -> complex",
            "(__pos_only, *dont_shadow_, do_rename: bytes = b'', dont_shadow: int, do_rename_: str, **do_rename__: frozenset) -> complex",
        ),
        (
            "(pos_only, some_name, /, dst_both_arg, *args, dst_kwarg=1, **dup_dict: complex)",
            "(dst_both_arg, dst_kwarg: int = 10, /, dup_dict: dict = (), **kwargs) -> None",
            "(dst_both_arg_, dst_kwarg_: int = 10, /, dup_dict: dict = (), *, dst_both_arg, dst_kwarg=1, **dup_dict_: complex) -> None",
        ),
    ],
)
def test_forwards_kwargs(dst_sig: str, outer_sig: str, want: str):
    code = textwrap.dedent(
        f"""
        def dst{dst_sig}: ...
        @forwards_kwargs_to(dst)
        def outer{outer_sig}: ...
        """
    )
    variables = dict(forwards_kwargs_to=_funcs.forwards_kwargs_to)
    exec(code, variables)
    got_sig = inspect.signature(variables["outer"])
    assert str(got_sig) == want


def test_forwards_kwargs_exclude():
    def dst(one, two, three): ...

    @_funcs.forwards_kwargs_to(dst, exclude=("one", "three"))
    def wrapped(one, four, five, **kwargs) -> None: ...

    got_sig = inspect.signature(wrapped)
    assert str(got_sig) == "(one, four, five, *, two) -> None"
