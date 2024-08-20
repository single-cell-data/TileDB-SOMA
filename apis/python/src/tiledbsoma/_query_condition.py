# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""A high level wrapper around the Pybind11 query_condition.cc implementation for
filtering query results on attribute values.
"""
import ast
from typing import Any, Callable, List, Optional, Tuple, Union

import attrs
import numpy as np
import pyarrow as pa

from . import pytiledbsoma as clib
from ._exception import SOMAError
from ._util import pa_types_is_string_or_bytes

QueryConditionNodeElem = Union[ast.Name, ast.Constant, ast.NameConstant, ast.Call]


@attrs.define
class QueryCondition:
    """Class representing a TileDB query condition object for attribute filtering
    pushdown.

    A query condition is set with a string representing an expression
    as defined by the grammar below. A more straight forward example of usage is
    given beneath.

    When querying a sparse array, only the values that satisfy the given
    condition are returned (coupled with their associated coordinates). An example
    may be found in `examples/query_condition_sparse.py`.

    For dense arrays, the given shape of the query matches the shape of the output
    array. Values that DO NOT satisfy the given condition are filled with the
    TileDB default fill value. Different attribute types have different default
    fill values as outlined here
    (https://docs.tiledb.com/main/background/internal-mechanics/writing#default-fill-values).
    An example may be found in `examples/query_condition_dense.py`.

    **BNF:**

    A query condition is made up of one or more Boolean expressions. Multiple
    Boolean expressions are chained together with Boolean operators. The ``or_op``
    Boolean operators are given lower presedence than ``and_op``.

        ``query_cond ::= bool_term | query_cond or_op bool_term``

        ``bool_term ::= bool_expr | bool_term and_op bool_expr``

    Logical ``and`` and bitwise ``&`` Boolean operators are given equal precedence.

        ``and_op ::= and | &``

    Likewise, ``or`` and ``|`` are given equal precedence.

        ``or_op ::= or | |``

    We intend to support ``not`` in future releases.

    A Boolean expression may either be a comparison expression or membership
    expression.

        ``bool_expr ::= compare_expr | member_expr``

    A comparison expression contains a comparison operator. The operator works on a
    TileDB attribute name and value.

        ``compare_expr ::= attr compare_op val
            | val compare_op attr
            | val compare_op attr compare_op val``

    All comparison operators are supported.

        ``compare_op ::= < | > | <= | >= | == | !=``

    A memership expression contains the membership operator, ``in``. The operator
    works on a TileDB attribute and list of values.

        ``member_expr ::= attr in <list>``

    TileDB attribute names are Python valid variables or a ``attr()`` casted string.

        ``attr ::= <variable> | attr(<str>)``

    Values are any Python-valid number or string. datetime64 values should first be
    cast to UNIX seconds. Values may also be casted with ``val()``.

        ``val ::= <num> | <str> | val(val)``

    Example:
        with tiledbsoma.open(df_uri) as dataframe:
            # Select cells where the attribute values for `foo` are less than 5
            # and `bar` equal to string "asdf".
            # Note precedence is equivalent to:
            # "foo > 5 or ('asdf' == attr('b a r') and baz <= val(1.0))"
            foo_bar_baz = dataframe.read(
                value_filter="foo > 5 or 'asdf' == attr('b a r') and baz <= val(1.0)")

            # Select cells where the attribute values for `foo` are equal to
            # 1, 2, or 3.
            # Note this is equivalent to:
            # "foo == 1 or foo == 2 or foo == 3"
            foo_123 = dataframe.read(value_filter="foo in [1, 2, 3]")
    """

    expression: str
    tree: ast.Expression = attrs.field(init=False, repr=False)
    c_obj: clib.PyQueryCondition = attrs.field(init=False, repr=False)

    def __attrs_post_init__(self):
        try:
            self.tree = ast.parse(self.expression, mode="eval")
        except Exception as pex:
            raise SOMAError(
                "Could not parse the given QueryCondition statement: "
                f"{self.expression}"
            ) from pex

        if not self.tree:
            raise SOMAError(
                "The query condition statement could not be parsed properly. "
                "(Is this an empty expression?)"
            )

    def init_query_condition(
        self,
        schema: pa.Schema,
        query_attrs: Optional[List[str]],
    ):
        try:
            qctree = QueryConditionTree(schema, query_attrs)
            self.c_obj = qctree.visit(self.tree.body)
        except Exception as pex:
            raise SOMAError(pex)

        if not isinstance(self.c_obj, clib.PyQueryCondition):
            raise SOMAError(
                "Malformed query condition statement. A query condition must "
                "be made up of one or more boolean expressions."
            )

        return query_attrs


@attrs.define
class QueryConditionTree(ast.NodeVisitor):
    schema: pa.Schema
    query_attrs: List[str]

    def visit_BitOr(self, node):
        return clib.TILEDB_OR

    def visit_Or(self, node):
        return clib.TILEDB_OR

    def visit_BitAnd(self, node):
        return clib.TILEDB_AND

    def visit_And(self, node):
        return clib.TILEDB_AND

    def visit_Gt(self, node):
        return clib.TILEDB_GT

    def visit_GtE(self, node):
        return clib.TILEDB_GE

    def visit_Lt(self, node):
        return clib.TILEDB_LT

    def visit_LtE(self, node):
        return clib.TILEDB_LE

    def visit_Eq(self, node):
        return clib.TILEDB_EQ

    def visit_NotEq(self, node):
        return clib.TILEDB_NE

    def visit_In(self, node):
        return node

    def visit_NotIn(self, node):
        return node

    def visit_Is(self, node):
        raise SOMAError("the `is` operator is not supported")

    def visit_IsNot(self, node):
        raise SOMAError("the `is not` operator is not supported")

    def visit_List(self, node):
        return list(node.elts)

    def visit_Attribute(self, node) -> clib.PyQueryCondition:
        raise SOMAError(
            f"Unhandled dot operator in {ast.dump(node)} -- if your attribute name "
            'has a dot in it, e.g. `orig.ident`, please wrap it with `attr("...")`, '
            'e.g. `attr("orig.ident")`'
        )

    def visit_Compare(self, node: ast.Compare) -> clib.PyQueryCondition:
        operator = self.visit(node.ops[0])

        if operator in (
            clib.TILEDB_GT,
            clib.TILEDB_GE,
            clib.TILEDB_LT,
            clib.TILEDB_LE,
            clib.TILEDB_EQ,
            clib.TILEDB_NE,
        ):
            result = self.aux_visit_Compare(
                self.visit(node.left),
                operator,
                self.visit(node.comparators[0]),
            )

            # Handling cases val < attr < val
            for lhs, op, rhs in zip(
                node.comparators[:-1], node.ops[1:], node.comparators[1:]
            ):
                value = self.aux_visit_Compare(
                    self.visit(lhs), self.visit(op), self.visit(rhs)
                )
                result = result.combine(value, clib.TILEDB_AND)
        elif isinstance(operator, (ast.In, ast.NotIn)):
            rhs = node.comparators[0]
            if not isinstance(rhs, ast.List):
                raise SOMAError(
                    "`in` operator syntax must be written as `attr in ['l', 'i', 's', 't']`"
                )

            variable = node.left.id
            values = [self.get_val_from_node(val) for val in self.visit(rhs)]
            if len(values) == 0:
                raise SOMAError(
                    "At least one value must be provided to the set membership"
                )

            dt = self.schema.field(variable).type
            if pa.types.is_dictionary(dt):
                dt = dt.value_type

            if pa_types_is_string_or_bytes(dt):
                dtype = "string"
            else:
                dtype = np.dtype(dt.to_pandas_dtype()).name

            # sdf.read(column_names=["foo"], value_filter='bar == 999') should
            # result in bar being added to the column names. See also
            # https://github.com/single-cell-data/TileDB-SOMA/issues/755
            att = self.get_att_from_node(node.left)
            if att not in self.query_attrs:
                self.query_attrs.append(att)

            op = clib.TILEDB_IN if isinstance(operator, ast.In) else clib.TILEDB_NOT_IN
            result = self.create_pyqc(dtype)(node.left.id, values, op)

        else:
            raise SOMAError(f"unrecognized operator in <<{ast.dump(node)}>>")

        return result

    def aux_visit_Compare(
        self,
        lhs: QueryConditionNodeElem,
        op_node: clib.tiledb_query_condition_op_t,
        rhs: QueryConditionNodeElem,
    ) -> clib.PyQueryCondition:
        att, val, op = self.order_nodes(lhs, rhs, op_node)

        att = self.get_att_from_node(att)
        val = self.get_val_from_node(val)

        dt = self.schema.field(att).type
        if pa.types.is_dictionary(dt):
            dt = dt.value_type

        if pa_types_is_string_or_bytes(dt):
            dtype = "string"
        else:
            dtype = np.dtype(dt.to_pandas_dtype()).name
        val = self.cast_val_to_dtype(val, dtype)

        pyqc = clib.PyQueryCondition()
        self.init_pyqc(pyqc, dtype)(att, val, op)

        return pyqc

    def is_att_node(self, att: QueryConditionNodeElem) -> bool:
        if isinstance(att, ast.Call):
            if not isinstance(att.func, ast.Name):
                raise SOMAError(f"Unrecognized expression {att.func}.")

            if att.func.id != "attr":
                return False

            return isinstance(att.args[0], (ast.Constant, ast.NameConstant))

        return isinstance(att, ast.Name)

    def order_nodes(
        self,
        att: QueryConditionNodeElem,
        val: QueryConditionNodeElem,
        op: clib.tiledb_query_condition_op_t,
    ) -> Tuple[
        QueryConditionNodeElem,
        QueryConditionNodeElem,
        clib.tiledb_query_condition_op_t,
    ]:
        if not self.is_att_node(att):
            REVERSE_OP = {
                clib.TILEDB_GT: clib.TILEDB_LT,
                clib.TILEDB_GE: clib.TILEDB_LE,
                clib.TILEDB_LT: clib.TILEDB_GT,
                clib.TILEDB_LE: clib.TILEDB_GE,
                clib.TILEDB_EQ: clib.TILEDB_EQ,
                clib.TILEDB_NE: clib.TILEDB_NE,
            }

            op = REVERSE_OP[op]
            att, val = val, att

        return att, val, op

    def get_att_from_node(self, node: QueryConditionNodeElem) -> Any:
        if self.is_att_node(node):
            att_node = node

            if isinstance(att_node, ast.Call):
                if not isinstance(att_node.func, ast.Name):
                    raise SOMAError(f"Unrecognized expression {att_node.func}.")
                att_node = att_node.args[0]

            if isinstance(att_node, ast.Name):
                att = str(att_node.id)
            elif isinstance(att_node, ast.Constant):
                att = str(att_node.value)
            else:
                raise SOMAError(
                    f"Incorrect type for attribute name: {ast.dump(att_node)}"
                )
        else:
            raise SOMAError(f"Incorrect type for attribute name: {ast.dump(node)}")

        if not att not in self.schema:
            raise SOMAError(f"`{att}` not found in schema.")

        # sdf.read(column_names=["foo"], value_filter='bar == 999') should
        # result in bar being added to the column names. See also
        # https://github.com/single-cell-data/TileDB-SOMA/issues/755
        if att not in self.query_attrs:
            self.query_attrs.append(att)

        return att

    def get_val_from_node(self, node: QueryConditionNodeElem) -> Any:
        val_node = node

        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name):
                raise SOMAError(f"Unrecognized expression {node.func}.")

            if node.func.id == "val":
                val_node = node.args[0]
            else:
                raise SOMAError(f"Incorrect type for cast value: {node.func.id}")

        if isinstance(val_node, ast.Constant) or isinstance(val_node, ast.NameConstant):
            val = val_node.value
        else:
            raise SOMAError(
                f"Incorrect type for comparison value: {ast.dump(val_node)}: right-hand sides must be constant"
                " expressions, not variables -- did you mean to quote the right-hand side as a string?"
            )

        return val

    def cast_val_to_dtype(
        self,
        val: Union[str, int, float, bytes, np.int32, np.int64, np.float32],
        dtype: str,
    ) -> Union[str, int, float, bytes, np.int32, np.int64, np.float32]:
        if dtype != "string":
            try:
                # this prevents numeric strings ("1", '123.32') from getting
                # casted to numeric types
                if isinstance(val, str):
                    raise SOMAError(f"Cannot cast `{val}` to {dtype}.")
                if np.issubdtype(dtype, np.datetime64):
                    cast = getattr(np, "int64")
                # silence DeprecationWarning: `np.bool`
                elif dtype == "bool":
                    cast = bool
                else:
                    cast = getattr(np, dtype)
                val = cast(val)
            except ValueError:
                raise SOMAError(f"Cannot cast `{val}` to {dtype}.")

        return val

    def init_pyqc(self, pyqc: clib.PyQueryCondition, dtype: str) -> Callable:
        if dtype != "string" and np.issubdtype(dtype, np.datetime64):
            dtype = "int64"

        init_fn_name = f"init_{dtype}"

        if not hasattr(pyqc, init_fn_name):
            raise SOMAError(f"PyQueryCondition.{init_fn_name}() not found.")

        return getattr(pyqc, init_fn_name)

    def create_pyqc(self, dtype: str) -> Callable:
        if dtype != "string":
            if np.issubdtype(dtype, np.datetime64):
                dtype = "int64"
            elif np.issubdtype(dtype, bool):
                dtype = "uint8"

        create_fn_name = f"create_{dtype}"

        try:
            return getattr(clib.PyQueryCondition, create_fn_name)
        except AttributeError as ae:
            raise SOMAError(f"PyQueryCondition.{create_fn_name}() not found.") from ae

    def visit_BinOp(self, node: ast.BinOp) -> clib.PyQueryCondition:
        try:
            op = self.visit(node.op)
        except KeyError:
            raise SOMAError(
                f"Unsupported binary operator: {ast.dump(node.op)}. Only & is currently supported."
            )

        result = self.visit(node.left)
        rhs = node.right[1:] if isinstance(node.right, list) else [node.right]
        for value in rhs:
            visited = self.visit(value)
            if not isinstance(result, clib.PyQueryCondition):
                raise Exception(
                    f"Unable to parse expression component {ast.dump(node)} -- did you mean to quote it as a string?"
                )
            result = result.combine(visited, op)

        return result

    def visit_BoolOp(self, node: ast.BoolOp) -> clib.PyQueryCondition:
        try:
            op = self.visit(node.op)
        except KeyError:
            raise SOMAError(f"Unsupported Boolean operator: {ast.dump(node.op)}.")

        result = self.visit(node.values[0])
        for value in node.values[1:]:
            result = result.combine(self.visit(value), op)

        return result

    def visit_Call(self, node: ast.Call) -> ast.Call:
        if not isinstance(node.func, ast.Name):
            raise SOMAError(f"Unrecognized expression {node.func}.")

        if node.func.id not in ["attr", "val"]:
            raise SOMAError("Valid casts are attr() or val().")

        if len(node.args) != 1:
            raise SOMAError(
                f"Exactly one argument must be provided to {node.func.id}()."
            )

        return node

    def visit_Name(self, node: ast.Name) -> ast.Name:
        return node

    def visit_Constant(self, node: ast.Constant) -> ast.Constant:
        return node

    def visit_NameConstant(self, node: ast.NameConstant) -> ast.NameConstant:
        return node

    def visit_UnaryOp(self, node: ast.UnaryOp, sign: int = 1):
        if isinstance(node.op, ast.UAdd):
            sign *= 1
        elif isinstance(node.op, ast.USub):
            sign *= -1
        else:
            raise SOMAError(f"Unsupported UnaryOp type. Saw {ast.dump(node)}.")

        if isinstance(node.operand, ast.UnaryOp):
            return self.visit_UnaryOp(node.operand, sign)
        else:
            if isinstance(node.operand, ast.Constant) or isinstance(
                node.operand, ast.NameConstant
            ):
                node.operand.value *= sign
            elif isinstance(node.operand, ast.Num):
                node.operand.n *= sign
            else:
                raise SOMAError(
                    f"Unexpected node type following UnaryOp. Saw {ast.dump(node)}."
                )

            return node.operand
