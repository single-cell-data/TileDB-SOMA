from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple, Union

import tiledb

from .soma import SOMA
from .soma_options import SOMAOptions
from .soma_slice import SOMASlice
from .tiledb_group import TileDBGroup
from .types import Ids


class SOMACollection(TileDBGroup):
    """
    Implements a collection of ``SOMA`` objects.
    """

    # This is a cache to avoid the overhead of calling the SOMA constructor repeatedly.  That
    # constructor isn't particuarly expensive, except that for tiledb-cloud URIs it needs to
    # interrogate the server repeatedly for recursive group-member URIs, which has web-request
    # latency.
    _somas: Dict[str, SOMA]

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        *,
        name: str = "soco",
        soma_options: Optional[SOMAOptions] = None,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        Create a new ``SOMACollection`` object. The existing group is opened at the
        specified ``uri`` if one is present, otherwise a new group will be created upon ingest.

        :param uri: URI of the TileDB group
        """

        # People can (and should) call by name. However, it's easy to forget. For example,
        # if someone does 'tiledbsoma.SOMACollection("myuri", ctx)' instead of 'tiledbsoma.SOMA("myury", ctx)',
        # behavior will not be what they expect, and we should let them know sooner than later.
        if name is not None:
            assert isinstance(name, str)
        if soma_options is not None:
            assert isinstance(soma_options, SOMAOptions)
        if config is not None:
            assert isinstance(config, tiledb.Config)
        if ctx is not None:
            assert isinstance(ctx, tiledb.Ctx)
        if parent is not None:
            assert isinstance(parent, TileDBGroup)

        if ctx is None and config is not None:
            ctx = tiledb.Ctx(config)
        if soma_options is None:
            soma_options = SOMAOptions()  # Use default values from the constructor
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            soma_options=soma_options,
            ctx=ctx,
        )

        self._somas = {}

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of SOMACollection.
        """
        lines = [
            "URI:        " + self.uri,
        ]
        if self.exists():
            lines.append(f"SOMA count: {len(self)}")
        else:
            lines.append("Unpopulated")
        return "\n".join(lines)

    # ----------------------------------------------------------------
    def __len__(self) -> int:
        """
        Implements ``len(soco)``. Returns the number of elements in the collection.
        """
        return len(self.get_member_names())

    # ----------------------------------------------------------------
    def add(self, soma: SOMA, relative: Optional[bool] = None) -> None:
        """
        Adds a ``SOMA`` to the ``SOMACollection``.

        * If ``relative`` is not supplied, it's taken from the ``soma_options`` the collection was instantiated with.

        * If ``relative`` is ``False``, either via the ``relative`` argument or via ``soma_options.member_uris_are_relative``, then the collection will have the absolute path of the SOMA. For populating SOMA elements within a SOMACollection on local disk, this can be useful if you want to be able to move the SOMACollection storage around and have it remember the (unmoved) locations of SOMA objects elsewhere, i.e.  if the SOMACollection is in one place while its members are in other places.  If the SOMAs in the collection are contained within the SOMACollection directory, you probably want ``relative=True``.

        * If ``relative`` is ``True``, either via the ``relative`` argument or via ``soma_options.member_uris_are_relative``, then the group will have the relative path of the member. For TileDB Cloud, this is never the right thing to do. For local-disk or S3 storage, this is essential if you want to move a SOMA to another directory and have it remember the locations of the members within it.  In this case the SOMA storage must be located as a direct component of the collection storage.  Example: ``s3://mybucket/soco`` and ``s3://mybucket/soco/soma1``.

        * If ``relative`` is ``None``, either via the ``relative`` argument or via ``soma_options.member_uris_are_relative``, then we select ``relative=False`` if the URI starts with ``tiledb://``, else we select ``relative=True``. This is the default.
        """
        self._add_object(soma, relative=relative, check_is_direct_child=True)

    # ----------------------------------------------------------------
    def remove(self, soma: Union[SOMA, str]) -> None:
        """
        Removes a ``SOMA`` from the ``SOMACollection``, when invoked as ``soco.remove("namegoeshere")``.
        """
        if isinstance(soma, str):
            self._remove_object_by_name(soma)
        else:
            self._remove_object(soma)

    def __delattr__(self, matrix_name: str) -> None:
        """
        Removes a ``SOMA`` from the ``SOMACollection``, when invoked as ``del soco.namegoeshere``.
        """
        self.remove(matrix_name)

    def __delitem__(self, matrix_name: str) -> None:
        """
        Removes a ``SOMA`` from the ``SOMACollection``, when invoked as ``del soco["namegoeshere"]``.
        """
        self.remove(matrix_name)

    # ----------------------------------------------------------------
    def keys(self) -> Sequence[str]:
        """
        Returns the names of the SOMAs in the collection.
        """
        return self.get_member_names()

    # ----------------------------------------------------------------
    def __iter__(self) -> Iterator[SOMA]:
        """
        Implements ``for soma in soco: ...``
        """
        for _, soma in self._get_all_somas().items():
            yield soma

    def _get_all_somas(self) -> Dict[str, SOMA]:
        self._populate_all()  # Parallelized cache pre-fill
        return self._somas

    def _populate_all(self) -> None:
        """
        Run constructors for all SOMAs in the collection, parallelized.

        Note that we intentionally do not run this in the ``SOMACollection`` constructor --
        collections may have very many SOMAs in them, and when an user instantiates a collection in
        a notebook, we want the response to be rapid. If their user experience is that they only
        want to examine a single SOMA via ``soco.keys()`` and ``soma = soco['namehere']``, we should not
        run hundreds of constructors (parallelized or not) -- this is not a good use of their
        resources or time.

        However, if the user is indeed iterating over the collection, we do need to run all the
        constructors, and to do so quickly.

        Thus, the nominal use-case for this method is in our ``__iter__``.

        The reason this is a good candidate for thread-pooling is member name-to-URI resolution in
        the constructors. These read from storage using TileDB's core C++ storage engine, which
        releases the GIL -- moreover, for TileDB Cloud storage, the member name-to-URI resolution in
        core uses HTTP requests which are a fine candidate for parallelization.
        """
        futures = []
        with ThreadPoolExecutor(
            max_workers=self._soma_options.max_thread_pool_workers
        ) as executor:
            for name, uri in self.get_member_names_to_uris().items():
                if name not in self._somas:
                    future = executor.submit(self._populate_aux, name=name, uri=uri)
                    futures.append(future)

        for future in futures:
            name, soma = future.result()
            self._somas[name] = soma

    def _populate_aux(self, name: str, uri: str) -> Tuple[str, SOMA]:
        """
        Helper method for ``_populate``.``
        """
        return (name, SOMA(uri=uri, name=name, parent=self, ctx=self._ctx))

    # ----------------------------------------------------------------
    def __contains__(self, name: str) -> bool:
        """
        Implements ``name in soco``
        """
        if not self.exists():
            return False
        with self._open("r") as G:
            return name in G

    # At the tiledb-py API level, *all* groups are name-indexable.  But here at the tiledbsoma-py
    # level, we implement name-indexing only for some groups:
    #
    # * Most soma member references are done using Python's dot syntax. For example, rather than
    #   soma['X'], we have simply soma.X, and likewise, soma.raw.X.  Likewise soma.obs and soma.var.
    #
    # * Index references are supported for obsm, varm, obsp, varp, and uns. E.g.
    #   soma.obsm['X_pca'] or soma.uns['neighbors']['params']['method']
    #
    # * Overloading the ``[]`` operator at the TileDBGroup level isn't necessary -- e.g. we don't need
    #   soma['X'] when we have soma.X -- but also it causes circular-import issues in Python.
    #
    # * Rather than doing a TileDBIndexableGroup which overloads the ``[]`` operator, we overload
    #   the ``[]`` operator separately in the various classes which need indexing. This is again to
    #   avoid circular-import issues, and means that [] on ``AnnotationMatrixGroup`` will return an
    #   ``AnnotationMatrix, [] on ``UnsGroup`` will return ``UnsArray`` or ``UnsGroup``, etc.
    def __getitem__(self, name: str) -> Optional[SOMA]:
        """
        Returns a ``SOMA`` element at the given name within the group, or ``None`` if no such
        member exists.  Overloads the ``[...]`` operator.
        """
        if name in self._somas:
            # SOMA-constructor cache
            return self._somas[name]

        with self._open("r") as G:
            try:
                obj = G[name]  # This returns a tiledb.object.Object.
            except tiledb.TileDBError:
                return None

            if obj.type != tiledb.group.Group:
                raise Exception(
                    f"Internal error: found element which is not a subgroup: type is {obj.type}"
                )

            return SOMA(uri=obj.uri, name=name, parent=self)

    # ----------------------------------------------------------------
    def query(
        self,
        *,
        obs_attrs: Optional[Sequence[str]] = None,
        obs_query_string: Optional[str] = None,
        obs_ids: Optional[Ids] = None,
        var_attrs: Optional[Sequence[str]] = None,
        var_query_string: Optional[str] = None,
        var_ids: Optional[Ids] = None,
        X_layer_names: Optional[Sequence[str]] = None,
        return_arrow: bool = False,
    ) -> List[SOMASlice]:
        """
        Subselects the ``obs``, ``var``, and ``X/data`` using the specified queries on ``obs`` and ``var``,
        concatenating across SOMAs in the collection.  Queries use the TileDB-Py ``QueryCondition``
        API.

        If ``obs_query_string`` is ``None``, the ``obs`` dimension is not filtered and all of ``obs`` is
        used; similiarly for ``var``. Return value of ``None`` indicates an empty slice.  If ``obs_ids``
        or ``var_ids`` are not ``None``, they are effectively ANDed into the query.  For example, you
        can pass in a known list of ``obs_ids``, then use ``obs_query_string`` to further restrict the
        query.

        If ``X_layer_names`` is `None`, they are all returned; otherwise you can specify which layer(s)
        you want to be operated on.

        If ``obs_attrs`` or ``var_attrs`` are unspecified, slices will take all ``obs``/``var`` attributes
        from their source SOMAs; if they are specified, slices will take the specified ``obs``/``var``
        attributes.  If all SOMAs in the collection have the same ``obs``/``var`` attributes, then you
        needn't specify these; if they don't, you must.
        """

        return SOMA.queries(
            list(self._get_all_somas().values()),
            obs_attrs=obs_attrs,
            obs_query_string=obs_query_string,
            obs_ids=obs_ids,
            var_attrs=var_attrs,
            var_query_string=var_query_string,
            var_ids=var_ids,
            X_layer_names=X_layer_names,
            return_arrow=return_arrow,
        )

    # ----------------------------------------------------------------
    def find_unique_obs_values(self, obs_label: str) -> Set[str]:
        """
        Given an ``obs`` label such as ``cell_type`` or ``tissue``, returns a set of unique
        values for that label among all SOMAs in the collection.
        """
        return self._find_unique_obs_or_var_values(obs_label, True)

    def find_unique_var_values(self, var_label: str) -> Set[str]:
        """
        Given an ``var`` label such as ``feature_name``, returns a set of unique values for
        that label among all SOMAs in the collection.
        """
        return self._find_unique_obs_or_var_values(var_label, False)

    def _find_unique_obs_or_var_values(
        self, obs_or_var_label: str, use_obs: bool
    ) -> Set[str]:
        """
        Helper method for ``find_unique_obs_values`` and ``find_unique_var_values``.
        """
        unique_values_in_soco = set()

        for soma in self:
            annotation_matrix = soma.obs if use_obs else soma.var
            if obs_or_var_label not in annotation_matrix.keys():
                continue
            unique_values_in_soco.update(annotation_matrix.df()[obs_or_var_label])

        return unique_values_in_soco
