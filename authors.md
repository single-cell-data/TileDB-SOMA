# Authors and Citation

## Authors

- **Paul Hoffman**. Maintainer, author.
  [](https://orcid.org/0000-0002-7693-8957)

- **Aaron Wolen**. Author. [](https://orcid.org/0000-0003-2542-2202)

- **Julia Dark**. Contributor.

- **John Kerl**. Contributor.

- **TileDB, Inc.**. Copyright holder, funder.

## Citation

Source:
[`DESCRIPTION`](https://github.com/single-cell-data/TileDB-SOMA/blob/2.3.0rc0/DESCRIPTION)

Hoffman P, Wolen A (2026). *tiledbsoma: 'TileDB' Stack of Matrices,
Annotated ('SOMA')*. R package version 2.3.0,
<https://single-cell-data.github.io/TileDB-SOMA/>.

    @Manual{,
      title = {tiledbsoma: 'TileDB' Stack of Matrices, Annotated ('SOMA')},
      author = {Paul Hoffman and Aaron Wolen},
      year = {2026},
      note = {R package version 2.3.0},
      url = {https://single-cell-data.github.io/TileDB-SOMA/},
    }

## Additional details

    The 'tiledbsoma' R package has been written by a team comprised of
    members from Chan Zuckerberg Initiative (CZI) and TileDB, Inc., as
    part of the 'SOMA' ("Stack of Matrices, Annotated") initative.

    The repository provides detailed commit statistics across the C++,
    Python, and R components of the implementation:
      https://github.com/single-cell-data/TileDB-SOMA/graphs/contributors

    The package also includes code written by other contributors as
    detailed below:


    -- libtiledbsoma/src/external/{src,include}/thread_pool/ ----------------------------------------
    The thread_pool implementation is from TileDB, Inc., and part of TileDB Embedded released at
    https://github.com/tiledb-inc/tiledb

     * @copyright Copyright (c) 2017-2022 TileDB, Inc.
     *            Copyright (c) 2011 The LevelDB Authors. All rights reserved.
     *
     * Redistribution and use in source and binary forms, with or without
     * modification, are permitted provided that the following conditions are met:

     * Redistributions of source code must retain the above copyright notice, this
     * list of conditions and the following disclaimer.
     *
     * Redistributions in binary form must reproduce the above copyright notice,
     * this list of conditions and the following disclaimer in the documentation
     * and/or other materials provided with the distribution.
     *
     * Neither the name of Google Inc. nor the names of its contributors may be used
     * to endorse or promote products derived from this software without specific
     * prior written permission.
     *
     * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
     * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
     * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
     * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
     * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
     * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
     * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
     * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
     * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


    -- apis/r/src/nanoarrow.{h,cc}/ -------------------------------------------------------------
    The nanoarrow package was written by Dewey Dunnington and other and released at
    https://github.com/apache/arrow-nanoarrow

    // Licensed to the Apache Software Foundation (ASF) under one
    // or more contributor license agreements.  See the NOTICE file
    // distributed with this work for additional information
    // regarding copyright ownership.  The ASF licenses this file
    // to you under the Apache License, Version 2.0 (the
    // "License"); you may not use this file except in compliance
    // with the License.  You may obtain a copy of the License at
    //
    //   http://www.apache.org/licenses/LICENSE-2.0
    //
    // Unless required by applicable law or agreed to in writing,
    // software distributed under the License is distributed on an
    // "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
    // KIND, either express or implied.  See the License for the
    // specific language governing permissions and limitations
    // under the License.
