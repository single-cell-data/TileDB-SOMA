# Reading from SOMA objects

## Overview

In this tutorial we’ll learn how to read data from various SOMA objects.
We will assume familiarity with SOMA objects, so it is recommended to go
through the [`vignette("soma-objects")`](../articles/soma-objects.md)
first.

A core feature of SOMA is the ability to read *subsets* of data from
disk into memory as slices. SOMA uses [Apache
Arrow](https://arrow.apache.org/) as an intermediate in-memory storage.
From here, the slices can be further converted into native R objects,
like data frames and matrices.

``` r
library(tiledbsoma)
```

## Example data

Load the bundled `SOMAExperiment` containing a subsetted version of the
10X genomics [PBMC
dataset](https://satijalab.github.io/seurat-object/reference/pbmc_small.html)
provided by SeuratObject. This will return a `SOMAExperiment` object.
This is a small dataset that easily fits into memory, but we’ll focus on
operations that can easily scale to larger datasets as well.

``` r
experiment <- load_dataset("soma-exp-pbmc-small")
```

## SOMA DataFrame

We’ll start with the `obs` dataframe. Simply calling the
`read()$concat()` method will load all of the data in memory as an
[Arrow Table](https://arrow.apache.org/docs/r/reference/table.html).

``` r
obs <- experiment$obs
obs$read()$concat()
```

    ## Table
    ## 80 rows x 9 columns
    ## $soma_joinid <int64 not null>
    ## $orig.ident <dictionary<values=string, indices=int8>>
    ## $nCount_RNA <double>
    ## $nFeature_RNA <int32>
    ## $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
    ## $letter.idents <dictionary<values=string, indices=int8>>
    ## $groups <large_string>
    ## $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
    ## $obs_id <large_string>

This is easily converted into a `data.frame` using Arrow’s methods:

``` r
obs$read()$concat()$to_data_frame()
```

    ##    soma_joinid    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
    ## 1            0 SeuratProject         70           47               0
    ## 2            1 SeuratProject         85           52               0
    ## 3            2 SeuratProject         87           50               1
    ## 4            3 SeuratProject        127           56               0
    ## 5            4 SeuratProject        173           53               0
    ## 6            5 SeuratProject         70           48               0
    ## 7            6 SeuratProject         64           36               0
    ## 8            7 SeuratProject         72           45               0
    ## 9            8 SeuratProject         52           36               0
    ## 10           9 SeuratProject        100           41               0
    ## 11          10 SeuratProject         62           31               0
    ## 12          11 SeuratProject        191           61               0
    ## 13          12 SeuratProject        101           41               0
    ## 14          13 SeuratProject         51           26               0
    ## 15          14 SeuratProject         99           45               0
    ## 16          15 SeuratProject        168           44               0
    ## 17          16 SeuratProject         67           33               0
    ## 18          17 SeuratProject        135           45               0
    ## 19          18 SeuratProject         79           43               0
    ## 20          19 SeuratProject        109           41               0
    ## 21          20 SeuratProject        298           65               1
    ## 22          21 SeuratProject        406           74               1
    ## 23          22 SeuratProject        213           48               1
    ## 24          23 SeuratProject        231           49               1
    ## 25          24 SeuratProject        463           77               1
    ## 26          25 SeuratProject         87           42               1
    ## 27          26 SeuratProject        327           62               1
    ## 28          27 SeuratProject        224           50               1
    ## 29          28 SeuratProject        361           76               1
    ## 30          29 SeuratProject        353           80               1
    ## 31          30 SeuratProject        246           59               0
    ## 32          31 SeuratProject        115           51               0
    ## 33          32 SeuratProject        189           53               0
    ## 34          33 SeuratProject        187           61               0
    ## 35          34 SeuratProject        156           48               0
    ## 36          35 SeuratProject        164           47               0
    ## 37          36 SeuratProject        221           67               0
    ## 38          37 SeuratProject        151           59               0
    ## 39          38 SeuratProject        126           53               0
    ## 40          39 SeuratProject        316           65               0
    ## 41          40 SeuratProject        156           60               0
    ## 42          41 SeuratProject        139           61               0
    ## 43          42 SeuratProject        108           44               0
    ## 44          43 SeuratProject         41           32               0
    ## 45          44 SeuratProject        146           47               0
    ## 46          45 SeuratProject        104           40               0
    ## 47          46 SeuratProject        126           48               0
    ## 48          47 SeuratProject         94           55               0
    ## 49          48 SeuratProject        204           52               0
    ## 50          49 SeuratProject         99           45               0
    ## 51          50 SeuratProject        371           75               1
    ## 52          51 SeuratProject        387           83               1
    ## 53          52 SeuratProject        139           50               1
    ## 54          53 SeuratProject         99           42               1
    ## 55          54 SeuratProject        443           77               1
    ## 56          55 SeuratProject        417           75               0
    ## 57          56 SeuratProject        502           81               1
    ## 58          57 SeuratProject        324           76               1
    ## 59          58 SeuratProject        292           71               1
    ## 60          59 SeuratProject        443           81               0
    ## 61          60 SeuratProject        787           88               0
    ## 62          61 SeuratProject        612           69               1
    ## 63          62 SeuratProject        286           68               0
    ## 64          63 SeuratProject        462           86               1
    ## 65          64 SeuratProject        872           96               1
    ## 66          65 SeuratProject        709           94               1
    ## 67          66 SeuratProject        745           84               1
    ## 68          67 SeuratProject        328           72               1
    ## 69          68 SeuratProject        389           73               1
    ## 70          69 SeuratProject        754           83               0
    ## 71          70 SeuratProject        212           38               0
    ## 72          71 SeuratProject        172           29               0
    ## 73          72 SeuratProject        168           37               0
    ## 74          73 SeuratProject        210           33               0
    ## 75          74 SeuratProject        228           39               0
    ## 76          75 SeuratProject        527           47               0
    ## 77          76 SeuratProject        202           30               0
    ## 78          77 SeuratProject        157           29               0
    ## 79          78 SeuratProject        150           30               0
    ## 80          79 SeuratProject        233           76               1
    ##    letter.idents groups RNA_snn_res.1         obs_id
    ## 1              A     g2             0 ATGCCAGAACGACT
    ## 2              A     g1             0 CATGGCCTGTGCAT
    ## 3              B     g2             0 GAACCTGATGAACC
    ## 4              A     g2             0 TGACTGGATTCTCA
    ## 5              A     g2             0 AGTCAGACTGCACA
    ## 6              A     g1             0 TCTGATACACGTGT
    ## 7              A     g1             0 TGGTATCTAAACAG
    ## 8              A     g1             0 GCAGCTCTGTTTCT
    ## 9              A     g1             0 GATATAACACGCAT
    ## 10             A     g1             0 AATGTTGACAGTCA
    ## 11             A     g2             2 AGGTCATGAGTGTC
    ## 12             A     g1             2 AGAGATGATCTCGC
    ## 13             A     g2             2 GGGTAACTCTAGTG
    ## 14             A     g2             2 CATGAGACACGGGA
    ## 15             A     g2             2 TACGCCACTCCGAA
    ## 16             A     g1             2 CTAAACCTGTGCAT
    ## 17             A     g2             2 GTAAGCACTCATTC
    ## 18             A     g1             2 TTGGTACTGAATCC
    ## 19             A     g1             2 CATCATACGGAGCA
    ## 20             A     g2             2 TACATCACGCTAAC
    ## 21             B     g1             1 TTACCATGAATCGC
    ## 22             B     g1             1 ATAGGAGAAACAGA
    ## 23             B     g2             1 GCGCACGACTTTAC
    ## 24             B     g2             1 ACTCGCACGAAAGT
    ## 25             B     g1             1 ATTACCTGCCTTAT
    ## 26             B     g2             1 CCCAACTGCAATCG
    ## 27             B     g2             1 AAATTCGAATCACG
    ## 28             B     g2             1 CCATCCGATTCGCC
    ## 29             B     g2             1 TCCACTCTGAGCTT
    ## 30             B     g1             1 CATCAGGATGCACA
    ## 31             A     g1             0 CTAAACCTCTGACA
    ## 32             A     g1             2 GATAGAGAAGGGTG
    ## 33             A     g1             0 CTAACGGAACCGAT
    ## 34             A     g2             0 AGATATACCCGTAA
    ## 35             A     g1             0 TACTCTGAATCGAC
    ## 36             A     g1             0 GCGCATCTTGCTCC
    ## 37             A     g2             0 GTTGACGATATCGG
    ## 38             A     g1             0 ACAGGTACTGGTGT
    ## 39             A     g1             0 GGCATATGCTTATC
    ## 40             A     g2             0 CATTACACCAACTG
    ## 41             A     g1             0 TAGGGACTGAACTC
    ## 42             A     g2             2 GCTCCATGAGAAGT
    ## 43             A     g2             0 TACAATGATGCTAG
    ## 44             A     g2             0 CTTCATGACCGAAT
    ## 45             A     g1             2 CTGCCAACAGGAGC
    ## 46             A     g2             2 TTGCATTGAGCTAC
    ## 47             A     g1             0 AAGCAAGAGCTTAG
    ## 48             A     g2             0 CGGCACGAACTCAG
    ## 49             A     g1             0 GGTGGAGATTACTC
    ## 50             A     g2             0 GGCCGATGTACTCT
    ## 51             B     g1             1 CGTAGCCTGTATGC
    ## 52             B     g2             1 TGAGCTGAATGCTG
    ## 53             B     g2             2 CCTATAACGAGACG
    ## 54             B     g2             1 ATAAGTTGGTACGT
    ## 55             B     g1             1 AAGCGACTTTGACG
    ## 56             A     g1             1 ACCAGTGAATACCG
    ## 57             B     g1             1 ATTGCACTTGCTTT
    ## 58             B     g1             1 CTAGGTGATGGTTG
    ## 59             B     g2             1 GCACTAGACCTTTA
    ## 60             A     g1             0 CATGCGCTAGTCAC
    ## 61             A     g1             2 TTGAGGACTACGCA
    ## 62             B     g1             1 ATACCACTCTAAGC
    ## 63             A     g1             2 CATATAGACTAAGC
    ## 64             B     g1             1 TTTAGCTGTACTCT
    ## 65             B     g1             2 GACATTCTCCACCT
    ## 66             B     g2             1 ACGTGATGCCATGA
    ## 67             B     g2             1 ATTGTAGATTCCCG
    ## 68             B     g1             1 GATAGAGATCACGA
    ## 69             B     g1             1 AATGCGTGGACGGA
    ## 70             A     g1             2 GCGTAAACACGGTT
    ## 71             A     g2             0 ATTCAGCTCATTGG
    ## 72             A     g1             0 GGCATATGGGGAGT
    ## 73             A     g2             0 ATCATCTGACACCA
    ## 74             A     g2             0 GTCATACTTCGCCT
    ## 75             A     g1             0 TTACGTACGTTCAG
    ## 76             A     g1             0 GAGTTGTGGTAGCT
    ## 77             A     g2             0 GACGCTCTCTCTCG
    ## 78             A     g1             0 AGTCTTACTTCGGA
    ## 79             A     g2             0 GGAACACTTCAGAC
    ## 80             B     g1             1 CTTGATTGATCTTC

### Slicing

Slices of data can be read by passing coordinates to the `read()`
method. Before we do that, let’s take a look at the schema of `obs`:

``` r
obs$schema()
```

    ## Schema
    ## soma_joinid: int64 not null
    ## orig.ident: dictionary<values=string, indices=int8>
    ## nCount_RNA: double
    ## nFeature_RNA: int32
    ## RNA_snn_res.0.8: dictionary<values=string, indices=int8>
    ## letter.idents: dictionary<values=string, indices=int8>
    ## groups: large_string
    ## RNA_snn_res.1: dictionary<values=string, indices=int8>
    ## obs_id: large_string

With any SOMA object, you can only slice across an indexed column (a
“dimension” in TileDB parlance). You can use
[`dimnames()`](https://rdrr.io/r/base/dimnames.html) to retrieve the
names of any SOMA object’s indexed dimensions:

``` r
obs$dimnames()
```

    ## [1] "soma_joinid"

In this case, there is a single dimension called `soma_joinid`. From the
schema above we can see this contains integers.

Let’s look at a few ways to slice the dataframe.

Select a single row:

``` r
obs$read(coords = 0)$concat()
```

    ## Table
    ## 1 rows x 9 columns
    ## $soma_joinid <int64 not null>
    ## $orig.ident <dictionary<values=string, indices=int8>>
    ## $nCount_RNA <double>
    ## $nFeature_RNA <int32>
    ## $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
    ## $letter.idents <dictionary<values=string, indices=int8>>
    ## $groups <large_string>
    ## $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
    ## $obs_id <large_string>

Select multiple, non-contiguous rows:

``` r
obs$read(coords = c(0, 2))$concat()
```

    ## Table
    ## 2 rows x 9 columns
    ## $soma_joinid <int64 not null>
    ## $orig.ident <dictionary<values=string, indices=int8>>
    ## $nCount_RNA <double>
    ## $nFeature_RNA <int32>
    ## $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
    ## $letter.idents <dictionary<values=string, indices=int8>>
    ## $groups <large_string>
    ## $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
    ## $obs_id <large_string>

Select multiple, contiguous rows:

``` r
obs$read(coords = 0:4)$concat()
```

    ## Table
    ## 5 rows x 9 columns
    ## $soma_joinid <int64 not null>
    ## $orig.ident <dictionary<values=string, indices=int8>>
    ## $nCount_RNA <double>
    ## $nFeature_RNA <int32>
    ## $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
    ## $letter.idents <dictionary<values=string, indices=int8>>
    ## $groups <large_string>
    ## $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
    ## $obs_id <large_string>

### Selecting columns

As TileDB is a columnar format, it is possible to select a subset of
columns to read by using the `column_names` argument:

``` r
obs$read(coords = 0:4, column_names = c("obs_id", "groups"))$concat()
```

    ## Table
    ## 5 rows x 2 columns
    ## $obs_id <large_string>
    ## $groups <large_string>

### Filtering

In addition to slicing by coordinates you can also apply filters to the
data using the `value_filter` argument. These expressions are pushed
down to the TileDB engine and efficiently applied to the data on disk.
Here are a few examples.

Identify all cells in the `"g1"` group:

``` r
obs$read(value_filter = "groups == 'g1'")$concat()$to_data_frame()
```

    ##    soma_joinid    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
    ## 1            1 SeuratProject         85           52               0
    ## 2            5 SeuratProject         70           48               0
    ## 3            6 SeuratProject         64           36               0
    ## 4            7 SeuratProject         72           45               0
    ## 5            8 SeuratProject         52           36               0
    ## 6            9 SeuratProject        100           41               0
    ## 7           11 SeuratProject        191           61               0
    ## 8           15 SeuratProject        168           44               0
    ## 9           17 SeuratProject        135           45               0
    ## 10          18 SeuratProject         79           43               0
    ## 11          20 SeuratProject        298           65               1
    ## 12          21 SeuratProject        406           74               1
    ## 13          24 SeuratProject        463           77               1
    ## 14          29 SeuratProject        353           80               1
    ## 15          30 SeuratProject        246           59               0
    ## 16          31 SeuratProject        115           51               0
    ## 17          32 SeuratProject        189           53               0
    ## 18          34 SeuratProject        156           48               0
    ## 19          35 SeuratProject        164           47               0
    ## 20          37 SeuratProject        151           59               0
    ## 21          38 SeuratProject        126           53               0
    ## 22          40 SeuratProject        156           60               0
    ## 23          44 SeuratProject        146           47               0
    ## 24          46 SeuratProject        126           48               0
    ## 25          48 SeuratProject        204           52               0
    ## 26          50 SeuratProject        371           75               1
    ## 27          54 SeuratProject        443           77               1
    ## 28          55 SeuratProject        417           75               0
    ## 29          56 SeuratProject        502           81               1
    ## 30          57 SeuratProject        324           76               1
    ## 31          59 SeuratProject        443           81               0
    ## 32          60 SeuratProject        787           88               0
    ## 33          61 SeuratProject        612           69               1
    ## 34          62 SeuratProject        286           68               0
    ## 35          63 SeuratProject        462           86               1
    ## 36          64 SeuratProject        872           96               1
    ## 37          67 SeuratProject        328           72               1
    ## 38          68 SeuratProject        389           73               1
    ## 39          69 SeuratProject        754           83               0
    ## 40          71 SeuratProject        172           29               0
    ## 41          74 SeuratProject        228           39               0
    ## 42          75 SeuratProject        527           47               0
    ## 43          77 SeuratProject        157           29               0
    ## 44          79 SeuratProject        233           76               1
    ##    letter.idents groups RNA_snn_res.1         obs_id
    ## 1              A     g1             0 CATGGCCTGTGCAT
    ## 2              A     g1             0 TCTGATACACGTGT
    ## 3              A     g1             0 TGGTATCTAAACAG
    ## 4              A     g1             0 GCAGCTCTGTTTCT
    ## 5              A     g1             0 GATATAACACGCAT
    ## 6              A     g1             0 AATGTTGACAGTCA
    ## 7              A     g1             2 AGAGATGATCTCGC
    ## 8              A     g1             2 CTAAACCTGTGCAT
    ## 9              A     g1             2 TTGGTACTGAATCC
    ## 10             A     g1             2 CATCATACGGAGCA
    ## 11             B     g1             1 TTACCATGAATCGC
    ## 12             B     g1             1 ATAGGAGAAACAGA
    ## 13             B     g1             1 ATTACCTGCCTTAT
    ## 14             B     g1             1 CATCAGGATGCACA
    ## 15             A     g1             0 CTAAACCTCTGACA
    ## 16             A     g1             2 GATAGAGAAGGGTG
    ## 17             A     g1             0 CTAACGGAACCGAT
    ## 18             A     g1             0 TACTCTGAATCGAC
    ## 19             A     g1             0 GCGCATCTTGCTCC
    ## 20             A     g1             0 ACAGGTACTGGTGT
    ## 21             A     g1             0 GGCATATGCTTATC
    ## 22             A     g1             0 TAGGGACTGAACTC
    ## 23             A     g1             2 CTGCCAACAGGAGC
    ## 24             A     g1             0 AAGCAAGAGCTTAG
    ## 25             A     g1             0 GGTGGAGATTACTC
    ## 26             B     g1             1 CGTAGCCTGTATGC
    ## 27             B     g1             1 AAGCGACTTTGACG
    ## 28             A     g1             1 ACCAGTGAATACCG
    ## 29             B     g1             1 ATTGCACTTGCTTT
    ## 30             B     g1             1 CTAGGTGATGGTTG
    ## 31             A     g1             0 CATGCGCTAGTCAC
    ## 32             A     g1             2 TTGAGGACTACGCA
    ## 33             B     g1             1 ATACCACTCTAAGC
    ## 34             A     g1             2 CATATAGACTAAGC
    ## 35             B     g1             1 TTTAGCTGTACTCT
    ## 36             B     g1             2 GACATTCTCCACCT
    ## 37             B     g1             1 GATAGAGATCACGA
    ## 38             B     g1             1 AATGCGTGGACGGA
    ## 39             A     g1             2 GCGTAAACACGGTT
    ## 40             A     g1             0 GGCATATGGGGAGT
    ## 41             A     g1             0 TTACGTACGTTCAG
    ## 42             A     g1             0 GAGTTGTGGTAGCT
    ## 43             A     g1             0 AGTCTTACTTCGGA
    ## 44             B     g1             1 CTTGATTGATCTTC

Identify all cells in the `"g1"` or `"g2"` group:

``` r
obs$read(value_filter = "groups == 'g1' | groups == 'g2'")$concat()$to_data_frame()
```

    ##    soma_joinid    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
    ## 1            0 SeuratProject         70           47               0
    ## 2            1 SeuratProject         85           52               0
    ## 3            2 SeuratProject         87           50               1
    ## 4            3 SeuratProject        127           56               0
    ## 5            4 SeuratProject        173           53               0
    ## 6            5 SeuratProject         70           48               0
    ## 7            6 SeuratProject         64           36               0
    ## 8            7 SeuratProject         72           45               0
    ## 9            8 SeuratProject         52           36               0
    ## 10           9 SeuratProject        100           41               0
    ## 11          10 SeuratProject         62           31               0
    ## 12          11 SeuratProject        191           61               0
    ## 13          12 SeuratProject        101           41               0
    ## 14          13 SeuratProject         51           26               0
    ## 15          14 SeuratProject         99           45               0
    ## 16          15 SeuratProject        168           44               0
    ## 17          16 SeuratProject         67           33               0
    ## 18          17 SeuratProject        135           45               0
    ## 19          18 SeuratProject         79           43               0
    ## 20          19 SeuratProject        109           41               0
    ## 21          20 SeuratProject        298           65               1
    ## 22          21 SeuratProject        406           74               1
    ## 23          22 SeuratProject        213           48               1
    ## 24          23 SeuratProject        231           49               1
    ## 25          24 SeuratProject        463           77               1
    ## 26          25 SeuratProject         87           42               1
    ## 27          26 SeuratProject        327           62               1
    ## 28          27 SeuratProject        224           50               1
    ## 29          28 SeuratProject        361           76               1
    ## 30          29 SeuratProject        353           80               1
    ## 31          30 SeuratProject        246           59               0
    ## 32          31 SeuratProject        115           51               0
    ## 33          32 SeuratProject        189           53               0
    ## 34          33 SeuratProject        187           61               0
    ## 35          34 SeuratProject        156           48               0
    ## 36          35 SeuratProject        164           47               0
    ## 37          36 SeuratProject        221           67               0
    ## 38          37 SeuratProject        151           59               0
    ## 39          38 SeuratProject        126           53               0
    ## 40          39 SeuratProject        316           65               0
    ## 41          40 SeuratProject        156           60               0
    ## 42          41 SeuratProject        139           61               0
    ## 43          42 SeuratProject        108           44               0
    ## 44          43 SeuratProject         41           32               0
    ## 45          44 SeuratProject        146           47               0
    ## 46          45 SeuratProject        104           40               0
    ## 47          46 SeuratProject        126           48               0
    ## 48          47 SeuratProject         94           55               0
    ## 49          48 SeuratProject        204           52               0
    ## 50          49 SeuratProject         99           45               0
    ## 51          50 SeuratProject        371           75               1
    ## 52          51 SeuratProject        387           83               1
    ## 53          52 SeuratProject        139           50               1
    ## 54          53 SeuratProject         99           42               1
    ## 55          54 SeuratProject        443           77               1
    ## 56          55 SeuratProject        417           75               0
    ## 57          56 SeuratProject        502           81               1
    ## 58          57 SeuratProject        324           76               1
    ## 59          58 SeuratProject        292           71               1
    ## 60          59 SeuratProject        443           81               0
    ## 61          60 SeuratProject        787           88               0
    ## 62          61 SeuratProject        612           69               1
    ## 63          62 SeuratProject        286           68               0
    ## 64          63 SeuratProject        462           86               1
    ## 65          64 SeuratProject        872           96               1
    ## 66          65 SeuratProject        709           94               1
    ## 67          66 SeuratProject        745           84               1
    ## 68          67 SeuratProject        328           72               1
    ## 69          68 SeuratProject        389           73               1
    ## 70          69 SeuratProject        754           83               0
    ## 71          70 SeuratProject        212           38               0
    ## 72          71 SeuratProject        172           29               0
    ## 73          72 SeuratProject        168           37               0
    ## 74          73 SeuratProject        210           33               0
    ## 75          74 SeuratProject        228           39               0
    ## 76          75 SeuratProject        527           47               0
    ## 77          76 SeuratProject        202           30               0
    ## 78          77 SeuratProject        157           29               0
    ## 79          78 SeuratProject        150           30               0
    ## 80          79 SeuratProject        233           76               1
    ##    letter.idents groups RNA_snn_res.1         obs_id
    ## 1              A     g2             0 ATGCCAGAACGACT
    ## 2              A     g1             0 CATGGCCTGTGCAT
    ## 3              B     g2             0 GAACCTGATGAACC
    ## 4              A     g2             0 TGACTGGATTCTCA
    ## 5              A     g2             0 AGTCAGACTGCACA
    ## 6              A     g1             0 TCTGATACACGTGT
    ## 7              A     g1             0 TGGTATCTAAACAG
    ## 8              A     g1             0 GCAGCTCTGTTTCT
    ## 9              A     g1             0 GATATAACACGCAT
    ## 10             A     g1             0 AATGTTGACAGTCA
    ## 11             A     g2             2 AGGTCATGAGTGTC
    ## 12             A     g1             2 AGAGATGATCTCGC
    ## 13             A     g2             2 GGGTAACTCTAGTG
    ## 14             A     g2             2 CATGAGACACGGGA
    ## 15             A     g2             2 TACGCCACTCCGAA
    ## 16             A     g1             2 CTAAACCTGTGCAT
    ## 17             A     g2             2 GTAAGCACTCATTC
    ## 18             A     g1             2 TTGGTACTGAATCC
    ## 19             A     g1             2 CATCATACGGAGCA
    ## 20             A     g2             2 TACATCACGCTAAC
    ## 21             B     g1             1 TTACCATGAATCGC
    ## 22             B     g1             1 ATAGGAGAAACAGA
    ## 23             B     g2             1 GCGCACGACTTTAC
    ## 24             B     g2             1 ACTCGCACGAAAGT
    ## 25             B     g1             1 ATTACCTGCCTTAT
    ## 26             B     g2             1 CCCAACTGCAATCG
    ## 27             B     g2             1 AAATTCGAATCACG
    ## 28             B     g2             1 CCATCCGATTCGCC
    ## 29             B     g2             1 TCCACTCTGAGCTT
    ## 30             B     g1             1 CATCAGGATGCACA
    ## 31             A     g1             0 CTAAACCTCTGACA
    ## 32             A     g1             2 GATAGAGAAGGGTG
    ## 33             A     g1             0 CTAACGGAACCGAT
    ## 34             A     g2             0 AGATATACCCGTAA
    ## 35             A     g1             0 TACTCTGAATCGAC
    ## 36             A     g1             0 GCGCATCTTGCTCC
    ## 37             A     g2             0 GTTGACGATATCGG
    ## 38             A     g1             0 ACAGGTACTGGTGT
    ## 39             A     g1             0 GGCATATGCTTATC
    ## 40             A     g2             0 CATTACACCAACTG
    ## 41             A     g1             0 TAGGGACTGAACTC
    ## 42             A     g2             2 GCTCCATGAGAAGT
    ## 43             A     g2             0 TACAATGATGCTAG
    ## 44             A     g2             0 CTTCATGACCGAAT
    ## 45             A     g1             2 CTGCCAACAGGAGC
    ## 46             A     g2             2 TTGCATTGAGCTAC
    ## 47             A     g1             0 AAGCAAGAGCTTAG
    ## 48             A     g2             0 CGGCACGAACTCAG
    ## 49             A     g1             0 GGTGGAGATTACTC
    ## 50             A     g2             0 GGCCGATGTACTCT
    ## 51             B     g1             1 CGTAGCCTGTATGC
    ## 52             B     g2             1 TGAGCTGAATGCTG
    ## 53             B     g2             2 CCTATAACGAGACG
    ## 54             B     g2             1 ATAAGTTGGTACGT
    ## 55             B     g1             1 AAGCGACTTTGACG
    ## 56             A     g1             1 ACCAGTGAATACCG
    ## 57             B     g1             1 ATTGCACTTGCTTT
    ## 58             B     g1             1 CTAGGTGATGGTTG
    ## 59             B     g2             1 GCACTAGACCTTTA
    ## 60             A     g1             0 CATGCGCTAGTCAC
    ## 61             A     g1             2 TTGAGGACTACGCA
    ## 62             B     g1             1 ATACCACTCTAAGC
    ## 63             A     g1             2 CATATAGACTAAGC
    ## 64             B     g1             1 TTTAGCTGTACTCT
    ## 65             B     g1             2 GACATTCTCCACCT
    ## 66             B     g2             1 ACGTGATGCCATGA
    ## 67             B     g2             1 ATTGTAGATTCCCG
    ## 68             B     g1             1 GATAGAGATCACGA
    ## 69             B     g1             1 AATGCGTGGACGGA
    ## 70             A     g1             2 GCGTAAACACGGTT
    ## 71             A     g2             0 ATTCAGCTCATTGG
    ## 72             A     g1             0 GGCATATGGGGAGT
    ## 73             A     g2             0 ATCATCTGACACCA
    ## 74             A     g2             0 GTCATACTTCGCCT
    ## 75             A     g1             0 TTACGTACGTTCAG
    ## 76             A     g1             0 GAGTTGTGGTAGCT
    ## 77             A     g2             0 GACGCTCTCTCTCG
    ## 78             A     g1             0 AGTCTTACTTCGGA
    ## 79             A     g2             0 GGAACACTTCAGAC
    ## 80             B     g1             1 CTTGATTGATCTTC

Altenratively, you can use the `%in%` operator:

``` r
obs$read(value_filter = "groups %in% c('g1', 'g2')")$concat()$to_data_frame()
```

    ##    soma_joinid    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
    ## 1            0 SeuratProject         70           47               0
    ## 2            1 SeuratProject         85           52               0
    ## 3            2 SeuratProject         87           50               1
    ## 4            3 SeuratProject        127           56               0
    ## 5            4 SeuratProject        173           53               0
    ## 6            5 SeuratProject         70           48               0
    ## 7            6 SeuratProject         64           36               0
    ## 8            7 SeuratProject         72           45               0
    ## 9            8 SeuratProject         52           36               0
    ## 10           9 SeuratProject        100           41               0
    ## 11          10 SeuratProject         62           31               0
    ## 12          11 SeuratProject        191           61               0
    ## 13          12 SeuratProject        101           41               0
    ## 14          13 SeuratProject         51           26               0
    ## 15          14 SeuratProject         99           45               0
    ## 16          15 SeuratProject        168           44               0
    ## 17          16 SeuratProject         67           33               0
    ## 18          17 SeuratProject        135           45               0
    ## 19          18 SeuratProject         79           43               0
    ## 20          19 SeuratProject        109           41               0
    ## 21          20 SeuratProject        298           65               1
    ## 22          21 SeuratProject        406           74               1
    ## 23          22 SeuratProject        213           48               1
    ## 24          23 SeuratProject        231           49               1
    ## 25          24 SeuratProject        463           77               1
    ## 26          25 SeuratProject         87           42               1
    ## 27          26 SeuratProject        327           62               1
    ## 28          27 SeuratProject        224           50               1
    ## 29          28 SeuratProject        361           76               1
    ## 30          29 SeuratProject        353           80               1
    ## 31          30 SeuratProject        246           59               0
    ## 32          31 SeuratProject        115           51               0
    ## 33          32 SeuratProject        189           53               0
    ## 34          33 SeuratProject        187           61               0
    ## 35          34 SeuratProject        156           48               0
    ## 36          35 SeuratProject        164           47               0
    ## 37          36 SeuratProject        221           67               0
    ## 38          37 SeuratProject        151           59               0
    ## 39          38 SeuratProject        126           53               0
    ## 40          39 SeuratProject        316           65               0
    ## 41          40 SeuratProject        156           60               0
    ## 42          41 SeuratProject        139           61               0
    ## 43          42 SeuratProject        108           44               0
    ## 44          43 SeuratProject         41           32               0
    ## 45          44 SeuratProject        146           47               0
    ## 46          45 SeuratProject        104           40               0
    ## 47          46 SeuratProject        126           48               0
    ## 48          47 SeuratProject         94           55               0
    ## 49          48 SeuratProject        204           52               0
    ## 50          49 SeuratProject         99           45               0
    ## 51          50 SeuratProject        371           75               1
    ## 52          51 SeuratProject        387           83               1
    ## 53          52 SeuratProject        139           50               1
    ## 54          53 SeuratProject         99           42               1
    ## 55          54 SeuratProject        443           77               1
    ## 56          55 SeuratProject        417           75               0
    ## 57          56 SeuratProject        502           81               1
    ## 58          57 SeuratProject        324           76               1
    ## 59          58 SeuratProject        292           71               1
    ## 60          59 SeuratProject        443           81               0
    ## 61          60 SeuratProject        787           88               0
    ## 62          61 SeuratProject        612           69               1
    ## 63          62 SeuratProject        286           68               0
    ## 64          63 SeuratProject        462           86               1
    ## 65          64 SeuratProject        872           96               1
    ## 66          65 SeuratProject        709           94               1
    ## 67          66 SeuratProject        745           84               1
    ## 68          67 SeuratProject        328           72               1
    ## 69          68 SeuratProject        389           73               1
    ## 70          69 SeuratProject        754           83               0
    ## 71          70 SeuratProject        212           38               0
    ## 72          71 SeuratProject        172           29               0
    ## 73          72 SeuratProject        168           37               0
    ## 74          73 SeuratProject        210           33               0
    ## 75          74 SeuratProject        228           39               0
    ## 76          75 SeuratProject        527           47               0
    ## 77          76 SeuratProject        202           30               0
    ## 78          77 SeuratProject        157           29               0
    ## 79          78 SeuratProject        150           30               0
    ## 80          79 SeuratProject        233           76               1
    ##    letter.idents groups RNA_snn_res.1         obs_id
    ## 1              A     g2             0 ATGCCAGAACGACT
    ## 2              A     g1             0 CATGGCCTGTGCAT
    ## 3              B     g2             0 GAACCTGATGAACC
    ## 4              A     g2             0 TGACTGGATTCTCA
    ## 5              A     g2             0 AGTCAGACTGCACA
    ## 6              A     g1             0 TCTGATACACGTGT
    ## 7              A     g1             0 TGGTATCTAAACAG
    ## 8              A     g1             0 GCAGCTCTGTTTCT
    ## 9              A     g1             0 GATATAACACGCAT
    ## 10             A     g1             0 AATGTTGACAGTCA
    ## 11             A     g2             2 AGGTCATGAGTGTC
    ## 12             A     g1             2 AGAGATGATCTCGC
    ## 13             A     g2             2 GGGTAACTCTAGTG
    ## 14             A     g2             2 CATGAGACACGGGA
    ## 15             A     g2             2 TACGCCACTCCGAA
    ## 16             A     g1             2 CTAAACCTGTGCAT
    ## 17             A     g2             2 GTAAGCACTCATTC
    ## 18             A     g1             2 TTGGTACTGAATCC
    ## 19             A     g1             2 CATCATACGGAGCA
    ## 20             A     g2             2 TACATCACGCTAAC
    ## 21             B     g1             1 TTACCATGAATCGC
    ## 22             B     g1             1 ATAGGAGAAACAGA
    ## 23             B     g2             1 GCGCACGACTTTAC
    ## 24             B     g2             1 ACTCGCACGAAAGT
    ## 25             B     g1             1 ATTACCTGCCTTAT
    ## 26             B     g2             1 CCCAACTGCAATCG
    ## 27             B     g2             1 AAATTCGAATCACG
    ## 28             B     g2             1 CCATCCGATTCGCC
    ## 29             B     g2             1 TCCACTCTGAGCTT
    ## 30             B     g1             1 CATCAGGATGCACA
    ## 31             A     g1             0 CTAAACCTCTGACA
    ## 32             A     g1             2 GATAGAGAAGGGTG
    ## 33             A     g1             0 CTAACGGAACCGAT
    ## 34             A     g2             0 AGATATACCCGTAA
    ## 35             A     g1             0 TACTCTGAATCGAC
    ## 36             A     g1             0 GCGCATCTTGCTCC
    ## 37             A     g2             0 GTTGACGATATCGG
    ## 38             A     g1             0 ACAGGTACTGGTGT
    ## 39             A     g1             0 GGCATATGCTTATC
    ## 40             A     g2             0 CATTACACCAACTG
    ## 41             A     g1             0 TAGGGACTGAACTC
    ## 42             A     g2             2 GCTCCATGAGAAGT
    ## 43             A     g2             0 TACAATGATGCTAG
    ## 44             A     g2             0 CTTCATGACCGAAT
    ## 45             A     g1             2 CTGCCAACAGGAGC
    ## 46             A     g2             2 TTGCATTGAGCTAC
    ## 47             A     g1             0 AAGCAAGAGCTTAG
    ## 48             A     g2             0 CGGCACGAACTCAG
    ## 49             A     g1             0 GGTGGAGATTACTC
    ## 50             A     g2             0 GGCCGATGTACTCT
    ## 51             B     g1             1 CGTAGCCTGTATGC
    ## 52             B     g2             1 TGAGCTGAATGCTG
    ## 53             B     g2             2 CCTATAACGAGACG
    ## 54             B     g2             1 ATAAGTTGGTACGT
    ## 55             B     g1             1 AAGCGACTTTGACG
    ## 56             A     g1             1 ACCAGTGAATACCG
    ## 57             B     g1             1 ATTGCACTTGCTTT
    ## 58             B     g1             1 CTAGGTGATGGTTG
    ## 59             B     g2             1 GCACTAGACCTTTA
    ## 60             A     g1             0 CATGCGCTAGTCAC
    ## 61             A     g1             2 TTGAGGACTACGCA
    ## 62             B     g1             1 ATACCACTCTAAGC
    ## 63             A     g1             2 CATATAGACTAAGC
    ## 64             B     g1             1 TTTAGCTGTACTCT
    ## 65             B     g1             2 GACATTCTCCACCT
    ## 66             B     g2             1 ACGTGATGCCATGA
    ## 67             B     g2             1 ATTGTAGATTCCCG
    ## 68             B     g1             1 GATAGAGATCACGA
    ## 69             B     g1             1 AATGCGTGGACGGA
    ## 70             A     g1             2 GCGTAAACACGGTT
    ## 71             A     g2             0 ATTCAGCTCATTGG
    ## 72             A     g1             0 GGCATATGGGGAGT
    ## 73             A     g2             0 ATCATCTGACACCA
    ## 74             A     g2             0 GTCATACTTCGCCT
    ## 75             A     g1             0 TTACGTACGTTCAG
    ## 76             A     g1             0 GAGTTGTGGTAGCT
    ## 77             A     g2             0 GACGCTCTCTCTCG
    ## 78             A     g1             0 AGTCTTACTTCGGA
    ## 79             A     g2             0 GGAACACTTCAGAC
    ## 80             B     g1             1 CTTGATTGATCTTC

Identify all cells in the `"g1"` group with more than more than 60
features:

``` r
obs$read(value_filter = "groups == 'g1' & nFeature_RNA > 60")$concat()$to_data_frame()
```

    ##    soma_joinid    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
    ## 1           11 SeuratProject        191           61               0
    ## 2           20 SeuratProject        298           65               1
    ## 3           21 SeuratProject        406           74               1
    ## 4           24 SeuratProject        463           77               1
    ## 5           29 SeuratProject        353           80               1
    ## 6           50 SeuratProject        371           75               1
    ## 7           54 SeuratProject        443           77               1
    ## 8           55 SeuratProject        417           75               0
    ## 9           56 SeuratProject        502           81               1
    ## 10          57 SeuratProject        324           76               1
    ## 11          59 SeuratProject        443           81               0
    ## 12          60 SeuratProject        787           88               0
    ## 13          61 SeuratProject        612           69               1
    ## 14          62 SeuratProject        286           68               0
    ## 15          63 SeuratProject        462           86               1
    ## 16          64 SeuratProject        872           96               1
    ## 17          67 SeuratProject        328           72               1
    ## 18          68 SeuratProject        389           73               1
    ## 19          69 SeuratProject        754           83               0
    ## 20          79 SeuratProject        233           76               1
    ##    letter.idents groups RNA_snn_res.1         obs_id
    ## 1              A     g1             2 AGAGATGATCTCGC
    ## 2              B     g1             1 TTACCATGAATCGC
    ## 3              B     g1             1 ATAGGAGAAACAGA
    ## 4              B     g1             1 ATTACCTGCCTTAT
    ## 5              B     g1             1 CATCAGGATGCACA
    ## 6              B     g1             1 CGTAGCCTGTATGC
    ## 7              B     g1             1 AAGCGACTTTGACG
    ## 8              A     g1             1 ACCAGTGAATACCG
    ## 9              B     g1             1 ATTGCACTTGCTTT
    ## 10             B     g1             1 CTAGGTGATGGTTG
    ## 11             A     g1             0 CATGCGCTAGTCAC
    ## 12             A     g1             2 TTGAGGACTACGCA
    ## 13             B     g1             1 ATACCACTCTAAGC
    ## 14             A     g1             2 CATATAGACTAAGC
    ## 15             B     g1             1 TTTAGCTGTACTCT
    ## 16             B     g1             2 GACATTCTCCACCT
    ## 17             B     g1             1 GATAGAGATCACGA
    ## 18             B     g1             1 AATGCGTGGACGGA
    ## 19             A     g1             2 GCGTAAACACGGTT
    ## 20             B     g1             1 CTTGATTGATCTTC

## SOMA SparseNDArray

For `SOMASparseNDArray`, let’s consider the `X` layer containing the
`"counts"` data:

``` r
counts <- experiment$ms$get("RNA")$X$get("counts")
counts
```

    ## <SOMASparseNDArray>
    ##   uri: file:///tmp/RtmpPSFqT3/soma-exp-pbmc-small/ms/RNA/X/counts
    ##   dimensions: soma_dim_0, soma_dim_1 
    ##   attributes: soma_data

Similar to `SOMADataFrame`, we can load the data into memory as an Arrow
Table:

``` r
counts$read()$tables()$concat()
```

    ## Table
    ## 4456 rows x 3 columns
    ## $soma_dim_0 <int64 not null>
    ## $soma_dim_1 <int64 not null>
    ## $soma_data <double not null>

Or as a \$\$\`Matrix::sparseMatrix()\`\$\$:

``` r
counts$read()$sparse_matrix()$concat()
```

    ## 80 x 230 sparse Matrix of class "dgTMatrix"
    ##                                                                                
    ##  [1,] . 1  .   . .  1 . .  3 . . 1 . . . . . . . . . .  1 . . . .  . . .  4 . .
    ##  [2,] . .  .   1 .  . . .  7 . . . . . . . . . . . . 1  1 . 2 . 1  . . .  4 3 1
    ##  [3,] . .  .   . .  . . . 11 . . 1 . . . . . . . . . .  . 1 . . .  . . .  4 2 .
    ##  [4,] . .  .   . .  . . . 13 . . 1 . . . . . . . . . .  6 . . . .  . . .  5 2 1
    ##  [5,] . .  .   1 .  . . .  3 . . . . . . . . . . . . .  . . . . .  . . .  4 3 .
    ##  [6,] . .  .   1 .  . . .  4 . . . . 1 . . . . . . . .  2 1 . . .  . . .  4 1 1
    ##  [7,] . .  .   . .  . . .  6 . . . . . . . . . . . . .  4 . . . .  . . .  3 1 1
    ##  [8,] . .  .   1 .  . . .  4 . . . . . . . . . . . . .  1 1 . . .  . . .  2 3 .
    ##  [9,] . .  .   . .  . . .  2 . . . . . . . . . . . . .  . . . . .  . . .  2 2 .
    ## [10,] . 1  .   . .  . . . 21 . . 1 . . . . . . . . . .  4 . 1 . .  . . .  2 1 1
    ## [11,] 2 2  .  14 3  1 3 .  2 . . . 1 . 3 . . . . 1 1 1  2 2 . 2 .  1 1 1  . . .
    ## [12,] 2 4  5  28 .  6 1 4  9 2 1 3 . 1 . 3 1 1 . . . .  . . 1 2 4  . . .  . . .
    ## [13,] 4 3  2  18 2  2 . 1  2 . 1 2 1 . 1 . . . . . 1 .  4 . 1 1 . 15 . .  . . .
    ## [14,] 4 3  2   7 4  2 . 1  4 1 . . . 1 . 1 2 . . . 1 .  1 . . . .  . 1 .  . . .
    ## [15,] 2 2  5  15 .  2 2 2  4 1 . 1 . . 1 . . 1 . . . 2  . 2 2 . 1  . . .  . . .
    ## [16,] 3 3  8  28 .  8 . 2  . 1 . 2 . 3 1 . 1 . 1 1 . .  . . 3 . 1  . . 1  . 2 .
    ## [17,] 3 1  1   7 3  2 2 1  3 . 1 2 1 . . . . . 2 1 . .  4 . . 1 .  . . 1  . . .
    ## [18,] 4 2  5  26 3  2 1 2  6 . 1 . 2 . . 1 . 1 1 1 . 2  2 . 1 . 1 23 1 .  . . .
    ## [19,] 2 2  5  10 3  1 1 .  5 1 . 1 2 1 . . 1 . . . 1 1  6 1 . 1 .  . . .  . . .
    ## [20,] 3 5 12  16 2  2 2 1  7 . 1 . . . . 2 . 3 1 . . .  2 1 1 1 .  . . 1  . . .
    ## [21,] . .  .   7 .  . 1 .  1 . . . . . . . . . . . . .  3 . . . .  . . .  . . .
    ## [22,] . .  .  22 .  3 . 1  . . . . . . . . . . . . . .  . . . . .  . . .  . 2 .
    ## [23,] . .  1   . .  . 1 .  . . . 1 1 . . . . . . . . .  . . . . .  . . .  . . .
    ## [24,] . .  .  10 .  . . 1  1 . . . . . . . . . . . . .  . 1 . . .  . . .  . . .
    ## [25,] 1 .  .   6 .  . . .  1 . . . 1 . . . . . . . . .  . . . . .  . . .  1 1 .
    ## [26,] . .  .   . .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [27,] . .  .   4 .  1 . .  . . . . . . . . . . . . . .  1 . . . .  . . .  . . .
    ## [28,] . .  1   3 .  . . .  . . . . . . . . . . . . . .  . 1 . . .  . . .  . . .
    ## [29,] . .  .   7 .  1 . 2  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [30,] . 1  .  13 .  . . .  1 . . . 1 . . . . . . . . .  . . 1 . .  . . .  . . .
    ## [31,] . 1  .   . .  . . .  1 . . . . . . . . . . . . 1  . . . . .  . . .  . 1 .
    ## [32,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . 1 . . .  . . .  7 . .
    ## [33,] . 2  .   . .  . . .  . . . . . . . . . . . . . .  . . . 1 .  . . .  . 1 .
    ## [34,] . .  .   . .  . . .  1 . . . . . . . . . . . . .  . 1 . . .  . . .  . 1 .
    ## [35,] . .  .   1 .  . . .  . . . 1 . . . . . . . . . .  . . . 1 .  . . .  . . .
    ## [36,] . .  .   . .  . . .  . . . 1 . 1 . . . . . . . .  1 . . . .  . . .  . . .
    ## [37,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  1 . .
    ## [38,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . . 1 . .  . . .  . 2 .
    ## [39,] . .  .   . .  . . .  . . . . . 1 . . . . . . . .  4 . . . .  . . .  1 . .
    ## [40,] . .  .   . .  . . .  . . . . . . . . . . . . . .  . . . 1 .  . . .  . . .
    ## [41,] . .  .   . .  . . .  . . . . . . . . . . . . . .  7 . . . .  . . .  . 1 .
    ## [42,] . .  .   . .  . . .  1 . . . . . . . . . . . . .  1 . . . .  . . .  2 . .
    ## [43,] . .  .   . .  . . .  1 . . . . . . . . . . . . .  3 . . . .  . . .  3 . .
    ## [44,] . .  .   . .  . . .  1 . . . . . . . . . . . . 1  . . . . .  . . .  . . .
    ## [45,] . .  .   . .  . . .  7 . . 1 . . . . . . . . . .  6 1 1 . .  . . .  3 . .
    ## [46,] . .  .   1 .  . . .  1 . . . . . . . . . . . . .  1 . . 3 .  . . . 15 . .
    ## [47,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  1 1 .
    ## [48,] . .  .   1 .  . . .  1 . . . . 1 . . . . . . . .  1 . . . .  . . .  3 . .
    ## [49,] . .  .   . .  . . .  5 . . . . . . . . . . . . .  . . . . .  . . .  6 . .
    ## [50,] . .  .   . .  . . .  3 . . . . . . . . . . . . .  1 . . . .  . . .  4 1 .
    ## [51,] . 1  .  10 .  . . .  1 . . . . . 1 . . . . . . .  . 1 . . .  . . .  . . .
    ## [52,] . .  .  10 .  1 . .  2 . . . . . . . . . . . . .  1 . . . .  . . .  . 2 .
    ## [53,] . 1  .   4 .  1 . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [54,] . 1  .   1 .  . . 1  . . . . . . . . . . . . . .  1 . . . .  . . .  . . .
    ## [55,] . 2  .   6 .  2 . 1  1 . . 1 . . 1 . . . . . . 1  . . . . .  . . .  . . .
    ## [56,] . 2  .  28 .  . 1 .  1 . . . . . . . . . . . . 1  . . . . .  . . .  . 1 .
    ## [57,] . .  .  10 .  . . .  1 . . 1 . . . . . . . . . 1  . . . . .  . . .  . 1 .
    ## [58,] . .  .  13 .  1 . 1  1 . . . 1 . . . . . . . . 1  . 1 . . .  . . .  . . .
    ## [59,] . 3  .   5 .  1 1 .  2 . . . . . . . . . . . . .  1 . . . .  . . .  . . .
    ## [60,] . .  .   8 .  . . 1  1 . . . . . . . . . . . . .  2 . 1 . .  . . .  . . .
    ## [61,] . .  . 108 . 21 . 3  . . . 1 . . 1 . . . . . . 2 12 . . 1 .  . . .  . 1 .
    ## [62,] . .  .  93 . 21 . 2  1 . . . . . . . . . . . . .  3 . 2 . .  . . .  . 1 .
    ## [63,] . .  .  41 .  3 . 1  . . . . . . . . . . . . . .  1 . 1 1 .  . . .  . . .
    ## [64,] . 4  8  42 4  5 . 4  5 . . 3 . . . 1 . . . . 1 .  3 . 1 1 .  . . 2  . . .
    ## [65,] . 1  . 138 . 11 . 5  . . . 1 . . 1 . . . . . . .  . . . 1 .  . . .  . . .
    ## [66,] . .  .  77 . 11 . 2  . . . . . . . . . . . . . 1  1 . . . .  . . .  . 1 .
    ## [67,] . .  .  76 . 10 . 1  . . . . . . 1 . . . . . . .  2 . . . .  . . .  . . .
    ## [68,] 1 .  .  15 .  1 . 1  . . . . . . . . . . . . . .  . . . . .  . . .  . 1 .
    ## [69,] . .  1  19 .  2 . 1  1 . . . . . . . . . . . . .  . . . 1 1  . . .  . . .
    ## [70,] . .  . 104 . 11 . 5  4 . . 1 . . . . . . . . . .  2 . . . .  . . .  1 1 .
    ## [71,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [72,] . .  .   . .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [73,] . .  .   . .  . . .  1 . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [74,] . .  .   . .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [75,] . .  .   2 .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [76,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . . 1 . 1  . . .  . . .
    ## [77,] . .  .   1 .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [78,] . .  .   . .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [79,] . .  .   2 .  . . .  . . . . . . . . . . . . . .  . . . . .  . . .  . . .
    ## [80,] . .  .   7 .  1 . .  . . . . 1 . . . . . . . . .  . . . 1 .  . . .  . . .
    ##                                                                               
    ##  [1,] 1  5 . 1  . 1 . . . 3 3 . 1 3 . 1 1 1 1  2 . 1 2 1 2 1  .  . .  . 1 .  .
    ##  [2,] .  2 . 1  1 1 5 2 2 2 3 . . 1 2 . . 2 .  . . 3 2 1 . .  .  . .  . . .  1
    ##  [3,] 2  1 1 .  1 . . 1 . 1 . . 1 2 1 . 2 . .  1 1 2 2 1 . .  .  . .  . 9 .  1
    ##  [4,] 2  2 . .  1 1 . 1 2 6 4 2 . 5 4 . . . .  1 2 3 3 1 2 1  .  . .  . 8 .  1
    ##  [5,] .  2 . . 36 . 2 . 1 5 2 . . 2 3 . 2 1 . 54 . 2 2 . . .  .  . .  1 1 .  3
    ##  [6,] 1  . 2 .  . . . 2 2 3 1 1 1 4 1 1 . . .  2 1 . 1 . . 1  .  . .  . . 1  .
    ##  [7,] .  1 3 1  . . 1 1 . 4 1 . . 3 3 1 . . .  1 2 . . . 2 1  .  . .  . 3 .  .
    ##  [8,] 1 12 2 1  . 1 1 . 2 . 2 1 1 2 4 1 . . 1  1 2 . . 1 1 .  .  . .  . 3 .  .
    ##  [9,] 2  . 3 1  . . . 1 3 1 1 . 1 3 2 . 1 . .  1 . 1 3 2 1 .  .  . .  . . .  .
    ## [10,] 1  9 . 1  . 1 1 . 1 6 1 . 1 . . . . . 1  3 1 . . . 1 .  .  . .  . 3 .  .
    ## [11,] .  . . .  . . . . 1 . . . . . . . . . .  1 . . . . . .  .  . .  . . .  .
    ## [12,] .  . . .  . . . . 1 1 1 . . 2 . . . . .  2 1 . 1 . . .  .  . .  . . .  .
    ## [13,] .  . 1 .  . . . . . . 1 . . 1 . . . . .  . . . . . . .  .  . .  . . 1  .
    ## [14,] .  . . .  . . . . . . . . . 1 . . . . .  1 . . . . . .  .  . .  2 . .  .
    ## [15,] .  1 . .  . . . . . . . . 1 5 . . . . .  . . 1 . . . .  1  . .  . . .  .
    ## [16,] .  . 1 .  . . . . . . . . 1 . . . . . .  . . 1 . . . .  1  . .  . . .  .
    ## [17,] .  . . .  . . . . . 2 1 . . . . . . . .  . . . . 1 . .  .  . .  . . .  2
    ## [18,] .  . . .  . . . . . . 3 . . 1 1 . . . .  . . . . . . .  1  . .  1 . 1  .
    ## [19,] .  1 . .  . . . . . 1 1 . 1 . 1 . . . .  . . 1 . . . .  .  1 1  . . .  .
    ## [20,] .  . 1 .  . . 1 . . . 1 . . . 1 . . . .  . . 1 . . . .  .  . .  . . .  .
    ## [21,] .  1 . .  . . . . . 1 . . . . 1 1 . . .  . . . . . . .  .  . .  2 . .  .
    ## [22,] .  . . .  . . . . . 2 . . . . . . . . .  2 . . . . . .  .  . .  1 . .  .
    ## [23,] .  . . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [24,] .  . . .  . . . . . 2 . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [25,] .  . . .  . . . 1 . 1 3 . . 1 . . . . .  1 . 1 . 1 . .  1  . .  . . 1  .
    ## [26,] .  1 . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [27,] .  . . 1  . . . . . 1 . . . 1 . . . . .  . . . . 1 . .  .  . .  1 . .  .
    ## [28,] .  . . .  . . . . . . . . . 1 . . . . .  . . 1 . 1 . .  .  . .  . . .  .
    ## [29,] .  . . .  1 . . . 1 . . . . . . . . . .  . . . . . 1 .  1  . .  . . .  .
    ## [30,] .  1 1 .  . . . . . . 1 . . 1 . . 1 . .  3 . . . . . .  .  . .  1 . .  .
    ## [31,] .  2 . .  1 . . . . . . . 1 2 . . . . .  1 . . 3 . . .  1  1 . 35 . .  .
    ## [32,] 2  . . .  . . . . . . . . . 2 1 . . . .  . . . . . . .  .  3 . 14 2 .  .
    ## [33,] .  . . .  . . . . . . . . . . . . . . .  1 . . 1 . . .  .  2 . 12 . .  .
    ## [34,] .  . . .  . . . . . 2 . . . 2 1 . . . . 15 1 1 1 . 1 .  .  3 . 30 5 .  .
    ## [35,] .  . . .  1 . . . . 2 . . . . . . . . .  . . . 1 . . .  .  2 . 20 4 .  .
    ## [36,] .  . . .  . . . 1 . . . . . . . . 1 . .  . . . 3 . . .  .  4 . 27 . .  2
    ## [37,] .  . . .  . . . . . 1 . . . . . . 1 . .  . . 1 4 . 2 .  1  8 . 28 . .  .
    ## [38,] .  . . 1  . . . . . . . . . 1 . . . . .  . 1 . 2 . . .  1  6 . 10 . .  .
    ## [39,] .  . . .  . . . . . . . . . 3 1 . . . .  . . . 1 . . .  1  1 . 25 . .  1
    ## [40,] 2  . . .  . . . . . . . . 1 1 . . . . .  1 1 . 1 . . .  . 11 . 27 7 .  1
    ## [41,] .  . . .  . . . . . 2 2 . . 1 1 1 . . .  2 . . 2 . . .  .  1 1 31 8 .  .
    ## [42,] .  1 . .  . . . . . 1 2 . . . . . . . .  1 . . 1 . 1 .  1  4 . 22 5 1  1
    ## [43,] .  3 . .  . . . . . 4 . . . . 2 . . . .  3 1 1 4 . . .  1  1 2  7 5 .  1
    ## [44,] .  1 . .  . . . . . . 1 . . 1 . . . . .  1 . . . . . .  1  2 1  2 . .  2
    ## [45,] .  1 . .  . . . 1 . 4 . 1 . 2 1 . . . .  . . . 2 . . .  .  1 2  4 7 . 47
    ## [46,] .  1 . .  . . . . 1 4 . . . . 2 . . . .  1 . . . 2 . .  .  2 . 14 1 .  .
    ## [47,] 2  . . .  . . . . 1 . . . . 2 . . . . .  1 . . . . . .  1  2 . 16 6 1  1
    ## [48,] .  2 1 .  . . . 1 . . 1 . . 3 1 . . . .  1 . . . . . .  .  1 2  4 7 1  1
    ## [49,] .  . . .  2 . . . . . . . . . 5 . . . .  . . . . . . . 39  5 . 29 6 1  1
    ## [50,] .  2 1 .  . . . . . 2 . . . 1 2 . . . .  1 3 . 1 . 1 1  .  1 3  8 1 1  1
    ## [51,] .  . . 1  . . . 1 . . 1 . . . . 1 . . .  . . 1 . . 1 .  .  . .  5 . .  .
    ## [52,] .  . . .  . 1 . . . . 1 . . 2 . . 1 . .  1 . . . 1 1 .  .  . .  3 . .  2
    ## [53,] .  . . 1  . . . . . 1 . . . . . . . . .  . . . . . 1 .  .  . .  . . .  .
    ## [54,] .  . . .  . . . . . . 1 . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [55,] .  . 1 .  . . . . . . 1 . . 2 . . . . .  . . . . . . .  .  . .  . . .  .
    ## [56,] .  . . .  . . . . 2 2 . . 2 . . . . . .  1 . 1 . . . .  .  . .  . . .  .
    ## [57,] .  . . .  1 . . . 3 . . . . 1 . . 1 . .  . . . . . 1 .  1  . .  5 1 1  2
    ## [58,] .  . . .  . . . . . 1 . . . . . 1 . . .  3 . . . . . .  .  . .  . 1 .  .
    ## [59,] .  . . .  1 . . . . . . . . 2 . . . . .  1 . . . . . .  .  . .  . . .  1
    ## [60,] .  . . .  . . . . 1 1 1 . 1 1 . . . . .  . . 1 . . 1 .  2  . .  . . 1  .
    ## [61,] .  1 . .  1 . . 1 . 2 . 1 . 4 1 . . . .  . . . . . . .  3  . .  . . .  .
    ## [62,] .  1 . .  . . . . . . . . . . . 1 . . .  2 . 1 . 1 2 .  .  . .  1 . .  .
    ## [63,] .  1 . .  . . . . 1 . . . . . . . . . .  1 . 1 . . . .  .  1 .  . . .  .
    ## [64,] .  . . .  . . . . 1 5 2 . . 4 . . . . .  1 . . . . 2 .  2  . .  . . 1  .
    ## [65,] 2  1 . .  1 . . . 1 2 3 . . 2 . . . . .  3 1 4 . . . .  .  1 .  1 1 .  .
    ## [66,] 1  . . .  . . . . 1 2 1 . . 4 . . 1 . .  . . . . . 1 .  .  . .  3 . 1  1
    ## [67,] .  . . .  1 . . . . . . . . 1 . . . . .  1 . 2 . . . .  1  . .  . . .  .
    ## [68,] .  1 . .  . . . . 1 1 . . . . . . . . .  5 . . . . . .  .  . .  1 . .  .
    ## [69,] .  . . .  . . . . . 2 . . . . . . . . . 13 . . . . 1 .  .  1 .  . . .  .
    ## [70,] .  . . .  2 . . . . . . . . 1 1 . . . .  2 . . 1 . . .  .  . 1  1 . .  .
    ## [71,] .  . . .  . . . . . 1 . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [72,] .  . . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [73,] .  . . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  1
    ## [74,] .  . . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [75,] .  . . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [76,] .  1 . .  . . . . . 1 . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [77,] .  . . .  . . . . . . . . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [78,] .  . . .  . . . . . . 1 . . . . . . . .  . . . . . . .  .  . .  . . .  .
    ## [79,] .  . . .  . . . . . . . . . . 1 . . . .  . . . . . . .  .  . .  . . .  .
    ## [80,] .  . . .  . . . . 1 1 1 . . . . . . . .  . . . . . . 1  .  . .  . 1 .  .
    ##                                                                               
    ##  [1,] . .  . . . .  . . . . . 1 . . . 1  . . . . . .  .  .   1 .  1  . .  .  .
    ##  [2,] . 3  . . 2 .  2 . . 1 . 2 . . . .  . . . . . .  .  1   1 .  1  . .  .  .
    ##  [3,] . 2  . . 2 .  . . . . . 2 . . . .  . . . . 1 .  .  1   1 .  .  . .  1  .
    ##  [4,] . .  2 . 4 .  . . 1 . . . 1 . . .  . . . 1 . 1  .  .   . .  .  2 .  1  .
    ##  [5,] 1 1  . . 1 .  . 1 . . . . . . . .  . . . . 1 2  .  .   . .  .  . .  .  1
    ##  [6,] . 1  . . 1 .  . . . . . 1 1 . . .  . . . . . .  .  .   1 .  .  . .  .  .
    ##  [7,] . 2  1 . . .  . . . . . . . . . .  . . . . . .  .  .   . .  .  1 .  .  .
    ##  [8,] . .  . 1 . .  1 . . . . . 1 . . .  1 . . . . .  .  1   . .  .  . .  .  .
    ##  [9,] . .  . . . .  . . . . . 2 2 . 1 .  . . . . 1 .  .  .   1 .  1  1 .  .  .
    ## [10,] . 2  . . 6 .  2 . . . . 1 . . . .  . . . . . .  .  .   . .  .  . .  4  .
    ## [11,] . .  . . . .  . . . . . . . . . .  . 1 . . . 1  .  .   1 .  1  . .  .  .
    ## [12,] 1 .  . . 1 .  . 2 . . . . 2 . . 1  . . . . 1 1  .  .   4 .  .  . .  1  .
    ## [13,] . .  . . . .  . . . . . . 1 1 . .  . . . . . .  .  .   . .  .  . .  1  .
    ## [14,] . .  . . . .  . . . . . . . . . .  . . . . . .  .  .   1 .  .  . .  .  .
    ## [15,] . .  . . 1 .  . . . . . . 1 . . 1  . . . . . .  .  .   . .  .  . .  1  .
    ## [16,] . .  . 1 1 .  . . . . 1 . . . . .  . . . . . .  .  .   . .  .  . .  .  .
    ## [17,] . .  . . . .  . . . . . . . . . .  . . . . . .  .  .   . .  .  . .  .  .
    ## [18,] . .  . . 1 .  . . . . . . 1 . . .  . . . . . .  .  .   1 .  .  . .  .  .
    ## [19,] . .  . . . .  . . . . . . . . . .  . . . . . 1  .  .   1 .  .  . .  1  .
    ## [20,] . .  . . . .  . 1 . . . 2 1 . . .  . . . . . .  .  .   . .  .  . .  1  .
    ## [21,] 1 .  . . . .  . . . . . . 1 . . .  . . . . . 2 18 30  50 1 10 14 3  3  4
    ## [22,] . .  1 . . .  1 . 1 . . . . 1 . .  . . . . . .  5 12  29 2  6 13 2 13  7
    ## [23,] . .  . . . .  1 1 . . . . . . . .  . . . . . . 25 51  25 2  5  3 .  5  1
    ## [24,] . .  . . . .  1 . . . . . . . . .  . . . . . .  5 22  49 4  9 10 .  .  6
    ## [25,] . .  . . 3 .  1 . . . . . 1 . . .  . . . . . . 25 85  98 1  7 16 1 11  5
    ## [26,] . .  . . . .  . . . . . . . . . .  . . . . . .  6  3  11 .  1  4 .  .  1
    ## [27,] . .  . . 1 .  . . . . . . . 1 . .  . . . . . . 24 54  59 1  1 13 1  2  6
    ## [28,] . .  . . . .  . . . . . . . . . .  . . . . . . 40 55  28 1  2 12 .  3  4
    ## [29,] . .  . . 1 .  . . . . 1 1 . . . .  . . . . . . 16 35  34 3  8 19 1  5  5
    ## [30,] . .  1 . . .  . . . . . . 1 . . 1  . . . 1 . . 11 17  16 .  7 12 . 10  1
    ## [31,] 4 .  . . . .  . . . . . 1 . . . 1  . . . . . .  1  .   . .  .  3 .  .  .
    ## [32,] 4 1  2 . 2 .  . . . . . 1 . . . 1  . . 1 . . .  .  .   . .  .  . .  1  .
    ## [33,] 2 .  5 . . .  1 . . 1 . . . 1 . 1  . . 1 . . .  .  .   1 .  .  4 .  .  .
    ## [34,] 7 2 14 1 . .  1 . . . . 1 . . . .  . . 1 . . 2  .  1   . .  .  3 .  1  2
    ## [35,] 2 1  . . 1 .  . . . . . 1 . . 2 .  . . . . . 2  .  1   2 .  .  6 .  .  .
    ## [36,] 4 1 29 . . .  2 1 . . . . 1 . . .  . . . . . 1  1  .   . .  .  7 .  1  .
    ## [37,] 3 1  1 . 1 .  . . . . . . 1 . . .  . . . 2 . 1  .  1   . .  .  3 .  .  .
    ## [38,] 3 2  7 . . .  1 . . . . 1 1 . 3 .  . . . . . 1  .  .   1 .  1  4 .  .  .
    ## [39,] 2 1  5 . 1 .  1 . . . . . 2 . . .  . . . . . 1  .  .   . .  .  5 .  .  .
    ## [40,] 5 . 25 . 2 .  1 1 . . . 1 1 . 2 2  . . . . 1 .  .  .   . .  . 15 .  .  1
    ## [41,] 2 1  . 1 1 1  1 1 1 . . . 1 . . 1  1 . 1 . . 2  .  .   . .  .  2 .  5  .
    ## [42,] 3 1 14 . 2 .  . 2 . . 1 1 2 1 1 1  1 . . 1 . 1  .  .   1 .  .  . .  1  .
    ## [43,] 1 2 27 1 1 .  2 . . 3 . 3 1 . . .  1 1 . . 1 2  .  .   1 .  .  1 .  1  .
    ## [44,] 1 .  3 . 1 .  1 1 . . . . 1 1 . 1  . 1 . 1 1 .  .  .   . .  .  1 .  .  .
    ## [45,] . 1 13 . 1 1  . 1 1 1 . 1 1 . . 1  . 1 2 1 2 2  .  .   . .  .  . .  .  .
    ## [46,] 2 2 17 1 2 .  . 1 . 3 . . . 1 . 1  . 2 1 . . 1  .  .   . .  .  . .  .  .
    ## [47,] 8 1  7 . 4 .  1 . 1 3 . . 1 . 1 1  . . . . . .  .  .   . .  1  . .  1  .
    ## [48,] 4 1  3 1 1 .  1 1 . . 1 7 2 1 4 .  . . 1 1 . 1  .  .   1 .  .  . .  .  .
    ## [49,] 5 1 16 1 2 1 17 1 . . 1 . . 1 1 1  . . . 1 . 2  .  .   . .  .  . .  1  .
    ## [50,] 2 1 12 . 4 1  . . . . 2 1 2 . 1 . 13 . . . 1 .  .  .   . .  2  . .  .  .
    ## [51,] . .  3 . 2 .  . 1 . . . . . . . .  . . . . . .  2 20  41 . 13 11 .  2  6
    ## [52,] . .  1 . 1 .  . . . . . . 1 . . .  . . 1 . . 1  2  6   4 .  7 21 .  2  5
    ## [53,] . .  . . . .  1 . . . . . . . . .  . . . . . .  .  1   3 .  5  2 .  2  1
    ## [54,] . .  . . . .  . . . . . . . . . .  . . . . . .  4  .   3 .  1  5 .  .  1
    ## [55,] . .  . . 2 .  1 . . . . . . . . .  . . . . . .  3 10  14 .  4 21 .  2  6
    ## [56,] . .  . . . .  . . . . . . . . . .  . . 1 . . 2  .  4  17 1  3 13 .  1  4
    ## [57,] . 1  . . 3 .  1 . . . . . 1 . . .  . . . . . 1  1  8   7 .  1 16 .  1  3
    ## [58,] . .  . . 1 .  . . . . . . 1 . . .  . . . . . 1  1  6   6 .  1  9 .  1  2
    ## [59,] . .  . . 3 .  2 . . . . . . . . 2  . . . . . .  2  .   9 .  2 16 .  2  4
    ## [60,] . .  . . 1 .  . . . . . . 1 . . .  . . . . 1 1  .  .   6 .  . 17 .  9  5
    ## [61,] . .  . 2 . .  1 1 . . . . . 1 . 1  1 . . 1 1 1  .  .  76 .  .  2 .  2  1
    ## [62,] . .  1 . 2 .  . . . . . . . . . .  . . . . . .  2  .  20 .  .  8 .  2  3
    ## [63,] . .  . . 3 .  . . . . . 1 1 . . 1  . . . . . .  .  1  24 .  3  6 .  .  2
    ## [64,] . .  . . . .  1 1 . . . . 1 . . .  . . . . . .  2 10  79 2  1  9 1  1  5
    ## [65,] 1 .  1 1 1 .  1 2 . . . . . . . 1  . . 1 . . 1  1  .  53 2  2 11 .  1 14
    ## [66,] 1 .  1 . 3 .  . 2 . . . . . 1 . .  . . . . . 4  9 41  53 1  4 14 1  6 11
    ## [67,] . .  . . 2 .  . . . . . . 2 1 . .  . . . . 1 .  1 11  87 1  6 10 .  1  3
    ## [68,] . .  . . . .  1 . . . . . . . . .  . . . . . 1 23 32  76 1  1 10 .  3  4
    ## [69,] . .  . . 4 .  . . . . . . . . . .  . . . . . .  4 17  42 .  .  6 1  2  8
    ## [70,] 1 .  . . 2 .  . . . 1 . 2 . . . .  . . . . . 2  .  . 114 .  .  7 1  4  4
    ## [71,] . .  8 . . .  . . . . . . . . . .  . . . . . .  .  3   3 .  .  . .  .  .
    ## [72,] . .  5 . . .  . . . . . . . . . .  . . . . . .  1  .   1 .  .  . .  .  .
    ## [73,] . .  4 . . .  . . . . . . . . . 1  . . . . . .  .  .   1 .  .  . .  .  .
    ## [74,] . . 10 . . .  . . . . . . . . . .  . . . 1 . .  .  .   . .  .  . .  .  .
    ## [75,] . . 11 . . .  . . . . . . . . . .  . . . . . .  .  .   1 .  .  . .  .  .
    ## [76,] . . 30 . . .  . . . . . . . . . .  . . . . . .  1  .   . .  .  . .  .  .
    ## [77,] . .  8 . . .  . . . . . . . . . .  . . . . . .  .  .   . .  .  . .  .  .
    ## [78,] . .  5 . . .  . . . . . . . . . .  . . . . . .  .  .   . .  .  . .  .  .
    ## [79,] . .  9 . . .  . . . . . . . . . .  . . . . . .  .  .   . .  .  . .  .  .
    ## [80,] . .  2 . . .  . . . . . . . . . .  . . . . . .  2  7  22 .  . 14 .  6  2
    ##                                                                               
    ##  [1,]  1  . . .  1 . . . .  1 .  . .  . .  . . . . . 1  .  . . . .   .  . .  .
    ##  [2,]  1  . . .  . . . . .  . .  1 .  . .  1 1 1 . . .  .  . . . .   .  . .  .
    ##  [3,]  .  . . .  1 . . 1 .  . .  . .  . .  . . . . . .  .  . . 1 .   .  . .  .
    ##  [4,]  .  . . .  2 . . . .  . .  . 1  . .  . . . . . .  .  . . . .   .  1 .  .
    ##  [5,]  1  1 . .  1 . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . 1  .
    ##  [6,]  2  1 . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ##  [7,]  .  1 . .  . . . . .  . .  . .  1 .  . . . . . .  .  . . . 1   .  . .  1
    ##  [8,]  1  . . .  . . . . .  . .  . .  . .  . . . . . 1  .  . . . .   .  . .  .
    ##  [9,]  1  . . .  . . . . .  . .  . 1  . .  . . . . . .  .  . . . .   .  . .  .
    ## [10,]  1  1 . .  . . . . .  . .  . .  . .  1 . . . . .  .  . . . .   .  . .  .
    ## [11,]  1  . . .  . . . . .  . .  . .  . .  . . . . . .  3  . . . .   4  . .  1
    ## [12,]  .  1 . .  . . . 1 .  . .  . 1  . .  . . . . . .  8  . . 1 .  10  4 .  4
    ## [13,]  1  . . .  . . . . 1  . .  . .  . .  . . . . . .  2  . . . .   4  4 .  3
    ## [14,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  2  . . . .   4  1 .  .
    ## [15,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  5  . . . .   8  . .  4
    ## [16,]  .  1 . .  . . . . .  . .  . .  . .  . . . . . .  9  . . 4 .  23  8 .  8
    ## [17,]  1  1 . .  . . . 1 .  . .  . .  . .  . . . . . .  .  . . 1 .   7  1 .  1
    ## [18,]  2  . . .  . . . . .  . .  . .  . .  . . . . . .  5  . . 1 .   .  5 .  2
    ## [19,]  .  . 1 .  . . . . .  . .  . .  . .  . . . . . .  1  . . . .   4  . .  2
    ## [20,]  2  . . .  . . . . .  . .  . .  . .  . . . . . .  5  . . 1 .   6  1 .  4
    ## [21,] 15  1 . .  2 1 . 1 1  2 1  3 5 12 .  2 . . . 1 5  .  . . . .   .  1 .  .
    ## [22,]  9  2 . 1 14 1 1 3 1 27 1  4 1  6 .  . 4 . 1 . 3 13  . . 4 .  18  5 1  8
    ## [23,]  1  6 1 . 10 1 . 1 .  . .  2 .  2 1  . . . 1 . .  2  . . 1 .   1  . .  1
    ## [24,]  5  . . 4  8 . . 2 1  1 2  1 .  1 1  1 4 . . . .  1  . . 1 .   2  . .  1
    ## [25,]  7 36 2 1 11 1 1 . 1  1 1  1 1  6 2 14 4 1 . 1 3  .  . . . .   .  . 2  .
    ## [26,]  3  1 1 1  4 1 . 1 1  1 .  . 1  . 2  1 1 1 . 1 .  1  . . . .   3  . .  .
    ## [27,]  4  5 4 7  6 1 1 . .  . 1  . .  . 1  2 3 3 . 1 1  .  . . . .   .  . .  .
    ## [28,]  4  . 1 1  7 1 2 1 .  . .  2 .  . 3  . . 2 2 . 5  .  . . . .   1  . 1  .
    ## [29,] 11  3 1 1 22 . . 1 2  1 2 15 2  5 .  1 2 1 . . .  7  . . 1 .   7  . 1  4
    ## [30,]  7  5 . 2 37 . 1 3 1  1 1  2 1  2 1  2 . . 5 . 4  6  . . 2 .   7  1 .  1
    ## [31,]  .  1 . .  3 . . . .  . .  . .  . .  . . . . . .  .  . . . .   2  . .  .
    ## [32,]  .  . . .  4 . . . .  . .  . .  . .  . . . . . .  1  . . . .   4  . .  1
    ## [33,]  1  . . .  9 . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [34,]  .  . . .  6 . . 1 .  . .  . .  . .  . . . . . .  2  . . . .   .  . .  .
    ## [35,]  .  . . .  1 . . . .  . .  . .  . .  . . . . . 1  .  . . 1 .   .  . .  .
    ## [36,]  .  . . .  3 . . . .  . .  . .  . .  . 1 . . . .  .  . . . .   .  . .  .
    ## [37,]  1  1 . 1 14 . . . .  . .  . .  . .  . . . . . 1  1  . . . .   .  . .  .
    ## [38,]  .  1 . .  2 . . . .  1 .  . .  . .  . . . . . 1  .  . . . .   .  . .  .
    ## [39,]  .  1 . .  1 . . . .  . .  . .  . .  . . . . . 2  .  . . . .   .  . .  .
    ## [40,]  4  . . .  4 . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [41,]  1  1 . .  1 . . . .  . .  . .  . .  . . . . . .  .  1 . . .   .  . .  .
    ## [42,]  1  1 . .  3 . . . .  . .  1 .  . .  . . . . . .  .  . . . .   4  1 .  .
    ## [43,]  .  . . .  . . . 1 .  . .  . .  1 .  . . . . . .  .  . . . .   .  . .  .
    ## [44,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [45,]  1  2 . .  . . . . .  . .  . .  . .  . . . . . .  1  . . . .   1  . .  .
    ## [46,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  3  . . . .   2  . .  1
    ## [47,]  .  1 . .  1 . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [48,]  .  . . .  . . . 1 .  . .  . .  . .  . 1 . . . .  .  . . . .   .  . .  .
    ## [49,]  .  . . 1  1 . . . .  . .  . .  . .  . . . . . 1  1  . . . .   .  . .  .
    ## [50,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [51,]  8  2 1 4  5 . . . .  1 .  4 .  3 .  3 2 . 1 1 . 12  . . 1 .   8  . .  4
    ## [52,]  8  4 . . 12 . . . .  1 1  5 .  . .  2 2 . . . 2  4  . . . .   3  2 .  5
    ## [53,]  7  . . 1  4 . . . .  1 1  2 .  . 1  . . . . . .  2  . . . .   5  . .  .
    ## [54,]  3  1 . .  2 . . . .  . .  . .  . .  1 . . . . 1  1  . . . .   2  . 1  .
    ## [55,] 10  2 . 1 16 2 1 2 .  3 1  . 1  . .  3 4 2 1 . 1  5  . . . .   3  . .  3
    ## [56,] 15  3 . . 10 . . 2 .  6 2  5 .  . .  1 2 . 1 . 3  5  . . 4 .   7  1 .  3
    ## [57,] 18  6 . 3  6 1 . 5 1  1 1  2 .  . .  1 . 1 . . 1  7  . . . .   6  1 1  6
    ## [58,] 19  4 . .  2 . . 1 .  2 1  3 1  . .  . 2 2 . . 3 14  . . . .   5  2 1  3
    ## [59,]  4  2 . . 12 2 . 2 .  2 .  2 3  . .  1 1 . 1 . .  5  . . . .   9  . 2  6
    ## [60,] 17  5 . 2 16 1 1 . .  4 .  3 1  . .  1 . . 2 . 2 11  . . . .   4  . 2  2
    ## [61,]  5  1 . 1  8 . . 1 .  2 .  . 2  3 4  3 . . . . . 75 16 . 6 3 102 25 2 11
    ## [62,]  3  . . 1 13 . . . .  2 .  . 1 10 1  1 5 . . . 4 52  1 5 6 3  78 39 2 26
    ## [63,]  1  . . . 21 . . 1 1  2 .  . 1  1 .  . 1 1 2 . . 11  2 2 5 1  23  5 .  5
    ## [64,]  5  4 . 2  9 . . . .  2 .  . .  2 7  2 . 1 . . 2 19  4 4 4 3  25  2 1  2
    ## [65,]  .  2 . . 20 . . . 2  . .  1 3  3 7  2 6 . 1 . 6 54  8 2 6 .  69 16 1 31
    ## [66,]  3  5 . 2 10 2 . . 1  . .  2 1  4 .  3 2 1 . . 2 23  5 3 5 1  24  6 1 21
    ## [67,]  6 10 . 2 23 1 . 1 .  3 .  . 1  4 2  3 7 . . . 5 45  8 6 6 3  43 11 6 21
    ## [68,]  2  6 . 1  5 . . . 1  1 1  . 1  1 1  1 2 1 . . 1 10  4 4 5 .   8  3 3  2
    ## [69,]  .  4 . . 28 . . . .  . .  . 1  3 2  4 2 . . . . 23  7 2 3 1  10  4 5  3
    ## [70,]  3  2 . . 13 . . . .  2 .  . 1  6 .  3 . . . . . 37  . 1 5 2  50  9 1 10
    ## [71,]  1  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   1  . .  .
    ## [72,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [73,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [74,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [75,]  .  . . 1  1 . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [76,]  .  2 . .  . . . . 1  . .  . .  . .  . 1 . . . .  .  . . . .   .  . .  .
    ## [77,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [78,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [79,]  .  . . .  . . . . .  . .  . .  . .  . . . . . .  .  . . . .   .  . .  .
    ## [80,]  3  3 . 1 10 . . . .  . .  3 1  3 1  2 1 1 . . 4  5  . 1 1 .   5  . .  1
    ##                                                                              
    ##  [1,]  .   . .  . . .  .  1  . .  .  .  . .  .  2 .  .  2  .  .  .  .  . .  .
    ##  [2,]  .   . .  . . .  .  1  . .  .  .  . 1  .  . .  .  .  .  .  .  1  . .  .
    ##  [3,]  .   . .  . . .  1  1  . .  1  .  . .  .  1 .  .  1  2  .  2  .  1 .  1
    ##  [4,]  .   . .  . . .  .  1  3 .  .  .  . .  .  2 .  3  .  .  1  1  .  1 .  .
    ##  [5,]  .   . .  1 . .  1  2  . .  1  .  1 .  .  3 .  2  .  3  1  .  .  . .  .
    ##  [6,]  .   1 .  . . .  .  .  . .  1  .  1 .  .  2 .  1  .  2  .  1  1  . .  .
    ##  [7,]  .   . .  . . .  .  .  . .  .  .  . .  .  2 .  .  2  .  .  .  1  . .  .
    ##  [8,]  .   . .  . . .  .  1  . .  .  .  . .  .  1 .  .  1  .  1  .  1  . .  .
    ##  [9,]  .   . .  . . .  .  .  1 .  .  .  . .  1  . .  .  .  .  .  .  1  . .  .
    ## [10,]  .   . .  . . .  .  2  . .  .  .  1 .  .  3 1  .  .  .  1  .  2  . .  .
    ## [11,]  2   . .  . . .  .  .  . 1  .  .  . .  .  . .  .  .  .  .  .  .  . .  .
    ## [12,] 10   . .  . 2 .  .  1  . 1  .  .  1 .  .  . .  .  .  .  .  1  1  . .  .
    ## [13,]  6   . .  . . .  .  .  1 .  .  .  1 .  .  1 .  .  .  .  2  1  .  . .  .
    ## [14,]  1   . .  . . .  .  .  . 1  .  .  . .  .  . .  .  .  .  .  .  .  . .  .
    ## [15,]  5   . .  . 1 .  .  .  2 .  1  .  . .  .  . .  .  1  2  1  1  2  . .  .
    ## [16,] 16   . .  . . .  3  2  . .  .  .  1 .  .  2 .  .  .  .  .  .  5  . 1  .
    ## [17,]  5   . .  . 1 .  .  .  . .  .  .  . .  .  4 .  .  .  .  .  .  .  . .  .
    ## [18,] 11   . .  . 1 .  1  .  . 2  .  .  . .  .  . .  .  .  .  1  .  .  . .  .
    ## [19,]  5   1 .  . . .  .  .  . 1  .  .  . .  .  . .  1  .  .  1  .  .  . .  .
    ## [20,]  8   1 .  . . .  .  .  . 1  1  .  . .  .  1 .  .  1  .  .  .  .  . .  .
    ## [21,]  2  13 4  . . .  2  .  . .  .  .  1 4  .  . .  3  5  6  1  .  4  3 6  1
    ## [22,] 12  28 .  . 2 .  3  .  . .  1  .  3 1  7  . .  6  7  5  .  1 15  . 4  5
    ## [23,]  1  15 .  1 . .  1  1  . .  .  .  . 1  1  . .  1  6  1  .  1  8  . .  .
    ## [24,]  5  11 .  . . .  6  .  1 .  .  .  3 .  1  . .  4  5  5  1  .  5  . 2  .
    ## [25,]  1  13 .  . . .  5  2  2 2  .  .  1 1  1  . .  8  4  3  2  .  4  3 .  .
    ## [26,]  .   7 1  . . .  1  .  . .  .  .  1 .  .  . .  3  3  2  .  .  2  . .  2
    ## [27,]  3  37 .  . . .  3  .  1 .  .  .  1 5  2  . .  5  1  1  .  .  8  1 1  4
    ## [28,]  .   5 .  . . .  1  1  1 .  .  .  . 1  .  . .  .  2  1  1  .  2  1 .  1
    ## [29,]  5  20 2  1 . .  4  2  1 1  .  .  2 .  2  3 .  7 10  6  2  2 11  1 3  2
    ## [30,]  3  18 3  . . .  2  1  . .  .  .  3 2  1  3 . 13 12  4  .  . 18  1 3  7
    ## [31,]  .   1 .  . . .  1  3  . .  .  .  1 .  .  1 .  .  1  .  2  1  3  . .  .
    ## [32,]  2   . .  . . .  2  .  . .  .  1  . .  .  4 .  .  .  1  .  .  .  1 .  .
    ## [33,]  .   . .  . . .  .  1  1 .  .  .  . .  .  5 .  .  .  2  1  1  .  . .  .
    ## [34,]  1   . .  . . .  1  .  1 .  .  .  4 .  .  1 .  .  .  .  1  .  .  . .  .
    ## [35,]  .   . .  . . .  2  3  1 .  .  .  1 .  1  1 .  .  .  1  .  .  .  . .  .
    ## [36,]  .   . .  . . .  .  .  . .  .  .  . .  .  . .  1  .  .  1  .  .  . .  .
    ## [37,]  .   1 .  . . .  1  .  1 .  .  .  4 .  .  1 .  .  .  1  .  .  1  . 1  .
    ## [38,]  .   1 .  . . .  2  1  1 .  .  .  1 .  .  . .  .  .  1  1  .  .  . .  .
    ## [39,]  1   . .  . . .  .  .  1 .  .  .  . .  .  . .  .  .  .  .  1  1  . .  .
    ## [40,]  .   . .  1 . .  1  3  1 .  .  .  1 .  .  . .  .  .  .  2  1  1  . .  .
    ## [41,]  .   . .  . . .  2  1  . .  .  1  . .  .  . .  .  .  3  .  1  .  1 .  .
    ## [42,]  .   . .  . 1 .  .  1  . .  1  .  . .  .  . .  .  .  1  1  1  1  . .  .
    ## [43,]  .   . .  1 . .  .  .  . .  .  .  1 .  .  . .  .  .  .  1  .  3  1 .  .
    ## [44,]  .   . .  . . .  .  .  . .  .  .  . .  .  2 .  .  .  .  .  .  .  . .  .
    ## [45,]  .   . .  . . .  .  .  . .  .  .  . .  .  . .  .  2  1  .  1  .  . .  .
    ## [46,]  3   . .  . . .  .  .  . .  .  .  . .  .  1 .  .  .  .  1  1  .  . .  .
    ## [47,]  .   . .  . . .  .  2  . .  .  .  1 .  .  . .  .  .  .  2  .  .  1 .  .
    ## [48,]  .   1 .  . . .  .  1  . .  1  .  1 .  .  1 .  .  .  1  2  .  .  . .  .
    ## [49,]  .   . .  . . .  1  1  . .  .  .  . .  .  3 .  .  .  2  1  1  2  . .  .
    ## [50,]  .   . .  . . .  .  1  . .  .  .  2 .  .  . .  .  .  1  1  .  1  . .  .
    ## [51,]  8  16 .  . . .  4  3  1 .  .  .  9 2  2  3 . 15  7  8  2  . 21  2 2  5
    ## [52,]  4  32 .  . . .  1  3  1 .  .  .  3 .  1  1 . 17 12  8  2  1 25  3 1  3
    ## [53,]  .   7 .  . . .  2  .  . .  .  .  1 1  3  . .  8  7  6  1  .  6 16 1  4
    ## [54,]  .   9 .  . . .  .  1  1 .  .  .  . 1  .  . . 11  6  2  .  . 10  . 1  1
    ## [55,]  7  11 .  . . .  1  2  . .  .  .  4 2  .  . . 18 32  9 50  3 26  1 3 11
    ## [56,]  7  17 .  . . .  5  2  . .  1  .  2 3  1  1 . 13 33  9  1  3 26 11 4  9
    ## [57,] 13  33 .  2 1 .  .  .  . .  .  .  3 .  1  . . 36 12 10  1  1 16  3 5  2
    ## [58,]  6  10 .  1 2 .  .  1  2 2  .  .  2 1  .  . . 17 19  8  1 27 15  5 5  5
    ## [59,]  6  15 1  . . .  1  2  3 .  .  .  . 1  3  1 . 12 18  5  3  1 11  4 3  7
    ## [60,]  4  25 .  . . .  1  .  . .  .  .  6 3  2  1 . 27 29 10  1  1 22  6 6 10
    ## [61,] 50  61 1  . 7 2  9  1  4 2  .  2  5 6  4  3 . 12  6  1  5  . 10  8 1  .
    ## [62,] 53  31 8  1 9 5  4  .  3 .  .  .  1 1  .  1 4  7  7  2  .  .  5  1 1 12
    ## [63,] 10  25 .  1 1 .  5  1  . 3  .  .  5 .  2  . 1  7  1  1  .  .  5  1 .  2
    ## [64,]  9  14 3 33 . .  7  2  1 2  1  .  1 2  .  2 .  4  3  6  .  . 16  . 3  1
    ## [65,] 68  58 1  . 6 3  2  1  . 3  1  . 22 5  1 64 1  8 11  6  2  .  2  . .  3
    ## [66,] 36 112 2  . 1 3  5  .  1 1  .  . 10 4 39  2 . 10  7  4  5  .  3  . 1  4
    ## [67,] 49  37 3  . 4 . 12  1 38 2 21  .  9 8  1  3 .  4  9  4  4  2 16  1 .  4
    ## [68,]  3  18 6  . 1 .  7  1  . .  .  1  1 2  3  1 1  2  4  2  .  2  3  1 .  .
    ## [69,]  9  29 1  1 . . 10  3  . 1  1 32  3 4  5  1 .  6  1  2  1  1  4  1 .  .
    ## [70,] 26 125 .  3 5 1 18 43  1 8  .  .  3 5  1  1 .  6  4  7  3  2  5  2 1  1
    ## [71,]  .   5 .  . . .  .  .  . 1  .  .  . .  .  . .  .  .  .  .  .  3  . .  .
    ## [72,]  .   1 .  . . .  .  .  . .  .  .  . 1  .  . .  .  .  .  .  .  4  . .  .
    ## [73,]  .   . 1  . . .  .  .  1 .  .  .  . .  .  . .  .  .  1  .  .  2  . .  .
    ## [74,]  .   . .  . . .  .  .  . .  .  .  . .  .  . .  .  .  1  .  .  6  . .  .
    ## [75,]  .   5 .  . . .  1  .  . .  .  .  . .  .  . .  .  1  1  .  .  3  . .  .
    ## [76,]  .   1 .  . . .  .  .  2 .  .  .  . .  .  . .  .  .  .  .  . 17  . .  .
    ## [77,]  .   3 .  . . .  .  .  . .  .  .  . .  .  . .  .  .  .  .  .  3  . .  .
    ## [78,]  .   . .  . . .  .  .  . .  .  .  . .  .  . .  .  .  .  .  .  6  . .  .
    ## [79,]  .   . .  . . .  .  .  . .  .  .  . .  .  . .  .  .  .  .  .  4  . .  .
    ## [80,]  4  16 5  . . .  4  3  2 .  .  .  4 .  1  . .  7  5  1  1  1  3  . 2  1
    ##                                                                              
    ##  [1,]  . .  2  . .  3 1 . . . .  . . . . . 1 .  .  .  .  .  .  .  . . .  .  .
    ##  [2,]  . .  .  . .  . . . . . .  . . . . . 1 .  .  .  .  .  .  .  . . .  .  .
    ##  [3,]  . .  1  1 .  3 . . . . .  4 . . 1 . . .  1  .  1  .  .  .  . . .  .  .
    ##  [4,]  . .  2  . .  3 . . 1 . .  2 1 . . 1 . .  .  1  .  1  .  .  . . .  .  .
    ##  [5,]  1 .  1  . .  1 . . . . .  1 . . . 1 . .  .  1  .  1  1  .  . . .  .  .
    ##  [6,]  . .  .  . .  3 . . . . .  2 . . . . . .  .  .  .  1  .  .  . . .  .  .
    ##  [7,]  . .  .  . .  3 . . . . .  . . . . . . .  1  .  .  .  .  .  . . .  .  .
    ##  [8,]  1 .  1  . .  1 . . . . .  1 . . . . . .  .  .  .  1  1  .  . . .  .  .
    ##  [9,]  . .  1  . .  . . . . . .  1 . . . . . .  1  .  .  .  .  .  . . .  .  .
    ## [10,]  . .  .  . .  3 . . . . .  3 . . . . . .  .  .  .  1  .  .  . . .  .  .
    ## [11,]  . .  .  . .  . . . . . .  . . . . . . .  .  .  .  .  .  .  . . .  .  .
    ## [12,]  1 .  1  . 1  1 . . 2 . .  2 . . . . . .  .  1  1  1  .  .  . . .  .  .
    ## [13,]  . .  .  . .  2 . . . . .  . . . . . . .  .  .  .  1  .  .  . . .  .  .
    ## [14,]  . .  .  . .  . . . . . .  . . . . . . 1  .  1  .  .  .  .  . . .  .  .
    ## [15,]  . .  .  . .  1 . . . . .  . . . . . . .  .  .  .  1  .  .  . . .  .  .
    ## [16,]  . .  .  1 .  . . . . . .  1 . . . 1 . .  .  .  .  .  .  .  . . .  .  .
    ## [17,]  . .  .  . .  . . . . . .  . . . . . . .  .  .  .  1  .  .  . . .  .  .
    ## [18,]  . .  .  . .  . . . . . .  1 . . 1 . . .  .  .  .  .  1  .  . . .  .  .
    ## [19,]  . .  .  . .  2 . . . . .  . . . . . . .  .  1  .  .  .  .  . . .  .  .
    ## [20,]  . .  .  . .  3 . . . . .  . . . . . . .  .  .  .  .  .  .  . . .  .  .
    ## [21,]  1 .  2  6 .  6 . . . 1 .  6 1 . . . 2 1  .  .  .  4  .  .  . 1 .  .  .
    ## [22,]  1 . 10  3 1  4 2 . . 1 1 15 1 2 . 1 1 .  .  1  .  5  .  .  . . .  .  .
    ## [23,]  . .  4  4 .  . . . . . .  2 . 4 . . 1 .  .  1  .  3  .  .  . . .  .  .
    ## [24,]  . .  2  6 .  . . . . . .  4 . . . . 1 .  .  .  .  5  .  .  . . .  .  .
    ## [25,]  . .  2  1 1  1 . 1 2 . .  7 . . . . 1 .  .  .  1 12  1  .  . . .  .  .
    ## [26,]  . .  2  2 .  1 . . . . .  3 . . . . 1 .  .  .  .  1  .  .  . . .  .  .
    ## [27,]  . .  1  4 .  . . . . . 1  6 . 3 . . 1 .  .  1  . 15  1  .  . . .  .  .
    ## [28,]  . .  6  4 .  1 . 2 . 1 .  . . . . . . 1  .  .  .  2  .  .  . . .  .  .
    ## [29,]  . 1  5  9 .  3 1 . . 1 3  4 . 1 2 . . .  .  .  .  3  .  .  . . .  .  .
    ## [30,]  1 .  6  8 .  6 . 2 . . . 20 . 3 1 1 2 .  .  1  .  1  .  .  . 2 1  .  .
    ## [31,]  6 .  6  8 .  8 . . . . .  . . . . 2 . .  .  .  .  .  2  .  . . .  .  .
    ## [32,]  2 .  .  . .  2 . . . . .  1 . 1 . . 2 .  .  .  .  1  1  .  . . .  .  .
    ## [33,]  2 .  .  3 .  3 . . . . .  1 . . . . 1 .  1  .  .  .  .  .  . . .  .  .
    ## [34,]  1 .  .  1 .  5 . . . . .  . 1 . . 1 . .  .  .  .  .  .  .  . . .  .  .
    ## [35,]  . .  .  2 .  2 . . . . .  . 1 . . . . .  .  .  .  .  .  .  . . .  .  .
    ## [36,]  1 .  .  5 .  1 . . . . .  . . . . 1 . .  .  .  .  .  .  .  . 1 .  .  .
    ## [37,]  2 .  1  6 1  5 . . . . .  1 . . . 3 . .  1  .  .  2  1  .  . . .  .  .
    ## [38,]  1 .  1  6 .  1 . . . . .  . . . . 1 . .  .  .  .  .  1  .  . . .  .  .
    ## [39,]  2 .  .  1 1  3 . . . . .  . . . . 1 . .  .  .  .  1  1  .  . . .  .  .
    ## [40,]  6 .  .  6 .  2 . . . . .  . . . . 2 . .  3  .  .  .  .  .  . . .  .  .
    ## [41,]  2 .  1  3 .  7 . . . . .  1 . . . . . .  .  .  .  .  .  .  . . .  .  .
    ## [42,]  . .  .  . .  4 . . . . .  2 . . . . . .  1  .  .  .  .  .  . . .  .  .
    ## [43,]  . .  .  . .  2 . . . . .  . . . . . . .  .  .  .  .  .  .  . . .  .  .
    ## [44,]  . .  .  . .  2 . . . . .  . . . . . . .  .  1  .  1  .  .  . . .  .  .
    ## [45,]  . .  1  . .  5 1 . . . .  5 . . . . . .  .  .  .  1  .  .  . . .  .  .
    ## [46,]  . .  1  . .  1 . . . . .  1 . . . . 1 .  .  .  .  .  .  .  . . .  .  .
    ## [47,]  1 .  .  . .  1 . . . . .  . . . . 1 1 .  .  .  .  1  .  .  . . .  .  .
    ## [48,]  . .  1  . .  4 . . . . .  . . . . 1 . .  .  .  .  2  .  .  . . .  .  .
    ## [49,]  . .  .  . .  1 . . . . .  2 . . . . . .  1  .  .  .  .  .  . . .  .  .
    ## [50,]  . .  .  . .  2 1 . . . .  . . . . . . .  .  .  .  .  .  .  . . .  .  .
    ## [51,]  . 1 17 12 .  5 . 3 . 2 1  9 2 6 . 1 2 1  .  .  .  5  2  .  . . .  .  .
    ## [52,]  5 1 13 12 2 10 2 3 5 4 2 20 3 4 5 6 2 . 26  .  .  3  .  .  . . .  .  .
    ## [53,]  1 3  1  2 1  1 2 . 1 2 .  9 . 1 . . . 1  .  1  .  .  .  .  . . .  .  .
    ## [54,]  2 .  2  4 .  4 1 1 . 1 1  3 . . . 1 2 .  .  1  .  .  .  .  . . .  .  .
    ## [55,] 14 .  9 35 3 17 1 2 2 1 .  6 6 4 4 1 2 1  2  2  .  1  .  .  . . .  .  .
    ## [56,]  4 1 12 16 4  8 . 1 4 3 3  9 . 3 . 2 1 6  2 25  .  1  .  .  . . .  .  .
    ## [57,] 18 1 14 24 1 33 3 6 3 . 1 91 . . 5 7 . 4  1  .  .  .  .  .  . . .  .  .
    ## [58,]  9 .  8  9 2  8 1 1 1 1 2 11 3 4 2 2 6 1  .  1  .  1  .  .  . . .  .  .
    ## [59,]  5 1  7  9 1 14 1 5 . . 3 18 . 2 1 6 3 .  .  .  .  1  .  .  . . .  .  .
    ## [60,] 11 1 13 30 1 19 2 6 2 . 5 18 3 8 4 3 6 2  1  .  .  2  1  .  . . .  .  .
    ## [61,]  . .  5  8 .  4 . 6 . . . 18 . . . 2 1 .  .  .  .  6  1  .  . . .  .  .
    ## [62,]  . .  4  8 .  7 . . . . 1  2 . . . . . 1  .  .  .  7  .  .  . . .  .  .
    ## [63,]  . .  5  3 .  4 . 3 1 . .  9 1 . . 1 . 1  .  .  .  2  .  1  . . .  .  .
    ## [64,]  . .  3  3 .  3 . 1 . 1 . 11 1 . . . 1 .  1  .  .  6  .  .  . . .  .  .
    ## [65,]  . . 11 13 .  2 . 1 . . 1 12 . 1 . 2 2 .  1  .  . 24  1  .  . . .  .  .
    ## [66,]  . .  9  8 1  2 2 3 . . 1 11 1 1 . . 1 .  1  .  . 16  2  .  . 2 .  .  .
    ## [67,]  1 .  9  7 .  . 1 3 . . .  7 . . . . . 1  1  .  . 28  1  .  . . .  .  .
    ## [68,]  . .  4  5 .  1 . . . . .  5 1 . . . 1 1  .  .  .  3  .  .  . . .  .  .
    ## [69,]  1 .  5  8 .  6 1 1 1 . 1  4 1 . . . . .  1  1  .  6  1  .  . . .  .  .
    ## [70,]  . .  2  3 2  4 . 4 1 . 4 25 1 1 . . 1 2  .  .  .  3  1  .  1 . .  .  .
    ## [71,]  . .  .  . .  . . . . . .  1 . . . . . .  .  . 43 18  4 14 11 1 8  6 14
    ## [72,]  . .  .  . .  . . . . . .  2 . . . . . .  .  . 41  8  4 11  3 5 3  5  5
    ## [73,]  . .  .  . .  . . . . . .  . . . . . . .  .  . 36 12  2 14 13 3 2  9  8
    ## [74,]  . .  .  1 .  . . . . . .  3 . . . . . .  .  . 55 18  2 18  8 3 2 10 11
    ## [75,]  . .  .  . .  1 . . . . .  . . 1 . . . .  .  . 58 18  2 23  8 2 3  7 15
    ## [76,]  . .  .  . .  1 . . . . .  2 . . . . . .  .  . 54 28 15 62 29 7 9 23  6
    ## [77,]  . .  .  . .  . . . . . .  3 . . . . . .  .  . 66 11  2  9  3 3 3 12  4
    ## [78,]  . .  .  1 .  . . . . . .  . . 1 . 1 . .  .  . 34 13  1 14  6 1 3  6  3
    ## [79,]  . .  .  . .  . . . . . .  4 . 1 . . . .  .  . 30 16  3  6  5 1 4 11  5
    ## [80,]  . .  1  4 .  1 . 1 . . 1  7 . 1 . . . 1  1  .  6  9  2  .  2 2 2  1  2
    ##                                                                              
    ##  [1,]  .  .  . .  1  . .  .  1 . . .  . . . . .  .  1  .  .  .  .  . . . .  .
    ##  [2,]  .  1  . 1  .  . .  1  2 . . .  . . . . .  .  .  .  .  .  .  . . . .  1
    ##  [3,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ##  [4,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  1  .  .  1 . . .  .
    ##  [5,]  .  .  . .  .  . .  .  2 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ##  [6,]  .  .  . .  1  . .  .  . . 1 .  . . . . .  .  .  .  .  .  .  . . . .  .
    ##  [7,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  1  .  1  . . . .  .
    ##  [8,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ##  [9,]  .  1  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [10,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  1  .  .  . . . .  .
    ## [11,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [12,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  . 29  .  .  . . . .  .
    ## [13,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . 1 . .  .
    ## [14,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [15,]  .  .  . .  1  . 1  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [16,]  .  .  . .  1  . .  .  2 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [17,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [18,]  1  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [19,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [20,]  .  .  . .  .  . 1  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [21,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  .  1  .  . . . 1  .
    ## [22,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [23,]  .  .  . .  .  . .  .  2 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [24,]  .  .  . .  .  3 .  .  3 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [25,]  .  .  . .  .  . .  .  5 . . .  . . . . .  .  .  .  1  1  .  . . . .  .
    ## [26,]  .  1  . .  .  . .  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [27,]  .  .  . .  .  . .  .  2 . . .  . . . . 1  .  .  .  .  .  .  . 1 . .  .
    ## [28,]  .  .  . 1  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [29,]  .  1  . 1  .  . .  .  3 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [30,]  .  .  . .  .  . .  .  2 . . .  . 1 . . 1  .  .  .  1  .  .  . . . .  .
    ## [31,]  .  .  . 1  .  . .  .  . . . 1  . . . . .  .  .  .  . 27  2 35 5 7 5 14
    ## [32,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  2  5  . 3 . 3  1
    ## [33,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  1  3 15 9 1 1  4
    ## [34,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  . 10  4  3 2 1 .  9
    ## [35,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  8 10 29 6 . 3  7
    ## [36,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  5  8 11 3 1 1 10
    ## [37,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  2 10 12 22 6 5 1 10
    ## [38,]  .  .  . 1  .  . .  .  . . . .  . . . . .  .  .  .  .  7 10 15 8 4 2  2
    ## [39,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  4  3 18 2 1 1  4
    ## [40,]  .  .  . .  .  . 1  .  2 . . .  . . . . .  .  .  .  1 11 13 18 5 1 1  7
    ## [41,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  1  3  1 10 4 . 3  6
    ## [42,]  .  .  . 1  .  . .  .  . . . .  . . . . 1  .  .  .  .  .  8  . 1 . . 13
    ## [43,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  2  . . . .  .
    ## [44,]  .  .  . 1  .  . .  .  . . . .  . . . . .  .  .  .  .  .  1  3 . . .  .
    ## [45,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [46,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . 2 . .  .
    ## [47,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  6  .  4 9 . .  6
    ## [48,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  3  1 . . .  .
    ## [49,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  1  2  3  3 3 . .  5
    ## [50,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  2  . . . .  3
    ## [51,]  .  .  . 1  1  . .  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [52,]  .  .  . 1  .  . 1  .  2 . . .  . . . . .  .  .  .  1  1  .  . . . .  .
    ## [53,]  .  1  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [54,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . 1 . .  .
    ## [55,]  .  .  . 1  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [56,]  .  .  . .  .  . .  .  2 . . .  . . . . 1  .  .  .  .  .  .  . . . .  .
    ## [57,]  .  1  . .  1  . .  .  2 . . .  . . . . .  .  .  .  .  .  .  1 . . .  .
    ## [58,]  .  .  . .  .  . .  .  2 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [59,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  .  .  1  . . . .  1
    ## [60,]  .  2  . 2  .  . .  .  2 . . .  . . . . .  .  .  .  2  .  .  1 . . .  .
    ## [61,]  .  1  . .  .  . 1  .  1 . . .  . . . . .  .  .  .  1  .  .  . . . .  .
    ## [62,]  .  1  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [63,]  .  .  . .  .  . .  .  . . . .  . . . . .  .  .  .  .  .  .  1 . . .  1
    ## [64,]  .  .  . 2  .  . .  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [65,]  .  1  . .  2  . .  .  3 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [66,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [67,]  .  .  . 1  1  . .  .  3 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [68,]  .  .  . .  .  . .  .  1 . . .  1 . . . .  .  .  .  .  .  .  . . . .  .
    ## [69,]  .  .  . .  .  . .  .  1 . . .  . . . . .  .  .  .  .  .  .  . . . .  .
    ## [70,]  .  .  . 1  .  . .  .  2 . . .  . . . . .  .  .  .  .  .  .  1 1 . .  .
    ## [71,]  5  8  1 2  3  6 3  4  2 3 4 4  3 5 4 . 5  2  1  .  .  .  .  . . . .  .
    ## [72,]  3  2  3 5  .  4 3  3  . . 1 1  1 1 . . 4  .  .  . 41  .  .  . . . .  .
    ## [73,]  5  .  3 4  1  4 2  5  1 2 2 .  4 4 3 3 1  4  1  6  .  .  .  . . . .  .
    ## [74,]  5 12  2 4  2  3 3  2  2 7 . 4  1 2 4 1 4  .  .  1  .  .  .  . . . .  .
    ## [75,]  2  8  3 1  1  4 2 14  1 4 2 2  3 4 8 1 .  .  .  .  2  .  .  . . . .  1
    ## [76,] 42  7 11 6 14  3 9 32 10 . 3 6  8 1 1 2 . 20 25 26  1  .  .  . . . .  .
    ## [77,]  2  3  6 .  2  4 .  2 37 1 1 2  . 4 2 1 1  2  .  1  .  .  .  . . . .  .
    ## [78,]  1  2  5 4  . 20 .  .  . 3 2 2 13 1 . 2 .  2  3  .  1  .  .  . . . .  .
    ## [79,]  2  6  3 .  4  5 1  8  2 5 4 .  2 . . 3 .  1  1  .  1  .  .  . . . .  .
    ## [80,]  1  .  . 1  1  . .  .  3 2 . .  . . 1 . .  .  1  1  1  .  .  . . . .  .
    ##                                                   
    ##  [1,]  . 1 .  . . . .  . . . . . 1 . . 1  . 1 .  .
    ##  [2,]  . . 1  . . . .  1 . . 1 1 . . . .  . . 1  .
    ##  [3,]  . 1 .  . . . .  . . . . 1 . . . .  . . .  .
    ##  [4,]  . 1 .  . . . .  . . . . . . . . 3  . . .  .
    ##  [5,]  . 2 1  . . . .  . . . . 2 1 . . .  . . .  .
    ##  [6,]  . . 1  . . . .  . . . 1 1 . . . 1  . . .  .
    ##  [7,]  . . .  . . . .  . . . . . . . . 1  . 1 .  .
    ##  [8,]  . 1 .  . . . .  . . . . . 1 . . .  . . .  .
    ##  [9,]  . . .  . . . .  . . . . . . . . 1  . . .  .
    ## [10,]  . . .  . . . .  . . . 1 . . . . .  . . .  .
    ## [11,]  . . .  . . . .  . . . . . 1 . . .  . . .  .
    ## [12,]  . 1 .  . . . .  . . . . . . . . 2  . . .  .
    ## [13,]  . . .  . . . .  . . 1 1 . . . . .  . . .  .
    ## [14,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [15,]  . 1 .  . . . .  . . . 1 . . . . 1  . . .  .
    ## [16,]  . . .  . . . .  1 . . . . 2 . . .  . . .  .
    ## [17,]  . . .  . . . .  . . . . . 1 . . .  . . .  .
    ## [18,]  . . .  . . . .  . 1 . . . . . . .  . . .  .
    ## [19,]  . . .  . . . .  1 . . 1 . 1 . . .  . . .  .
    ## [20,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [21,]  . 1 .  . . . .  . . . . . . . . .  . . .  .
    ## [22,]  . . .  . . . .  . . . 2 . . . . .  . . .  .
    ## [23,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [24,]  . . .  . . . .  . . . . . . . . .  1 . .  .
    ## [25,]  . . .  . . . .  . . . . . . . . .  . . 1  .
    ## [26,]  . . .  . . . .  . . . . . . . . .  1 . .  .
    ## [27,]  . . .  . 1 . .  1 . . . . . . . 1  . . .  .
    ## [28,]  . . .  . . . .  1 . . 1 . 1 . . .  . . .  .
    ## [29,]  . . .  . 1 . 1  1 . . . . . . . .  . . .  .
    ## [30,]  . . .  . . . .  . . . . 1 . . . .  2 . .  .
    ## [31,]  . 1 .  . 4 2 1  1 1 3 2 3 3 . 4 7  1 . 3  .
    ## [32,]  3 1 4  1 . 1 .  1 . . 2 1 . 1 1 3  . . 5  .
    ## [33,]  5 2 2  1 . 1 1  1 . 4 1 1 2 3 . 2 58 1 1  .
    ## [34,]  7 2 1  1 3 1 1  1 . . 1 3 1 2 1 .  . 2 3  .
    ## [35,]  1 4 1  . . . .  1 2 1 1 . 1 . . 1  . 1 . 10
    ## [36,]  . 1 2  . 1 . 1  . 1 2 1 2 3 . . 3  1 3 .  .
    ## [37,]  3 . 3  1 7 1 2  1 1 3 3 2 1 . 4 3  . 2 1  .
    ## [38,]  1 2 2  . 4 . 2  . 1 . . . 2 . 2 5  2 1 2  1
    ## [39,]  . 1 2  1 . 1 1  1 . 1 3 1 1 1 1 .  1 2 .  .
    ## [40,]  2 3 6 51 3 1 . 25 3 . 1 1 2 2 2 1  . 1 2  1
    ## [41,]  6 1 2  . 1 1 1  1 . . . 2 . . . .  . . 3  .
    ## [42,]  . . 1  . . . .  . . . 1 1 2 1 1 2  1 1 .  .
    ## [43,]  . . .  . . . 1  . . . 2 . . . . 1  . . .  .
    ## [44,]  . . .  . . . .  . . . . . . . . 1  . 1 .  .
    ## [45,]  . . 1  1 . . .  . . . . . . . . .  . . .  .
    ## [46,]  . . .  . . . .  . . . 1 . . . . 2  . . .  .
    ## [47,] 10 1 3  . . . .  . . . . 1 . . . 2  1 1 3  .
    ## [48,]  . 1 2  . . . .  . 1 . . . 1 . . .  . 2 1  .
    ## [49,]  9 . .  . . . 1  . . . . 1 . 1 . 1  . . 3  .
    ## [50,]  . . 1  . . . .  . . . . . . . 1 1  . 1 .  .
    ## [51,]  . . .  . 1 . .  2 . . . . . . . .  . . .  .
    ## [52,]  . . .  . . . .  . . 1 . . . . . .  1 1 .  .
    ## [53,]  . . .  . . . .  1 . . . . . . . .  . . .  .
    ## [54,]  . . .  . 1 . .  . . . . . . . . .  . . .  .
    ## [55,]  . . .  . . 1 .  . . . . . 1 . . 2  2 . .  .
    ## [56,]  . . .  . . . .  3 . . . . . . . .  . . .  .
    ## [57,]  . . .  . . . .  2 . . . . . . . .  1 . .  .
    ## [58,]  . . .  . . . .  1 . . . . . . . .  . . .  .
    ## [59,]  . . .  . 1 . .  . . . . . . . . .  . . .  .
    ## [60,]  . . .  . 1 . .  . . . . . . . . .  . 3 .  .
    ## [61,]  . 1 1  . . . .  1 . . . . . . 1 .  2 . 1  .
    ## [62,]  . . .  . . . .  . . . . . . . . 1  . . 1  .
    ## [63,]  . . .  . . . .  1 . . . . . . . .  . . .  .
    ## [64,]  . 1 .  . 1 . .  1 . . . . . . . .  . . .  .
    ## [65,]  . . .  . 3 . .  1 . . 1 . 1 . . 2  4 . .  .
    ## [66,]  . . .  . 2 . .  2 . . 1 . 1 . . 1  2 . .  .
    ## [67,]  . 2 .  . . . .  4 . . . . . . 1 1  1 1 .  .
    ## [68,]  . . .  . . . .  1 . . . . . . . .  . . 1  .
    ## [69,]  . . .  . 1 . .  . . . . . 1 . . .  3 . .  .
    ## [70,]  . . .  . . . .  . . . 1 . . . . .  1 . .  .
    ## [71,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [72,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [73,]  . . .  . . . .  . . . . . . . . 1  . . .  .
    ## [74,]  . . .  . . . .  . . . . . . . 2 .  . . .  .
    ## [75,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [76,]  . . .  1 . . .  1 . . . . . . . 1  . 1 .  .
    ## [77,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [78,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [79,]  . . .  . . . .  . . . . . . . . .  . . .  .
    ## [80,]  . . .  . . . .  . . . 1 . . . . .  . . .  .

### Slicing

Just as with a `SOMADataFrame`, we can also retrieve subsets of the data
from a `SOMASparseNDArray` that can fit in memory.

Unlike `SOMADataFrame`s, `SOMASparseNDArray`s are always indexed using a
zero-based offset integer on each dimension, named `soma_dim_N`.
Therefore, if the array is `N`-dimensional, the `read()` method can
accept a list of length `N` that specifies how to slice the array.

`SOMASparseNDArray` dimensions are always named `soma_dim_N` where `N`
is the dimension number. As before you could use `schema()` or
[`dimnames()`](https://rdrr.io/r/base/dimnames.html) to retrieve the
dimension names.

``` r
counts$schema()
```

    ## Schema
    ## soma_dim_0: int64 not null
    ## soma_dim_1: int64 not null
    ## soma_data: double not null

For example, here’s how to fetch the first 5 rows of the matrix:

``` r
counts$read(coords = list(soma_dim_0 = 0:4))$tables()$concat()
```

    ## Table
    ## 258 rows x 3 columns
    ## $soma_dim_0 <int64 not null>
    ## $soma_dim_1 <int64 not null>
    ## $soma_data <double not null>
