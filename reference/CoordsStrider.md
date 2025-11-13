# The Coordinate Strider

The `CoordsStrider` allows creating coordinate slices in an iterated
manner. Alternatively, it can chunk an existing vector of coordinates

## Usage

``` r
# S3 method for class 'CoordsStrider'
as.list(x, ...)

# S3 method for class 'CoordsStrider'
length(x)

# S3 method for class 'CoordsStrider'
nextElem(obj, ...)

# S3 method for class 'CoordsStrider'
hasNext(obj, ...)
```

## Note

The `CoordsStrider` operates using [64-bit
integer](https://rdrr.io/pkg/bit64/man/bit64-package.html) objects; as
such, accessing fields, such as `strider$start` or `strider$stride` will
return an `integer64` object, which functions differently than a regular
`integer`. Use with caution and convert back to integers or numerics as
necessary

## Active bindings

- `coords`:

  If set, the coordinates to iterate over

- `start`:

  If set, the starting point of the iterated coordinates; otherwise the
  minimum value of `self$coords`

- `end`:

  If set, the end point of the iterated coordinates; otherwise the
  maximum value of `self$coords`

- `stride`:

  The stride, or how many coordinates to generate per iteration; note:
  this field is settable, which will reset the iterator

## Methods

### Public methods

- [`CoordsStrider$new()`](#method-CoordsStrider-new)

- [`CoordsStrider$print()`](#method-CoordsStrider-print)

- [`CoordsStrider$length()`](#method-CoordsStrider-length)

- [`CoordsStrider$has_next()`](#method-CoordsStrider-has_next)

- [`CoordsStrider$next_element()`](#method-CoordsStrider-next_element)

------------------------------------------------------------------------

### Method `new()`

Create a coordinate strider

#### Usage

    CoordsStrider$new(coords, ..., stride = NULL, start = NULL, end = NULL)

#### Arguments

- `coords`:

  An integer vector of coordinates

- `...`:

  Ignored

- `stride`:

  The stride of how many coordinates to yield per iteration; by default,
  will try to yield all coordinates per iteration

- `start`:

  If `coords` is missing, the starting coordinate to generate

- `end`:

  If `coords` is missing, the ending coordinate to generate

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the coordinate strider to the screen

#### Usage

    CoordsStrider$print()

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

Get the length (span) of the coordinates within the strider

#### Usage

    CoordsStrider$length()

#### Returns

The length (span) of the coordinate strider

------------------------------------------------------------------------

### Method `has_next()`

Determine if there are more coordinates to yield

#### Usage

    CoordsStrider$has_next()

#### Returns

`TRUE` if there are more coordinates to yield or `FALSE` if otherwise

------------------------------------------------------------------------

### Method `next_element()`

Generate the next set of coordinates to yield. If there are no more
coordinates to yield, raises a `stopIteration` error

#### Usage

    CoordsStrider$next_element()

#### Returns

An integer vector of the next set of coordinates

## Examples

``` r
(strider <- CoordsStrider$new(start = 1L, end = 200L, stride = 60L))
#> <CoordsStrider>
#>   start: 1 
#>   end: 200 
#>   stride: 60 
while (strider$has_next()) {
  str(strider$next_element())
}
#> integer64 [1:60] 1 2 3 4 5 6 7 8 ... 
#> integer64 [1:60] 61 62 63 64 65 66 67 68 ... 
#> integer64 [1:60] 121 122 123 124 125 126 127 128 ... 
#> integer64 [1:20] 181 182 183 184 185 186 187 188 ... 

(strider <- CoordsStrider$new(start = 1L, end = 200L, stride = 60L))
#> <CoordsStrider>
#>   start: 1 
#>   end: 200 
#>   stride: 60 
as.list(strider)
#> [[1]]
#> integer64
#>  [1] 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
#> [51] 51 52 53 54 55 56 57 58 59 60
#> 
#> [[2]]
#> integer64
#>  [1] 61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79 
#> [20] 80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98 
#> [39] 99  100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117
#> [58] 118 119 120
#> 
#> [[3]]
#> integer64
#>  [1] 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139
#> [20] 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158
#> [39] 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177
#> [58] 178 179 180
#> 
#> [[4]]
#> integer64
#>  [1] 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199
#> [20] 200
#> 

length(strider)
#> [1] 199

(strider <- CoordsStrider$new(start = 1L, end = 200L, stride = 60L))
#> <CoordsStrider>
#>   start: 1 
#>   end: 200 
#>   stride: 60 
while (itertools::hasNext(strider)) {
  str(iterators::nextElem(strider))
}
#> integer64 [1:60] 1 2 3 4 5 6 7 8 ... 
#> integer64 [1:60] 61 62 63 64 65 66 67 68 ... 
#> integer64 [1:60] 121 122 123 124 125 126 127 128 ... 
#> integer64 [1:20] 181 182 183 184 185 186 187 188 ... 
```
