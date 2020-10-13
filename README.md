# Snap-Analyze

Analyzes the bonding structure of an [ESPResSo](https://github.com/espressomd/espresso) MPI-IO snapshot.

Snap-analyze supports two modes:

1. `--print-all-to-files` Outputs all agglomerates to individual files. These can be read in e.g. by `numpy.loadtxt('...')`.
1. `--df` Calculates the radii of gyration and the fractal dimension of each individual agglomerate and prints them to stdout.

## Build

```sh
git clone https://github.com/hirschsn/snap-analyze.git
cd snap-analyze
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## Usage

```
Usage: ./sana-bond [OPTIONS...] SNAP-PREFIX
Options:
 --sigma SIGMA | -s SIGMA   Set sigma for Df calculation
 --box_l BOX | -b BOX       Set box size for Df calculation
Modi:
 --print-all-to-files | -p  Print all agglomerates to files named
                              <NPART>_POS_<IDX>
 --df | -d                  Calculate fractal dimensions
```

For large numbers of particles and bonds, have sufficiently RAM available.

Note: The code can be used as is for ESPResSo < v4.0.0. For newer versions, a binary file offset has to be changed. For recent versions of ESPResSo, see the comment in `snapshot.hpp:53ff`.

## License (ISC)

Copyright 2018-2020 Steffen Hirschmann

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

