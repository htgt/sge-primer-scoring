# SGE Primer Scoring

## Description
Tool to score primers based on specificity using output from [Exonerate iPCRess](https://www.ebi.ac.uk/about/vertebrate-genomics/software/ipcress-manual).

For each primer pair, occurrences of primer A and B in PCR products with each number of mismatches are counted. The number of PCR products with each total number of mismatches is also counted. The total counts are then used to score each primer pair, with a penalty added for each hit. The penalty increases exponentially as the number of mismatches decreases from 8. Hits with 0 or 1 mismatches are given the highest penalty, although the first hit with 0 mismatches is exempt as it should represent the on-target hit. An error is raised if this is not found.

Primer pairs are ranked by score, individually for each targeton if a CSV mapping targetons and primer pairs is provided. The lowest score shows highest specificity.

## Installation
```
git clone https://gitlab.internal.sanger.ac.uk/sci/sge-primer-scoring.git
cd sge-primer-scoring
pip3 install --upgrade pip
pip3 install -r requirements.txt
```

## Usage
```
usage: score_primers.py [-h] [--targeton_csv TARGETON_CSV] [--version]
                        ipcress_file mismatch output_tsv

Tool to score primer pairs using output from Exonerate iPCRess

positional arguments:
  ipcress_file          File containing output from Exonerate iPCRess
  mismatch              Mismatch number used for Exonerate iPCRess
  output_tsv            Path for output TSV file

optional arguments:
  -h, --help            show this help message and exit
  --targeton_csv TARGETON_CSV
                        CSV of primer pairs and corresponding targetons - adds
                        targeton column to output
  --version             show program's version number and exit
```

Example command: ```./score_primers.py examples/example_ipcress_file.txt 4 example_output.tsv```

Example command with targeton CSV (to be used if looking at multiple targetons):  
```./score_primers.py examples/example_ipcress_file.txt 4 example_targeton_output.tsv --targeton_csv examples/example_targetons.csv```

CSV format:
```
primer_pair_1,targeton_1
primer_pair_2,targeton_1
primer_pair_3,targeton_2
```

The mismatch number provided dictates the number of mismatch columns in the output TSV, so please use the same value as used with iPCRess or results could be misleading. The mismatch number used with iPCRess affects the score, so bear this in mind if comparing results.

Parent directories in the output path are created if required. Example output files can be found in the examples folder along with the input files.

**Raises:**
- ArgumentTypeError if an input file does not exist
- ArgumentTypeError if an input file is empty
- ArgumentTypeError if mismatch number is negative
- ArgumentTypeError if output file is a directory
- ArgumentTypeError if output file already exists
- ScoringError if an input file format is invalid
- ScoringError if mismatch number is not negative but still too low for ipcress file provided
- ScoringError if there is no data in the ipcress file
- ScoringError if no on-target hit is found in the ipcress file
- ScoringError if a primer pair in the targeton csv appears again with a different targeton

## Contributing
Linting:
`pycodestyle`

Testing:
`python3 -m unittest`

Tests are written using the Arrange Act Assert pattern.

## License
```
Copyright (c) 2022, 2023 Genome Research Ltd.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
```
