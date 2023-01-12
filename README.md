# SGE Primer Scoring

## Description
Tool to score SGE primers based on output from [Exonerate iPCRess](https://www.ebi.ac.uk/about/vertebrate-genomics/software/ipcress-manual).

For each primer pair, occurrences of primer A and B in PCR products with each number of mismatches are counted. The number of PCR products with each total number of mismatches is also counted. The total counts are then used to score each primer pair, with a penalty added for each hit. The penalty increases exponentially as the number of mismatches decreases from 8. Hits with 0 or 1 mismatches are given the highest penalty, although the first hit with 0 mismatches is exempt as it should represent the on-target hit. An error is raised if this is not found. Primer pairs are ranked by score, individually for each targeton if a CSV mapping targetons and primer pairs is provided. The lowest score shows highest specificity.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

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
                        CSV of primer pairs and corresponding targetons- adds
                        targeton column to output
  --version             show program's version number and exit
```

Example command: ```./score_primers.py examples/example_input.txt 4 examples/example_output.tsv```

The mismatch number provided dictates the number of mismatch columns in the output TSV, so please use the same value as used with iPCRess or results could be misleading. The mismatch number used with iPCRess affects the score, so bear this in mind if comparing results.

Parent directories in the output path are created if required.

Example command with targeton CSV (to be used if looking at multiple targetons):  
```./score_primers.py examples/example_input.txt 4 example_targeton_output.tsv --targeton_csv examples/example_targetons.csv```  

CSV format:
```
SMARCA4_exon24_1,SMARCA4_exon24
SMARCA4_exon24_2,SMARCA4_exon24
BRCA1_exon1_1,BRCA1_exon1
```

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

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
Linting:
`pycodestyle`

Testing:
`python3 -m unittest`

Tests are written using the Arrange Act Assert pattern.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
