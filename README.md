# SGE Primer Scoring

## Description
Tool to score SGE primers based on output from [Exonerate iPCRess](https://www.ebi.ac.uk/about/vertebrate-genomics/software/ipcress-manual).

For each primer pair, occurrences of primer A and B in PCR products with each number of mismatches are counted. The number of PCR products with each total number of mismatches is also counted.

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
An example input file and TSV of the formatted mismatches can be found in the examples folder.

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
