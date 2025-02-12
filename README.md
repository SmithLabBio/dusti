# Quartet Inference

A Python package for quartet inference from gene alignments.

## Installation

Install the package with:

```bash
pip install git+https://github.com/SmithLabBio/dusti.git
```

Get the wQFM program from: https://github.com/Mahim1997/wQFM-2020

## Running

For a simple test run:
```bash
dusti --input example/alignments --map example/s_map.txt -o example_results/results --qfm ~/Documents/programs/wQFM-2020/wQFM-v1.4.jar --svd --parsimony
```

For more information, view the manual at ./manual/manual.pdf
