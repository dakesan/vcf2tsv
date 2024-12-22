# vcf2tsv.py

A simple standalone bioinformatics tool for converting VCF (Variant Calling Format) files into TSV files.

`vcf2tsv.py` is a modified version of the `vcf2tsv.py` script from **vcf-kit**. Since installing vcf-kit can be a bit challenging for now and its `vcf2tsv.py` script may not work seamlessly in Snakemake workflows, this standalone version was refactored from the original code to simplify usage.

## Features

- Supports only the `--wide` mode of the original `vcf2tsv.py`.
- Retains the `--print-header` option for extracting column headers.

## Installation

Clone this repository and hit command below.

```bash
pip instlal .
```
## Usage

To use this script, run the following command:

```bash
vcf2tsv --print-header variants.vcf.gz --output output.tsv
```

## Options

--print-header: Prints the header of the VCF file to the output.
--output: Specifies the output file path.

## License

This script inherits the MIT license from the original vcf-kit project.
