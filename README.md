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

If you want to use this code in snakemake, you can build conda env be this requirements.yml.

```
name: vcf2tsv
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - bwa>0.7.17
  - samtools>=1.10
  - bcftools>=1.10
  - blast>=2.2.31
  - muscle>=3.8.31
  - primer3>=2.5.0
  - openblas
  - pip
  - setuptools>=64.0.0
  - wheel
  - pip:
    - git+https://github.com/dakesan/vcf2tsv
```

## Options

--print-header: Prints the header of the VCF file to the output.
--output: Specifies the output file path.

## License

This script inherits the MIT license from the original vcf-kit project.
