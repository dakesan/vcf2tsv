import gzip
import os
import re
import sys
import tempfile
from pathlib import Path
from signal import SIG_DFL, SIGPIPE, signal
from subprocess import PIPE, Popen

import typer
from clint.textui import colored, puts_err

signal(SIGPIPE, SIG_DFL)

app = typer.Typer()


# Utility Functions
def command_out(message: str):
    typer.secho(message, fg=typer.colors.GREEN)


def message(message: str, color: str = "blue"):
    color_map = {
        "blue": typer.colors.BLUE,
        "red": typer.colors.RED,
    }
    typer.secho(message, fg=color_map.get(color, typer.colors.WHITE))


def parse_region(region):
    return re.split("[:-]+", region)


def which(program: str) -> Path | None:
    """
    Mimics the Unix 'which' command to locate a program in the system PATH.
    Uses pathlib instead of os.
    """
    program_path = Path(program)

    # If program includes a path and is executable, return it
    if program_path.is_file() and program_path.stat().st_mode & 0o111:
        return program_path

    # Search through the system PATH
    for path in map(Path, os.environ["PATH"].split(os.pathsep)):
        exe_file = path / program
        if exe_file.is_file() and exe_file.stat().st_mode & 0o111:
            return exe_file

    return None


def check_program_exists(program):
    if which(program) is None:
        exit(puts_err(colored.red("\n\t" + program + " not installed or on PATH.\n")))


def get_genome_directory_file() -> Path:
    return Path(__file__).resolve().parent / ".genome_directory"


def get_genome_directory() -> Path:
    genome_directory_file = get_genome_directory_file()
    if genome_directory_file.exists():
        genome_directory = Path(genome_directory_file.read_text().strip()).expanduser()
    else:
        genome_directory = Path("~/.genome").expanduser()
        genome_directory_file.write_text(str(genome_directory))
    genome_directory.mkdir(parents=True, exist_ok=True)
    return genome_directory


def get_genome_list() -> list[str]:
    return [p.name for p in get_genome_directory().iterdir() if p.is_dir()]


def output_genome_list():
    typer.secho("\n".join(get_genome_list()), fg=typer.colors.BLUE)


def resolve_reference_genome(loc: str) -> Path:
    genome_dir = get_genome_directory()
    if not loc:
        typer.secho("You must specify a genome:", fg=typer.colors.RED)
        output_genome_list()
        raise typer.Exit(code=1)

    loc_path = Path(loc)
    if loc_path.exists():
        return loc_path
    elif loc in get_genome_list():
        reference_location = genome_dir / loc / f"{loc}.fa.gz"
        typer.secho(
            f"Using reference located at {reference_location}", fg=typer.colors.GREEN
        )
        return reference_location
    else:
        typer.secho(f"Genome '{loc}' does not exist", fg=typer.colors.RED)
        raise typer.Exit(code=1)


class Vcf:
    def __init__(self, filename: Path, reference=None):
        filename = Path(filename)

        if not filename.is_file() and str(filename) != "-":
            typer.secho(f"Error: {filename} does not exist", fg=typer.colors.RED)
            raise typer.Exit(code=1)

        self.filename = filename
        self.tempfile = None

        # Handle .gz files
        if self.filename.suffix == ".gz":
            self.filename = self._decompress_gzip(self.filename)

        self.raw_header = self._read_header()
        self.samples = self._read_samples()

    def _decompress_gzip(self, gz_file: Path) -> Path:
        """
        Decompress a .gz file and return the path to the decompressed file.
        """
        temp_dir = tempfile.TemporaryDirectory()
        decompressed_file = Path(temp_dir.name) / gz_file.stem

        with gzip.open(gz_file, "rt") as gz, open(decompressed_file, "w") as out_file:
            out_file.write(gz.read())

        # Store the temporary directory to ensure it's cleaned up later
        self.tempfile = temp_dir
        return decompressed_file

    def _read_header(self):
        with open(self.filename) as f:
            return "".join(line for line in f if line.startswith("##"))

    def _read_samples(self):
        with open(self.filename) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    return line.strip().split("\t")[9:]
        return []

    def __del__(self):
        """
        Clean up temporary files when the object is destroyed.
        """
        if self.tempfile:
            self.tempfile.cleanup()


# Info regex
r_info = re.compile(
    r"""\#\#INFO=<
  ID=(?P<id>[^,]+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>Integer|Float|Flag|Character|String),
  Description="(?P<desc>[^"]*)".*
  >""",
    re.VERBOSE,
)

# Format regex
r_format = re.compile(
    r"""\#\#FORMAT=<
  ID=(?P<id>.+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>.+),
  Description="(?P<desc>.*)".*
  >""",
    re.VERBOSE,
)

# ANN header for snpEff
ANN_header = [
    "allele",
    "effect",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "exon_intron_rank",
    "nt_change",
    "aa_change",
    "cDNA_position/cDNA_len",
    "protein_position",
    "distance_to_feature",
    "error",
]


@app.command()
def vcf2tsv(
    vcf: Path = typer.Argument(
        ..., exists=True, readable=True, help="Path to the VCF file to process."
    ),
    output: Path = typer.Option(
        None, "--output", "-o", help="Path to the output TSV file."
    ),
    print_header: bool = typer.Option(
        False, "--print-header", help="Print header line"
    ),
    ann: bool = typer.Option(
        False, "--ANN", help="Include ANN fields from snpEff annotations."
    ),
):
    """
    Convert a VCF file to a TSV format and save to a file or print to stdout.
    """
    v = Vcf(str(vcf))

    # Extract INFO and FORMAT fields from the VCF header
    info = [m.groupdict()["id"] for m in r_info.finditer(v.raw_header)]
    format = [m.groupdict()["id"] for m in r_format.finditer(v.raw_header)]

    # Construct Query String for bcftools
    query_start = repr(
        "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t"
        + "\t".join(["%INFO/" + x for x in info])
        + "[\t%SAMPLE\t"
        + "\t".join(["%" + x for x in format])
        + "]\n"
    ).strip("'")

    comm = [
        "bcftools",
        "query",
        "--print-header" if print_header else "",
        "-f",
        query_start,
        v.filename,
    ]

    # Filter empty arguments
    comm = [arg for arg in comm if arg]

    process = Popen(comm, stdout=PIPE, stderr=PIPE)
    output_lines = []
    if ann and "ANN" in info:
        ANN_loc = info.index("ANN") + 7

    for n, line in enumerate(process.stdout):
        line = line.decode("utf-8")
        if n == 0 and print_header:
            # Add snpEff headers if requested
            if "ANN" in info and ann:
                line = line.split("\t")
                line = line[: ANN_loc - 1] + ANN_header + line[ANN_loc + 1 :]
                line = "\t".join(line)
            output_lines.append(
                re.sub(r"\[[0-9]+\]", "", line).strip("#\n ").replace(":", "_")
            )
        else:
            # Parse annotations if ANN field exists
            if "ANN" in info and ann:
                line = line.split("\t")
                for var_effect in line[ANN_loc].split(","):
                    out_line = (
                        line[: ANN_loc - 1]
                        + var_effect.split("|")
                        + line[ANN_loc + 2 :]
                    )
                    output_lines.append("\t".join(out_line).strip("\n"))
            else:
                output_lines.append(line.strip("\n"))

    # Write to file or print to stdout
    if output:
        with open(output, "w") as f:
            f.write("\n".join(output_lines) + "\n")
        print(f"Output written to {output}")
    else:
        for line in output_lines:
            print(line)


if __name__ == "__main__":
    app(prog_name="vcf2tsv")
