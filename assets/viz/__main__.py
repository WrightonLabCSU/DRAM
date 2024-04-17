import click
from pathlib import Path

from .definitions import DEFAULT_GROUPBY_COLUMN
from .make_product import main as make_product


@click.command()
@click.option(
    "--annotations",
    "-a",
    type=Path,
    exists=True,
    help="Path to the annotations tsv file"
)
@click.option(
    "--groupby_column",
    "-g",
    type=str,
    default=DEFAULT_GROUPBY_COLUMN,
    help="Column to group by"
)
@click.option(
    "--output_dir",
    "-o",
    type=Path,
    help="Path to the output directory",
    default=Path.cwd()
)
@click.option(
    "--dashboard",
    "-d",
    is_flag=True,
)
def main(annotations, groupby_column, output_dir, dashboard):
    """
    A visualization tool for DRAM2
    """
    make_product(
        annotations_tsv_path=annotations,
        groupby_column=groupby_column,
        output_dir=output_dir,
        dashboard=dashboard
    )


if __name__ == "__main__":
    main()
