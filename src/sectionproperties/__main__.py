"""Command-line interface."""

import click


@click.command()
@click.version_option()
def main() -> None:
    """sectionproperties."""


if __name__ == "__main__":
    main(prog_name="sectionproperties")  # pragma: no cover
