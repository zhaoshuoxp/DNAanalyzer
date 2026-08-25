#!/usr/bin/env python3
"""Generate the static restriction-enzyme table used by the Rust build.

This is a maintainer tool only; the application has no Python dependency.
"""

from Bio import __version__ as biopython_version
from Bio.Restriction import AllEnzymes


def value(number):
    return "-" if number is None else str(number)


print(f"# Generated from Biopython {biopython_version} Bio.Restriction")
print("# name\tsite\tfst5\tfst3\tscd5\tscd3")
for enzyme in sorted(AllEnzymes, key=lambda item: item.__name__.upper()):
    print(
        "\t".join(
            [
                enzyme.__name__,
                enzyme.site,
                value(enzyme.fst5),
                value(enzyme.fst3),
                value(enzyme.scd5),
                value(enzyme.scd3),
            ]
        )
    )
