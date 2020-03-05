import click

from checseqtools import DyadPosition, PrepareGenomeCoverage


@click.group()
def chectools():
    pass


chectools.add_command(DyadPosition.dyadposition)
chectools.add_command(PrepareGenomeCoverage.prepgenomecov)

if __name__ == '__main__':
   chectools()
