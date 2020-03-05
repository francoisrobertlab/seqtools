import click

from mnaseseqtools import DyadCoverage, FitDoubleGaussian, FitGaussian, FirstDyadPosition, PrepareGenomeCoverage


@click.group()
def mnasetools():
    pass


mnasetools.add_command(DyadCoverage.dyadcov)
mnasetools.add_command(FitGaussian.fitgaussian)
mnasetools.add_command(FitDoubleGaussian.fitdoublegaussian)
mnasetools.add_command(FirstDyadPosition.firstdyadposition)
mnasetools.add_command(PrepareGenomeCoverage.prepgenomecov)

if __name__ == '__main__':
   mnasetools()
