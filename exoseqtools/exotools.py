import click

from exoseqtools import MoveAnnnotations, RemoveSecondMateBam


@click.group()
def exotools():
    pass


exotools.add_command(MoveAnnnotations.moveannotations)
exotools.add_command(RemoveSecondMateBam.removesecondmate)

if __name__ == '__main__':
   exotools()
