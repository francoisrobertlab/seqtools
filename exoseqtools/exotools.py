import click

from exoseqtools import RemoveSecondMateBam


@click.group()
def exotools():
    pass


exotools.add_command(RemoveSecondMateBam.removesecondmate)

if __name__ == '__main__':
   exotools()
