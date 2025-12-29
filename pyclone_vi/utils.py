import click


def print_command_header(command_title):
    click.echo()
    click.echo("#" * 100)
    click.secho(f"PyClone-VI: {command_title}")
    click.echo("#" * 100)
    click.echo()
