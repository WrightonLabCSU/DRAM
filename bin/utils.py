import logging
import shlex
import subprocess
from pathlib import Path

from bin.definitions import IS_WINDOWS

logger = logging.getLogger(__name__)


def split_shell_command(cmd: str):
    """
    split shell command for passing to python subproccess.
    This should correctly split commands like "echo 'Hello, World!'"
    to ['echo', 'Hello, World!'] (2 items) and not ['echo', "'Hello,", "World!'"] (3 items)

    It also works for posix and windows systems appropriately
    """
    return shlex.split(cmd, posix=not IS_WINDOWS)


def run_process(command, shell:bool=False, stop_on_error:bool=True) -> str:
    """
    Standardization of parameters for using subprocess.run, provides verbose mode and option to run via shell
    """
    # TODO just remove check
    
    logger.debug(f'Running command: {command}')
    results = subprocess.run(split_shell_command(command), shell=shell,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if results.returncode != 0:
        logger.critical(f'The subcommand {command} experienced an error: {results.stderr}')
        logging.debug(results.stdout)
        if stop_on_error:
           raise subprocess.SubprocessError(f"The subcommand {' '.join(command)} experienced an error, see the log for more info.")


def get_fasta_sample_name(fasta: Path) -> str:
    return fasta.stem.replace(".", "_")
