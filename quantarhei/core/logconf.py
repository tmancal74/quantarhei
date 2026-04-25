from __future__ import annotations

from queue import LifoQueue
from typing import IO


class LogConf:
    """Holds information of the configuration for Quantarhei logging system

    When initialized, it holds hard default configuration of the logging
    system.

    """

    def __init__(self) -> None:

        self.verbosity: int = 5
        self.fverbosity: int = 5
        self.verbose: bool = True
        self.log_on_screen: bool = True
        self.log_indent: int = 0
        self.log_to_file: bool = False
        self.log_file_opened: bool = False
        self.log_file_name: str = "qrhei.log"
        self.log_file_appendix: str = ""
        self.log_file: IO[str] | None = None
        self.initialized: bool = False
        self.time_stamp: LifoQueue[float] = LifoQueue()
        self.is_serial: bool = True

    def __del__(self) -> None:
        """Destructor to be used on garbage collection

        It closes openned log_file if it was not closed by something before

        """
        if self.log_file_opened and self.log_file is not None:
            self.log_file.close()
