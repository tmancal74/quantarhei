"""Module saveable

Defines the class Saveable. Saveable objects can save themselves to and
load themselves from a file.


"""

from __future__ import annotations

import copy
import os

# import time
# import hashlib
import uuid
from tempfile import TemporaryDirectory
from typing import IO, Any

from ..exceptions import QuantarheiError
from .parcel import Parcel, load_parcel


class Saveable:
    """Base class for objects that can save and load themselves to/from files.

    Subclasses gain the ability to persist their state to a ``.qrp`` file
    (a pickled :class:`Parcel`) and restore it later, supporting both named
    files and file-like objects.
    """

    hashes: dict = {}

    def save(
        self, filename: str | IO[bytes], comment: str | None = None, test: bool = False
    ) -> None:
        """Save the object and all its content to a file.

        Parameters
        ----------
        filename : str or file-like
            Path to the output file or an open binary file object.
        comment : str, optional
            A comment stored alongside the object content.
        test : bool, optional
            If ``True`` and a file object is given, seek back to the
            start after writing so the file is ready for immediate reading.
            Default is ``False``.
        """
        p = Parcel()
        p.set_content(self)
        p.set_comment(comment)

        p.save(filename)

        if test:
            if not isinstance(filename, str):
                filename.seek(0)

    def load(self, filename: str | IO[bytes], test: bool = False) -> Any:
        """Load an object from a file and return it.

        Parameters
        ----------
        filename : str or file-like
            Path to the file or an open binary file object from which
            the object should be loaded.
        test : bool, optional
            If ``True`` and a file object is given, seek to the start
            before reading. Default is ``False``.

        Returns
        -------
        object
            The object stored in the file.
        """
        if test:
            if not isinstance(filename, str):
                filename.seek(0)

        return load_parcel(filename)

    def _get_fname(self) -> str:

        # hashstr = str(hash(time.time()))
        # str40 = str(hashlib.sha1(hashstr.encode()).hexdigest())
        str40 = str(uuid.uuid4())
        return str40

    def savedir(
        self,
        dirname: str,
        tag: Any = None,
        comment: str | None = None,
        test: bool = False,
    ) -> None:
        """Saves an object into directory containing a file with unique name"""
        hfile = os.path.join(dirname, "_hashes_.qrp")
        try:
            os.makedirs(dirname)
            self.hashes = {}
        except FileExistsError:
            self.hashes = load_parcel(hfile)

        if tag is None:
            try:
                last = list(self.hashes.keys())[-1]
            except IndexError:
                last = 0
            tag = last + 1

        # get a unique name for the file
        str40 = self._get_fname()
        fname = os.path.join(dirname, str40 + ".qrp")

        # we try once more if the file already exists
        if os.path.isfile(fname):
            str40 = self._get_fname()
            fname = os.path.join(dirname, str40 + ".qrp")
            if os.path.isfile(fname):
                raise QuantarheiError("File already exists")

        self.save(fname)

        self.hashes[tag] = str40
        # print(tag, str40)
        p = Parcel()
        p.set_content(self.hashes)
        p.set_comment("Hashes")
        p.save(hfile)

    def loaddir(self, dirname: str) -> dict[Any, Any]:
        """Returns a directory of objects saved into a directory"""
        out: dict[Any, Any] = {}
        hfile = os.path.join(dirname, "_hashes_.qrp")
        hashes = load_parcel(hfile)

        for tag in hashes:
            fname = hashes[tag] + ".qrp"
            fdname = os.path.join(dirname, fname)
            obj = load_parcel(fdname)
            out[tag] = obj

        return out

    def scopy(self) -> Any:
        """Creates a copy of the object by saving and loading it"""
        with TemporaryDirectory() as td:
            fname = os.path.join(td, "ssave.qrp")
            self.save(fname)

            no = load_parcel(fname)

        return no

    def deepcopy(self) -> Any:
        """Returns a deep copy of the self"""
        return copy.deepcopy(self)

    def copy(self) -> Any:
        """Returns a shallow copy of the self"""
        return copy.copy(self)
