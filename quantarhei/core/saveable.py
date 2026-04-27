"""Module saveable

Defines the class Saveable. Saveable objects can save themselves to and
load themselves from a file.


"""
from __future__ import annotations

import copy
import os

#import time
#import hashlib
import uuid
from tempfile import TemporaryDirectory
from typing import IO, Any

from .parcel import Parcel, load_parcel


class Saveable:
    """Defines a class of objects that can save and load themselves


    """

    hashes: dict = {}

    def save(self, filename: str | IO[bytes], comment: str | None = None,
             test: bool = False) -> None:
        """Saves the object with all its content into a file


        Parameters
        ----------
        filename : str or File
            Name of the file or a file object to which the content of
            the object will be saved

        comment : str
            A comment which will be saved together with the content of
            the object

        test : bool
            If test is True, and file descriptor is submitted, we save into
            the file and move to its start for subsequent reading


        """
        p = Parcel()
        p.set_content(self)
        p.set_comment(comment)

        p.save(filename)

        if test:
            if not isinstance(filename, str):
                filename.seek(0)


    def load(self, filename: str | IO[bytes], test: bool = False) -> Any:
        """Loads an object from a file and returns it

        Parameters
        ----------
        filename : str or File
            Filename of the file or file descriptor of the file from which
            and object should be loaded.

        """
        if test:
            if not isinstance(filename, str):
                filename.seek(0)

        return load_parcel(filename)


    def _get_fname(self) -> str:

        #hashstr = str(hash(time.time()))
        #str40 = str(hashlib.sha1(hashstr.encode()).hexdigest())
        str40 = str(uuid.uuid4())
        return str40


    def savedir(self, dirname: str, tag: Any = None,
                comment: str | None = None, test: bool = False) -> None:
        """Saves an object into directory containing a file with unique name


        """
        hfile = os.path.join(dirname,"_hashes_.qrp")
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
        fname = os.path.join(dirname, str40+".qrp")

        # we try once more if the file already exists
        if os.path.isfile(fname):
            str40 = self._get_fname()
            fname = os.path.join(dirname, str40+".qrp")
            if os.path.isfile(fname):
                raise Exception("File already exists")

        self.save(fname)

        self.hashes[tag] = str40
        #print(tag, str40)
        p = Parcel()
        p.set_content(self.hashes)
        p.set_comment("Hashes")
        p.save(hfile)


    def loaddir(self, dirname: str) -> dict[Any, Any]:
        """Returns a directory of objects saved into a directory

        """
        out: dict[Any, Any] = {}
        hfile = os.path.join(dirname,"_hashes_.qrp")
        hashes = load_parcel(hfile)

        for tag in hashes:
            fname = hashes[tag]+".qrp"
            fdname = os.path.join(dirname,fname)
            obj = load_parcel(fdname)
            out[tag] = obj

        return out


    def scopy(self) -> Any:
        """Creates a copy of the object by saving and loading it

        """
        with TemporaryDirectory() as td:
            fname = os.path.join(td, "ssave.qrp")
            self.save(fname)

            no = load_parcel(fname)

        return no


    def deepcopy(self) -> Any:
        """Returns a deep copy of the self


        """
        return copy.deepcopy(self)


    def copy(self) -> Any:
        """Returns a shallow copy of the self


        """
        return copy.copy(self)
