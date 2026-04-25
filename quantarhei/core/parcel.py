from __future__ import annotations

from typing import IO, Any

import dill as pickle

from .managers import Manager


class Parcel:
    def set_content(self, obj: Any) -> None:
        """Set the content of the parcel"""
        self.content = obj
        self.class_name = f"{obj.__class__.__module__}.{obj.__class__.__name__}"
        self.qrversion = Manager().version
        self.comment = ""

    def set_comment(self, comm: str | None) -> None:
        """Sets a string value to a comment saved togethet with the object"""
        if comm is not None:
            self.comment = comm

    def save(self, filename: str | IO[bytes]) -> None:
        """Saves the parcel to a file

        Parameters
        ----------
        filename : str or File
            Name of the file or a file object to which the content of
            the object will be saved


        """
        if isinstance(filename, str):
            with open(filename, "wb") as f:
                pickle.dump(self, f)
        else:
            pickle.dump(self, filename)


def save_parcel(
    obj: Any, filename: str | IO[bytes], comment: str | None = None
) -> None:
    """Saves a given object as a parcel

    Parameters
    ----------
    filename : str or File
        Name of the file or a file object to which the content of
        the object will be saved

    comment : str
        A comment which will be saved together with the content of
        the object

    """
    p = Parcel()
    p.set_content(obj)
    p.set_comment(comment)

    p.save(filename)


def load_parcel(filename: str | IO[bytes]) -> Any:
    """Loads the object saved as parcel

    Parameters
    ----------
    filename : str or File
        Filename of the file or file descriptor of the file from which
        and object should be loaded.

    """
    if isinstance(filename, str):
        with open(filename, "rb") as f:
            obj = pickle.load(f)
    else:
        obj = pickle.load(filename)

    if isinstance(obj, Parcel):
        return obj.content
    raise Exception("Only Quantarhei Parcels can be loaded")


def check_parcel(filename: str | IO[bytes]) -> dict[str, Any]:
    """Checks the content of a Quantarhei parcel

    Parameters
    ----------
    filename : str or File
        Filename of the file or file descriptor of the file from which
        and object should be loaded.

    """
    if isinstance(filename, str):
        with open(filename, "rb") as f:
            obj = pickle.load(f)
    else:
        obj = pickle.load(filename)

    if isinstance(obj, Parcel):
        return dict(
            class_name=obj.class_name, qrversion=obj.qrversion, comment=obj.comment
        )
    raise Exception("The file does not represent a Quantarhei parcel")
