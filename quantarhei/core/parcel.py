from __future__ import annotations

import os
import tempfile
import warnings
from typing import IO, Any

import dill as pickle

from ..exceptions import QuantarheiError
from .managers import Manager


class DeserializationWarning(UserWarning):
    """Warning emitted when loading a .qrp file without trusted=True."""


_DESERIALIZATION_WARNING = (
    "Loading a .qrp file deserializes Python objects using dill, which can "
    "execute arbitrary code. Only load files from trusted sources. "
    "Pass trusted=True to suppress this warning. "
    "See https://github.com/tmancal74/quantarhei/issues/208 for details."
)


class Parcel:
    """Container that serialises a single Python object to a ``.qrp`` file.

    Content is set via :meth:`set_content` and persisted via :meth:`save`.
    Use the module-level :func:`~quantarhei.save_parcel` and
    :func:`~quantarhei.load_parcel` helpers for the common save/load workflow.
    """

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
            dest_dir = os.path.dirname(os.path.abspath(filename))
            tmp_fd, tmp_path = tempfile.mkstemp(dir=dest_dir)
            try:
                with os.fdopen(tmp_fd, "wb") as f:
                    pickle.dump(self, f)
                os.replace(tmp_path, filename)
            except Exception:
                os.unlink(tmp_path)
                raise
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


def load_parcel(filename: str | IO[bytes], *, trusted: bool = False) -> Any:
    """Loads the object saved as parcel

    Parameters
    ----------
    filename : str or File
        Filename of the file or file descriptor of the file from which
        and object should be loaded.

    trusted : bool
        When *False* (default), a :class:`DeserializationWarning` is emitted to
        inform the caller that deserialization may execute arbitrary code.
        Set to *True* only when you are certain the file comes from a
        trusted source.

    """
    if not trusted:
        warnings.warn(_DESERIALIZATION_WARNING, DeserializationWarning, stacklevel=2)

    if isinstance(filename, str):
        with open(filename, "rb") as f:
            obj = pickle.load(f)
    else:
        obj = pickle.load(filename)

    if isinstance(obj, Parcel):
        return obj.content
    raise QuantarheiError("Only Quantarhei Parcels can be loaded")


def check_parcel(filename: str | IO[bytes], *, trusted: bool = False) -> dict[str, Any]:
    """Checks the content of a Quantarhei parcel

    Parameters
    ----------
    filename : str or File
        Filename of the file or file descriptor of the file from which
        and object should be loaded.

    trusted : bool
        When *False* (default), a :class:`DeserializationWarning` is emitted to
        inform the caller that deserialization may execute arbitrary code.
        Set to *True* only when you are certain the file comes from a
        trusted source.

    """
    if not trusted:
        warnings.warn(_DESERIALIZATION_WARNING, DeserializationWarning, stacklevel=2)

    if isinstance(filename, str):
        with open(filename, "rb") as f:
            obj = pickle.load(f)
    else:
        obj = pickle.load(filename)

    if isinstance(obj, Parcel):
        return dict(
            class_name=obj.class_name, qrversion=obj.qrversion, comment=obj.comment
        )
    raise QuantarheiError("The file does not represent a Quantarhei parcel")
