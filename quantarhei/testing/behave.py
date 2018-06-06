# -*- coding: utf-8 -*-
"""
    Support module for tests using `behave` package



"""
# standard imports
import os
import tempfile
import re
from subprocess import check_output

# non-standard imports
import pkg_resources


def quantarhei_installed(context, version=None):
    """Tests if quantarhei is installed and sets version to context if it is


    Parameters
    ----------

    context :
        context object of behave

    version : str, None
        version string, if not none, presence of a specific version is asserted

    """
    import quantarhei as qr
    context.version = qr.Manager().version
    if version is None:
        assert isinstance(context.version, str)
    else:
        assert context.version == version


def shell_command(context, cmd, err_msg="Shell command error"):
    """Runs a shell command in current directory


    Parameters
    ----------

    context :
        context object of behave

    cmd : str
        shell command to run

    err_msg : str
        Error message to put into Exception if the command fails

    """
    try:

        context.last_cmd = cmd
        output = check_output(cmd, shell=True, cwd=os.getcwd())
        context.output = output

    except:
        raise Exception(err_msg)


def check_output_contains(context, text, err_msg):
    """Checks that message contains certain text


    Parameters
    ----------

    context :
        context object of behave

    text : str
        text which is to be found in the context.output string

    err_msg: str
        message to be used for Exception if the text is not found

    """
    res = re.search(text, context.output.decode('utf-8'))
    if res is None:
        print(context.output.decode('utf-8'))
        raise Exception(err_msg)


def secure_temp_dir(context):
    """Creates temporary directory and stores its info into context

    """
    tmpd = tempfile.TemporaryDirectory()
    context.tempdir = tmpd


def cleanup_temp_dir(context):
    """Cleans up temporary directory


    """

    try:
        os.chdir(context.cwd)
    except:
        print("Current working file record does not exist")

    try:
        context.tempdir.cleanup()
    except:
        print("Temporary directory cannot be cleaned up - does it exist?")


def fetch_test_feature_file(context, filename):
    """Fetches the file with a given name from the storage of test files

    """
    resource_package = "quantarhei"
    resource_path = '/'.join(('testing', 'resources', 'behave', filename))

    content = pkg_resources.resource_string(resource_package, resource_path)

    with open(filename, "w") as file:
        file.write(content.decode("utf-8"))

    context.output = ""


class testdir():
    """Context manager for test directory

    With this context manager we enter temporary directory which was
    prepared in the startup phase of the test. The context manager makes sure
    that when tests fail the temporary directory is properly cleaned up.

    """

    def __init__(self, context):

        self.context = context
        try:
            tempdir = context.tempdir
            if tempdir is None:
                raise Exception()
        except:
            raise Exception("Context does not contain info about tempdir")


    def __enter__(self):
        self.context.cwd = os.getcwd()
        os.chdir(self.context.tempdir.name)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.context.cwd)
        if exc_type is not None:
            cleanup_temp_dir(self.context)
