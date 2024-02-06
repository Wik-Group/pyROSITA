"""
Utility functions for basic functionality of the py:module:`cluster_generator` package.
"""
import logging
import os
import pathlib as pt
import sys

import yaml
from unyt import unyt_array

# -- configuration directory -- #
_config_directory = os.path.join(pt.Path(__file__).parents[0], "bin", "config.yaml")


# defining the custom yaml loader for unit-ed objects
def _yaml_unit_constructor(loader: yaml.FullLoader, node: yaml.nodes.MappingNode):
    kw = loader.construct_mapping(node)
    i_s = kw["input_scalar"]
    del kw["input_scalar"]
    return unyt_array(i_s, **kw)


def _yaml_lambda_loader(loader: yaml.FullLoader, node: yaml.nodes.ScalarNode):
    return eval(loader.construct_scalar(node))


def _get_loader():
    loader = yaml.FullLoader
    loader.add_constructor("!unyt", _yaml_unit_constructor)
    loader.add_constructor("!lambda", _yaml_lambda_loader)
    return loader


try:
    with open(_config_directory, "r+") as config_file:
        cgparams = yaml.load(config_file, _get_loader())

except FileNotFoundError as er:
    raise FileNotFoundError(
        f"Couldn't find the configuration file! Is it at {_config_directory}? Error = {er.__repr__()}"
    )
except yaml.YAMLError as er:
    raise yaml.YAMLError(
        f"The configuration file is corrupted! Error = {er.__repr__()}"
    )


stream = (
    sys.stdout
    if cgparams["system"]["logging"]["main"]["stream"] in ["STDOUT", "stdout"]
    else sys.stderr
)
cgLogger = logging.getLogger("pyROSITA")

cg_sh = logging.StreamHandler(stream=stream)

# create formatter and add it to the handlers
formatter = logging.Formatter(cgparams["system"]["logging"]["main"]["format"])
cg_sh.setFormatter(formatter)
# add the handler to the logger
cgLogger.addHandler(cg_sh)
cgLogger.setLevel(cgparams["system"]["logging"]["main"]["level"])
cgLogger.propagate = False

mylog = cgLogger

# -- Setting up the developer debugger -- #
devLogger = logging.getLogger("development_logger")

if cgparams["system"]["logging"]["developer"][
    "enabled"
]:  # --> We do want to use the development logger.
    # -- checking if the user has specified a directory -- #
    if cgparams["system"]["logging"]["developer"]["output_directory"] is not None:
        from datetime import datetime

        dv_fh = logging.FileHandler(
            os.path.join(
                cgparams["system"]["logging"]["developer"]["output_directory"],
                f"{datetime.now().strftime('%m-%d-%y_%H-%M-%S')}.log",
            )
        )

        # adding the formatter
        dv_formatter = logging.Formatter(
            cgparams["system"]["logging"]["main"]["format"]
        )

        dv_fh.setFormatter(dv_formatter)
        devLogger.addHandler(dv_fh)
        devLogger.setLevel("DEBUG")
        devLogger.propagate = False

    else:
        mylog.warning(
            "User enabled development logger but did not specify output directory. Dev logger will not be used."
        )
else:
    devLogger.propagate = False
    devLogger.disabled = True


class EHalo:
    def __new__(cls, *args, **kwargs):
        if "text" in kwargs:
            kwargs["text"] = "[pyROSITA] " + kwargs["text"]
        kwargs["stream"] = sys.stderr

        try:
            from halo import Halo

            return Halo(*args, **kwargs)
        except ImportError:
            return cls(*args, **kwargs)

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, exc_tb):
        pass


def split(a, n):
    k, m = divmod(len(a), n)
    return [a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n)]
