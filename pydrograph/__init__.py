from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .nwis import Nwis
from .attributes import gw_attributes, streamflow_attributes
from .baseflow import get_upstream_area, IHmethod, WI_statewide_eqn