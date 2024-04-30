import importlib.metadata as _meta

from packaging.version import parse as _parse

version = _meta.version("streamtracer")
_version = _parse(version)
major, minor, bugfix = [*_version.release, 0][:3]
release = not _version.is_devrelease
