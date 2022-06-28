from .build import BaseBuild as Build

Build(library=True, f2cmap=None).run()

__all__ = ['mapping']

from ._interface import mapping_ as mapping
