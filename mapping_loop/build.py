"""
Package to build f2py modules

typical use may be

create derived class Build(BaseBuild)

Build().run()

from _interface import *
"""

import os
import subprocess
import importlib
import sys
import shutil
import platform
import re
import multiprocessing
from itertools import zip_longest, chain
from shutil import rmtree

from pathlib import Path

from importlib.machinery import EXTENSION_SUFFIXES

import numpy as np

# a sample file

SAMPLE_MAKEFILE = """
$(SOURCE) ?= .

.DEFAULT_GOAL := library.a

LIBRARY_OBJECTS = library.o

library.o: ${SOURCE}/library.f90
	gfortran -c -Ofast -fPIC -o library.o ${SOURCE}/library.f90

library.a: $(LIBRARY_OBJECTS)
	rm -f library.a
	ar cvr $@ $(LIBRARY_OBJECTS)

.PHONY:	clean

clean:
	-rm -f *.o *.a *.mod *.smod *~ \#*\# .*~ .\#*
"""
F2CMAP="""
{
   'real': {'real32': 'float', 'real64': 'double', 'real128': 'long_double'},
   'integer': {'int8': 'signed_char', 'int16': 'short', 'int32': 'int', 'int64': 'long', 'int128': 'long_long'},
   'complex': {'comp32': 'complex_float', 'comp64': 'complex_double', 'comp128': 'complex_long_double'},
   'character': {'char8' : 'char'},
}
"""

class BaseBuild():
    """
    Class to build module binaries.

    Having this in a class and not storing instance keeps namespace clean.
    """

    f2py_options = (
         f'f2py{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}',
         f'f2py{sys.version_info.major}.{sys.version_info.minor}',
         f'f2py{sys.version_info.major}',
         f'f2py'
         )
    for f in f2py_options:
        if shutil.which(f):
            f2py_exec = f
            break
    else:
        raise Exception('f2py not found.')

    f2pycomp = (
        f2py_exec,
        '--verbose',
        # '--debug-capi',
        )
    f2pycomp1 = ()
    f2pycomp2 = ()

    fcomp = 'gfortran -v'

    def _defaults(self, defaults):

        # we may want to add some checks ...
        self.__dict__.update(defaults)

        # set defaults where not provided
        if not hasattr(self, 'path'):
            path = Path(sys.modules[self.__module__].__file__).resolve()
            self.path = path.parent

        if not hasattr(self, 'package'): self.package = 'interface'
        if not hasattr(self, 'parent'): self.parent = self.__module__.rsplit('.', 1)[0]

        # TODO - better config for sources
        if not hasattr(self, 'macos'): self.macros = dict()
        if not hasattr(self, 'sources'): self.sources = (f'{self.package}.f90',)
        if not hasattr(self, 'processed'): self.processed = (None,)
        if not hasattr(self, 'objects'): self.objects= (f'{self.package}.o',)
        if not hasattr(self, 'intermediate_path'): self.intermediate_path = None

        if not hasattr(self, 'include_libraries'): self.include_libraries = () # e.g., 'uuid'
        if not hasattr(self, 'include_paths'): self.include_paths = ()
        if not hasattr(self, 'module'): self.module = f'_{self.package}'
        if not hasattr(self, 'signature_file'): self.signature_file = f'{self.module}.pyf'
        if not hasattr(self, 'compile_flags'): self.compile_flags = (
            '-fPIC',
            '-O3',
            '-funroll-loops',
            '-fno-second-underscore',
            '-fconvert=big-endian',
            )
        if not hasattr(self, 'f2cmap'): self.f2cmap = Path(self.path) / '.f2py_f2cmap'

        if not hasattr(self, 'build_path'): self.build_path = Path(self.path) / '_build'

        if not hasattr(self, 'executable'): self.executable = False
        if not hasattr(self, 'executable_file'): self.executable_file = '{self.package}.exe'
        if not hasattr(self, 'executable_link_flags'): self.executable_link_flags = ()

        # todo - multi-libraries may need more features for native Kepler-like support
        if not hasattr(self, 'library'): self.library = False
        if not hasattr(self, 'libraries'): self.libraries = (
                dict(
                    update = True,
                    build_dir = '_library',
                    source_path = self.path,
                    name = 'library.a',
                    files = ('library.f90',),
                    makefile_path = self.path,
                    makefile = 'Makefile',
                    makeflags = [],
                    ),
                )

        if not hasattr(self, 'clean_build_path'): self.clean_build_path = False

    # does not yet work
    # from numpy.f2py import f2py2e
    # def f2py(self, args):
    #     result = self.f2py2e.run_main(args)
    #     print(result)

    def __init__(self, debug=True, ncpu=None, **defaults):
        """
        init routine, likely to be called before doing own initialisations
        """
        self.debug = debug
        self.ncpu = ncpu

        self._defaults(defaults)

    def f2py(self, args, debug=None):
        assert np.all([isinstance(a, str) for a in args]), f' [{self.__class__.__name__}] all arguments need to be of type str: {args}'
        assert np.all([len(a) > 0 for a in args]), f' [{self.__class__.__name__}] all arguments need to be of length > 0: {args}'
        args = [*self.f2pycomp] + list(args)
        if debug is None:
            debug = self.debug
        if debug:
            print(' [DEBUG][f2py] ' + ' '.join(args))
        result = subprocess.run(
            args,
            check = True,
            shell = False)
        if debug:
            print(f' [DEBUG][f2py] Result: {result}')

    def f2bin(self, args):
        args = [*(self.fcomp.split())] + list(args)
        print(' [DEBUG][f2bin] ' + ' '.join(args))
        result = subprocess.run(
            args,
            shell = True,
            check = True)
        print(f' [DEBUG][f2bin] Result: {result}')

    def run(self, debug=None, ncpu=None):
        """
        execute tests and build
        """
        if debug is None:
            debug = self.debug
        if ncpu is None:
            ncpu = self.ncpu
        kw = dict(debug=debug, ncpu=ncpu)
        if self.build_library_check(**kw) or self.build_check(**kw):
            self.build_module(**kw)
        if self.clean_build_path:
            if self.build_path is not None:
                rmtree(self.build_path)

    def test_executable(self):
        """
        CUSTOM - test exectuable, raise error if there is a problem
        """
        try:
            result = subprocess.run(Path(self.path) / self.executable_file, check = True)
        except subprocess.CalledProcessError:
            raise Exception("module executable failed")
        print(result)
        # TODO - a real test

    def make_executable(self):
        """
        build an executable a test case
        """
        extra_flags = tuple()
        if self.f2cmap is not None:
            if not Path(self.f2cmap).exists():
                Path(self.f2cmap).write_text(F2CMAP)
            extra_flags = ('--f2cmap', str(self.f2cmap),)
        libraries = ()
        if self.library:
            path = self.path
            for l in self.libraries:
                if l['build_dir'] is not None:
                    p = path / l['build_dir']
                else:
                    p = path
                libraries += [p / l['name']]
        try:
            for s,o in zip(self.sources, self.objects):
                self.f2bin([
                    '-c', s, '-o', o,
                    # '-DF2PY_USE_PYTHON_TLS',
                    *self.compile_flags,
                    *self.extra_flags,
                    *chain(('-I', str(p)) for p in self.include_paths),
                    ])
            self.f2bin([
                *self.objects,
                *libraries,
                *self.include_libraries,
                *self.executable_link_flags,
                '-o', self.executable_file,
                ])
        except subprocess.CalledProcessError:
            raise Exception("executable compilation failed")

    def clean_executable(self):
        """
        remove executable file
        """
        try:
            os.remove(self.executable_file)
        except FileNotFoundError:
            pass


    def build_module(self, debug=None, ncpu=None):
        """
        Build python module binary library.

        We also do a test of the executable version
        """

        if ncpu is None:
            ncpu = self.ncpu
        if debug is None:
            debug = self.debug

        cwd = os.getcwd()
        os.chdir(self.path)

        if self.executable:
            # build executable
            self.make_executable()

            # test executable
            self.test_executable()

            # remove executable
            self.clean_executable()

        extra_flags = tuple()
        if self.f2cmap is not None:
            if not Path(self.f2cmap).exists():
                Path(self.f2cmap).write_text(F2CMAP)
            extra_flags = ('--f2cmap', str(self.f2cmap),)

        sources = []
        for s,i in zip_longest(self.sources, self.processed):
            if i is None:
                sources.append(s)
                continue
            if self.intermediate_path is not None:
                i = Path(self.intermediate_path) / i
            else:
                i = Path(self.path) / i
            self.process_macros(s, i, self.macros)
            self.process_includes(i, i)
            sources.append(i)

        try:
            args = [
                '-m', str(self.module),
                ]
            if self.build_path is not None:
                path = Path(self.build_path)
                if not path.exists():
                    path.mkdir(parents=True)
                    if debug:
                        print(f' [DEBUG] creating directory {path}')
            else:
                path = Path(self.path)
            args += [
                '-h', str(path / self.signature_file),
                ]
            if len(self.include_paths) > 0:
                args += [
                    '--include-paths', ':'.join(self.include_paths),
                    ]
            args += [
                *extra_flags,
                *self.f2pycomp1,
                *sources,
                '--overwrite-signature',
               ]
            self.f2py(args)
        except subprocess.CalledProcessError:
            raise Exception("creating f2py signature failed")
        libraries = list()
        if self.library:
            path = Path(self.path)
            for l in self.libraries:
                if l['build_dir'] is not None:
                    p = path / l['build_dir']
                else:
                    p = path
                libraries += [p / l['name']]
        try:
            args = []
            if self.build_path is not None:
                args = [
                    # '--debug',
                    '--build-dir', str(self.build_path),
                ]
            fflags = [
                *chain(
                    self.compile_flags,
                    *chain(('-I', str(p)) for p in self.include_paths),
                    )
                ]
            # include the paths form related libraries
            if self.library:
                path = self.path
                for l in self.libraries:
                    if l['build_dir'] is not None:
                        p = path / l['build_dir']
                    else:
                        p = path
                    fflags += ['-I', str(p)]
            fflags = ' '.join(fflags)
            if len(fflags) > 0:
                args += [
                    f'--f90flags={fflags}',
                    f'--f77flags={fflags}',
                    ]
            if len(self.include_paths) > 0:
                args += [
                    '--include-paths',
                    ':'.join(str(p) for p in self.include_paths),
                    ]
            args += [
                *(f'-l{l}' for l in self.include_libraries),
                *extra_flags,
                *self.f2pycomp2,
                '-c',
                '-m', str(self.module),
                # '-DF2PY_USE_PYTHON_TLS',
                ]
            if self.library:
                path = self.path
                for l in self.libraries:
                    if l['build_dir'] is not None:
                        p = path / l['build_dir']
                    else:
                        p = path
                    args += [str(p / l['name'])]
            args += [
                *sources,
                ]
            if self.build_path is not None:
                path = Path(self.build_path)
            else:
                path = Path(self.path)
            args += [
                str(path / self.signature_file),
                ]
            self.f2py(args)
        except subprocess.CalledProcessError:
            raise Exception("creating module failed")
        os.chdir(cwd)

    def build_library_check(self, debug=None, ncpu=None):
        """
        CUSTOM check whether required libraries are up to date

        return value is whether library needs to be built

        """
        if debug is None:
            sebug = self.debug

        if self.library == False:
            return False
        library_time = 0
        update = False
        for l in self.libraries:
            path = Path(self.path)
            if l['build_dir'] is not None:
                path = path / l['build_dir']
            try:
                filename = path / l['name']
                library_time =  os.stat(filename).st_mtime
                if l['update'] == False:
                    continue
            except FileNotFoundError:
                library_time = 0

            makefile = Path(l['makefile_path']) / l['makefile']

            if not makefile.exists():
                makefile.write_text(SAMPLE_MAKEFILE)

            last_time = os.stat(makefile).st_mtime
            for f in l['files']:
                f = Path(l['source_path']) / f
                last_time = max(last_time, os.stat(f).st_mtime)
            if last_time > library_time:
                cwd = os.getcwd()
                try:
                    os.mkdir(path)
                    if debug:
                        print(f' [DEBUG] creating directory {path}')
                except FileExistsError:
                    pass
                os.chdir(path)
                if ncpu is None:
                    ncpu = multiprocessing.cpu_count()
                if ncpu is Ellipsis:
                    ncpu = 2**10
                cmd = ['make']
                if ncpu > 1:
                    cmd += ['-j', f'{ncpu:d}']
                cmd += ['-f', str(makefile)] + l.get('makeflags', [])
                if l.get('source_path', None) is not None:
                    cmd += [f'SOURCE={l["source_path"]!s}']
                if debug:
                    print(f' [DEBUG] call: {cmd}')
                subprocess.run(cmd, shell=False, check=True)
                os.chdir(cwd)
                update = True
        return update

    def build_check(self, debug=None, ncpu=None):
        """
        check whether build is OK
        """
        if debug is None:
            debug = self.debug

        so_file_base = Path(self.path) / self.module

        for extension in EXTENSION_SUFFIXES:
            so_file = so_file_base.with_suffix(extension)
            if so_file.exists():
                break
        else:
            if debug:
                print(' [DEBUG][build_check] so file does not exist.')
            return True

        source_files = [Path(self.path) / s for s in self.sources]
        so_file_date = os.stat(so_file).st_mtime
        for f in source_files:
            if (so_file_date < os.stat(f).st_mtime):
                if debug:
                    print(f' [DEBUG][build_check] {f} newer than {so_file}.')
                return True
        for i in self.processed:
            if i is None:
                continue
            if self.intermediate_path is not None:
                i = Path(self.intermediate_path) /  i
            else:
                i = Path(self.path) / i
            if not i.exists():
                return True
        try:
            module = importlib.import_module('.' + self.module, self.parent)
        except ImportError as e:
            if debug:
                print(f' [DEBUG][build_check] Import Error: {e}')
                # sys.exit()
            return True
        # check f2py numpy version (since 1.20)
        try:
            np_f2py_version = module.__f2py_numpy_version__
            print(f' [DEBUG] f2py version {np_f2py_version}')
            np_version = np.__version__
            if np_version != np_f2py_version:
                if debug:
                    print(f' [DEBUG][build_check] Library/numpy version mismatch:')
                    print(f' [DEBUG][build_check] Library Version {np_f2py_version}')
                    print(f' [DEBUG][build_check]   Numpy Version {np_version}')
                return True
        except Exception as e:
            print(f' [DEBUG][build_check] module f2py used < 1.20? {e}')

        # check for changed compiler version
        # (works on Fedora 18+)
        # other updates are welcome!
        # seems to have changed with gcc 9.0
        # entirely different in MacOSX - here we require Homebrew
        try:
            if platform.system() == 'Darwin':
                result = subprocess.check_output("gfortran --version", shell=True).decode('ASCII', errors='ignore')
                compiler_version = (result.splitlines()[0]).split(' ', 5)[-1]
                result = subprocess.check_output("strings - {so_file!s} | grep 'GCC version'", shell = True).decode('ASCII', errors='ignore')
                library_version = []
                library_version.append(result.splitlines()[0].split(' ', 2)[-1]) # Homebrew GCC 10.2.0_4
            elif platform.system() == 'Linux':
                result = subprocess.check_output("gcc --version", shell=True).decode('ASCII', errors='ignore')
                compiler_version = (result.splitlines()[0]).split(' ', 2)[-1]
                result = subprocess.check_output(f"strings - {so_file!s} | grep GCC:", shell = True).decode('ASCII', errors='ignore')
                library_version = []
                library_version.append(result.splitlines()[0].split(' ', 2)[2]) # pre 9.0
                library_version.append(result.splitlines()[0].split('GCC:')[1].split(' ', 2)[2]) # 9.1
            if not compiler_version in library_version:
                if debug:
                    print(f' [DEBUG][build_check] Compiler/library version mismatch:')
                    print(f' [DEBUG][build_check] Compiler Version {compiler_version}')
                    print(f' [DEBUG][build_check]  Library Version {library_version}')
                return True
        except:
            if debug:
                print(" [DEBUG][build_check] Compiler comparison failed.")
            return True

        return False


    def process_macros(self, infile, outfile, macros):
        """
        use macro {DATA}
        """
        source = Path(infile).expanduser().read_text()
        for k,v in macros.items():
            source = source.replace(f'{{{k}}}', str(v))
        Path(outfile).expanduser().write_text(source)

    def process_includes(self, infile, outfile):
        """
        use include files

        {include [SOURCE]/filename.f90:section}

        to be replaced by

        !$PY:BEGIN:section
        <stuff to repleace>
        !$PY:END:section

        If no section is provided, insert entire file.

        This allows multiple replacements from multiple files.
        """
        source = Path(infile).expanduser().read_text()
        includes = re.findall(r'(?m)^( *\{insert\s[^\}]+\} *\n)', source)
        for i in includes:
            filesection = re.findall(r'{insert\s+(\S+)\}$', i)[0]
            x = filesection.split(':')
            if len(x) == 2:
                filename, section = x
            else:
                filename = x[0]
                section = None
            if filename.count('[SOURCE]') > 0:
                for inc in self.include_paths:
                    fn = filename.replace('[SOURCE]', inc)
                    if fn.exists():
                        filename = fn
                        break
            include = Path(filename).expanduser().read_text()
            if section is not None:
                include = re.findall(f'(?ms)^\s*\!\$PY:BEGIN:{section} *\n(.*\n) *\!\$PY:END:{section}', include)[0]
            source = source.replace(i, include)
        Path(outfile).expanduser().write_text(source)
