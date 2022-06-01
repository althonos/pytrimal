import configparser
import glob
import os
import platform
import re
import setuptools
import subprocess
import sys
from distutils import log
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants -----------------------------------------------------------------

SETUP_FOLDER = os.path.realpath(os.path.join(__file__, os.pardir))
INCLUDE_FOLDER = os.path.join(SETUP_FOLDER, "vendor", "trimal", "include")

# --- Utils ------------------------------------------------------------------

def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

# --- Commands ------------------------------------------------------------------

class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    # --- Compatibility with `setuptools.Command`

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler

    # --- Build code ---

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "with", self.compiler.compiler_type, "compiler")

        # add debug symbols if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Z7")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # add C++11 flags
        if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
            ext.extra_compile_args.append("-std=c++11")
        elif self.compiler.compiler_type == "msvc":
            ext.extra_compile_args.append("/std:c11")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            ext.define_macros.append(("WIN32", 1))

        # update link and include directories
        for name in ext.libraries:
            lib = self._clib_cmd.get_library(name)
            libfile = self.compiler.library_filename(
                lib.name, output_dir=self._clib_cmd.build_clib
            )
            ext.depends.append(libfile)
            ext.extra_objects.append(libfile)

        # build the rest of the extension as normal
        ext._needs_stub = False
        _build_ext.build_extension(self, ext)

    def build_extensions(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {
                "cdivision": True,
                "nonecheck": False,
            },
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
            }
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["cdivision_warnings"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False

        # compile the C library
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # add the include dirs
        for ext in self.extensions:
            ext.include_dirs.append(self._clib_cmd.build_clib)

        # cythonize the extensions
        self.extensions = cythonize(self.extensions, **cython_args)

        # build the extensions as normal
        _build_ext.build_extensions(self)


class build_clib(_build_clib):
    """A custom `build_clib` that makes all C++ class attributes public.
    """

    # --- Autotools-like helpers ---

    def _publicize(self, input, output):
        with open(input, "rb") as src:
            with open(output, "wb") as dst:
                for line in src:
                    dst.write(line.replace(b"private:", b"public:"))

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    def get_library(self, name):
        return next(lib for lib in self.libraries if lib.name == name)

    # --- Build code ---

    def build_libraries(self, libraries):
        # build each library only if the sources are outdated
        self.mkpath(self.build_clib)
        for library in libraries:
            libname = self.compiler.library_filename(library.name, output_dir=self.build_clib)
            self.make_file(library.sources, libname, self.build_library, (library,))

    def build_library(self, library):
        # show the compiler being used
        _eprint("building", library.name, "with", self.compiler.compiler_type, "compiler")

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Z7")

        # add C++11 flags
        if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
            library.extra_compile_args.append("-std=c++11")
        elif self.compiler.compiler_type == "msvc":
            library.extra_compile_args.append("/std:c11")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            library.define_macros.append(("WIN32", 1))

        # expose all private members and copy headers to build directory
        for header in library.depends:
            output = os.path.join(self.build_clib, os.path.relpath(header, INCLUDE_FOLDER))
            self.mkpath(os.path.dirname(output))
            self.make_file(
                [header],
                output,
                self._publicize,
                (header, output)
            )

        # copy sources to build directory
        sources = [
            os.path.join(self.build_temp, os.path.basename(source))
            for source in library.sources
        ]
        for source, source_copy in zip(library.sources, sources):
            self.make_file(
                [source],
                source_copy,
                self.copy_file,
                (source, source_copy)
            )

        # store compile args
        compile_args = (
            library.define_macros,
            library.include_dirs + [self.build_clib],
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )
        # manually prepare sources and get the names of object files
        objects = [
            s.replace(".cpp", self.compiler.obj_extension)
            for s in sources
        ]
        # only compile outdated files
        for source, object in zip(sources, objects):
            self.make_file(
                [source],
                object,
                self.compiler.compile,
                ([source], None, *compile_args),
            )

        # link into a static library
        libfile = self.compiler.library_filename(
            library.name,
            output_dir=self.build_clib,
        )
        self.make_file(
            objects,
            libfile,
            self.compiler.create_static_lib,
            (objects, library.name, self.build_clib, None, self.debug)
        )


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython.
    """

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pymuscle5")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c", "*.cpp"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                log.info("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)

# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    libraries=[
        Library(
            "trimal",
            language="c++",
            libraries=["m"],
            sources=[
                # commonFiles
                os.path.join("vendor", "trimal", "source", "Cleaner.cpp"),
                os.path.join("vendor", "trimal", "source", "Alignment", "Alignment.cpp"),
                os.path.join("vendor", "trimal", "source", "Alignment", "sequencesMatrix.cpp"),
                os.path.join("vendor", "trimal", "source", "Statistics", "similarityMatrix.cpp"),
                # statisticsFiles
                os.path.join("vendor", "trimal", "source", "Statistics", "Mold.cpp"),
                os.path.join("vendor", "trimal", "source", "Statistics", "Gaps.cpp"),
                os.path.join("vendor", "trimal", "source", "Statistics", "Manager.cpp"),
                os.path.join("vendor", "trimal", "source", "Statistics", "Similarity.cpp"),
                os.path.join("vendor", "trimal", "source", "Statistics", "Consistency.cpp"),
                # reportSystemFiles
                # os.path.join("vendor", "trimal", "source", "reportsystem"),
                os.path.join("vendor", "trimal", "source", "reportMessages", "errorMessages.cpp"),
                os.path.join("vendor", "trimal", "source", "reportMessages", "infoMessages.cpp"),
                os.path.join("vendor", "trimal", "source", "reportMessages", "warningMessages.cpp"),
                # formatHandlerFiles
                *glob.iglob(os.path.join("vendor", "trimal", "source", "FormatHandling", "*_state.cpp")),
                # formatHandler
                os.path.join("vendor", "trimal", "source", "FormatHandling", "BaseFormatHandler.cpp"),
                # utils
                os.path.join("vendor", "trimal", "source", "utils.cpp"),
                # Internal Benchmarker
                os.path.join("vendor", "trimal", "source", "InternalBenchmarker.cpp"),
                # trimAl manager
                os.path.join("vendor", "trimal", "source", "trimalManager.cpp"),
                os.path.join("vendor", "trimal", "source", "VCFHandler.cpp"),
            ],
            depends=[
                os.path.join(d, file)
                for d, dirs, files in os.walk(INCLUDE_FOLDER)
                for file in files
                if file.endswith((".h", ".txt"))
            ],
            include_dirs=[
                os.path.join("pytrimal", "patch"),
                "pytrimal",
                "include",
            ],
        ),
    ],
    ext_modules=[
        Extension(
            "pytrimal._trimal",
            language="c++",
            sources=[
                os.path.join("pytrimal", "patch", "reportsystem.cpp"),
                os.path.join("pytrimal", "_trimal.pyx"),
            ],
            include_dirs=[
                os.path.join("pytrimal", "patch"),
                "pytrimal",
                "include",
            ],
            libraries=[
                "trimal",
            ],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean
    }
)
