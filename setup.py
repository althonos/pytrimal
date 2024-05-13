import configparser
import functools
import glob
import itertools
import io
import multiprocessing.pool
import os
import platform
import re
import setuptools
import setuptools.extension
import subprocess
import sys
import sysconfig
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants -----------------------------------------------------------------

SETUP_FOLDER = os.path.realpath(os.path.join(__file__, os.pardir))
INCLUDE_FOLDER = os.path.join(SETUP_FOLDER, "vendor", "trimal", "include")

# --- Utils ------------------------------------------------------------------

_HEADER_PATTERN = re.compile("^@@ -(\d+),?(\d+)? \+(\d+),?(\d+)? @@$")


def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def _patch_osx_compiler(compiler):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of SSE2.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next(
            (
                i
                for i in range(1, len(flags))
                if flags[i - 1] == "-arch" and flags[i] != platform.machine()
            ),
            None,
        )
        if i is not None:
            flags.pop(i)
            flags.pop(i - 1)


def _detect_target_machine(platform):
    if platform == "win32":
        return "x86"
    return platform.rsplit("-", 1)[-1]


def _detect_target_cpu(platform):
    machine = _detect_target_machine(platform)
    if re.match("^mips", machine):
        return "mips"
    elif re.match("^(aarch64|arm64)$", machine):
        return "aarch64"
    elif re.match("^arm", machine):
        return "arm"
    elif re.match("(x86_64)|AMD64|amd64", machine):
        return "x86_64"
    elif re.match("(x86)|(^i.86$)", machine):
        return "x86"
    elif re.match("^(powerpc|ppc)", machine):
        return "ppc"
    return None


def _detect_target_system(platform):
    if platform.startswith("win"):
        return "windows"
    elif platform.startswith("macos"):
        return "macos"
    elif platform.startswith("linux"):
        return "linux_or_android"
    elif platform.startswith("freebsd"):
        return "freebsd"
    return None


def _apply_patch(s, patch, revert=False):
    # see https://stackoverflow.com/a/40967337
    s = s.splitlines(keepends=True)
    p = patch.splitlines(keepends=True)
    t = []
    i = 0
    sl = 0
    midx, sign = (1, "+") if not revert else (3, "-")
    while i < len(p) and p[i].startswith(("---", "+++")):
        i += 1  # skip header lines

    while i < len(p):
        match = _HEADER_PATTERN.match(p[i])
        if not match:
            raise ValueError("Invalid line in patch: {!r}".format(p[i]))
        i += 1
        l = int(match.group(midx)) - 1 + (match.group(midx + 1) == "0")
        t.extend(s[sl:l])
        sl = l
        while i < len(p) and p[i][0] != "@":
            if i + 1 < len(p) and p[i + 1][0] == "\\":
                line = p[i][:-1]
                i += 2
            else:
                line = p[i]
                i += 1
            if len(line) > 0:
                if line[0] == sign or line[0] == " ":
                    t += line[1:]
                sl += line[0] != sign

    t.extend(s[sl:])
    return "".join(t)


# --- Extension with SIMD support --------------------------------------------


class Library(setuptools.extension.Library):
    def __init__(self, *args, **kwargs):
        self.platform_sources = kwargs.pop("platform_sources", {})
        super().__init__(*args, **kwargs)


# --- Commands ------------------------------------------------------------------


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly."""

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", "build-backend", '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode."""

    # --- Compatibility with `setuptools.Command`

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # check platform
        if self.plat_name is None:
            self.plat_name = sysconfig.get_platform()
        # detect platform options
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler
        self._clib_cmd.parallel = self.parallel
        self._clib_cmd.plat_name = self.plat_name
        self._clib_cmd.target_machine = self.target_machine
        self._clib_cmd.target_cpu = self.target_cpu
        self._clib_cmd.target_system = self.target_system

    # --- Autotools-like helpers ---

    def _check_getid(self):
        _eprint('checking whether `PyInterpreterState_GetID` is available')

        base = "have_getid"
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write("""
            #include <stdint.h>
            #include <stdlib.h>
            #include <Python.h>

            int main(int argc, char *argv[]) {{
                PyInterpreterState_GetID(NULL);
                return 0;
            }}
            """)

        if self.compiler.compiler_type == "msvc":
            flags = ["/WX"]
        else:
            flags = ["-Werror=implicit-function-declaration"]

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_postargs=flags)
        except CompileError:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

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
            ext.extra_compile_args.append("-funroll-loops")
            ext.extra_link_args.append("-Wno-alloc-size-larger-than")
        elif self.compiler.compiler_type == "msvc":
            ext.extra_compile_args.append("/std:c11")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            ext.define_macros.append(("WIN32", 1))

        # check `PyInterpreterState_GetID` is present
        if self._check_getid():
            ext.define_macros.append(("HAS_PYINTERPRETERSTATE_GETID", 1))

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
            raise RuntimeError(
                "Cython is required to run `build_ext` command"
            ) from cythonize

        # remove universal compilation flags for OSX
        if platform.system() == "Darwin":
            _patch_osx_compiler(self.compiler)

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
                "DEFAULT_BUFFER_SIZE": io.DEFAULT_BUFFER_SIZE,
                "TARGET_CPU": self.target_cpu,
                "TARGET_SYSTEM": self.target_system,
                "AVX2_BUILD_SUPPORT": False,
                "NEON_BUILD_SUPPORT": False,
                "SSE2_BUILD_SUPPORT": False,
            },
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
        
        # check which SIMD features are supported
        for simd in self._clib_cmd._simd_supported:
            if self._clib_cmd._simd_supported[simd]:
                cython_args["compile_time_env"][f"{simd}_BUILD_SUPPORT"] = True

        # add the include dirs
        for ext in self.extensions:
            ext.include_dirs.append(self._clib_cmd.build_clib)

        # cythonize the extensions (retaining platform-specific sources)
        self.extensions = cythonize(self.extensions, **cython_args)

        # build the extensions as normal
        _build_ext.build_extensions(self)


class build_clib(_build_clib):
    """A custom `build_clib` that makes all C++ class attributes public."""

    # --- Compatibility with `setuptools.Command`

    user_options = _build_clib.user_options + [
        ("parallel", "j", "number of parallel build jobs"),
        (
            "disable-avx2",
            None,
            "Force compiling the extension without AVX2 instructions",
        ),
        (
            "disable-mmx",
            None,
            "Force compiling the extension without MMX instructions",
        ),
        (
            "disable-sse2",
            None,
            "Force compiling the extension without SSE2 instructions",
        ),
        (
            "disable-neon",
            None,
            "Force compiling the extension without NEON instructions",
        ),
    ]

    def initialize_options(self):
        _build_clib.initialize_options(self)
        self.parallel = None
        self.disable_avx2 = False
        self.disable_sse2 = False
        self.disable_neon = False
        self.target_machine = None
        self.target_system = None
        self.target_cpu = None
        self.plat_name = None

    def finalize_options(self):
        _build_clib.finalize_options(self)
        if self.parallel is not None:
            self.parallel = int(self.parallel)
        # check platform
        if self.plat_name is None:
            self.plat_name = sysconfig.get_platform()
        # detect platform options
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
         # record SIMD-specific options
        self._simd_supported = dict(AVX2=False, SSE2=False, NEON=False)
        self._simd_defines = dict(AVX2=[], SSE2=[], NEON=[])
        self._simd_flags = dict(AVX2=[], SSE2=[], NEON=[])
        self._simd_disabled = {
            "AVX2": self.disable_avx2,
            "SSE2": self.disable_sse2,
            "NEON": self.disable_neon,
        }

    # --- Autotools-like helpers ---

    def _patch_file(self, input, output):
        basename = os.path.basename(input)
        patchname = os.path.realpath(
            os.path.join(__file__, os.pardir, "patches", "{}.patch".format(basename))
        )
        if os.path.exists(patchname):
            _eprint(
                "patching", os.path.relpath(input), "with", os.path.relpath(patchname)
            )
            with open(patchname, "r") as patchfile:
                patch = patchfile.read()
            with open(input, "r") as src:
                srcdata = src.read()
            with open(output, "w") as dst:
                dst.write(_apply_patch(srcdata, patch))
        else:
            self.copy_file(input, output)

    def _check_function(self, funcname, header, args="()"):
        _eprint("checking whether function", repr(funcname), "is available", end="... ")

        self.mkpath(self.build_temp)
        base = "have_{}".format(funcname)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        with open(testfile, "w") as f:
            f.write(
                """
                #include <{}>
                int main() {{
                    {}{};
                    return 0;
                }}
            """.format(
                    header, funcname, args
                )
            )
        try:
            objects = self.compiler.compile([testfile], debug=self.debug)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
        except:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _check_simd_generic(self, name, flags, program):
        _eprint("checking whether compiler can build", name, "code", end="... ")

        base = "have_{}".format(name)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write(program)

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_preargs=flags)
            self.compiler.link_executable(objects, base, extra_preargs=flags, output_dir=self.build_temp)
            subprocess.run([binfile], check=True)
        except CompileError:
            _eprint("no")
            return False
        except subprocess.CalledProcessError:
            _eprint("yes, but cannot run code")
            return True  # assume we are cross-compiling, and still build
        else:
            if not flags:
                _eprint("yes")
            else:
                _eprint("yes, with {}".format(" ".join(flags)))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _avx2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:AVX2"]
        return ["-mavx", "-mavx2"]

    def _check_avx2(self):
        return self._check_simd_generic(
            "AVX2",
            self._avx2_flags(),
            program="""
                #include <immintrin.h>
                int main() {{
                    __m256i a = _mm256_set1_epi16(-1);
                            a = _mm256_abs_epi16(a);
                    short   x = _mm256_extract_epi16(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """,
        )

    def _sse2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:SSE2"]
        return ["-msse", "-msse2"]

    def _check_sse2(self):
        return self._check_simd_generic(
            "SSE2",
            self._sse2_flags(),
            program="""
                #include <emmintrin.h>
                int main() {{
                    __m128i a = _mm_set1_epi16(-1);
                            a = _mm_and_si128(a, a);
                    short   x = _mm_extract_epi16(a, 1);
                    return (x == -1) ? 0 : 1;
                }}
            """,
        )

    def _neon_flags(self):
        return ["-mfpu=neon"] if self.target_cpu == "arm" else []

    def _check_neon(self):
        return self._check_simd_generic(
            "NEON",
            self._neon_flags(),
            program="""
                #include <arm_neon.h>
                int main() {{
                    int16x8_t a = vdupq_n_s16(-1);
                              a = vabsq_s16(a);
                    short     x = vgetq_lane_s16(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """,
        )

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [lib.name for lib in self.libraries]

    def get_source_files(self):
        return [source for lib in self.libraries for source in lib.sources]

    def get_library(self, name):
        return next(lib for lib in self.libraries if lib.name == name)

    # --- Build code ---

    def build_simd_code(self, lib):
        # build platform-specific code
        for simd, sources in lib.platform_sources.items():
            if self._simd_supported[simd] and not self._simd_disabled[simd]:
                objects = [
                    os.path.join(
                        self.build_temp, s.replace(".cpp", self.compiler.obj_extension)
                    )
                    for s in sources
                ]
                for source, object in zip(sources, objects):
                    self.make_file(
                        [source],
                        object,
                        self.compiler.compile,
                        (
                            [source],
                            self.build_temp,
                            lib.define_macros + self._simd_defines[simd],
                            lib.include_dirs,
                            self.debug,
                            lib.extra_compile_args + self._simd_flags[simd],
                            None,
                            lib.depends,
                        ),
                    )
                lib.extra_objects.extend(objects)

    def build_libraries(self, libraries):
        # remove universal compilation flags for OSX
        if platform.system() == "Darwin":
            _patch_osx_compiler(self.compiler)

        # check if we can build platform-specific code
        if self.target_cpu == "x86" or self.target_cpu == "x86_64":
            if not self._simd_disabled["AVX2"] and self._check_avx2():
                self._simd_supported["AVX2"] = True
                self._simd_flags["AVX2"].extend(self._avx2_flags())
                self._simd_defines["AVX2"].append(("__AVX2__", 1))
            if not self._simd_disabled["SSE2"] and self._check_sse2():
                self._simd_supported["SSE2"] = True
                self._simd_flags["SSE2"].extend(self._sse2_flags())
                self._simd_defines["SSE2"].append(("__SSE2__", 1))
        elif self.target_cpu == "arm" or self.target_cpu == "aarch64":
            if not self._simd_disabled["NEON"] and self._check_neon():
                self._simd_supported["NEON"] = True
                self._simd_flags["NEON"].extend(self._neon_flags())
                self._simd_defines["NEON"].append(("__ARM_NEON", 1))

        # build each library only if the sources are outdated
        self.mkpath(self.build_clib)
        for library in libraries:
            libname = self.compiler.library_filename(
                library.name, output_dir=self.build_clib
            )
            self.make_file(library.sources, libname, self.build_library, (library,))

    def build_library(self, library):
        # show the compiler being used
        _eprint(
            "building", library.name, "with", self.compiler.compiler_type, "compiler"
        )

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Z7")

        # add C++11 flags
        if library.language == "c++":
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-std=c++11")
                library.extra_link_args.append("-Wno-alloc-size-larger-than")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/std:c11")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            library.define_macros.append(("WIN32", 1))

        # add SIMD flags
        if library.name == "trimal":
            for simd in self._simd_supported:
                if self._simd_supported[simd]:
                    library.define_macros.append((f"HAVE_{simd}", 1))

        # expose all private members and copy headers to build directory
        for header in library.depends:
            output = os.path.join(
                self.build_clib, os.path.relpath(header, INCLUDE_FOLDER)
            )
            self.mkpath(os.path.dirname(output))
            self.make_file([header], output, self._patch_file, (header, output))

        # copy sources to build directory
        sources = [
            os.path.join(self.build_temp, os.path.basename(source))
            for source in library.sources
        ]
        for source, source_copy in zip(library.sources, sources):
            self.make_file(
                [source], source_copy, self._patch_file, (source, source_copy)
            )

        # store compile args
        compile_args = (
            None,
            library.define_macros,
            library.include_dirs + [self.build_clib],
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )
        # manually prepare sources and get the names of object files
        objects = [
            re.sub(r"(.cpp|.c)$", self.compiler.obj_extension, s) for s in sources
        ]
        # only compile outdated files
        with multiprocessing.pool.ThreadPool(self.parallel) as pool:
            pool.starmap(
                functools.partial(self._compile_file, compile_args=compile_args),
                zip(sources, objects),
            )

        # build SIMD code
        self.build_simd_code(library)
        objects.extend(library.extra_objects)

        # link into a static library
        libfile = self.compiler.library_filename(
            library.name,
            output_dir=self.build_clib,
        )
        self.make_file(
            objects,
            libfile,
            self.compiler.create_static_lib,
            (objects, library.name, self.build_clib, None, self.debug),
        )

    def _compile_file(self, source, object, compile_args):
        self.make_file(
            [source], object, self.compiler.compile, ([source], *compile_args)
        )


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython."""

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pytrimal")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c", "*.cpp"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                _eprint("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)


# --- Setup ---------------------------------------------------------------------

CPU_FEATURES = Library(
    "cpu_features",
    language="c",
    define_macros=[("STACK_LINE_READER_BUFFER_SIZE", 1024)],
    sources=[
        os.path.join("vendor", "trimal", "vendor", "cpu_features", "src", "filesystem.c"),
        os.path.join("vendor", "trimal", "vendor", "cpu_features", "src", "stack_line_reader.c"),
        os.path.join("vendor", "trimal", "vendor", "cpu_features", "src", "string_view.c"),
        # os.path.join("vendor", "trimal", "vendor", "cpu_features", "src", "hwcaps.c"),
        *glob.glob(os.path.join("vendor", "trimal", "vendor", "cpu_features", "src", "impl_*.c")),
    ],
    include_dirs=[
        os.path.join("vendor", "trimal", "vendor", "cpu_features", "include"),
        os.path.join("vendor", "trimal", "vendor", "cpu_features", "src"),
    ]
)

TRIMAL = Library(
    "trimal",
    language="c++",
    libraries=["m"],
    sources=[
        # commonFiles
        os.path.join("vendor", "trimal", "source", "Cleaner.cpp"),
        os.path.join(
            "vendor", "trimal", "source", "Alignment", "Alignment.cpp"
        ),
        os.path.join(
            "vendor", "trimal", "source", "Alignment", "sequencesMatrix.cpp"
        ),
        os.path.join(
            "vendor", "trimal", "source", "Statistics", "similarityMatrix.cpp"
        ),
        # statisticsFiles
        os.path.join("vendor", "trimal", "source", "Statistics", "Mold.cpp"),
        os.path.join("vendor", "trimal", "source", "Statistics", "Gaps.cpp"),
        os.path.join("vendor", "trimal", "source", "Statistics", "Manager.cpp"),
        os.path.join(
            "vendor", "trimal", "source", "Statistics", "Similarity.cpp"
        ),
        os.path.join(
            "vendor", "trimal", "source", "Statistics", "Consistency.cpp"
        ),
         os.path.join(
            "vendor", "trimal", "source", "Statistics", "Identity.cpp"
        ),
        os.path.join(
            "vendor", "trimal", "source", "Statistics", "Overlap.cpp"
        ),
        os.path.join(
            "vendor", "trimal", "source", "Statistics", "similarityMatrix.cpp"
        ),
        # reportSystemFiles
        # os.path.join("vendor", "trimal", "source", "reportsystem"),
        os.path.join(
            "vendor", "trimal", "source", "reportMessages", "errorMessages.cpp"
        ),
        os.path.join(
            "vendor", "trimal", "source", "reportMessages", "infoMessages.cpp"
        ),
        os.path.join(
            "vendor",
            "trimal",
            "source",
            "reportMessages",
            "warningMessages.cpp",
        ),
        # formatHandlerFiles
        *glob.iglob(
            os.path.join(
                "vendor", "trimal", "source", "FormatHandling", "*_state.cpp"
            )
        ),
        # formatHandler
        os.path.join(
            "vendor",
            "trimal",
            "source",
            "FormatHandling",
            "BaseFormatHandler.cpp",
        ),
        # utils
        os.path.join("vendor", "trimal", "source", "utils.cpp"),
        # Internal Benchmarker
        os.path.join("vendor", "trimal", "source", "InternalBenchmarker.cpp"),
        # trimAl manager
        os.path.join("vendor", "trimal", "source", "trimalManager.cpp"),
        os.path.join("vendor", "trimal", "source", "VCFHandler.cpp"),
    ],

    platform_sources={
        "SSE2": [os.path.join("vendor", "trimal", "source", "Platform", "x86", "SSE2.cpp")],
        "NEON": [os.path.join("vendor", "trimal", "source", "Platform", "Arm", "NEON.cpp")],
        "AVX2": [os.path.join("vendor", "trimal", "source", "Platform", "x86", "AVX2.cpp")],
    },

    depends=[
        os.path.join(d, file)
        for d, dirs, files in os.walk(INCLUDE_FOLDER)
        for file in files
        if file.endswith((".h", ".txt"))
    ],
    include_dirs=[
        os.path.join("vendor", "trimal", "vendor", "cpu_features", "include"),
        os.path.join("vendor", "trimal", "include"),
        os.path.join("pytrimal", "patch"),
        "pytrimal",
        "include",
    ],
)

setuptools.setup(
    libraries=[CPU_FEATURES, TRIMAL],
    ext_modules=[
        Extension(
            "pytrimal._trimal",
            language="c++",
            sources=[
                os.path.join("pytrimal", "fileobj", "pywritebuf.cpp"),
                os.path.join("pytrimal", "fileobj", "pyreadbuf.cpp"),
                os.path.join("pytrimal", "fileobj", "pyreadintobuf.cpp"),
                os.path.join("pytrimal", "patch", "reportsystem.cpp"),
                # os.path.join("pytrimal", "impl", "generic.cpp"),
                os.path.join("pytrimal", "_trimal.pyx"),
            ],
            include_dirs=[
                os.path.join("pytrimal", "patch"),
                os.path.join("pytrimal", "fileobj"),
                os.path.join("pytrimal", "impl"),
                os.path.join("vendor", "trimal", "vendor", "cpu_features", "include"),
                "pytrimal",
                "include",
            ],
            libraries=[
                "cpu_features",
                "trimal",
            ],
            depends=[
                os.path.join("pytrimal", "impl", "template.h"),
            ]
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean,
    },
)
