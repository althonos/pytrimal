import configparser
import glob
import itertools
import os
import platform
import re
import setuptools
import setuptools.extension
import subprocess
import sys
from distutils import log
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants -----------------------------------------------------------------

SETUP_FOLDER = os.path.realpath(os.path.join(__file__, os.pardir))
INCLUDE_FOLDER = os.path.join(SETUP_FOLDER, "vendor", "trimal", "include")

MACHINE = platform.machine()
if re.match("^mips", MACHINE):
    TARGET_CPU = "mips"
elif re.match("^(aarch64|arm64)$", MACHINE):
    TARGET_CPU = "aarch64"
elif re.match("^arm", MACHINE):
    TARGET_CPU = "arm"
elif re.match("(x86_64)|(AMD64|amd64)|(^i.86$)", MACHINE):
    TARGET_CPU = "x86"
elif re.match("^(powerpc|ppc)", MACHINE):
    TARGET_CPU = "ppc"
else:
    TARGET_CPU = None

SYSTEM  = platform.system()
if SYSTEM == "Linux" or SYSTEM == "Java":
    TARGET_SYSTEM = "linux_or_android"
elif SYSTEM.endswith("FreeBSD"):
    TARGET_SYSTEM = "freebsd"
elif SYSTEM == "Darwin":
    TARGET_SYSTEM = "macos"
elif SYSTEM.startswith(("Windows", "MSYS", "MINGW", "CYGWIN")):
    TARGET_SYSTEM = "windows"
else:
    TARGET_SYSTEM = None

# --- Utils ------------------------------------------------------------------

def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


# --- Extension with SIMD support --------------------------------------------

class Extension(setuptools.extension.Extension):

    def __init__(self, *args, **kwargs):
        self._needs_stub = False
        self.platform_sources = kwargs.pop("platform_sources", {})
        super().__init__(*args, **kwargs)


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

    user_options = _build_ext.user_options + [
        ("disable-avx2", None, "Force compiling the extension without AVX2 instructions"),
        ("disable-sse2", None, "Force compiling the extension without SSE2 instructions"),
        ("disable-neon", None, "Force compiling the extension without NEON instructions"),
    ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.disable_avx2 = False
        self.disable_sse2 = False
        self.disable_neon = False

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # record SIMD-specific options
        self._simd_supported = dict(AVX2=False, SSE2=False, NEON=False)
        self._simd_defines = dict(AVX2=[], SSE2=[], NEON=[])
        self._simd_flags = dict(AVX2=[], SSE2=[], NEON=[])
        self._simd_disabled = {
            "AVX2": self.disable_avx2,
            "SSE2": self.disable_sse2,
            "NEON": self.disable_neon,
        }
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler

    # --- Autotools-like helpers ---

    def _check_simd_generic(self, name, flags, header, vector, set, op, extract):
        _eprint('checking whether compiler can build', name, 'code', end="... ")

        base = "have_{}".format(name)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}      a = {}(1);
                            a = {}(a);
                    short   x = {}(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """.format(header, vector, set, op, extract))

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_preargs=flags)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
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
            header="immintrin.h",
            vector="__m256i",
            set="_mm256_set1_epi16",
            op="_mm256_abs_epi32",
            extract="_mm256_extract_epi16",
        )

    def _sse2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:SSE2"]
        return ["-msse", "-msse2"]

    def _check_sse2(self):
        return self._check_simd_generic(
            "SSE2",
            self._sse2_flags(),
            header="emmintrin.h",
            vector="__m128i",
            set="_mm_set1_epi16",
            op="_mm_move_epi64",
            extract="_mm_extract_epi16",
        )

    def _neon_flags(self):
        return ["-mfpu=neon"] if TARGET_CPU == "arm" else []

    def _check_neon(self):
        return self._check_simd_generic(
            "NEON",
            self._neon_flags(),
            header="arm_neon.h",
            vector="int16x8_t",
            set="vdupq_n_s16",
            op="vabsq_s16",
            extract="vgetq_lane_s16"
        )

    # --- Build code ---

    def build_simd_code(self, ext):
        # build platform-specific code
        for simd, sources in ext.platform_sources.items():
            if self._simd_supported[simd] and not self._simd_disabled[simd]:
                objects = [
                    os.path.join(self.build_temp, s.replace(".cpp", self.compiler.obj_extension))
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
                            ext.define_macros + self._simd_defines[simd],
                            ext.include_dirs,
                            self.debug,
                            ext.extra_compile_args + self._simd_flags[simd],
                            None,
                            ext.depends
                        )
                    )
                ext.extra_objects.extend(objects)

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
            ext.extra_link_args.append("-Wno-alloc-size-larger-than")
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

        # build platform-specific code
        self.build_simd_code(ext)

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
                "TARGET_CPU": TARGET_CPU,
                "TARGET_SYSTEM": TARGET_SYSTEM,
                "AVX2_BUILD_SUPPORT": False,
                "NEON_BUILD_SUPPORT": False,
                "SSE2_BUILD_SUPPORT": False,
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

        # check if we can build platform-specific code
        if TARGET_CPU == "x86":
            # if not self._simd_disabled["AVX2"] and self._check_avx2():
            #     cython_args["compile_time_env"]["AVX2_BUILD_SUPPORT"] = True
            #     self._simd_supported["AVX2"] = True
            #     self._simd_flags["AVX2"].extend(self._avx2_flags())
            #     self._simd_defines["AVX2"].append(("__AVX2__", 1))
            if not self._simd_disabled["SSE2"] and self._check_sse2():
                cython_args["compile_time_env"]["SSE2_BUILD_SUPPORT"] = True
                self._simd_supported["SSE2"] = True
                self._simd_flags["SSE2"].extend(self._sse2_flags())
                self._simd_defines["SSE2"].append(("__SSE2__", 1))
        # elif TARGET_CPU == "arm" or TARGET_CPU == "aarch64":
        #     if not self._simd_disabled["NEON"] and self._check_neon():
        #         cython_args["compile_time_env"]["NEON_BUILD_SUPPORT"] = True
        #         self._simd_supported["NEON"] = True
        #         self._simd_flags["NEON"].extend(self._neon_flags())
        #         self._simd_defines["NEON"].append(("__ARM_NEON__", 1))

        # add the platform sources as dependencies
        for ext in self.extensions:
            ext.depends.extend(itertools.chain.from_iterable(ext.platform_sources.values()))

        # cythonize the extensions (retaining platform-specific sources)
        platform_sources = [ext.platform_sources for ext in self.extensions]
        self.extensions = cythonize(self.extensions, **cython_args)
        for ext, plat_src in zip(self.extensions, platform_sources):
            ext.platform_sources = plat_src

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
                    # make the `Similarity::calculateMatrixIdentity` virtual
                    # so it we can override it with an SSE implementation
                    if line.strip() == b"void calculateMatrixIdentity();":
                        dst.write(b"virtual void calculateMatrixIdentity();\n")
                    # make the `Similarity::calculateVectors` virtual so it
                    # we can override it with an SSE implementation
                    elif line.strip() == b"bool calculateVectors(bool cutByGap = true);":
                        dst.write(b"virtual bool calculateVectors(bool cutByGap = true);\n")
                    # make the `Similarity::setSimilarityMatrix` virtual so it
                    # we can override it with an SSE implementation
                    elif line.strip() == b"bool setSimilarityMatrix(similarityMatrix * sm);":
                        dst.write(b"virtual bool setSimilarityMatrix(similarityMatrix * sm);\n")
                    # make the `Cleaner::calculateSeqIdentity` virtual so it
                    # we can override it with an SSE implementation
                    elif line.strip() == b"void calculateSeqIdentity();":
                        dst.write(b"virtual void calculateSeqIdentity();\n")
                    # add a virtual destructor to `Cleaner` so it can be
                    # subclassed safely
                    elif line.strip() == b"explicit Cleaner(Alignment *parent);":
                        dst.write(line)
                        dst.write(b"virtual ~Cleaner() {};\n")
                    # expose all attributes as public by adding a `public`
                    # qualifier right at the beginning of a class declaration
                    elif re.match(rb'\W*class\W*.*\W*\{', line):
                        dst.write(line)
                        dst.write(b"public:\n")
                    # expose all attributes as public by preventing any
                    # `private` qualifier
                    else:
                        dst.write(line.replace(b"private:", b"public:"))

    def _check_function(self, funcname, header, args="()"):
        _eprint('checking whether function', repr(funcname), 'is available', end="... ")

        self.mkpath(self.build_temp)
        base = "have_{}".format(funcname)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}{};
                    return 0;
                }}
            """.format(header, funcname, args))
        try:
            objects = self.compiler.compile([testfile], debug=self.debug)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
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
            if os.path.isfile(binfile):
                os.remove(binfile)

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
        # check for functions required for libcpu_features on OSX
        if SYSTEM == "Darwin":
            if self._check_function("sysctlbyname", "sys/sysctl.h", args="(NULL, NULL, 0, NULL, 0)"):
                self.compiler.define_macro("HAVE_SYSCTLBYNAME", 1)

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
        if library.language == "c++":
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-std=c++11")
                library.extra_link_args.append("-Wno-alloc-size-larger-than")
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
        # if library.name == "trimal":
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
        # else:
        #     sources = library.sources[:]

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
            re.sub(r'(.cpp|.c)$', self.compiler.obj_extension, s)
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
            "cpu_features",
            language="c",
            sources=[
                os.path.join("vendor", "cpu_features", "src", base)
                for base in [
                    "impl_{}_{}.c".format(TARGET_CPU, TARGET_SYSTEM),
                    "filesystem.c",
                    "stack_line_reader.c",
                    "string_view.c",
                    "copy.inl",
                    "define_introspection.inl",
                    "define_introspection_and_hwcaps.inl",
                    "equals.inl",
                    "impl_x86__base_implementation.inl",
                ]
            ],
            include_dirs=[os.path.join("vendor", "cpu_features", "include")],
            define_macros=[("STACK_LINE_READER_BUFFER_SIZE", 1024)]
        ),
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
                os.path.join("pytrimal", "fileobj", "pyfilebuf.cpp"),
                os.path.join("pytrimal", "patch", "reportsystem.cpp"),
                os.path.join("pytrimal", "_trimal.pyx"),
            ],
            platform_sources={
                "SSE2": [os.path.join("pytrimal", "impl", "sse.cpp")]
            },
            include_dirs=[
                os.path.join("pytrimal", "patch"),
                os.path.join("pytrimal", "fileobj"),
                os.path.join("pytrimal", "impl"),
                os.path.join("vendor", "cpu_features", "include"),
                "pytrimal",
                "include",
            ],
            libraries=[
                "cpu_features",
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
