# add_library(opal STATIC opal.h ScoreMatrix.cpp)
# target_include_directories(opal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Set common files for all executables
set(TRIMAL_OBJECTS
    $<TARGET_OBJECTS:CoreOBJLib>
    $<TARGET_OBJECTS:FormatsOBJLib>
    $<TARGET_OBJECTS:FormatHandlerOBJLib>
    $<TARGET_OBJECTS:ReportSystemOBJLib>
    $<TARGET_OBJECTS:StatisticOBJLib>
    $<TARGET_OBJECTS:UtilsOBJLib>
    $<TARGET_OBJECTS:InternalBenchmarkOBJLib>
)

# Scripts to detect various SIMD support for the current compiler
if (NOT DISABLE_AVX2)
    include("./scripts/CMake/FindAVX2.cmake")
endif()
if (NOT DISABLE_NEON)
    include("./scripts/CMake/FindNEON.cmake")
endif()
if (NOT DISABLE_SSE2)
    include("./scripts/CMake/FindSSE2.cmake")
endif()

# Add SIMD code if compiler supports it
if (HAVE_SSE2)
    list(APPEND TRIMAL_OBJECTS $<TARGET_OBJECTS:SSE2OBJLib>)
endif()
if (HAVE_AVX2)
    list(APPEND TRIMAL_OBJECTS $<TARGET_OBJECTS:AVX2OBJLib>)
endif()
if (HAVE_NEON)
    list(APPEND TRIMAL_OBJECTS $<TARGET_OBJECTS:NEONOBJLib>)
endif()

# Create a library with trimAl code
add_library(libtrimal source/trimalManager.cpp source/VCFHandler.cpp ${TRIMAL_OBJECTS})
target_include_directories(libtrimal PUBLIC include)

# Add `cpu_features` for detecting supported SIMD at compile time
add_subdirectory(vendor/cpu_features)
include_directories(vendor/cpu_features/include)
target_link_libraries(libtrimal cpu_features)

# Link the mathematical library to the targets
target_link_libraries(libtrimal m)

# Script that sets Release type to default and prints information on configuration
include("./scripts/CMake/BuildTypeWrapper.cmake")

# Script to handle static and dynamic compilation.
#   Static compilation is useful if compiling for containers without SO
include("./scripts/CMake/staticCompilationWrapper.cmake")

# Script to create intermediate OBJLib objects to include in targets
include("./scripts/CMake/OBJ-LIB-creator.cmake")
set_target_properties(CatchOBJLib PROPERTIES LINKER_LANGUAGE CXX)

# Script to configure the overwrite policy of the format manager:
#   - Overwrite and warn.
#   - Dont overwrite unless needed. If needed, overwrite original.
#   - Dont overwrite unless needed. If needed, overwrite last available.
include("./scripts/CMake/FormatHandlerOverwritePolicy.cmake")

# Script to create the Format Handler constructor, which incorporates
#   all found stats to the internal pool of the handler
include("./scripts/CMake/FormatHandlerHeaderCreator.cmake")

# Ensure `Python.h` can be imported while loading the report system,
# which we patched so it reports with the Python system for exceptions
# and warnings.
target_include_directories(ReportSystemOBJLib PUBLIC ${Python_INCLUDE_DIRS})
