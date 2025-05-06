file(COPY ../vendor/trimal DESTINATION .)
file(COPY_FILE ./trimal/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/trimal/CMakeLists.txt)
file(COPY_FILE ./trimal/source/reportsystem.cpp ${CMAKE_CURRENT_BINARY_DIR}/trimal/source/reportsystem.cpp)

set(TRIMAL_SOURCES
    "source/Cleaner.cpp"
    "source/FormatHandling/BaseFormatHandler.cpp"
    "include/Statistics/Manager.h"
    "include/Statistics/similarityMatrix.h"
    "include/trimalManager.h"
)

foreach(_file IN ITEMS ${TRIMAL_SOURCES})
    cmake_path(GET _file FILENAME _name)
    execute_process(
        COMMAND ${Python_EXECUTABLE} 
            ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
            --input ${PROJECT_SOURCE_DIR}/vendor/trimal/${_file}
            --patch ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
            --output ${CMAKE_CURRENT_BINARY_DIR}/trimal/${_file}
    )
endforeach()

add_subdirectory(pystreambuf)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/trimal ${CMAKE_CURRENT_BINARY_DIR}/trimal/build)
add_subdirectory(pytrimal)

