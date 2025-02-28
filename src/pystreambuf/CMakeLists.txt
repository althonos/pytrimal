set(
    PYSTREAMBUF_SOURCES
    pyreadbuf.cpp
    pyreadintobuf.cpp
    pywritebuf.cpp
    pyreadbuf.h 
    pyreadintobuf.h 
    pywritebuf.h
)

add_library(pystreambuf ${PYSTREAMBUF_SOURCES})
target_include_directories(pystreambuf PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(pystreambuf PUBLIC ${Python_INCLUDE_DIRS})