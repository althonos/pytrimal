#include <Python.h>

#include "pyfilebuf.h"

pyfilebuf::pyfilebuf(PyObject* handle): std::streambuf(), handle(handle) {
    Py_INCREF(handle);
    method = PyUnicode_FromString("write");
    mview = PyMemoryView_FromMemory(buffer, 1, PyBUF_READ);
}

pyfilebuf::~pyfilebuf() {
    Py_DECREF(handle);
    Py_DECREF(method);
    Py_DECREF(mview);
}

int pyfilebuf::overflow(int c) {
    if (c != EOF) {
        buffer[0] = c;
        if (PyObject_CallMethodOneArg(handle, method, mview) == nullptr)
            return EOF;
    }
    return c;
}
