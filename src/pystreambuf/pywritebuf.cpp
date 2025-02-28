#include "pywritebuf.h"

pywritebuf::pywritebuf(PyObject *handle) : std::streambuf(), handle(handle) {
  Py_INCREF(handle);
  method = PyUnicode_FromString("write");
  mview = PyMemoryView_FromMemory(buffer, 1, PyBUF_READ);
}

pywritebuf::~pywritebuf() {
  Py_DECREF(handle);
  Py_DECREF(method);
  Py_DECREF(mview);
}

int pywritebuf::overflow(int c) {
  PyObject *result;
  if (c != EOF) {
    buffer[0] = c;
    result = PyObject_CallMethodObjArgs(handle, method, mview, NULL);
    if (result == nullptr)
      return EOF;
    Py_DECREF(result);
  }
  return c;
}
