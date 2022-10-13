#include "pyreadbuf.h"

pyreadbuf::pyreadbuf(PyObject *handle) : std::streambuf(), handle(handle) {
  Py_INCREF(handle);
  if (PyObject_HasAttrString(handle, "read1")) {
    method = PyUnicode_FromString("read1");
  } else {
    method = PyUnicode_FromString("read");
  }
  bufsize_py = PyLong_FromLong(1);
  setbuf(buffer, 1);
}

pyreadbuf::~pyreadbuf() {
  Py_DECREF(handle);
  Py_DECREF(method);
  Py_DECREF(bufsize_py);
}

std::streampos pyreadbuf::seekpos(std::streampos sp,
                                  std::ios_base::openmode which) {
  PyObject *n = PyObject_CallMethod(handle, "seek", "i", sp);
  if (n == nullptr)
    return std::streampos(std::streamoff(-1));

  long l = PyLong_AsLong(n);
  Py_DECREF(n);

  setg(eback(), eback(), eback());
  return std::streampos(std::streamoff(l));
}

pyreadbuf *pyreadbuf::setbuf(char *s, std::streamsize n) {
  setg(s, &s[n], &s[n]);
  bufsize = n;
  Py_DECREF(bufsize_py);
  bufsize_py = PyLong_FromLongLong(n);
  return this;
}

int pyreadbuf::underflow() {
  PyObject *bytes =
      PyObject_CallMethodObjArgs(handle, method, bufsize_py, NULL);
  if (bytes == nullptr) {
    return EOF;
  }
  if (!PyBytes_Check(bytes)) {
    Py_DECREF(bytes);
    PyErr_SetString(PyExc_TypeError, "a bytes-like object is required");
    return EOF;
  }
  Py_ssize_t size = PyBytes_Size(bytes);
  if (size <= 0) {
    // EOF
    Py_DECREF(bytes);
    return EOF;
  } else if (size > bufsize) {
    // error
    Py_DECREF(bytes);
    PyErr_SetString(PyExc_BufferError,
                    "more data returned by `read` than can fit in buffer");
    return EOF;
  } else {
    // OK: copy bytes to the internal buffer and reset pointers.
    memcpy(eback(), PyBytes_AS_STRING(bytes), size * sizeof(char));
    setg(eback(), eback(), eback() + size);
    return (size == 0) ? EOF : *eback();
  }
}
