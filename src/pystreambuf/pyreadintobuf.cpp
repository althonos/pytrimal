#include "pyreadintobuf.h"

pyreadintobuf::pyreadintobuf(PyObject *handle)
    : std::streambuf(), handle(handle) {
  method = PyUnicode_FromString("readinto");
  mview = PyMemoryView_FromMemory(buffer, 1, PyBUF_READ);
  setbuf(buffer, 1);
}

pyreadintobuf::~pyreadintobuf() {
  Py_DECREF(handle);
  Py_DECREF(method);
  Py_DECREF(mview);
}

std::streampos pyreadintobuf::seekpos(std::streampos sp,
                                      std::ios_base::openmode which) {
  PyObject *n = PyObject_CallMethod(handle, "seek", "i", sp);
  if (n == nullptr)
    return std::streampos(std::streamoff(-1));

  long l = PyLong_AsLong(n);
  Py_DECREF(n);

  setg(eback(), eback(), eback());
  return std::streampos(std::streamoff(l));
}

pyreadintobuf *pyreadintobuf::setbuf(char *s, std::streamsize n) {
  setg(s, &s[n], &s[n]);
  bufsize = n;
  Py_DECREF(mview);
  mview = PyMemoryView_FromMemory(s, bufsize, PyBUF_WRITE);
  return this;
}

int pyreadintobuf::underflow() {
  PyObject *nread = PyObject_CallMethodObjArgs(handle, method, mview, NULL);
  if (nread == nullptr) {
    Py_DECREF(mview);
    return EOF;
  }

  long n = PyLong_AsLong(nread);
  int c = (n == 0) ? EOF : *eback();

  Py_DECREF(nread);

  setg(eback(), eback(), eback() + n);
  return c;
}
