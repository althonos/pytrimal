#include "pyreadintobuf.h"

pyreadintobuf::pyreadintobuf(PyObject* handle):
    std::streambuf(),
    handle(handle)
{
    method = PyUnicode_FromString("readinto");
    mview = PyMemoryView_FromMemory(buffer, 1, PyBUF_READ);
    setbuf(buffer, 1);
}

pyreadintobuf::~pyreadintobuf() {
    Py_DECREF(handle);
    Py_DECREF(method);
    Py_DECREF(mview);
}

std::streampos pyreadintobuf::seekpos(std::streampos sp, std::ios_base::openmode which) {
    PyObject* n = PyObject_CallMethod(handle, "seek", "i", sp);
    if (n == nullptr)
        return std::streampos(std::streamoff(-1));

    long l = PyLong_AsLong(n);
    Py_DECREF(n);

    setg(eback(), eback(), eback());
    return std::streampos(std::streamoff(l));
}

pyreadintobuf* pyreadintobuf::setbuf(char* s, std::streamsize n) {
    setg(s, &s[n], &s[n]);
    bufsize = n;
    return this;
}

int pyreadintobuf::underflow() {
    PyObject* mview = PyMemoryView_FromMemory(eback(), bufsize, PyBUF_WRITE);
    if (mview == nullptr) {
        return EOF;
    }

    PyObject* nread = PyObject_CallMethodObjArgs(handle, method, mview, NULL);
    if (nread == nullptr) {
        Py_DECREF(mview);
        return EOF;
    }

    long n = PyLong_AsLong(nread);
    char c = (n == 0) ? EOF : *eback();

    Py_DECREF(mview);
    Py_DECREF(nread);

    setg(eback(), eback(), eback() + n);
    return c;
}
