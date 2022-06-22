#include <streambuf>
#include <Python.h>

// Wrapper class that exposes a Python file-like object for writing operations.
class pywritebuf: public std::streambuf {
private:
    PyObject* handle;
    PyObject* method;
    PyObject* mview;
    char      buffer[1];
public:
    pywritebuf(PyObject* handle);
    ~pywritebuf();
protected:
    int overflow(int c) override;
};
