#include <streambuf>
#include <Python.h>

class pyfilebuf: public std::streambuf {
private:
    PyObject* handle;
    PyObject* method;
    PyObject* mview;
    char      buffer[1];
public:
    pyfilebuf(PyObject* handle);
    ~pyfilebuf();
protected:
    int overflow(int c = EOF) override;
};
