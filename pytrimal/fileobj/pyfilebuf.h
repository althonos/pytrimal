#include <streambuf>
#include <Python.h>

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

class pyreadbuf: public std::streambuf {
private:
    PyObject* handle;
    PyObject* method;
    PyObject* mview;
    char      buffer[1];

    std::streamsize bufsize;


public:
    pyreadbuf(PyObject* handle);
    ~pyreadbuf();
protected:
    // int uflow() override;
    int underflow() override;
    pyreadbuf* setbuf(char* s, std::streamsize n) override;
    std::streampos seekpos(std::streampos sp, std::ios_base::openmode which = std::ios_base::in | std::ios_base::out) override;
    // void bump();
    // int underflow() override;
};
