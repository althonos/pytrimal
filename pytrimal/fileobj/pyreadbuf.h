#include <Python.h>
#include <streambuf>

// Wrapper class that exposes a Python file-like object for reading operations.
// Internally uses `read1` if possible, and `read` otherwise.
class pyreadbuf : public std::streambuf {
protected:
  PyObject *handle;
  PyObject *method;
  char buffer[1];
  PyObject *bufsize_py;
  std::streamsize bufsize;

public:
  pyreadbuf(PyObject *handle);
  ~pyreadbuf();

protected:
  int underflow() override;
  pyreadbuf *setbuf(char *s, std::streamsize n) override;
  std::streampos
  seekpos(std::streampos sp,
          std::ios_base::openmode which = std::ios_base::in |
                                          std::ios_base::out) override;
};
