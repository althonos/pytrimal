#include <Python.h>
#include <streambuf>

// Wrapper class that exposes a Python file-like object for reading operations.
// Internally uses `readinto` for zero-copy read from Python to the C++ buffer.
class pyreadintobuf : public std::streambuf {
protected:
  PyObject *handle;
  PyObject *method;
  PyObject *mview;
  char buffer[1];
  std::streamsize bufsize;

public:
  pyreadintobuf(PyObject *handle);
  ~pyreadintobuf();

protected:
  int underflow() override;
  pyreadintobuf *setbuf(char *s, std::streamsize n) override;
  std::streampos
  seekpos(std::streampos sp,
          std::ios_base::openmode which = std::ios_base::in |
                                          std::ios_base::out) override;
};
