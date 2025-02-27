/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

/* *****************************************************************************

    This file contains a modified version of the trimAl report manager that
    emits exceptions and warnings with the Python C API. It requires functions
    that can raise exceptions to be declared as such in the Cython `.pxd`
    files.

***************************************************************************** */

#include <Python.h>

#include "InternalBenchmarker.h"
#include "reportsystem.h"

static PyObject *error_from_errorcode(ErrorCode code) {
  switch (code) {
  case UnknownCharacter:
  case UndefinedSymbol:
  case IncorrectSymbol:
    return PyExc_ValueError;
  default:
    return PyExc_RuntimeError;
  }
}

reporting::reportManager debug = reporting::reportManager();

void reporting::reportManager::PrintCodesAndMessages() {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void reporting::reportManager::PrintCodesAndMessages() ");
  switch (Level) {
  case VerboseLevel::NONE:
    std::cout << "[VerboseLevel] None" << std::endl;
    break;
  case VerboseLevel::INFO:
    std::cout << "[VerboseLevel] Info" << std::endl;
    break;
  case VerboseLevel::WARNING:
    std::cout << "[VerboseLevel] Warning" << std::endl;
    break;
  case VerboseLevel::ERROR:
    std::cout << "[VerboseLevel] Error" << std::endl;
    break;
  }

  for (int i = 1; i < InfoCode::__MAXINFO; i++) {
    report((InfoCode)i);
  }

  for (int i = 1; i < WarningCode::__MAXWARNING; i++) {
    report((WarningCode)i);
  }

  for (int i = 1; i < ErrorCode::__MAXERROR; i++) {
    report((ErrorCode)i);
  }
}

void reporting::reportManager::report(ErrorCode message, std::string *vars) {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void PythonReportManager::report(ErrorCode message, std::string "
              "*vars) ");

  std::string s(reporting::reportManager::ErrorMessages.at(message));

  if (vars != nullptr) {
    int counter = 0;
    std::size_t index;
    std::string FindWord = "[tag]";
    while ((index = s.find(FindWord)) != std::string::npos)
      s.replace(index, FindWord.length(), vars[counter++]);
    delete[] vars;
  }

  PyGILState_STATE state = PyGILState_Ensure();
  PyErr_SetString(error_from_errorcode(message), s.c_str());
  PyGILState_Release(state);
}

void reporting::reportManager::report(ErrorCode message, const char *vars) {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming(
      "void reporting::reportManager::report(ErrorCode message, char *vars) ");

  std::string s(reporting::reportManager::ErrorMessages.at(message));

  if (vars != nullptr) {
    std::size_t index;
    std::string FindWord = "[tag]";
    std::string Vars = vars;
    while ((index = s.find(FindWord)) != std::string::npos)
      s.replace(index, FindWord.length(), Vars);
  }

  PyGILState_STATE state = PyGILState_Ensure();
  PyErr_SetString(error_from_errorcode(message), s.c_str());
  PyGILState_Release(state);
}

void reporting::reportManager::report(WarningCode message, std::string *vars) {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void reporting::reportManager::report(WarningCode message, "
              "std::string *vars) ");

  std::string s(reporting::reportManager::WarningMessages.at(message));

  if (vars != nullptr) {
    int counter = 0;
    std::size_t index;
    std::string FindWord = "[tag]";
    while ((index = s.find(FindWord)) != std::string::npos)
      s.replace(index, FindWord.length(), vars[counter++]);
    delete[] vars;
  }

  PyGILState_STATE state = PyGILState_Ensure();
  PyErr_WarnEx(PyExc_RuntimeWarning, s.c_str(), 1);
  PyGILState_Release(state);
}

void reporting::reportManager::report(WarningCode message, const char *vars) {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void reporting::reportManager::report(WarningCode message, char "
              "*vars) ");

  std::string s(reporting::reportManager::WarningMessages.at(message));

  if (vars != nullptr) {
    std::size_t index;
    std::string FindWord = "[tag]";
    std::string Vars = vars;
    while ((index = s.find(FindWord)) != std::string::npos)
      s.replace(index, FindWord.length(), Vars);
  }

  PyGILState_STATE state = PyGILState_Ensure();
  PyErr_WarnEx(PyExc_RuntimeWarning, s.c_str(), 1);
  PyGILState_Release(state);
}

void reporting::reportManager::report(InfoCode message, std::string *vars) {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void reporting::reportManager::report(InfoCode message, "
              "std::string *vars) ");
  delete[] vars;
}

void reporting::reportManager::report(InfoCode message, const char *vars) {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming(
      "void reporting::reportManager::report(InfoCode message, char *vars) ");
}
