#pragma once
#include <stdexcept>

struct ParserError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct DimensionError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct MeshError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct SolverError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

#ifdef __GNUC__
/* GNU C specific stacktrace routine */

#include "stdio.h"
#include "execinfo.h"
#include "stdlib.h"
#include "unistd.h"

inline void PrintStacktrace() {
  constexpr int nb_calls = 10;
  void* calls[nb_calls];
  size_t size = backtrace(calls, nb_calls);
  backtrace_symbols_fd(calls, size, STDERR_FILENO);
  exit(1);
}

#endif
