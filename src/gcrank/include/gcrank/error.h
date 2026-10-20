#pragma once
//
// Typed exception hierarchy.
//
// Every exception thrown by the gcrank library derives from gcrank::Error.
// The CLI tool maps these types onto its documented exit codes.
//
#include <stdexcept>
#include <string>

namespace gcrank {

// Base class of all gcrank errors.
class Error : public std::runtime_error {
public:
  explicit Error(const std::string& what) : std::runtime_error(what) {}
};

// Low-level I/O failure (open/read/write/mmap).
class IoError : public Error {
public:
  using Error::Error;
};

// Malformed or unreadable input (FASTA data, input streams).
class InputError : public IoError {
public:
  using IoError::IoError;
};

// Failure while writing or replacing an index on disk.
class OutputError : public IoError {
public:
  using IoError::IoError;
};

// Malformed, corrupt, or unsupported index data (manifest, contigs, bitvectors).
// This is an on-disk-format contract violation.
class FormatError : public Error {
public:
  using Error::Error;
};

// Invalid query arguments (bad range, zero denominator).
class QueryError : public Error {
public:
  using Error::Error;
};

// Invalid API/CLI usage (bad option, invalid argument value).
class UsageError : public Error {
public:
  using Error::Error;
};

}  // namespace gcrank
