#pragma once
# include "Includes.h"
#include <string>
#include <charconv>
#include <unordered_map>

struct ParserError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

class Parser {
  using ini_t = std::unordered_map<std::string, std::unordered_map<std::string, std::string>>;
  ini_t m_ini;
public:
  Parser(const std::string& filename);
  std::string get(const std::string& section, const std::string& property) const;
  template <typename numtype>
  numtype get(const std::string& section, const std::string& property) const {
    std::string strEntry = get(section, property);
    numtype numEntry;
    if (auto [p, ec] = std::from_chars(strEntry.data(), strEntry.data() + strEntry.size(), numEntry);
        ec == std::errc())
      return numEntry;
    else {
      switch (ec) {
      case std::errc::invalid_argument:
        throw ParserError("Incorrect data type provided when parsing " + section + " : " + property);
      case std::errc::result_out_of_range:
        throw ParserError("Parsed value " + section + " : " + property + " is out of range for given type");
      default:
        throw ParserError("Unexpected error when parsing " + section + " : " + property);
      }
    }
  }
};