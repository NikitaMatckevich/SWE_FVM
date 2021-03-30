#pragma once
#include <Includes.h>
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
  double get(const std::string& section, const std::string& property) const;	
};
