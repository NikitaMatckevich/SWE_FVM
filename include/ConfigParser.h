#pragma once
#include <Includes.h>
#include <Exceptions.h>
#include <string>
#include <charconv>
#include <unordered_map>

class Parser {
  using ini_t = std::unordered_map<std::string, std::unordered_map<std::string, std::string>>;
  ini_t m_ini;
public:
  Parser(const std::string& filename);
  double get(const std::string& section, const std::string& property) const;	
};
