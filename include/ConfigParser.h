#pragma once
#include <Exceptions.h>
#include <Includes.h>
#include <charconv>
#include <string>
#include <unordered_map>

struct Parser {

  using HashMap = std::unordered_map<std::string, std::unordered_map<std::string, std::string>>;
  Parser(const std::string& filename);
  double Get(const std::string& section, const std::string& property) const;	

 private:
  HashMap m_ini; 
};
