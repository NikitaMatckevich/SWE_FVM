#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <exception>
#include <algorithm>
#include <cctype>

// removes whitespace from either side of a string
static inline void trim(std::string & s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) { return !std::isspace(ch); }));
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); }).base(), s.end());
}

class Parser {
  using ini_t = std::unordered_map<std::string, std::unordered_map<std::string, std::string>>;
private:
  ini_t m_ini;
public:
  Parser(std::string const & filename) {
    try {
      std::ifstream file(filename.c_str());
      if (!file.good()) {
        std::ostringstream message;
        message << "File '" << filename << "' does not exist";
        throw std::runtime_error(message.str());
      }
      m_ini.clear();
      unsigned int lineCtr = 0;
      std::string line;
      std::string section = "";
      const char assignment = '=', comment1 = '#', comment2 = ';';
      while (std::getline(file, line)) {
        ++lineCtr; // increment the line counter
        trim(line); // remove whitespace from the line
        // ignore blank lines and comments entirely
        if (line.size() == 0 || line.at(0) == comment1 || line.at(0) == comment2) continue;
        if (line.at(0) == '[') {
          if (line.back() == ']' && line.size() > 2) {
            section = line.substr(1, line.size() - 2);
          }
          else {
            std::ostringstream message;
            message << "Invalid Section Declaration on Line " << lineCtr;
            throw std::runtime_error(message.str());
          }
        }
        else {
          size_t pos = line.find_first_of(assignment);
          if (pos != std::string::npos) {
            std::string key(line.substr(0, pos));
            std::string note(line.substr(pos + 1));
            trim(key);
            trim(note);
            if ((key.size() == 0) || (note.size() == 0)) {
              std::ostringstream message;
              message << "Neither the Key nor Value can be Empty on Line " << lineCtr;
              throw std::runtime_error(message.str());
            }
            if (section.size() == 0) {
              throw std::runtime_error("No Section to Apply Settings");
            }
            m_ini[section][key] = note;
          }
          else {
            std::ostringstream message;
            message << "Unknown Settings on Line " << lineCtr;
            throw std::runtime_error(message.str());
          }
        }
      }
    }
    catch (std::runtime_error const & e) {
      std::cerr << e.what() << std::endl;
    }
  }
  std::string get(std::string const & section, std::string const & property) const {
    try {
    auto sectionEntry = m_ini.find(section);
    if (sectionEntry == m_ini.end()) {
      std::ostringstream message;
      message << "Section '" << section << "' Not Found";
      throw std::runtime_error(message.str());
    }
    auto properties = sectionEntry->second;
    auto propertyEntry = properties.find(property);
    if (propertyEntry == properties.end()) {
      std::ostringstream message;
      message << "Property '" << property << "' Not Found";
      throw std::runtime_error(message.str());
    }
    return propertyEntry->second;
    }
    catch (std::runtime_error const & e) {
      std::cerr << e.what() << std::endl;
    }
  }
  template <typename T> T get(std::string const & section, std::string const & property) const {
    try {
    T value;
    auto note = get(section, property);
    std::istringstream define_value(note);
      define_value >> value;
      if (define_value.fail()) {
        std::ostringstream message;
        message << "Invalid Value at '" << section << ":" << property << "'";
        throw std::runtime_error(message.str());
      }
      return value;
    }
    catch (std::runtime_error const & e) {
      std::cerr << e.what() << std::endl;
    }
  }
};