#include "ConfigParser.h"
#include <fstream>

using namespace std;

namespace {
  void trim(string& s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](int ch) { return !isspace(ch); }));
    s.erase(find_if(s.rbegin(), s.rend(), [](int ch) { return !isspace(ch); }).base(), s.end());
  }
  const char assignment = '=';
  bool is_comment(string_view line) {
    return line.size() == 0 || line.at(0) == '#' || line.at(0) == ';';
  }
  bool is_section(string_view line) {
    return line.at(0) == '[' && line.back() == ']' && line.size() > 2;
  }
}

Parser::Parser(const string& filename) {
  ifstream file(filename);
  if (!file.good())
    throw ParserError("no file with a name " + filename + " found by parser");
  unsigned int ctr = 0;
  string line;
  string section;
  while (getline(file, line)) {
    ++ctr;      // increment the line counter
    trim(line); // remove whitespace from the line
    if (is_comment(line))
      continue; // ignore blank lines and comments entirely
    if (is_section(line))
      section = line.substr(1, line.size() - 2); // no emty section names allowed
    else {
      size_t pos = line.find_first_of(assignment);
      if (pos == string::npos)
        throw ParserError("non-section line " + to_string(ctr) + " in file " + filename
          + " does not contain an assignement symbol");
      string key(line.substr(0, pos));
      string val(line.substr(pos + 1));
      trim(key);
      trim(val);
      if (key.size() == 0 || val.size() == 0)
        throw ParserError("empty entries in file " + filename + " in line " + to_string(ctr));
      m_ini[section][key] = val;
    }
  }
}
string Parser::get(const string& section, const string& property) const {
  auto sectionEntry = m_ini.find(section);
  if (sectionEntry == m_ini.end())
    throw ParserError("no section with name" + section + " found while calling parser get function");
  auto properties = sectionEntry->second;
  auto propertyEntry = properties.find(property);
  if (propertyEntry == properties.end())
    throw ParserError("no property with name" + property + " found while calling parser get function");
  return propertyEntry->second;
}