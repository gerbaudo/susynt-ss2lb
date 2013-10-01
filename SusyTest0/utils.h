#ifndef SUSYTEST0_UTILS_H
#define SUSYTEST0_UTILS_H
/*
  Generic utility functions

  davide.gerbaudo@gmail.com
  Aug 2013
 */
#include <fstream>
#include <string>
#include <vector>

bool dirExists(const std::string &dirname);
// mkdir if it is not already there. Return dir path; return empty string if there was a problem
std::string mkdirIfNeeded(const std::string &dirname);
std::string basedir(const std::string &path);
bool contains(const std::string &s, const std::string &sub);
bool endswith(const std::string &s, const std::string &end);
bool replace(std::string& str, const std::string& from, const std::string& to);
std::string rmLeadingTrailingWhitespaces(const std::string &str);
bool isInt(const std::string &s);
bool fileExists(const std::string &filename);
std::string getRootCoreDir(); //!< return empty string if env var not defined
std::string vdouble2str(const std::vector<double> &v);
std::string vfloat2str(const std::vector<float> &v);

#endif
