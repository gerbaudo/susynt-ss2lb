#ifndef SUSYTEST0_UTILS_H
#define SUSYTEST0_UTILS_H
/*
  Generic utility functions

  davide.gerbaudo@gmail.com
  Aug 2013
 */
#include <fstream>
#include <string>

// mkdir if it is not already there. Return dir path; return empty string if there was a problem
std::string mkdirIfNeeded(const std::string &dirname);

std::string basedir(const std::string &path);
bool endswith(const std::string &s, const std::string &end);
bool isInt(const std::string &s);
bool fileExists(const std::string &filename);

#endif
