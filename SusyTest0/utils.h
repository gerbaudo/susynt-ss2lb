#ifndef SUSYTEST0_UTILS_H
#define SUSYTEST0_UTILS_H
/*
  Generic utility functions

  davide.gerbaudo@gmail.com
  Aug 2013
 */

#include <algorithm>    // std::set_intersection, remove_if
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


// Build a vector that is the difference between two vectors.
// Caveat1: computing the difference triggers copies of the vectors
// Caveat2: the result is not guaranteed to be sorted
// Caveat3: a-b != b-a
// Based on: http://stackoverflow.com/questions/14175858/c-subtract-vectors
template <typename T>
std::vector<T> subtract_vector(std::vector<T>& a, const std::vector<T>& b)
{
  std::vector<T> difference;
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter( difference ));
  return difference;
}

// Filter a container with a predicate
// Lifted from:
// http://stackoverflow.com/questions/2797142/higher-order-function-filter-in-c
template <typename C, typename P>
  C filter(C const & container, P pred) {
  C filtered(container);
  filtered.erase(remove_if(filtered.begin(), filtered.end(), pred), filtered.end());
  return filtered;
}


#endif
