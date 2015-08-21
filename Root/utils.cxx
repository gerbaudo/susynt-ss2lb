#include "susynt-ss3l/utils.h"

#include <algorithm>    // std::copy
#include <cstdlib> // strtol
#include <iterator>     // std::ostream_iterator
#include <iostream>
#include <sstream>      // std::ostringstream
#include <sys/stat.h>
#include <sys/types.h>

namespace ss3l
{
//----------------------------------------------------------
std::string basedir(const std::string &path)
{
  size_t found = path.find_last_of('/');
  return path.substr(0, found+1);
}
//----------------------------------------------------------
bool dirExists(const std::string &dirname)
{
  bool doesExist(false);
  if(dirname.length()<1) return doesExist;
  typedef struct stat Stat; Stat st;
  doesExist = (0==stat(dirname.c_str(), &st));
  bool isDir(S_ISDIR(st.st_mode));
  return doesExist && isDir;
}
//----------------------------------------------------------
// from : http://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux
std::string mkdirIfNeeded(const std::string &dirname)
{
  std::string result;
  if(dirname.length()<1)      result = "";
  else if(dirExists(dirname)) result = dirname;
  else {
    int status(mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH));
    bool success(status==0);
    result = (success ? dirname : "");
  }
  return result;
}
//----------------------------------------------------------
bool contains(const std::string &s, const std::string &sub)
{
  return (s.find(sub) != std::string::npos);
}
//----------------------------------------------------------
// http://stackoverflow.com/questions/874134/find-if-string-endswith-another-string-in-c
bool endswith(const std::string &s, const std::string &end) {
  if(s.length()<end.length()) return false;
  else return (0==s.compare(s.length() - end.length(), end.length(), end));
}
//----------------------------------------------------------
// http://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
bool replace(std::string& str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos) return false;
  str.replace(start_pos, from.length(), to);
  return true;
}
//----------------------------------------------------------
std::string rmLeadingTrailingWhitespaces(const std::string &str)
{
  using std::string;
  size_t startpos = str.find_first_not_of(" \t");
  size_t endpos = str.find_last_not_of(" \t");
  if(( string::npos == startpos ) || ( string::npos == endpos)) return string("");
  else return str.substr(startpos, endpos-startpos+1);
}
//----------------------------------------------------------
// http://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
bool isInt(const std::string &s)
{
  std::string rs(rmLeadingTrailingWhitespaces(s));
  if(rs.empty() || ((!isdigit(rs[0])) && (rs[0] != '-') && (rs[0] != '+'))) return false ;
  char * p ;
  strtol(rs.c_str(), &p, 10) ;
  return (*p == 0) ;
}
//----------------------------------------------------------
bool fileExists(const std::string &filename)
{
  std::ifstream file(filename.c_str());
  bool doesExists = file;
  file.close();
  return doesExists;
}
//----------------------------------------------------------
std::string getRootCoreDir()
{
  using namespace std;
  string dir;
  char* rootcoredir = getenv("ROOTCOREDIR");
  bool envvarDefined(rootcoredir!=0);
  if (envvarDefined) { dir = rootcoredir; }
  else               { cout<<"getRootCoreDir: ROOTCOREDIR not defined"<<endl; }
  return dir;
}
//----------------------------------------------------------
std::string vdouble2str(const std::vector<double> &v)
{
  std::ostringstream oss;
  std::ostream_iterator<double> it (oss,", ");
  std::copy(v.begin(), v.end(), it);
  return oss.str();
}
//----------------------------------------------------------
std::string vfloat2str(const std::vector<float> &v)
{
  std::ostringstream oss;
  std::ostream_iterator<float> it (oss,", ");
  std::copy(v.begin(), v.end(), it);
  return oss.str();
}
//----------------------------------------------------------
} // ss3l
