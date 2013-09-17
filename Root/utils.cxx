#include "SusyTest0/utils.h"

#include <sys/stat.h>
#include <sys/types.h>

// from : http://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux
std::string mkdirIfNeeded(const std::string &dirname)
{
  using std::string;
  if(dirname.length()<1) return string("");
  typedef struct stat Stat;
  Stat st;
  int status;
  bool doesnotExist(stat(dirname.c_str(), &st) != 0);
  if(doesnotExist) {
    status = mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    bool success(status==0);
    return success ? dirname : string("");
  } else {
    bool isDir(S_ISDIR(st.st_mode));
    return isDir ? dirname : string("");
  }
}

