/* fileaux.cpp */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/


#include <sys/stat.h>
#include <string>
#include <experimental/filesystem>

#include "fileaux.hh"

using namespace std;
namespace fs = std::experimental::filesystem;

#define SYSTEM_VERBOSE_QUIET true

/*the aux_exists(), aux_ensure_directory() functions adapted from the ones
originally written by Pedro Garcia Freitas for JPEG Pleno activies*/

bool aux_exists(const char* filename) {
    fs::path file_path(filename);
    return fs::exists(file_path);
}

bool aux_exists(string filename) {
    return aux_exists(filename.c_str());
}

void aux_ensure_directory(string filename) {
    fs::path file_path(filename);
    fs::path parent_path = file_path.parent_path();
    if (!aux_exists(parent_path.u8string()))
        fs::create_directories(parent_path);
}

void aux_ensure_directory(const char* filename) {
    string fname(filename);
    aux_ensure_directory(fname);
}

int32_t system_1(char *str) {
  string sys_call_str(str);
if (!SYSTEM_VERBOSE_QUIET)
  printf("System command: %s\n", str);
#ifdef _WIN32
  sys_call_str.append(" > nul");
#endif
#ifdef __unix__
  sys_call_str.append(" > /dev/null");
#endif
  return system(sys_call_str.c_str());
}


long aux_GetFileSize(string filename) {
  struct stat stat_buf;
  int32_t rc = stat(filename.c_str(), &stat_buf);
  return rc == 0 ? stat_buf.st_size : -1;
}


long aux_GetFileSize(const char* filename) {
  string sbuffer(filename);
  return aux_GetFileSize(sbuffer);
}


long aux_GetFileSize(char* filename) {
  string sbuffer(filename);
  return aux_GetFileSize(sbuffer);
}

