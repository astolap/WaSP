/*BSD 2-Clause License
* Copyright(c) 2019, Pekka Astola
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
*
* 1. Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*     OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <sys/stat.h>
#include <string>
#include <experimental/filesystem>

#include "fileaux.hh"

using namespace std;
namespace fs = std::experimental::filesystem;

#define SYSTEM_VERBOSE_QUIET true

/*the aux_exists(), aux_ensure_directory() functions adapted from the ones
originally written by Pedro Garcia Freitas for JPEG Pleno activities */

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

