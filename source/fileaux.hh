/* fileaux.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef FILEAUX_HH
#define FILEAUX_HH

#ifdef __unix__
#define _popen popen
#define _pclose pclose
#endif

#include <string>

using namespace std;

/*the aux_exists(), aux_ensure_directory() functions adapted from the ones
originally written by Pedro Garcia Freitas for JPEG Pleno activies*/
bool aux_exists(const char* filename);
bool aux_exists(string filename);
void aux_ensure_directory(string filename);
void aux_ensure_directory(const char* filename);

int32_t system_1(char *str);

long aux_GetFileSize(const char* filename);
long aux_GetFileSize(char* filename);
long aux_GetFileSize(string filename);

#endif
