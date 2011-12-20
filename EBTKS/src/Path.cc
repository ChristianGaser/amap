/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: Path.cc,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2003-04-16 18:42:39 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Path.h"
#include <assert.h>
#include <dirent.h>
#include <fstream>		// (bert) changed from fstream.h
using namespace std;		// (bert) added
#include <pwd.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

const char *Path::_separator = ".";
const char *Path::_imageNumberWildCard = "#";

//
// Constructors
//

Path::Path(const MString& dir, const MString& file) 
: MString(dir.length() + file.length() + 1)
{
  if ((dir.length() > 0) && (dir.lastChar() != '/') && (file.firstChar() != '/')) {
    strcpy((char *) string(), dir.string());
    (*this)[dir.length()] = '/';
    strcpy((char *) string() + dir.length() + 1, file.string());
  }
  else if ((dir.lastChar() == '/') && (file.firstChar() == '/')) {
    strcpy((char *) string(), dir.string());
    (*this)[dir.length() - 1] = '/';
    strcpy((char *) string() + dir.length(), file.string());
  }
  else {
    strcpy((char *) string(), dir.string());
    strcpy((char *) string() + dir.length(), file.string());
  }
}

//
// Get functions
//

MString *
Path::dir() const
{
  int lastSlashIndex = this->indexOf('/', -1);

  if (lastSlashIndex >= 0) {
    MString *dirString = new MString((*this)(0, lastSlashIndex));
    assert(dirString);
    return(dirString);
  }

  MString *dirString = new MString;
  assert(dirString);
  return(dirString);
}

MString *
Path::file() const
{
  int lastSlashIndex = this->indexOf('/', -1);
  
  if ((lastSlashIndex >= 0) && (lastSlashIndex < (int) this->length() - 1)) {
    MString *fileString = new MString((*this)(lastSlashIndex + 1));
    assert(fileString);
    return(fileString);
  }

  MString *fileString = new MString(string());
  assert(fileString);
  return(fileString);
}

//
// Operators
//

//
// Special functions
//

Path
Path::expanded() const
{
  Path newPath(*this);

  if (firstChar() != '~') {
    MString *dirString = dir();
    if (dirString->isEmpty() || ((*this)[(unsigned int)0] == '.')) {
      char *pwdPtr = getenv("PWD");
      if (pwdPtr == NULL) {
	cerr << "Couldn't get PWD environment variable!" << endl;
	assert(0);
      }
      MString pwd(pwdPtr);

      if ((*this)(0, 2) == MString("..")) { // Expand .. to full path (bert)
	int slashIndex = pwd.indexOf('/', -1);
	if (slashIndex >= 0)
	  newPath = pwd(0, slashIndex) + (*this)(2);
      }
      else if ((*this)[(unsigned int)0] == '.')               // Expand . to full path
	newPath = pwd + (*this)(1);
      else
	newPath = Path(pwd, *this);
    }
    delete dirString;
  }
  else {
    if ((*this)[(unsigned int)1] == '/') {                    // Expand ~ to full path
      char *home;
      if ((home = getenv("HOME")) != NULL)
	newPath = Path(MString(home) + (*this)(1));
      else {
	cerr << "Couldn't get HOME environment variable!" << endl;
	assert(0);
      }
    }
    else {                                      // Expand ~user to full path
      int slashIndex = this->indexOf('/');
      MString userName((*this)(1, slashIndex - 1));
      struct passwd *passWordEntry;
      if ((passWordEntry = getpwnam(userName.string())) != NULL)
	newPath = Path(MString(passWordEntry->pw_dir) + (*this)(slashIndex));
      else {
	cerr << "Couldn't get password entry!" << endl;
	assert(0);
      }
    }
  }

  return(newPath);
}

Boolean
Path::exists() const
{
  return (ifstream((const char *) *this)) ? TRUE : FALSE;

  /* Changed because the readdir function fails on Sun (truncates two chars
     from filenames (??)
  MString *dirString  = dir();
  MString *fileString = file();

  if (dirString->length() == 0)
    *dirString += _separator;

  DIR *dirPtr = opendir(dirString->string());
  if (!dirPtr) {
    delete dirString;
    delete fileString;
    return FALSE;
  }
  struct dirent *dirEntryPtr;
  Boolean pathFound = FALSE;

  while (((dirEntryPtr = readdir(dirPtr)) != NULL) && !pathFound)
    pathFound = (strcmp(dirEntryPtr->d_name, fileString->string()) == 0);

  delete dirString;
  delete fileString;
  closedir(dirPtr);

  return pathFound;

  */
}

Boolean
Path::existsCompressed(MString *extension) const
{
  MString tempExtension;

  Path tempPath(*this);
  if (!tempPath.exists()) {
    tempExtension = ".z";
    tempPath += tempExtension;
    if (!tempPath.exists()) {
      tempExtension = ".Z";
      tempPath = *this + tempExtension;
      if (!tempPath.exists()) {
	tempExtension = ".gz";
	tempPath = *this + tempExtension;
	if (!tempPath.exists())
	  return FALSE;
      }
    }
  }

  if (extension)
    *extension = tempExtension;

  return TRUE;
}

Boolean
Path::hasCompressedExtension(MString *extension) const
{
  unsigned nChars = this->length();

  if (nChars < 3)
    return FALSE;

  MString tempExtension(".Z");
  if ((*this)(nChars - 2) != tempExtension) {
    tempExtension = ".z";
    if ((*this)(nChars - 2) != tempExtension) {
      tempExtension = ".gz";
      if ((*this)(nChars - 3) != tempExtension)
	return FALSE;
    }
  }

  if (extension)
    *extension = tempExtension;
  return TRUE;
}

Path&
Path::removeCompressedExtension()
{
  MString lastChars((*this)(MAX(0, int(length() - 3))));

  if (lastChars.contains(".gz"))
    chop(3);
  else if (lastChars.contains(".z") || lastChars.contains(".Z"))
    chop(2);

  return *this;
}

Boolean
Path::removeExtension(MString *extension)
{
  int period = indexOf('.', -1); // Locate last period in path

  MString ext;

  if (period >= 0)
    ext = chop(length() - period);

  if (extension)
    *extension = ext;

  return (period >= 0) ? TRUE : FALSE;
}

Boolean
Path::isWritable() const
{
  Path fullPath(this->expanded());
  Boolean exists = fullPath.exists();

  FILE *test = fopen(fullPath, "a");
  if (test) {
    fclose(test);
    if (!exists)
      unlink(fullPath);
    return TRUE;
  }

  return FALSE;
}

void
Path::applyTemplate(const Path& templatePath, const char *separator)
{
  MString *dirString  = dir();
  MString *fileString = file();
  MString *templateDirString  = templatePath.dir();
  MString *templateFileString = templatePath.file();

  Path *newPath = 
    new Path(dirString->applyTemplate(*templateDirString, separator),
	     fileString->applyTemplate(*templateFileString, separator));
  assert(newPath);

  delete dirString;
  delete fileString;
  delete templateDirString;
  delete templateFileString;
  
  *this = *newPath;
}

Boolean
Path::imageNumber(unsigned *value) const
{
  MString *fileString = file();

  if (!fileString)
    return FALSE;

  MStringIterator fileNameIt(*fileString, _separator);

  MString token(fileNameIt++);

  while (!token.isEmpty()) {
    if (token.isInteger((int *) value)) {
      delete fileString;
      return TRUE;
    }
    token = fileNameIt++;
  }

  delete fileString;
  return FALSE;
}

Boolean
Path::replaceImageNumber(unsigned newValue)
{
  char newValueText[500];
  sprintf(newValueText, "%d", newValue);

  return replaceImageNumber(newValueText);
}

Boolean
Path::replaceImageNumber(const char *newString)
{
  MString *fileString = file();

  if (!fileString)
    return FALSE;

  MStringIterator fileNameIt(*fileString, _separator);

  MString token(fileNameIt++);
  MString newFile(token);
  token = fileNameIt++;

  Boolean replaced = FALSE;

  while (!token.isEmpty()) {
    newFile += _separator;
    int value;
    if (!replaced && token.isInteger(&value)) {
      newFile += newString;
      replaced = TRUE;
    }
    else
      newFile += token;
    token = fileNameIt++;
  }

  if (replaced) {
    MString *dirString = dir();

    *this = Path(*dirString, newFile);

    delete dirString;
  }

  delete fileString;

  return replaced;
}

Boolean
Path::replaceImageNumberWildCard(unsigned newValue)
{
  char newValueText[500];
  sprintf(newValueText, "%d", newValue);

  return replaceImageNumberWildCard(newValueText);
}

Boolean
Path::replaceImageNumberWildCard(const char *newString)
{
  MString *fileString = file();

  if (!fileString)
    return FALSE;

  MStringIterator fileNameIt(*fileString, _separator);

  MString token(fileNameIt++);
  MString newFile(token);
  token = fileNameIt++;

  Boolean replaced = FALSE;

  while (!token.isEmpty()) {
    newFile += _separator;
    if (!replaced && (token == _imageNumberWildCard)) {
      newFile += newString;
      replaced = TRUE;
    }
    else
      newFile += token;
    token = fileNameIt++;
  }

  if (replaced) {
    MString *dirString = dir();

    *this = Path(*dirString, newFile);

    delete dirString;
  }

  delete fileString;

  return replaced;
}

//
// Private functions
//

Boolean
Path::_locateStringAfterSeparator(const char *string) const
{
  MString *fileString = file();

  char *filePtr = (char *) fileString->string();

  if (!strtok(filePtr, _separator)) 
    return FALSE;

  char *token = strtok(NULL, _separator);
  while (token) {
    if (strcmp(token, string) == 0) {
      delete fileString;
      return TRUE;
    }
    token = strtok(NULL, _separator);
  }
  
  delete fileString;
  return FALSE;
}
