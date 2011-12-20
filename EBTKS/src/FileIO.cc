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
$RCSfile: FileIO.cc,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2006-03-02 13:22:17 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "FileIO.h"
#include <assert.h>
#include <stdlib.h>
// #include <stream.h> - (bert) removed unnecessary header
#include <sys/types.h>
#include <unistd.h>
#ifdef HAVE_MATLAB
  #include "Dictionary.h"
#endif

using namespace std;		// (bert) added
//
// Constructors/destructor
//

InputFile::InputFile()
{
  _ipipe = 0;
}

InputFile::operator void *() const
{
  if (_ipipe)
    return *_ipipe;
  return 0;
}

//
// Set functions
//

Boolean
InputFile::attach(const Path& path)
{
  close();

  Path fullPath(path.expanded().removeCompressedExtension());
  Path testPath(fullPath);

  Boolean compressed = FALSE;

  if (!testPath.exists()) {
    testPath = fullPath + ".gz";
    if (!testPath.exists()) {
      testPath = fullPath + ".z";
      if (!testPath.exists()) {
	testPath = fullPath + ".Z";
	if (!testPath.exists())
	  return FALSE;
      }
    }
    
    compressed = TRUE;

    MString tmpPath(256);
    assert(get_temp_filename((char *)tmpPath));

    int status = system("gunzip -c " + testPath + " > " + tmpPath);

    if (status)
      return FALSE;

    testPath = tmpPath;
  }

  _ipipe = new ifstream((const char *) testPath);

  if (compressed)
    unlink(testPath);

  if (_ipipe && *_ipipe)
    return TRUE;
  
  return FALSE;
}

istream&
InputFile::skip(unsigned nBytes)
{
  // Dummy read function, as seekg doesn't seem to work on ipopen

  if (!nBytes || !*_ipipe)
    return *_ipipe;

  char *dummy = new char[nBytes];
  assert(dummy);

  _ipipe->read(dummy, nBytes);

  delete [] dummy;

  return *_ipipe;
}

Boolean
InputFile::close()
{
  if (_ipipe) {
    //    _ipipe->exit();
    delete _ipipe;
    _ipipe = 0;
    return TRUE;
  }

  return FALSE;
}

/******************
 * OutputFile class
 ******************/

const Boolean OutputFile::NO_COMPRESS = 0;
const Boolean OutputFile::COMPRESS    = 1;

OutputFile::OutputFile(const Path& path, int mode, int compress)
: _path(path.expanded().removeCompressedExtension()),
  _compress(compress)
{
  Boolean isCompressed = FALSE;

  Path tempPath(_path);
  if (!tempPath.exists()) {
    tempPath = _path + ".gz";
    if (!tempPath.exists()) {
      tempPath = _path + ".z";
      if (!tempPath.exists())
	tempPath = _path + ".Z";
    }

    isCompressed = TRUE;
  }

  if (tempPath.exists()) {
    if (isCompressed && ((mode == ios::app) || (mode == ios::ate))) {
      MString command("gunzip " + tempPath);
      system(command);
    }

    if (mode == ios::out) {
      MString command("mv " + tempPath + " " + tempPath + "~");
      system(command);
    }
  }
  
  open(_path, openmode(mode));	// (bert) made this standard(?)
}

OutputFile::~OutputFile()
{
  Boolean good = this->good();

  close();

  if (good && (_compress == COMPRESS)) {
    MString command("gzip -f " + _path); // Compress fast
    system(command);
  }
}

//
// Set functions
//

void
OutputFile::compress(Boolean status)
{
  _compress = (status) ? COMPRESS : NO_COMPRESS;    
}

/***************************
 * C file handling functions
 ***************************/

FILE *
openFile(const Path& path, const char *mode)
{
  Path expandedPath(path.expanded());
  Path tempPath(expandedPath);

  Boolean isCompressed = FALSE;

  if (!tempPath.exists()) {
    tempPath = expandedPath + ".gz";
    if (!tempPath.exists()) {
      tempPath = expandedPath + ".z";
      if (!tempPath.exists())
	tempPath = expandedPath + ".Z";
    }

    isCompressed = TRUE;
  }

  if (tempPath.exists()) {
    if (isCompressed) {
      if (strcmp(mode, "a") == 0) {
	MString command("gunzip " + expandedPath);
	system(command);
	FILE *file = fopen(expandedPath, mode);
	return file;
      }
      else if (strcmp(mode, "r") == 0) {
	MString command("gunzip -c " + expandedPath);
	FILE *file = popen(command, mode);
	return file;
      }
    }

    if (strcmp(mode, "w") == 0) {
      MString command("mv " + tempPath + " " + tempPath + "~");
      system(command);
      FILE *file = fopen(tempPath, mode);
      return file;
    }
  }

  FILE *file = fopen(expandedPath, mode);
  return file;
}

void
closeFile(FILE *file)
{
  if (file != NULL)
    if (pclose(file) == -1)
      fclose(file);
}


#ifdef HAVE_MATLAB
static Dictionary<MATFile *, MString> matFileDict;

MATFile *
openMatlabFile(const Path& path, const char *mode)
{
  Path expandedPath(path.expanded());
  Path tempPath(expandedPath);

  Boolean isCompressed = FALSE;

  if (!tempPath.exists()) {
    tempPath = expandedPath + ".gz";
    if (!tempPath.exists()) {
      tempPath = expandedPath + ".z";
      if (!tempPath.exists())
	tempPath = expandedPath + ".Z";
    }

    isCompressed = TRUE;
  }

  if (tempPath.exists()) {
    if (isCompressed) {
      if (strcmp(mode, "r") == 0) {
	MString  tempName("/tmp/midas.");
	tempName += (int) getpid();
	MString command("gunzip -c " + expandedPath + "> " + tempName);
	MATFile *file = matOpen(tempName, (char *) mode);
	if (file)
	  matFileDict[file] = tempName;
	return file;
      }
    }

    if (strcmp(mode, "w") == 0) {
      MString command("mv " + tempPath + " " + tempPath + "~");
      system(command);
      MATFile *file = matOpen(tempPath, (char *) mode);
      return file;
    }
  }

  MATFile *file = matOpen(expandedPath, (char *) mode);
  if (file) {
    char **dir;
    int    nDirEntries = 0;
    dir = matGetDir(file, &nDirEntries);
    matClose(file); // Close and reopen to make sure matGetNextMatrix can be used

    if (!dir || (nDirEntries <= 0))
      file = 0;
    else
      file = matOpen(expandedPath, (char *) mode);

    if (dir)
      mxFree(dir);
  }

  return file;
}

void
closeMatlabFile(MATFile *file)
{
  if (!file)
    return;

  matClose(file);

  if (matFileDict.containsKey(file))
    unlink(matFileDict[file]);
}
#endif

int
get_temp_filename(char *pathname)
{
#ifdef HAVE_MKSTEMP
  char *tmp_dir = getenv("TMPDIR");
  if (tmp_dir == NULL) {
#ifdef P_tmpdir
    tmp_dir = P_tmpdir;
#else
    tmp_dir = "/tmp";
#endif
  }
  sprintf(pathname, "%s/EBTKSXXXXXX", tmp_dir);
  int fd = mkstemp(pathname);
  if (fd < 0) {
    return (0);
  }
  close(fd);
#else
  char *tmp_ptr = tempnam(NULL, "EBTKS");
  if (tmp_ptr == NULL) {
    return (0);			// Error
  }
  strcpy(pathname, tmp_ptr);
  free(tmp_ptr);		// Free result from tempnam()
#endif
  return (1);
}

//
// FORCED INSTANTIATIONS
//
// Needed for SGI, but won't compile under GCC!!
//
#if !defined(__GNUC__) && defined(__sgi)

#include "Array.h"
#include "SimpleArray.h"
#include "ValueMap.h"
#include "CachedArray.h"

// From Array.h
template Array<double>& Array<double>::absorb(Array<double>& array);
template class Array<Path>;
template class Array<LinearMap>;
template unsigned Array<float>::debug(Boolean);
template Array<float>& Array<float>::append(float);

template int *Array<int>::asCarray(int *) const;
template double *Array<double>::asCarray(double *) const;

// From SimpleArray.h
template SimpleArray<double> SimpleArray<double>::operator() (unsigned) const;
template SimpleArray<double> SimpleArray<double>::operator() (unsigned, unsigned) const;
template SimpleArray<double>& map(SimpleArray<double>&,const Array<LinearMap>&);
template int SimpleArray<int>::min(unsigned *) const;
template SimpleArray<int>::SimpleArray(int, double, int);
template SimpleArray<float>& SimpleArray<float>::operator/=(const SimpleArray<float> &);
template SimpleArray<float>& SimpleArray<float>::prune();
template SimpleArray<float> SimpleArray<float>::log() const;
template unsigned SimpleArray<char>::occurrencesOf(char, unsigned, unsigned) const;
template double SimpleArray<float>::sum() const;

// From CachedArray.h
template class CachedArray<float>;

#endif // __sgi

