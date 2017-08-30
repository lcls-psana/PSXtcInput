//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: RandomAccess.cpp 7696 2017-08-18 00:40:59Z eslaught@SLAC.STANFORD.EDU $
//
// Description:
//	Class RandomAccess...
//
// Author List:
//      Elliott Slaughter
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "PSXtcInput/RandomAccess.h"

//-----------------
// C/C++ Headers --
//-----------------
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <string>
#include <iomanip>
#include <fcntl.h>
#include <stdlib.h>
#include <sstream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "MsgLogger/MsgLogger.h"
#include "PSXtcInput/Exceptions.h"
#include "pdsdata/index/IndexFileStruct.hh"
#include "pdsdata/index/IndexFileReader.hh"
#include "pdsdata/index/IndexList.hh"
#include "pdsdata/xtc/Sequence.hh"
#include "IData/Dataset.h"
#include "XtcInput/XtcFileName.h"
#include "pdsdata/xtc/XtcIterator.hh"
#include "pdsdata/psddl/epics.ddl.h"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

using namespace XtcInput;
using namespace std;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

namespace {
  const char* logger = "PSXtcInput::RandomAccess";
}

namespace PSXtcInput {

class RunMap {
public:
  std::vector<XtcInput::XtcFileName> files;
  typedef std::map<unsigned, std::vector<XtcInput::XtcFileName> > map;
  std::vector<unsigned> runs;
  map runFiles;

  RunMap(std::vector<std::string> &m_fileNames) {
    // Input can be a mixture of files and datasets.
    // Live mode is not supported. "one-stream mode"
    // is only supported if the users provides a list of
    // timestamps from one stream.

    typedef std::vector<std::string> FileList;

    // guess whether we have datasets or pure file names (or mixture)
    for (FileList::const_iterator it = m_fileNames.begin(); it != m_fileNames.end(); ++ it) {

      IData::Dataset ds(*it);
      if (ds.exists("live")) MsgLog(logger, fatal, "Live mode not supported with xtc indexing");

      if (ds.isFile()) {

        // must be file name
        files.push_back(XtcInput::XtcFileName(*it));

      } else {

        // Find files on disk and add to the list
        const IData::Dataset::NameList& strfiles = ds.files();
        if (strfiles.empty()) MsgLog(logger, fatal, "Empty file list");
        for (IData::Dataset::NameList::const_iterator it = strfiles.begin(); it != strfiles.end(); ++ it) {
          XtcInput::XtcFileName file(*it);
          files.push_back(file);
        }
      }
      // sort files to make sure we get a chunk0 first
      sort(files.begin(),files.end());

      // sort all files according run
      for (std::vector<XtcInput::XtcFileName>::const_iterator it = files.begin(); it != files.end(); ++ it) {
        runFiles[it->run()].push_back(*it);
      }
      for (map::const_iterator it = runFiles.begin(); it != runFiles.end(); ++ it) {
        runs.push_back(it->first);
      }
    }
  }
};

// class which manages xtc files, including "jump" function to do random access

class RandomAccessXtcReader {
public:
  RandomAccessXtcReader() {}

  bool has(const std::string& filename) {
    return _fd.count(filename);
  }

  int ensure(const std::string& filename) {
    if (has(filename)) return _fd[filename];

    int fd = ::open(filename.c_str(), O_RDONLY | O_LARGEFILE);
    if (fd==-1) MsgLog(logger, fatal,
                               "File " << filename.c_str() << " not found");
    _fd[filename] = fd;
    return fd;
  }

  ~RandomAccessXtcReader() {
    for  (std::map<std::string, int>::const_iterator it = _fd.begin(); it!= _fd.end(); it++)
      ::close(it->second);
  }

  Pds::Dgram* jump(const std::string& filename, int64_t offset) {
    int fd = ensure(filename);
    int64_t found = lseek64(fd,offset, SEEK_SET);
    if (found != offset) {
      stringstream ss;
      ss << "Jump to offset " << offset << " failed";
      MsgLog(logger, error, ss.str());
      throw RandomAccessSeekFailed(ERR_LOC);
    }
    Pds::Dgram dghdr;
    if (::read(fd, &dghdr, sizeof(dghdr))==0) {
      return 0;
    } else {
      if (dghdr.xtc.sizeofPayload()>MaxDgramSize)
        MsgLog(logger, fatal, "Datagram size exceeds sanity check. Size: " << dghdr.xtc.sizeofPayload() << " Limit: " << MaxDgramSize);
      Pds::Dgram* dg = (Pds::Dgram*)new char[sizeof(dghdr)+dghdr.xtc.sizeofPayload()];
      *dg = dghdr;
      ::read(fd, dg->xtc.payload(), dg->xtc.sizeofPayload());
      return dg;
    }
  }

private:
  enum {MaxDgramSize=0x2000000};
  std::map<std::string, int> _fd;
};

// this is the implementation of the per-run indexing.  shouldn't be too
// hard to make it work for for per-calibcycle indexing as well.

class RandomAccessRun {
private:
  // add a datagram with "event" data (versus nonEvent data, like epics)
  // to the vector of pieces (i.e. add another "piece")
  void _add(Pds::Dgram* dg, const std::string& filename) {
    _pieces.eventDg.push_back(XtcInput::Dgram(XtcInput::Dgram::make_ptr(dg),XtcFileName(filename)));
  }

  // copy the event-pieces onto the queue where the DgramSourceIndex object
  // can pick them up.
  void _post() {
    _queue.push(_pieces);
  }

  // add only one "event" datagram and post
  void _post(Pds::Dgram* dg, const std::string& filename) {
    _add(dg, filename);
    _post();
  }

  // post only this dg
  void _postOneDg(Pds::Dgram* dg, const std::string& filename) {
    _pieces.reset();
    if (dg) _post(dg, filename);
  }

  // look for configure in first 2 datagrams from the first file.  this will fail
  // if we don't get a chunk0 first in the list of files.  we have previously
  // sorted the files in RunMap to ensure this is the case.
  void _configure(const std::string &filename) {

    _pieces.reset();

    int64_t offset = 0;
    for (int i=0; i<2; i++) {
      Pds::Dgram* dg = _xtc.jump(filename, offset);
      if (dg->seq.service()==Pds::TransitionId::Configure) {
        _post(dg,filename);
        _beginrunOffset = dg->xtc.sizeofPayload()+sizeof(Pds::Dgram);
        return;
      }
      offset+=dg->xtc.sizeofPayload()+sizeof(Pds::Dgram);
    }
    MsgLog(logger, fatal, "Configure transition not found in first 2 datagrams");
  }

  // send beginrun from the first file
  void _beginrun(const std::string &filename) {
    Pds::Dgram* dg = _xtc.jump(filename, _beginrunOffset);
    if (dg->seq.service()!=Pds::TransitionId::BeginRun)
      MsgLog(logger, fatal, "BeginRun transition not found after configure transition");
    _postOneDg(dg,filename);
  }

public:

  RandomAccessRun(queue<DgramPieces>& queue, const vector<XtcFileName> &xtclist) :
    _xtc(), _queue(queue) {

    const std::string &filename = xtclist[0].path();

    _configure(filename);
    // send a beginrun transition
    _beginrun(filename);
  }

  ~RandomAccessRun() {}

  // jump to an event
  // can't be a const method because it changes the "pieces" object
  int jump(const std::vector<std::string>& filenames, const std::vector<int64_t> &offsets, const std::string &lastBeginCalibCycleDgram) {
    _pieces.reset();
    if (_beginCalibCycleDgram != lastBeginCalibCycleDgram) {
      _beginCalibCycleDgram = lastBeginCalibCycleDgram;

      const Pds::Dgram* dghdr = (const Pds::Dgram*)lastBeginCalibCycleDgram.c_str();
      Pds::Dgram* dg = (Pds::Dgram*)new char[sizeof(*dghdr)+dghdr->xtc.sizeofPayload()];
      memcpy(dg, dghdr, sizeof(*dghdr)+dghdr->xtc.sizeofPayload());
      _add(dg, "");
      _post();
    }

    _pieces.reset();

    bool accept = false;
    assert(filenames.size() == offsets.size());
    for (size_t i = 0; i < filenames.size(); i++) {
      const std::string &filename = filenames[i];
      int64_t offset = offsets[i];

      Pds::Dgram* dg=0;
      dg = _xtc.jump(filename, offset);
      _add(dg, filename);
      accept = accept || dg->seq.service()==Pds::TransitionId::L1Accept;
    }

    _post();
    return !accept; // zero means success
  }

private:
  RandomAccessXtcReader    _xtc;
  int64_t                  _beginrunOffset;
  std::string              _beginCalibCycleDgram;
  DgramPieces              _pieces;
  queue<DgramPieces>&      _queue;
};

// above is the "private" implementation (class RandomAccessRun), below this is the
// "public" implementation (class RandomAccess)

RandomAccess::RandomAccess(const std::string& name, std::queue<DgramPieces>& queue) : Configurable(name), _queue(queue),_raxrun(0),_run(-1) {
  _fileNames = configList("files");
  if ( _fileNames.empty() ) MsgLog(logger, fatal, "Empty file list");
  _rmap = new RunMap(_fileNames);
}

RandomAccess::~RandomAccess() {
  delete _raxrun;
  delete _rmap;
}

int RandomAccess::jump(const std::vector<std::string>& filenames, const std::vector<int64_t> &offsets, const std::string &lastBeginCalibCycleDgram) {
  return _raxrun->jump(filenames, offsets, lastBeginCalibCycleDgram);
}

void RandomAccess::setrun(int run) {
  // we can be called twice for the same run, because
  // at beginJob we "prefetch" the first configure transition
  // and then we will get another setrun from the run iterator
  if (run==_run) return;
  _run=run;
  delete _raxrun;
  _raxrun = new RandomAccessRun(_queue,_rmap->runFiles[run]);
}

const std::vector<unsigned>& RandomAccess::runs() {
  return _rmap->runs;
}

} // namespace PSXtcInput
