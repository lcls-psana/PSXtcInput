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

#include <pthread.h>
#include <legion.h>
#include <legion_c.h>
#include <legion_c_util.h>
#include <unistd.h>

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

  static bool has(const std::string& filename) {
    // WARNING: Only call this while holding _fd_mutex
    return _fd.count(filename);
  }

  static int ensure(const std::string& filename) {
    if (!pthread_mutex_lock(&_fd_mutex)) {
      assert(false && "pthread_mutex_lock failed\n");
    }

    if (has(filename)) {
      int result = _fd[filename];

      if (!pthread_mutex_unlock(&_fd_mutex)) {
        assert(false && "pthread_mutex_unlock failed\n");
      }

      return result;
    }

    int fd = ::open(filename.c_str(), O_RDONLY | O_LARGEFILE);
    if (fd==-1) MsgLog(logger, fatal,
                               "File " << filename.c_str() << " not found");
    _fd[filename] = fd;

    if (!pthread_mutex_unlock(&_fd_mutex)) {
      assert(false && "pthread_mutex_unlock failed\n");
    }

    return fd;
  }

  ~RandomAccessXtcReader() {
    // FIXME: Not safe to close fds in the destructor since the memory is now static
    // for  (std::map<std::string, int>::const_iterator it = _fd.begin(); it!= _fd.end(); it++)
    //   ::close(it->second);
  }

  static bool jump_internal(const std::string& filename, int64_t offset, Pds::Dgram* dg) {
    int fd = ensure(filename);
    if (::pread(fd, dg, sizeof(Pds::Dgram), offset)==0) {
      return false;
    } else {
      if (dg->xtc.sizeofPayload()>MaxDgramSize)
        MsgLog(logger, fatal, "Datagram size exceeds sanity check. Size: " << dg->xtc.sizeofPayload() << " Limit: " << MaxDgramSize);
      ::pread(fd, dg->xtc.payload(), dg->xtc.sizeofPayload(), offset+sizeof(Pds::Dgram));
      return true;
    }
  }

  static bool jump_task(const Legion::Task *task,
            const std::vector<Legion::PhysicalRegion> &regions,
            Legion::Context ctx, Legion::HighLevelRuntime *runtime) {
    // Unpack arguments.
    assert(task->arglen >= sizeof(Args));
    const Args &args = *(const Args *)task->args;
    const char *filename_start = (const char *)(((const char *)task->args) + sizeof(Args));
    std::string filename(filename_start, args.filename_size);
    int64_t offset = args.offset;

    // Fetch destination pointer out of region argument.
    LegionRuntime::Accessor::RegionAccessor<LegionRuntime::Accessor::AccessorType::Generic> accessor =
      regions[0].get_field_accessor(args.fid);
    LegionRuntime::Arrays::Rect<1> rect = runtime->get_index_space_domain(
      regions[0].get_logical_region().get_index_space()).get_rect<1>();
    LegionRuntime::Arrays::Rect<1> subrect;
    LegionRuntime::Accessor::ByteOffset stride;
    void *base_ptr = accessor.raw_rect_ptr<1>(rect, subrect, &stride);
    assert(base_ptr);
    assert(subrect == rect);
    assert(rect.lo == LegionRuntime::Arrays::Point<1>::ZEROES());
    assert(stride.offset == 1);

    // Call base jump.
    return jump_internal(filename, offset, (Pds::Dgram *)base_ptr);
  }

  static Legion::TaskID register_jump_task() {
    static const char * const task_name = "jump";
    Legion::TaskVariantRegistrar registrar(task_id, task_name);
    registrar.add_constraint(Legion::ProcessorConstraint(Legion::Processor::IO_PROC));
    Legion::Runtime::preregister_task_variant<bool, jump_task>(registrar, task_name);
    return task_id;
  }

  static Pds::Dgram *launch_jump_task(Legion::HighLevelRuntime *runtime,
                                      Legion::Context ctx,
                                      const std::string &filename,
                                      int64_t offset) {
    // Create destination region.
    Legion::Domain domain = Legion::Domain::from_rect<1>(
      LegionRuntime::Arrays::Rect<1>(LegionRuntime::Arrays::Point<1>(0), LegionRuntime::Arrays::Point<1>(MaxDgramSize-1)));
    Legion::IndexSpace ispace = runtime->create_index_space(ctx, domain);
    Legion::FieldSpace fspace = runtime->create_field_space(ctx);
    Legion::FieldID fid = 101;
    {
      Legion::FieldAllocator fsa(runtime->create_field_allocator(ctx, fspace));
      fsa.allocate_field(1, fid);
    }
    Legion::LogicalRegion region = runtime->create_logical_region(ctx, ispace, fspace);

    // Launch task.
    Legion::Future f;
    {
      Args args;
      args.filename_size = filename.size();
      args.offset = offset;
      args.fid = fid;
      size_t bufsize = sizeof(Args) + filename.size();
      char *buffer = (char *)malloc(bufsize);
      *(Args *)buffer = args;
      memcpy(buffer + sizeof(Args), filename.c_str(), filename.size());
      Legion::TaskArgument targs(buffer, bufsize);

      Legion::TaskLauncher task(task_id, targs);
      task.add_region_requirement(
        Legion::RegionRequirement(region, READ_WRITE, EXCLUSIVE, region).add_field(fid));
      f = runtime->execute_task(ctx, task);
      free(buffer);
    }

    // Map destination region.
    Pds::Dgram *dg = 0;
    {
      // Launch mapping first to avoid blocking analysis on task execution.
      Legion::InlineLauncher mapping(
        Legion::RegionRequirement(region, READ_WRITE, EXCLUSIVE, region).add_field(fid));
      Legion::PhysicalRegion physical = runtime->map_region(ctx, mapping);

      // Skip if task failed.
      if (f.get_result<bool>()) {
        physical.wait_until_valid();

        LegionRuntime::Accessor::RegionAccessor<LegionRuntime::Accessor::AccessorType::Generic> accessor =
          physical.get_field_accessor(fid);
        LegionRuntime::Arrays::Rect<1> rect = runtime->get_index_space_domain(
          physical.get_logical_region().get_index_space()).get_rect<1>();
        LegionRuntime::Arrays::Rect<1> subrect;
        LegionRuntime::Accessor::ByteOffset stride;
        void *base_ptr = accessor.raw_rect_ptr<1>(rect, subrect, &stride);
        assert(base_ptr);
        assert(subrect == rect);
        assert(rect.lo == LegionRuntime::Arrays::Point<1>::ZEROES());
        assert(stride.offset == 1);

        Pds::Dgram *dghdr = (Pds::Dgram *)base_ptr;
        if (dghdr->xtc.sizeofPayload()>MaxDgramSize)
          MsgLog(logger, fatal, "Datagram size exceeds sanity check. Size: " << dghdr->xtc.sizeofPayload() << " Limit: " << MaxDgramSize);
        Pds::Dgram* dg = (Pds::Dgram*)new char[sizeof(Pds::Dgram)+dghdr->xtc.sizeofPayload()];
        memcpy(dg, base_ptr, sizeof(Pds::Dgram)+dghdr->xtc.sizeofPayload());
      }
    }

    // Destroy temporary region.
    runtime->destroy_logical_region(ctx, region);
    runtime->destroy_index_space(ctx, ispace);
    runtime->destroy_field_space(ctx, fspace);

    return dg;
  }

  static const Legion::TaskID task_id = 501976; // chosen by fair dice roll
  struct Args {
    size_t filename_size;
    int64_t offset;
    Legion::FieldID fid;
  };

  Pds::Dgram* jump_blocking(const std::string& filename, int64_t offset) {
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

  Pds::Dgram* jump(const std::string& filename, int64_t offset, uintptr_t runtime_, uintptr_t ctx_) {
#define RANDOM_ACCESS_USE_JUMP_TASK 1
#if RANDOM_ACCESS_USE_JUMP_TASK
    ::legion_runtime_t c_runtime = *(::legion_runtime_t *)runtime_;
    ::legion_context_t c_ctx = *(::legion_context_t *)runtime_;
    Legion::Runtime *runtime = Legion::CObjectWrapper::unwrap(c_runtime);
    Legion::Context ctx = Legion::CObjectWrapper::unwrap(c_ctx)->context();
    return launch_jump_task(runtime, ctx, filename, offset);
#else
    return jump_blocking(filename, offset);
#endif
  }

private:
  enum {MaxDgramSize=0x2000000};
  static std::map<std::string, int> _fd;
  static pthread_mutex_t _fd_mutex;
};

static Legion::TaskID __attribute__((unused)) _force_jump_task_static_initialize =
  RandomAccessXtcReader::register_jump_task();

std::map<std::string, int> RandomAccessXtcReader::_fd;
pthread_mutex_t RandomAccessXtcReader::_fd_mutex = PTHREAD_MUTEX_INITIALIZER;

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
      Pds::Dgram* dg = _xtc.jump_blocking(filename, offset);
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
    Pds::Dgram* dg = _xtc.jump_blocking(filename, _beginrunOffset);
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
  int jump(const std::vector<std::string>& filenames, const std::vector<int64_t> &offsets, const std::string &lastBeginCalibCycleDgram, uintptr_t runtime, uintptr_t ctx) {
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
      dg = _xtc.jump(filename, offset, runtime, ctx);
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

int RandomAccess::jump(const std::vector<std::string>& filenames, const std::vector<int64_t> &offsets, const std::string &lastBeginCalibCycleDgram, uintptr_t runtime, uintptr_t ctx) {
  return _raxrun->jump(filenames, offsets, lastBeginCalibCycleDgram, runtime, ctx);
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
