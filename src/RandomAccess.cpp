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
#include <time.h>

#ifdef PSANA_USE_LEGION

#include <pthread.h>
#include <legion.h>
#include <legion_c.h>
#include <legion_c_util.h>
#include <unistd.h>

#endif

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
#ifdef PSANA_USE_LEGION
    if (pthread_mutex_lock(&_fd_mutex)) {
      assert(false && "pthread_mutex_lock failed\n");
    }
#endif

    if (has(filename)) {
      int result = _fd[filename];

#ifdef PSANA_USE_LEGION
      if (pthread_mutex_unlock(&_fd_mutex)) {
        assert(false && "pthread_mutex_unlock failed\n");
      }
#endif

      return result;
    }

    int fd = ::open(filename.c_str(), O_RDONLY | O_LARGEFILE);
    if (fd==-1) MsgLog(logger, fatal,
                               "File " << filename.c_str() << " not found");
    _fd[filename] = fd;

#ifdef PSANA_USE_LEGION
    if (pthread_mutex_unlock(&_fd_mutex)) {
      assert(false && "pthread_mutex_unlock failed\n");
    }
#endif

    return fd;
  }

  ~RandomAccessXtcReader() {
    // FIXME: Not safe to close fds in the destructor since the memory is now static
    // for  (std::map<std::string, int>::const_iterator it = _fd.begin(); it!= _fd.end(); it++)
    //   ::close(it->second);
  }

#ifdef PSANA_USE_LEGION
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
    size_t count;
    std::vector<int64_t> offsets;
    std::vector<std::string> filenames;
    {
      const char *current = (const char *)task->args;

      size_t total_bytes = *(size_t *)current; current += sizeof(size_t);
      count = *(size_t *)current; current += sizeof(size_t);

      for (size_t i = 0; i < count; i++) {
        offsets.push_back(*(int64_t *)current); current += sizeof(int64_t);
      }

      std::vector<size_t> filename_sizes;
      for (size_t i = 0; i < count; i++) {
        filename_sizes.push_back(*(size_t *)current); current += sizeof(size_t);
      }

      for (size_t i = 0; i < count; i++) {
        filenames.push_back(std::string(current, filename_sizes[i]));
        current += filename_sizes[i];
      }

      assert(current == ((const char *)task->args) + total_bytes);
    }

    // Fetch destination pointer out of region argument.
    LegionRuntime::Arrays::Rect<1> rect = runtime->get_index_space_domain(
      regions[0].get_logical_region().get_index_space()).get_rect<1>();
    LegionRuntime::Arrays::Rect<1> subrect;
    LegionRuntime::Accessor::ByteOffset stride;

    bool ok = true;
    for (size_t i = 0; i < count; i++) {
      LegionRuntime::Accessor::RegionAccessor<LegionRuntime::Accessor::AccessorType::Generic> accessor =
        regions[0].get_field_accessor(fid_base+i);
      void *base_ptr = accessor.raw_rect_ptr<1>(rect, subrect, &stride);
      assert(base_ptr);
      assert(subrect == rect);
      assert(rect.lo == LegionRuntime::Arrays::Point<1>::ZEROES());
      assert(stride.offset == 1);

      // Call base jump.
      // FIXME: Maybe not the right error behavior
      ok = ok && jump_internal(filenames[i], offsets[i], (Pds::Dgram *)base_ptr);
    }
    return ok;
  }

  static Legion::TaskID register_jump_task() {
    static const char * const task_name = "jump";
    Legion::TaskVariantRegistrar registrar(task_id, task_name, false /* global */);
    registrar.add_constraint(Legion::ProcessorConstraint(Legion::Processor::IO_PROC));
    Legion::Runtime *runtime = Legion::Runtime::get_runtime();
    runtime->register_task_variant<bool, jump_task>(registrar);
    runtime->attach_name(task_id, task_name);
    return task_id;
  }

  static std::vector<Pds::Dgram *> launch_jump_task(Legion::HighLevelRuntime *runtime,
                                                    Legion::Context ctx,
                                                    const std::vector<std::string> &filenames,
                                                    const std::vector<int64_t> &offsets) {
    size_t count = filenames.size();

    // Create destination region.
    Legion::Domain domain = Legion::Domain::from_rect<1>(
      LegionRuntime::Arrays::Rect<1>(LegionRuntime::Arrays::Point<1>(0), LegionRuntime::Arrays::Point<1>(MaxDgramSize-1)));
    Legion::IndexSpace ispace = runtime->create_index_space(ctx, domain);
    Legion::FieldSpace fspace = runtime->create_field_space(ctx);
    {
      Legion::FieldAllocator fsa(runtime->create_field_allocator(ctx, fspace));
      for (size_t i = 0; i < count; i++) {
        fsa.allocate_field(1, fid_base+i);
      }
    }
    Legion::LogicalRegion region = runtime->create_logical_region(ctx, ispace, fspace);

    // Launch task.
    Legion::Future f;
    {
      size_t total_filenames_size = 0;
      for (size_t i = 0; i < count; i++) {
        total_filenames_size += filenames[i].size();
      }

      size_t bufsize = sizeof(Args) + count*sizeof(int64_t) + count*sizeof(size_t) + total_filenames_size;
      Args args;
      args.total_bytes = bufsize;
      args.count = count;
      char *buffer = new char[bufsize];
      {
        char *current = buffer;

        *(Args *)current = args; current += sizeof(Args);

        for (size_t i = 0; i < count; i++) {
          *(int64_t *)current = offsets[i]; current += sizeof(int64_t);
        }

        for (size_t i = 0; i < count; i++) {
          *(size_t *)current = filenames[i].size(); current += sizeof(size_t);
        }

        for (size_t i = 0; i < count; i++) {
          memcpy(current, filenames[i].c_str(), filenames[i].size());
          current += filenames[i].size();
        }

        assert(current == buffer + bufsize);
      }
      Legion::TaskArgument targs(buffer, bufsize);

      Legion::TaskLauncher task(task_id, targs);
      task.add_region_requirement(
        Legion::RegionRequirement(region, READ_WRITE, EXCLUSIVE, region));
      for (size_t i = 0; i < count; i++) {
        task.add_field(0, fid_base+i);
      }
      f = runtime->execute_task(ctx, task);
      delete[] buffer;
    }

    // Map destination region.
    std::vector<Pds::Dgram *> dgs;
    {
      // Launch mapping first to avoid blocking analysis on task execution.
      Legion::InlineLauncher mapping(
        Legion::RegionRequirement(region, READ_WRITE, EXCLUSIVE, region));
      for (size_t i = 0; i < count; i++) {
        mapping.add_field(fid_base+i);
      }
      Legion::PhysicalRegion physical = runtime->map_region(ctx, mapping);

      // Skip if task failed.
      if (f.get_result<bool>()) {
        physical.wait_until_valid();

        LegionRuntime::Arrays::Rect<1> rect = runtime->get_index_space_domain(
          physical.get_logical_region().get_index_space()).get_rect<1>();
        LegionRuntime::Arrays::Rect<1> subrect;
        LegionRuntime::Accessor::ByteOffset stride;
        std::vector<void *> base_ptr;
        for (size_t i = 0; i < count; i++) {
          LegionRuntime::Accessor::RegionAccessor<LegionRuntime::Accessor::AccessorType::Generic> accessor =
            physical.get_field_accessor(fid_base+i);

          void * base_ptr = accessor.raw_rect_ptr<1>(rect, subrect, &stride);
          assert(base_ptr);
          assert(subrect == rect);
          assert(rect.lo == LegionRuntime::Arrays::Point<1>::ZEROES());
          assert(stride.offset == 1);


          Pds::Dgram *dghdr = (Pds::Dgram *)base_ptr;
          if (dghdr->xtc.sizeofPayload()>MaxDgramSize)
            MsgLog(logger, fatal, "Datagram size exceeds sanity check. Size: " << dghdr->xtc.sizeofPayload() << " Limit: " << MaxDgramSize);
          Pds::Dgram *dg = (Pds::Dgram*)new char[sizeof(Pds::Dgram)+dghdr->xtc.sizeofPayload()];
          memcpy(dg, base_ptr, sizeof(Pds::Dgram)+dghdr->xtc.sizeofPayload());
          dgs.push_back(dg);
        }
      }
    }

    // Destroy temporary region.
    runtime->destroy_logical_region(ctx, region);
    runtime->destroy_index_space(ctx, ispace);
    runtime->destroy_field_space(ctx, fspace);

    return dgs;
  }

  static const Legion::TaskID task_id = 501976; // chosen by fair dice roll
  static const Legion::FieldID fid_base = 100;
  struct Args {
    size_t total_bytes;
    size_t count;
    // std::vector<int64_t> offsets;
    // std::vector<size_t> filename_sizes;
    // std::vector<std::string content> filenames;
  };
#endif

  // retry with sleep to support live mode
  off64_t lseek64_retry(int fd, off64_t offset, int whence) {
    timespec req,rem;
    req.tv_sec = 1; req.tv_nsec = 0;
    int64_t found;
    for (int i=0; i<3; i++) {
      found = lseek64(fd,offset, SEEK_SET);
      if (found!=-1) return found;
      nanosleep(&req, &rem); // wait for the live-mode data
    }
    return found;
  };

  // retry with sleep to support live mode
  // interface tries to behave (somewhat) like ::read
  ssize_t read_retry(int fd, char* buf, size_t count) {
    timespec req,rem;
    req.tv_sec = 1; req.tv_nsec = 0;
    ssize_t bytes_read, remaining;
    remaining = count;
    for (int i=0; i<3; i++) {
      bytes_read = ::read(fd, buf, remaining);
      if (bytes_read<0) return -1; // abject failure
      remaining -= bytes_read;
      if (remaining==0) {
        return count;    // success
      } else {
        buf+=bytes_read; // partial success
      }
      nanosleep(&req, &rem); // wait for the live-mode data
    }
    return -1; // even more failure
  };

  Pds::Dgram* jump_blocking(const std::string& filename, int64_t offset) {
    int fd = ensure(filename);
    int64_t found = lseek64_retry(fd,offset, SEEK_SET);
    if (found != offset) {
      stringstream ss;
      ss << "Jump to offset " << offset << " failed";
      MsgLog(logger, error, ss.str());
      throw RandomAccessSeekFailed(ERR_LOC);
    }
    Pds::Dgram dghdr;
    if (read_retry(fd, (char*)&dghdr, sizeof(dghdr))<0) {
      return 0;
    } else {
      if (dghdr.xtc.sizeofPayload()>MaxDgramSize)
        MsgLog(logger, fatal, "Datagram size exceeds sanity check. Size: " << dghdr.xtc.sizeofPayload() << " Limit: " << MaxDgramSize);
      Pds::Dgram* dg = (Pds::Dgram*)new char[sizeof(dghdr)+dghdr.xtc.sizeofPayload()];
      *dg = dghdr;
      if (read_retry(fd, dg->xtc.payload(), dg->xtc.sizeofPayload())<0) {
        return 0;
      }
      return dg;
    }
  }

  std::vector<Pds::Dgram *> jump_async(const std::vector<std::string> &filenames, const std::vector<int64_t> &offsets, uintptr_t runtime_, uintptr_t ctx_) {
#ifdef PSANA_USE_LEGION
    ::legion_runtime_t c_runtime = *(::legion_runtime_t *)runtime_;
    ::legion_context_t c_ctx = *(::legion_context_t *)ctx_;
    Legion::Runtime *runtime = Legion::CObjectWrapper::unwrap(c_runtime);
    Legion::Context ctx = Legion::CObjectWrapper::unwrap(c_ctx)->context();
    return launch_jump_task(runtime, ctx, filenames, offsets);
#else
    std::vector<Pds::Dgram *> dgs;
    for (size_t i = 0; i < filenames.size(); i++) {
      const std::string &filename = filenames[i];
      int64_t offset = offsets[i];
      dgs.push_back(jump_blocking(filename, offset));
    }
    return dgs;
#endif
  }

private:
  enum {MaxDgramSize=0x2000000};
  static std::map<std::string, int> _fd;
#ifdef PSANA_USE_LEGION
  static pthread_mutex_t _fd_mutex;
#endif
};

#ifdef PSANA_USE_LEGION
static Legion::TaskID __attribute__((unused)) _force_jump_task_static_initialize =
  RandomAccessXtcReader::register_jump_task();
#endif

std::map<std::string, int> RandomAccessXtcReader::_fd;
#ifdef PSANA_USE_LEGION
pthread_mutex_t RandomAccessXtcReader::_fd_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

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
    bool accept = false;
    assert(filenames.size() == offsets.size());
    std::vector<Pds::Dgram*> dgs = _xtc.jump_async(filenames, offsets, runtime, ctx);

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

    for (size_t i = 0; i < dgs.size(); i++) {
      const std::string& filename = filenames[i];
      Pds::Dgram* dg = dgs[i];
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
