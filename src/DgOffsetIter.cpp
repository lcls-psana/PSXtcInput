#include "PSXtcInput/DgOffsetIter.h"

#include "pdsdata/xtc/XtcIterator.hh"
#include "pdsdata/xtc/Sequence.hh"
#include "pdsdata/xtc/Xtc.hh"

DgOffsetIter::DgOffsetIter(Pds::Xtc* xtc) :
  XtcIterator(xtc), payload(-1LL, 0) {}

int DgOffsetIter::process(Pds::Xtc* xtc) {
  Pds::Level::Type   level     = xtc->src.level();
  if (level < 0 || level >= Pds::Level::NumberOfLevels )
  {
      return Continue;
  }
  switch (xtc->contains.id()) {
  case (Pds::TypeId::Id_Xtc) : {
    iterate(xtc);
    break;
  }
  case (Pds::TypeId::Id_SmlDataOrigDgramOffset) : {
    payload = *reinterpret_cast<OffsetPayload *>(xtc->payload());
    return Stop;
  }
  default:
    break;
  }
  return Continue;
}
