#ifndef __UF_H__
#define __UF_H__
#include <boost/fusion/container/map.hpp>
#include <ginac/ginac.h>
#include <string>
using namespace boost;
using namespace GiNaC;
namespace UFX
{
  struct U;
  struct F;
  struct xlst;
}

typedef fusion::map<fusion::pair< UFX::U, ex>, fusion::pair< UFX::F, ex>, fusion::pair< UFX::xlst,lst> > UFXmap;
typedef fusion::map< fusion::pair< UFX::F, ex>, fusion::pair< UFX::xlst,lst> > FXmap;

/// Factory for symbol creation
const symbol & get_symbol(const std::string &);

UFXmap UF(lst k_lst,lst p_lst,lst subs_lst, unsigned int displacement = 0);

#endif // __UF_H__
