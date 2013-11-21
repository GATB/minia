/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _GATB_TOOLS_UTILS_HPP_
#define _GATB_TOOLS_UTILS_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

struct PATH : public std::vector<Edge>  {};

inline bool operator< (const PATH& a, const PATH& b)
{
    size_t N = std::min(a.size(),b.size());
    for (size_t i=0; i<N; i++)
    {
             if (ascii(a[i].nt) < ascii(b[i].nt)) { return true;  }
        else if (ascii(a[i].nt) > ascii(b[i].nt)) { return false; }
    }
    return a.size() < b.size();
}

inline std::ostream& operator<< (std::ostream& s, const PATH& p)
{
    for (size_t i=0; i<p.size(); i++)  { s << ascii(p[i].nt); }
    return s;
}

/********************************************************************************/

float needleman_wunch (const PATH& a, const PATH& b);

/********************************************************************************/

#endif /* _GATB_TOOLS_UTILS_HPP_ */

