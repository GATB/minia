/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _MINIA_HPP_
#define _MINIA_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Tool.hpp>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>

#include <string>

/********************************************************************************/

/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

/********************************************************************************/
/********************************************************************************/

class Minia : public misc::impl::Tool
{
private:

    size_t          _kmerSize;
    std::string     _solidFile;

public:

    /** */
    Minia ();

    /** */
    static const char* STR_URI_DEBLOOM;

private:

    /** */
    void execute ();

    /** */
    virtual collections::impl::Bloom<kmer_type>* createBloom ();
};

/********************************************************************************/

#endif /* _MINIA_HPP_ */

