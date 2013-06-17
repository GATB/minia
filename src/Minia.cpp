/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>

#include <omptl/omptl_numeric>
#include <omptl/omptl_algorithm>

#include <omp.h>

#include <Minia.hpp>
#include <DSK.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  printf a

/********************************************************************************/

const char* Minia::STR_URI_DEBLOOM = "-debloom";

/********************************************************************************/

class BuildKmerExtension : public IteratorFunctor
{
public:
    void operator() (const kmer_type& kmer)
    {
        /** We configure the neighbor kmers iterator for a given source kmer. */
        _itNeighbors.setSource (kmer);

        /** We iterate all neighbors. */
        for (_itNeighbors.first(); !_itNeighbors.isDone(); _itNeighbors.next())
        {
            /** If the bloom contains the current neighbor, we add it to the debloom file. */
            if (_bloom->contains (*_itNeighbors))
            {
                _extendBag.insert (*_itNeighbors);
            }
        }
    }

    BuildKmerExtension (KmerModel& model, Bloom<kmer_type>* bloom, Bag<kmer_type>* extendBag)
        : _bloom(bloom), _extendBag(extendBag, 5*1000, getSynchro()), _itNeighbors(model)  { }

    ~BuildKmerExtension ()  {  _extendBag.flush();  }

    Bloom<kmer_type>*   _bloom;
    BagCache<kmer_type> _extendBag;
    ModelAbstract<kmer_type>::KmerNeighborIterator _itNeighbors;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Minia::Minia () : Tool("minia"), _kmerSize(0)
{
    /** We add options specific to this tool. */
    _parser->add (new OptionOneParam (DSK::STR_KMER_SIZE,       "size of a kmer",   true                ));
    _parser->add (new OptionOneParam (DSK::STR_URI_SOLID_KMERS, "solid kmers file", false               ));
    _parser->add (new OptionOneParam (Minia::STR_URI_DEBLOOM,   "debloom file",     false,  "debloom"   ));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Minia::execute ()
{
    /*************************************************/
    // We set some attributes (shortcuts).
    /*************************************************/
    _kmerSize  = _input->getInt (DSK::STR_KMER_SIZE);
    _solidFile = _input->getStr (DSK::STR_URI_SOLID_KMERS);

    KmerModel model (_kmerSize);

    /*************************************************/
    /** We create a bloom with inserted solid kmers. */
    /*************************************************/
    Bloom<kmer_type>* bloom = createBloom ();
    LOCAL (bloom);

    /*************************************************/
    /** We create an iterator over the solid kmers.  */
    /*************************************************/
    Iterator<kmer_type>* itKmers = createIterator<kmer_type> (
        new IteratorFile<kmer_type> (_input->getStr (DSK::STR_URI_SOLID_KMERS)),
        System::file().getSize(_input->getStr (DSK::STR_URI_SOLID_KMERS)) / sizeof (kmer_type),
        "iterate solid kmers"
    );
    LOCAL (itKmers);

    /*************************************************/
    /** We create the debloom file.                 */
    /*************************************************/

    /** First, we delete the debloom file if already existing. */
    System::file().remove (getUriByKey(STR_URI_DEBLOOM));

    Bag<kmer_type>* debloomFile = new BagFile<kmer_type>(getUriByKey(STR_URI_DEBLOOM));
    LOCAL (debloomFile);

    /*************************************************/
    /** We fill the debloom file.                    */
    /*************************************************/
    {
        TIME_INFO (_timeInfo, "fill debloom file");

        _dispatcher->iterate (itKmers, BuildKmerExtension (model, bloom, debloomFile));
    }

    /** We make sure everything is put into the extension file. */
    debloomFile->flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bloom<kmer_type>* Minia::createBloom ()
{
    TIME_INFO (_timeInfo, "fill bloom filter");

    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

    u_int64_t solidFileSize = (System::file().getSize(_solidFile) / sizeof (kmer_type));

    u_int64_t estimatedBloomSize = solidFileSize * NBITS_PER_KMER;
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator<kmer_type>* itKmers = createIterator<kmer_type> (
        new IteratorFile<kmer_type> (_solidFile),
        solidFileSize,
        "fill bloom filter"
    );
    LOCAL (itKmers);

    /** We use a bloom builder. */
    BloomBuilder builder (itKmers, estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER), _input->getInt(Tool::STR_NB_CORES));

    /** We instantiate the bloom object. */
    IProperties* bloomProps = new Properties();
    Bloom<kmer_type>* bloom = builder.build (bloomProps);

    _info->add (1, bloomProps);

    /** We return the created bloom filter. */
    return bloom;
}
