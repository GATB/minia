/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _GATB_TOOLS_MINIA_HPP_
#define _GATB_TOOLS_MINIA_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

/**
 */
class Minia : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    Minia ();

private:

    /** \copydoc Tool::execute. */
    void  execute ();

    /** */
    void assemble (const Graph& graph);

    /** */
    void buildSequence (
        const Graph& graph,
        const Node& startingNode,
        size_t length,
        size_t nbContigs,
        const std::vector<Nucleotide>& consensusRight,
        const std::vector<Nucleotide>& consensusLeft,
        Sequence& seq
    );
};

/********************************************************************************/

#endif /* _GATB_TOOLS_MINIA_HPP_ */

