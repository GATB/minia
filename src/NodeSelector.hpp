/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _GATB_TOOLS_NODE_SELECTOR_HPP_
#define _GATB_TOOLS_NODE_SELECTOR_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>
#include <Terminator.hpp>

/********************************************************************************/

/** */
class INodeSelector
{
public:

    virtual ~INodeSelector()  {}

    virtual bool select (const Node& source, Node& result) = 0;
};

/********************************************************************************/

// initial k-mer selection function from original Minia release (up until summer 2013)
/* rationale:
 * it's established that it's not a good idea to assemble from a branching kmer
 *      -> what if it's a simplepathtraversal and it chooses to start in a deadend?
 *      -> what if it's a monumenttraversal and the kmer is a true branching: won't be traversed
 *  yet, branching kmers are the only indexed ones
 *  solution:
 *      detect a 2k+2 simple path (anything NOT deadend or snp) around the branching kmer and start to extend from it
 */
class NodeSelectorSimplePath : public INodeSelector
{
public:

    NodeSelectorSimplePath (const Graph& graph, Terminator& terminator)
        : _graph(graph), _terminator(terminator) {}

    bool select (const Node& source, Node& result);

private:

    const Graph&  _graph;
    Terminator&   _terminator;
};

/********************************************************************************/

#endif /* _GATB_TOOLS_NODE_SELECTOR_HPP_ */

