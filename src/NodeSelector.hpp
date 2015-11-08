/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _GATB_TOOLS_NODE_SELECTOR_HPP_
#define _GATB_TOOLS_NODE_SELECTOR_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

/** */
template <typename Node, typename Edge, typename GraphDataVariant>
class INodeSelector : public SmartPointer
{
public:

    virtual ~INodeSelector<Node,Edge,GraphDataVariant>()  {}

    virtual bool select (const Node& source, Node& result) = 0;

    virtual std::string getName () const = 0;
};

/********************************************************************************/

template <typename Node, typename Edge, typename GraphDataVariant>
class NodeSelectorFactory
{
public:

    static NodeSelectorFactory<Node,Edge,GraphDataVariant>& singleton()  { static NodeSelectorFactory<Node,Edge,GraphDataVariant> instance; return instance; }

    INodeSelector<Node,Edge,GraphDataVariant>* create (const std::string& type, const GraphTemplate<Node,Edge,GraphDataVariant>& graph, TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator);
};

/********************************************************************************/

/** Abstract factorization. */

template <typename Node, typename Edge, typename GraphDataVariant>
class NodeSelectorAbstract : public INodeSelector<Node,Edge,GraphDataVariant>
{
public:

    NodeSelectorAbstract (const GraphTemplate<Node,Edge,GraphDataVariant>& graph, TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator)
        : _graph(graph), _terminator(terminator) {}

protected:
    const GraphTemplate<Node,Edge,GraphDataVariant>&  _graph;
    TerminatorTemplate<Node,Edge,GraphDataVariant>&   _terminator;
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
template <typename Node, typename Edge, typename GraphDataVariant>
class NodeSelectorSimplePath : public NodeSelectorAbstract<Node,Edge,GraphDataVariant>
{
public:

    NodeSelectorSimplePath (const GraphTemplate<Node,Edge,GraphDataVariant>& graph, TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator)
        : NodeSelectorAbstract<Node,Edge,GraphDataVariant> (graph, terminator) {}

    bool select (const Node& source, Node& result);

    std::string getName () const  { return "simple"; }
};

/********************************************************************************/

template <typename Node, typename Edge, typename GraphDataVariant>
class NodeSelectorImproved : public NodeSelectorAbstract<Node,Edge,GraphDataVariant>
{
public:

    NodeSelectorImproved (const GraphTemplate<Node,Edge,GraphDataVariant>& graph, TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator)
        : NodeSelectorAbstract<Node,Edge,GraphDataVariant> (graph, terminator) {}

    bool select (const Node& source, Node& result);

    std::string getName () const  { return "improved"; }
};

/********************************************************************************/

template <typename Node, typename Edge, typename GraphDataVariant>
class NodeSelectorBest : public NodeSelectorAbstract<Node,Edge,GraphDataVariant>
{
public:

    NodeSelectorBest (const GraphTemplate<Node,Edge,GraphDataVariant>& graph, TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator)
        : NodeSelectorAbstract<Node,Edge,GraphDataVariant> (graph, terminator),  _firstSelector (graph, terminator) {}

    bool select (const Node& source, Node& result);

    std::string getName () const  { return "best"; }

private:
    NodeSelectorImproved<Node,Edge,GraphDataVariant> _firstSelector;
};

/********************************************************************************/

#endif /* _GATB_TOOLS_NODE_SELECTOR_HPP_ */

