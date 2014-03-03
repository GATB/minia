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
#include <Terminator.hpp>

/********************************************************************************/

/** */
class INodeSelector : public SmartPointer
{
public:

    virtual ~INodeSelector()  {}

    virtual bool select (const Node& source, Node& result) = 0;

    virtual std::string getName () const = 0;
};

/********************************************************************************/

class NodeSelectorFactory
{
public:

    static NodeSelectorFactory& singleton()  { static NodeSelectorFactory instance; return instance; }

    INodeSelector* create (const std::string& type, const Graph& graph, Terminator& terminator);
};

/********************************************************************************/

/** Abstract factorization. */

class NodeSelectorAbstract : public INodeSelector
{
public:

    NodeSelectorAbstract (const Graph& graph, Terminator& terminator)
        : _graph(graph), _terminator(terminator) {}

protected:
    const Graph&  _graph;
    Terminator&   _terminator;
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
class NodeSelectorSimplePath : public NodeSelectorAbstract
{
public:

    NodeSelectorSimplePath (const Graph& graph, Terminator& terminator)
        : NodeSelectorAbstract (graph, terminator) {}

    bool select (const Node& source, Node& result);

    std::string getName () const  { return "simple"; }
};

/********************************************************************************/

class NodeSelectorImproved : public NodeSelectorAbstract
{
public:

    NodeSelectorImproved (const Graph& graph, Terminator& terminator)
        : NodeSelectorAbstract (graph, terminator) {}

    bool select (const Node& source, Node& result);

    std::string getName () const  { return "improved"; }
};

/********************************************************************************/

class NodeSelectorBest : public NodeSelectorAbstract
{
public:

    NodeSelectorBest (const Graph& graph, Terminator& terminator)
        : NodeSelectorAbstract (graph, terminator),  _firstSelector (graph, terminator) {}

    bool select (const Node& source, Node& result);

    std::string getName () const  { return "best"; }

private:
    NodeSelectorImproved _firstSelector;
};

/********************************************************************************/

#endif /* _GATB_TOOLS_NODE_SELECTOR_HPP_ */

