/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _GATB_TOOLS_FRONTLINE_HPP_
#define _GATB_TOOLS_FRONTLINE_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>
#include <Terminator.hpp>
#include <set>
#include <queue>

/********************************************************************************/

// types using in advanced traversal functions
struct NodeNt
{
    Node       node;
    Nucleotide nt;

    NodeNt (const Node& node, const Nucleotide& nt)  : node(node), nt(nt)  {}

    bool operator< (const NodeNt &other) const
    {
        // need to define a strict weak ordering
        if (node.kmer   != other.node.kmer)    {  return (node.kmer   < other.node.kmer);    }
        if (node.strand != other.node.strand)  {  return (node.strand < other.node.strand);  }
        return (nt < other.nt);
    }
};

/********************************************************************************/

// auxiliary class that is used by MonumentTraversal and deblooming
class Frontline
{
public:

    /** Constructor. */
    Frontline (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode,
        const Node&       previousNode,
        std::set<Node>&   all_involved_extensions
    );

    /** Constructor. */
    Frontline (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode
    );

    /** */
    virtual ~Frontline() {}

    /** */
    bool go_next_depth();

    size_t size  () const  {  return _frontline.size();  }
    size_t depth () const  {  return _depth;             }

    NodeNt front () { return _frontline.front(); }

protected:

    virtual bool check (const Node& node)  { return true; }

    Direction _direction;

    const Graph& _graph;

    Terminator&  _terminator;

    typedef std::queue<NodeNt> queue_nodes;
    queue_nodes _frontline;

    int  _depth;

    std::set<Node>* _all_involved_extensions;

    std::set<Node::Type> _already_frontlined; // making it simpler now
};

/********************************************************************************/

// auxiliary class that is used by MonumentTraversal and deblooming
class FrontlineBranching : public Frontline
{
public:

    /** Constructor. */
    FrontlineBranching (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode,
        const Node&       previousNode,
        std::set<Node>&   all_involved_extensions
    )  : Frontline (direction,graph,terminator,startingNode,previousNode,all_involved_extensions) {}

    /** Constructor. */
    FrontlineBranching (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode
    ) : Frontline(direction,graph,terminator,startingNode) {}

private:

    bool check (const Node& node);
};

/********************************************************************************/

#endif /* _GATB_TOOLS_FRONTLINE_HPP_ */

