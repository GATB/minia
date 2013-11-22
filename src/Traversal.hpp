/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _GATB_TOOLS_TRAVERSAL_HPP_
#define _GATB_TOOLS_TRAVERSAL_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>
#include <Utils.hpp>
#include <Terminator.hpp>
#include <NodeSelector.hpp>
#include <set>

/********************************************************************************/
// semi-abstract class. implements traverse but not avance
class Traversal : public SmartPointer
{
public:

    /** */
    static Traversal* create (
        const std::string& type,
        const Graph& graph,
        Terminator& terminator
    );

    /** */
    virtual ~Traversal();

    /** */
    virtual std::string getName() const = 0;

    /** */
    int traverse (const Node& node, Direction dir, std::vector<Nucleotide>& resulting_sequence);

protected:

    /** */
    Traversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen,
        int max_depth,
        int max_breadth
    );

    const Graph& graph;
    Terminator&  terminator;

    int maxlen;
    int max_depth;
    int max_breadth;

    virtual char avance (const Node& node, Direction dir, bool first_extension, PATH& path, const Node& previousNode) = 0;

    void mark_extensions (std::set<Node>& extensions_to_mark);

    // record the start/end positions of traversed bubbles (only from the latest traverse() call)
    std::vector <std::pair<int, int> > bubbles_positions;
};

/********************************************************************************/

class SimplePathsTraversal: public Traversal
{
public:
    /** */
    SimplePathsTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = 1000000,
        int max_depth   = 500,
        int max_breadth = 20
    );

    std::string getName() const  { return std::string ("unitig"); }

private:

    char avance (const Node& node, Direction dir, bool first_extension, PATH& path, const Node& previousNode);
};

/********************************************************************************/

class MonumentTraversal: public Traversal
{
public:
    /** */
    MonumentTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = 1000000,
        int max_depth   = 500,
        int max_breadth = 20
    );

    std::string getName() const  { return std::string ("monument"); }

    bool explore_branching (
        const Node& node,
        Direction dir,
        PATH& consensus,
        const Node& previousNode,
        std::set<Node>& all_involved_extensions
    );

private:

    char avance (const Node& node, Direction dir, bool first_extension, PATH& path, const Node& previousNode);

    bool explore_branching (
        const Node& node,
        Direction dir,
        PATH& consensus,
        const Node& previousNode
    );


    int find_end_of_branching (
        Direction dir,
        const Node& startingNode,
        Node& endNode,
        const Node& previousNode,
        std::set<Node>& all_involved_extensions
    );

    std::set<PATH> all_consensuses_between (
        Direction    dir,
        const Node& startNode,
        const Node& endNode,
        int traversal_depth,
        std::set<Node> usedNode,
        PATH current_consensus,
        bool& success
    );

    std::set<PATH> all_consensuses_between (
        Direction    dir,
        const Node& startNode,
        const Node& endNode,
        int traversal_depth,
        bool &success
    );

    bool validate_consensuses (std::set<PATH>& consensuses, PATH& consensus);

    bool all_consensuses_almost_identical (std::set<PATH>& consensuses);

    void mark_extensions (std::set<Node>& extensions_to_mark);

    static const int consensuses_identity = 90; // traversing bubble if paths are all pair-wise identical by > 90%
};

/********************************************************************************/

#endif /* _GATB_TOOLS_TRAVERSAL_HPP_ */

