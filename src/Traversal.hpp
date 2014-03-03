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

#ifndef _GATB_TOOLS_TRAVERSAL_HPP_
#define _GATB_TOOLS_TRAVERSAL_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>
#include <Utils.hpp>
#include <NodeSelector.hpp>
#include <set>

// some stats
struct TraversalStats
{
    long ended_traversals;
    long couldnt_find_all_consensuses;
    long couldnt_validate_consensuses;
    long couldnt_traverse_bubble_breadth;
    long couldnt_traverse_bubble_depth;
    long couldnt_because_marked_kmer;
    long couldnt_inbranching_depth;
    long couldnt_inbranching_breadth;
    long couldnt_inbranching_other;
    long couldnt_find_extension;
};


/********************************************************************************/
// semi-abstract class. implements traverse but not avance
class Traversal : public SmartPointer
{
public:

    /** */
    static Traversal* create (
        const std::string&  type,
        const Graph&        graph,
        Terminator&         terminator,
        int                 max_len,
        int                 max_depth,
        int                 max_breadth
    );

    /** */
    virtual std::string getName() const = 0;

    /** */
    virtual int traverse (const Node& node, Direction dir, Path& resulting_sequence);

    /** */
    int getMaxDepth() const  { return max_depth; };

    /** */
    int getMaxBreadth () const  { return max_breadth; };

    /** */
    static const int defaultMaxLen     = 10*1000*1000;
    static const int defaultMaxDepth   = 500;
    static const int defaultMaxBreadth = 20;

    TraversalStats final_stats, stats;
    void commit_stats() { final_stats = stats; }; // save current stats into final stats
    void revert_stats() { stats = final_stats; }; // discard changes in stats (because contig was discarded)


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

    virtual char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode) = 0;

    void mark_extensions (std::set<Node>& extensions_to_mark);

    // record the start/end positions of traversed bubbles (only from the latest traverse() call)
    std::vector <std::pair<int, int> > bubbles_positions;
};

/********************************************************************************/

class NullTraversal: public Traversal
{
public:
    /** */
    NullTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = Traversal::defaultMaxLen,
        int max_depth   = Traversal::defaultMaxDepth,
        int max_breadth = Traversal::defaultMaxBreadth
    ) : Traversal (graph, terminator, maxlen, max_depth, max_breadth) {}

    std::string getName() const  { return std::string ("null"); }

private:

    char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode) { return 0; }
};

/********************************************************************************/

class SimplePathsTraversal: public Traversal
{
public:
    /** */
    SimplePathsTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = Traversal::defaultMaxLen,
        int max_depth   = Traversal::defaultMaxDepth,
        int max_breadth = Traversal::defaultMaxBreadth
    );

    std::string getName() const  { return std::string ("unitig"); }

private:

    char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode);
};

/********************************************************************************/

class MonumentTraversal: public Traversal
{
public:
    /** */
    MonumentTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = Traversal::defaultMaxLen,
        int max_depth   = Traversal::defaultMaxDepth,
        int max_breadth = Traversal::defaultMaxBreadth
    );

    std::string getName() const  { return std::string ("monument"); }

    bool explore_branching (
        const Node& node,
        Direction dir,
        Path& consensus,
        const Node& previousNode,
        std::set<Node>& all_involved_extensions
    );

private:

    char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode);

    bool explore_branching (
        const Node& node,
        Direction dir,
        Path& consensus,
        const Node& previousNode
    );


    int find_end_of_branching (
        Direction dir,
        const Node& startingNode,
        Node& endNode,
        const Node& previousNode,
        std::set<Node>& all_involved_extensions
    );

    std::set<Path> all_consensuses_between (
        Direction    dir,
        const Node& startNode,
        const Node& endNode,
        int traversal_depth,
        std::set<Node::Value> usedNode,
        Path current_consensus,
        bool& success
    );

    std::set<Path> all_consensuses_between (
        Direction    dir,
        const Node& startNode,
        const Node& endNode,
        int traversal_depth,
        bool &success
    );

    bool validate_consensuses (std::set<Path>& consensuses, Path& consensus);

    bool all_consensuses_almost_identical (std::set<Path>& consensuses);

    void mark_extensions (std::set<Node>& extensions_to_mark);

    static const int consensuses_identity = 90; // traversing bubble if paths are all pair-wise identical by > 90%
};

/********************************************************************************/

#endif /* _GATB_TOOLS_TRAVERSAL_HPP_ */

