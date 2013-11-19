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
#include <Terminator.hpp>

/********************************************************************************/
// semi-abstract class. implements traverse but not avance
class Traversal
{
public:

    /** */
    Traversal (const Graph& graph, Terminator& terminator, int maxlen, int max_depth, int max_breadth);

    /** */
    virtual ~Traversal();

    /** */
    virtual std::string getName() const = 0;

    /** */
    bool findStartingNode (const Node& from, Node& to);

    /** */
    int traverse (const Node& node, Direction dir, std::vector<Nucleotide>& resulting_sequence);

protected:
    const Graph& graph;
    Terminator&  terminator;

    int maxlen;
    int max_depth;
    int max_breadth;

    virtual char avance (const Node& node, Direction dir, bool first_extension, Edge* path, const Node& previousNode) = 0;

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

    std::string getName() const  { return std::string ("SimplePath"); }

private:

    char avance (const Node& node, Direction dir, bool first_extension, Edge* path, const Node& previousNode);
};

/********************************************************************************/


#endif /* _GATB_TOOLS_TRAVERSAL_HPP_ */

