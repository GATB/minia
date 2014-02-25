/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <Frontline.hpp>

using namespace std;

/********************************************************************************/

#define DEBUG(a)   //a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// a frontline is a set of nodes having equal depth in the BFS
Frontline::Frontline (
    Direction         direction,
    const Graph&      graph,
    Terminator&       terminator,
    const Node&       startingNode
) :
    _direction(direction), _graph(graph), _terminator(terminator), _depth(0),
    _all_involved_extensions(0)
{
    _already_frontlined.insert (startingNode.kmer);

    _frontline.push (NodeNt (startingNode, NUCL_UNKNOWN));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// a frontline is a set of nodes having equal depth in the BFS
Frontline::Frontline (
    Direction         direction,
    const Graph&      graph,
    Terminator&       terminator,
    const Node&       startingNode,
    const Node&       previousNode,
    std::set<Node>*   all_involved_extensions
) :
    _direction(direction), _graph(graph), _terminator(terminator), _depth(0),
    _all_involved_extensions(all_involved_extensions)
{
    _already_frontlined.insert (startingNode.kmer);
    _already_frontlined.insert (previousNode.kmer);

    _frontline.push (NodeNt (startingNode, NUCL_UNKNOWN));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Frontline::go_next_depth()
{
    // extend all nodes in this frontline simultaneously, creating a new frontline
    stopped_reason=NONE;
    queue_nodes new_frontline;

    while (!_frontline.empty())
    {
        /** We get the first item of the queue and remove it from the queue. */
        NodeNt current_node = _frontline.front();
        _frontline.pop();

        /** We check whether we use this node or not. */
        if (check(current_node.node) == false)  { return false; }

        /** We loop the neighbors edges of the current node. */
        Graph::Vector<Edge> edges = _graph.neighbors<Edge> (current_node.node, _direction);

        for (size_t i=0; i<edges.size(); i++)
        {
            /** Shortcuts. */
            const Edge& edge     = edges[i];
            const Node& neighbor = edge.to;

            // test if that node hasn't already been explored
            if (_already_frontlined.find (neighbor.kmer) != _already_frontlined.end())  { continue; }

            // if this bubble contains a marked (branching) kmer, stop everyone at once (to avoid redundancy)
            if (_terminator.isEnabled() && _terminator.is_branching (neighbor) &&  _terminator.is_marked_branching(neighbor))  
            {  
                stopped_reason=Frontline::MARKED;
                return false;  
            }

            // propagate information where this node comes from
            Nucleotide from_nt = (current_node.nt == NUCL_UNKNOWN) ? edge.nt : current_node.nt;

            /** We add the new node to the new front line. */
            new_frontline.push (NodeNt (neighbor, from_nt));

            /** We memorize the new node. */
            _already_frontlined.insert (neighbor.kmer);

            // since this extension is validated, insert into the list of involved ones
            if (_all_involved_extensions != 0)  {  _all_involved_extensions->insert (neighbor);  }
        }
    }

    _frontline = new_frontline;
    ++_depth;

    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
FrontlineBranching::FrontlineBranching (
    Direction         direction,
    const Graph&      graph,
    Terminator&       terminator,
    const Node&       startingNode,
    const Node&       previousNode,
    std::set<Node>*   all_involved_extensions
)  : Frontline (direction,graph,terminator,startingNode,previousNode,all_involved_extensions)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
FrontlineBranching::FrontlineBranching (
    Direction         direction,
    const Graph&      graph,
    Terminator&       terminator,
    const Node&       startingNode
) : Frontline(direction,graph,terminator,startingNode)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// new code, not in monument, to detect any in-branching longer than 3k
bool FrontlineBranching::check (const Node& node)
{
	/** We reverse the node for the inbranching path. */
    Node actual = _graph.reverse(node);

    /** We loop the neighbors nodes of the current node. */
    Graph::Vector<Node> neighbors = _graph.neighbors<Node> (actual, (_direction));

    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Node& neighbor = neighbors[i];

        // only check in-branching from kmers not already frontlined
        // which, for the first extension, includes the previously traversed kmer (previous_kmer)
        // btw due to avance() invariant, previous_kmer is always within a simple path
        if (_already_frontlined.find (neighbor.kmer) != _already_frontlined.end())  {   continue;  }

        // create a new frontline inside this frontline to check for large in-branching (i know, we need to go deeper, etc..)
        Frontline frontline (_direction, _graph, _terminator, neighbor, actual, _all_involved_extensions);

        do  {
            bool should_continue = frontline.go_next_depth();

            if (!should_continue)  
            {  
                stopped_reason=Frontline::IN_BRANCHING_OTHER;
                break;
            }

            // don't allow a depth > 3k
            if (frontline.depth() > 3 * _graph.getKmerSize())  
            {  
                stopped_reason=Frontline::IN_BRANCHING_DEPTH;
                break;
            }

            // don't allow a breadth too large
            if (frontline.size() > 10)  
            {  
                stopped_reason=Frontline::IN_BRANCHING_BREADTH;
                break;
            }

            // stopping condition: no more in-branching
            if (frontline.size() == 0)  {  break;  }
        }
        while (1);

        // found large in-branching
        if (frontline.size() > 0)  {  return false;  }
    }

    // didn't find any in-branching
    return true;
}
