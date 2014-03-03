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

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <NodeSelector.hpp>
#include <Frontline.hpp>
#include <Traversal.hpp>

using namespace std;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
INodeSelector* NodeSelectorFactory::create (const std::string& type, const Graph& graph, Terminator& terminator)
{
         if (type == "simple")  { return new NodeSelectorSimplePath (graph, terminator); }
    else if (type == "best")    { return new NodeSelectorBest       (graph, terminator); }
    else                        { return new NodeSelectorBest       (graph, terminator); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool NodeSelectorSimplePath::select (const Node& branchingNode, Node& startingNode)
{

    /** We retrieve the neighbors of the provided node. */
    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(branchingNode.kmer);

    /** We loop these neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** make sure this kmer isnt branching */
        if (_terminator.is_branching (neighbors[i].to))  {  continue;  }

        if (_terminator.is_marked (neighbors[i].to))  {  continue;  }

        /** only start from an unmarked nt/strand combo */
        if (_terminator.is_marked (neighbors[i]))  {  continue;  }

        /** We mark the current neighbor edge. */
        _terminator.mark (neighbors[i]);

        size_t len_extension = 0;

        /** We loop successive node on a simple path. */
        Graph::Iterator<Node> itNodes = _graph.simplePath<Node> (neighbors[i].to, neighbors[i].direction);

        for (itNodes.first(); !itNodes.isDone(); itNodes.next())
        {
            if (len_extension++ > 2*_graph.getKmerSize())
            {
                /** NOTE: By convention, the returned node is understood as the forward part of the node in the
                 * bi-directional De Bruijn graph.  */
                startingNode.strand = STRAND_FORWARD;

                /** Ok, we found a starting point. */
                return true;
            }

            startingNode = itNodes.item();
        }
    }

    return false;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool NodeSelectorImproved::select (const Node& source, Node& result)
{
    /** We retrieve the neighbors of the provided node. */
    Graph::Vector<Node> neighbors = _graph.neighbors<Node>(source.kmer);

    /** We loop these neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Node& node = neighbors[i];

        // alright let's use this convention now:
        // to select new kmers: mark non-branching neighbors of branching nodes
        // to mark as used in assembly: mark branching nodes

        // only start from a non-branching k-mer
        if (_terminator.is_branching(node))  { continue;  }

        if (_terminator.is_marked(node))     { continue;  }

        _terminator.mark(node);

        result = node;

        /** NOTE: By convention, the returned node is understood as the forward part of the node in the
         * bi-directional De Bruijn graph.  */
        result.strand = STRAND_FORWARD;

        return true;
    }

    return false;
}

/*********************************************************************
** METHOD  :
** PURPOSE : This is the NodeSelector to use when you're doing contigs traversal -- it avoids starting inside a bubble
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool NodeSelectorBest::select (const Node& source, Node& result)
{
    Direction dir = DIR_OUTCOMING;
    int sum_depths = 0;

    /** First selection. */
    if (_firstSelector.select (source, result) == false)  { return false; }

    /** We use a monument traversal. */
    MonumentTraversal monument (_graph, NullTerminator::singleton());

    /** We try each direction (outcoming and incoming). */
    Node nodes[] = { result, _graph.reverse(result) };

    for (size_t i=0; i<ARRAY_SIZE(nodes); i++)
    {
        Node current = nodes[i];
        Node previousNode;
        bool hasPreviousNode = false;

        // do a BFS to make sure we're not inside a bubble or tip
        Frontline frontline (dir, _graph, _terminator, current);

        do
        {
            bool should_continue = frontline.go_next_depth();
            if (!should_continue)  {  break;  }

            // put the same contraints as in a bubble
            if (frontline.depth() > monument.getMaxDepth() || frontline.size() > monument.getMaxBreadth())  {  break;  }

            // stopping condition: nothing more to explore
            if (frontline.size() == 0)  {  break;  }

            if (frontline.size() <= 1)
            {
                Node currentNode;
                if (frontline.size() == 1)  {  currentNode = frontline.front().node;  }

                if (hasPreviousNode && _terminator.is_branching(previousNode))
                {
                    /* the current situation is:
                     *
                     *    current_kmer  previous_kmer
                     *   -O-------------O------------- ... ---starting_kmer
                     *                  \_....
                     *
                     * or
                     *
                     *   [no extension] previous_kmer
                         *   X              O------------- ... ---starting_kmer
                     *                  \_....
                     *
                     *   so, by looking one k-mer ahead, we make sure that previous_kmer only branches to the right
                     *
                     */
                    set<Node> all_involved_extensions;
                    Path consensus;
                    
                    if (monument.explore_branching (_graph.reverse(previousNode), dir, consensus, currentNode, all_involved_extensions))
                    {
                        if (all_involved_extensions.find(current) != all_involved_extensions.end())
                        {
                            return false; // starting_kmer is in a tip/bubble starting from current_kmer
                        }
                    }
                }
            }

            // update previous_kmer
            if (frontline.size() == 1)
            {
                previousNode    = frontline.front().node;
                hasPreviousNode = true;
            }
            else
            {
                hasPreviousNode = false;
            }
        }
        while (1);

        sum_depths += frontline.depth();
    }

    // don't even assemble those regions which have no chance to yield a long contig
    if (sum_depths < (_graph.getKmerSize()+1))  {  return false;  }

    /** NOTE: By convention, the returned node is understood as the forward part of the node in the
     * bi-directional De Bruijn graph.  */
    result.strand = STRAND_FORWARD;

    return true;
}
