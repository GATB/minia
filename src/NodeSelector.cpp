/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <NodeSelector.hpp>

#undef NDEBUG
#include <cassert>

using namespace std;

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
    /** We try each direction (outcoming and incoming). */
    foreach_direction (dir)
    {
        /** We retrieve the neighbors of the provided node. */
        Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(branchingNode, dir);

        /** We loop these neighbors. */
        for (size_t i=0; i<neighbors.size(); i++)
        {
            if (_terminator.is_marked (neighbors[i].to))  {  continue;  }

            /** only start from an unmarked nt/strand combo */
            if (_terminator.is_marked (neighbors[i]))  { continue;  }

            /** We mark the current neighbor edge. */
            _terminator.mark (neighbors[i]);

            /** make sure this kmer isnt branching */
            if (_terminator.is_branching (neighbors[i].to))  { continue; }

            size_t len_extension = 0;

            /** We loop successive node on a simple path. */
            Graph::Iterator<Node> itNodes = GraphHelper(_graph).simplePathIterator<Node> (neighbors[i].to, dir);

            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                if (len_extension++ > 2*_graph.getKmerSize())
                {
                    startingNode.strand = STRAND_FORWARD;
                    return true;
                }

                startingNode = itNodes.item();
            }
        }
    }

    return false;
}
