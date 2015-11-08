// This is legacy code
// we don't use it anymore when MPHF is enabled
// I disabled its compilation because of the need to specialize the template for GraphFast. 
// Didn't want to customize CMakeFile in the same way as gatb-core. 
// Also, I'd like to get rid of this legacy code.

#if 0

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

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/


typedef boost::variant<GraphData<32>> GraphDataVariantT;
typedef GraphTemplate<Node_t<Kmer<32>::Type>, Edge_t<Node_t<Kmer<32>::Type>>, GraphDataVariantT> GraphT;
typedef Node_t<Kmer<32>::Type> NodeT;
typedef Edge_t<Node_t<Kmer<32>::Type>> EdgeT;
typedef BranchingNode_t<Node_t<Kmer<32>::Type>> BranchingNodeT;
typedef BranchingEdge_t<Node_t<Kmer<32>::Type>, Edge_t<Node_t<Kmer<32>::Type>>> BranchingEdgeT;

#include <gatb/debruijn/impl/Instantiations.hpp> // this might be a bit exotic to do it like that.. but it works.

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/


using namespace std;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
INodeSelector<Node,Edge,GraphDataVariant>* NodeSelectorFactory<Node,Edge,GraphDataVariant>::create (const std::string& type, const GraphTemplate<Node,Edge,GraphDataVariant>& graph, TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator)
{
         if (type == "simple")  { return new NodeSelectorSimplePath<Node,Edge,GraphDataVariant> (graph, terminator); }
    else if (type == "best")    { return new NodeSelectorBest<Node,Edge,GraphDataVariant>       (graph, terminator); }
    else                        { return new NodeSelectorBest<Node,Edge,GraphDataVariant>       (graph, terminator); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool NodeSelectorSimplePath<Node,Edge,GraphDataVariant>::select (const Node& branchingNode, Node& startingNode)
{

    /** We retrieve the neighbors of the provided node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Edge> neighbors = this->_graph.template neighbors<Edge>(branchingNode.kmer);

    /** We loop these neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** make sure this kmer isnt branching */
        if (this->_terminator.is_branching (neighbors[i].to))  {  continue;  }
        /* FIXME: if we skip unitigs that are just one branching kmer, it means that we don't return ALL unitigs */

        if (this->_terminator.is_marked (neighbors[i].to))  {  continue;  }

        /** only start from an unmarked nt/strand combo */
        if (this->_terminator.is_marked (neighbors[i]))  {  continue;  }

        /** We mark the current neighbor edge. */
        this->_terminator.mark (neighbors[i]);

        size_t len_extension = 0;

        /** We loop successive node on a simple path. */
        typename GraphTemplate<Node,Edge,GraphDataVariant>::template Iterator<Node> itNodes = this->_graph.template simplePath<Node> (neighbors[i].to, neighbors[i].direction);

        for (itNodes.first(); !itNodes.isDone(); itNodes.next())
        {
	    bool only_simple_path_longer_than_2k = false; 
	    /* no reason to discard small unitigs -- we generally want them all;
             * anyhow, minia has another length filter in output
		 */

            startingNode = itNodes.item();

            if (len_extension++ > 2*this->_graph.getKmerSize() || (!only_simple_path_longer_than_2k))
            {
                /** NOTE: By convention, the returned node is understood as the forward part of the node in the
                 * bi-directional De Bruijn graph.  */
                startingNode.strand = STRAND_FORWARD;

                /** Ok, we found a starting point. */
                return true;
            }
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
template <typename Node, typename Edge, typename GraphDataVariant>
bool NodeSelectorImproved<Node,Edge,GraphDataVariant>::select (const Node& source, Node& result)
{
    /** We retrieve the neighbors of the provided node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Node> neighbors = this->_graph.template neighbors<Node>(source.kmer);

    /** We loop these neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Node& node = neighbors[i];

        // alright let's use this convention now:
        // to select new kmers: mark non-branching neighbors of branching nodes
        // to mark as used in assembly: mark branching nodes

        // only start from a non-branching k-mer
        if (this->_terminator.is_branching(node))  { continue;  }

        if (this->_terminator.is_marked(node))     { continue;  }

        this->_terminator.mark(node);

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
template <typename Node, typename Edge, typename GraphDataVariant>
bool NodeSelectorBest<Node,Edge,GraphDataVariant>::select (const Node& source, Node& result)
{
    Direction dir = DIR_OUTCOMING;
    unsigned int sum_depths = 0;

    /** First selection. */
    if (_firstSelector.select (source, result) == false)  { return false; }

    /** We use a monument traversal. */
    MonumentTraversalTemplate<Node,Edge,GraphDataVariant> monument (this->_graph, NullTerminatorTemplate<Node,Edge,GraphDataVariant>::singleton());

    /** We try each direction (outcoming and incoming). */
    Node nodes[] = { result, this->_graph.reverse(result) };

    for (size_t i=0; i<ARRAY_SIZE(nodes); i++)
    {
        Node current = nodes[i];
        Node previousNode;
        bool hasPreviousNode = false;

        // do a BFS to make sure we're not inside a bubble or tip
        FrontlineTemplate<Node,Edge,GraphDataVariant> frontline (dir, this->_graph, this->_terminator, current);

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

                if (hasPreviousNode && this->_terminator.is_branching(previousNode))
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
                    Path_t<Node> consensus;
                   
                    Node rev_node = (this->_graph.reverse(previousNode)); 
                    if (monument.explore_branching (rev_node, dir, consensus, currentNode, all_involved_extensions))
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
    if (sum_depths < (this->_graph.getKmerSize()+1))  {  return false;  }

    /** NOTE: By convention, the returned node is understood as the forward part of the node in the
     * bi-directional De Bruijn graph.  */
    result.strand = STRAND_FORWARD;

    return true;
}

// TODO make that more clean or just remove NodeFactory altogether
//
template class NodeSelectorFactory<Node_t<Kmer<32>::Type>,Edge_t<Node_t<Kmer<32>::Type >>,  boost::variant<GraphData<32>> >; 
//template class NodeSelectorFactory<Node_t<Kmer<64>::Type>,Edge_t<Node_t<Kmer<64>::Type >>,  boost::variant<GraphData<64>> >; 
//template class NodeSelectorFactory<Node_t<Kmer<96>::Type>,Edge_t<Node_t<Kmer<96>::Type >>,  boost::variant<GraphData<96>> >; 
//template class NodeSelectorFactory<Node_t<Kmer<128>::Type>,Edge_t<Node_t<Kmer<128>::Type >>,  boost::variant<GraphData<128>> >; 
//
//
#endif
