/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <Traversal.hpp>
#include <NodeSelector.hpp>

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
Traversal::Traversal (
    const Graph& graph,
    Terminator& terminator,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : graph(graph), terminator(terminator), maxlen(1000000),max_depth(500),max_breadth(20)
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
Traversal::~Traversal ()
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
bool Traversal::findStartingNode (const Node& from, Node& to)
{
    NodeSelectorSimplePath startSelector (graph, terminator);

    return startSelector.select (from, to);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int Traversal::traverse (const Node& startingNode, Direction dir, std::vector<Nucleotide>& consensus)
{
    Node currentNode = startingNode;
    Node previousNode;

    int nnt = 0;

    bool looping = false;

    Edge path[max_depth+1];

    int bubble_start, bubble_end;
    bubbles_positions.clear();

    /** We reset the consensus to be filled. */
    consensus.clear ();

    while( (nnt = avance (currentNode, dir, consensus.size() == 0, path, previousNode)))
    {
        // found branching or marked kmers
        if (nnt < 0)  {   break;  }

        if (nnt > 1)  {  bubble_start = consensus.size();  }

        // keep re-walking the nucleotides we just discovered, to append to consensus and mark kmers as we go
        for (int i = 0; i < nnt; i++)
        {
            /** We add the current nucleotide into the contig. */
            consensus.push_back (path[i].nt);

            /** We update previous and current nodes. */
            previousNode = currentNode;
            currentNode  = path[i].to;

            /** We mark the node as used in assembly. */
            terminator.mark (currentNode);

            /** perfectly circular regions with no large branching can happen (rarely) */
            if (currentNode.kmer == startingNode.kmer)  {  looping = true;  }
        }

        if (nnt > 1)
        {
            bubble_end = consensus.size();
            bubbles_positions.push_back(std::make_pair(bubble_start,bubble_end));
        }

        if (looping)  {  break;  }

        if (consensus.size() > maxlen)  {  break;  }

    }  /* end of while( (nnt = avance ... */

    return consensus.size();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
SimplePathsTraversal::SimplePathsTraversal (
    const Graph& graph,
    Terminator& terminator,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : Traversal (graph, terminator, maxlen, max_depth, max_breadth)
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
char SimplePathsTraversal::avance (
    const Node& node,
    Direction dir,
    bool first_extension,
    Edge* path,
    const Node& previousNode
)
{
    return max (GraphHelper(graph).simplePathAvance (node, dir, path[0]),  0);
}
