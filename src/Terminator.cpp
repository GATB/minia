/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <Terminator.hpp>

#include <cassert>

using namespace std;

/********************************************************************************/

#define DEBUG(a)   //a

extern u_int64_t incomingTable[];
extern bool hack;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BranchingTerminator::BranchingTerminator (const Graph& graph)
    : Terminator (graph)
{
    /** We loop over the branching nodes. */
    Graph::Iterator<BranchingNode> itBranching = _graph.iterator<BranchingNode>();
    for (itBranching.first(); !itBranching.isDone(); itBranching.next())
    {
        /** We add the current branching node into the map. */
        branching_kmers.insert (itBranching.item().kmer);
    }

    /** We finalize the map. */
    branching_kmers.finalize (false);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BranchingTerminator::~BranchingTerminator()
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
bool BranchingTerminator::is_branching (const Node& node) const
{
    return _graph.isBranching (node);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BranchingTerminator::mark (const Edge& edge)
{
    // BranchingTerminator ignores non-branching kmers
    if (!is_indexed (edge.from))  {   return;  }

    Value val=0;
    branching_kmers.get (edge.from.kmer, val);

    Nucleotide nt = edge.nt;

	if (hack)
	{
	    if (edge.direction == DIR_INCOMING)
	    {
	        size_t span = _graph.getKmerSize();
	
	        if (edge.to.strand == STRAND_FORWARD)
	        {
	            nt = ((Nucleotide) (edge.to.kmer[span-1]));
	            nt = (Nucleotide) incomingTable[nt];
	        }
	        else
	        {
	            nt = reverse ((Nucleotide) (edge.to.kmer[0]));
	            nt = (Nucleotide) incomingTable[nt];
	        }
	    }
	}

    int delta = 0;

    // set a 1 at the right NT & strand position
         if (edge.direction == DIR_OUTCOMING && edge.from.strand == STRAND_FORWARD)  { delta = 0; }
    else if (edge.direction == DIR_OUTCOMING && edge.from.strand == STRAND_REVCOMP)  { delta = 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == STRAND_FORWARD)  { delta = 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == STRAND_REVCOMP)  { delta = 0; }
    else { throw "ERROR..."; }

    val |= 1 << (nt+delta);

    branching_kmers.set (edge.from.kmer,val); //was insert for Hash16

    assert (is_marked(edge) == true);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BranchingTerminator::is_marked (const Edge& edge)  const
{
    Value val = 0;
    int is_present = branching_kmers.get (edge.from.kmer, val);

    if (!is_present)  {   return false;  }

    int extension_nucleotide_marked;

    Nucleotide nt = edge.nt;

    if (hack)
    {
        size_t span = _graph.getKmerSize();
        if (edge.direction == DIR_INCOMING)
        {
            if (edge.to.strand == STRAND_FORWARD)
            {
                nt = ((Nucleotide) (edge.to.kmer[span-1]));
                nt = (Nucleotide) incomingTable[nt];
            }
            else
            {
                nt = reverse ((Nucleotide) (edge.to.kmer[0]));
                nt = (Nucleotide) incomingTable[nt];
            }
        }
    }

    int delta = 0;
         if (edge.direction == DIR_OUTCOMING && edge.from.strand == STRAND_FORWARD)  { delta = 0; }
    else if (edge.direction == DIR_OUTCOMING && edge.from.strand == STRAND_REVCOMP)  { delta = 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == STRAND_FORWARD)  { delta = 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == STRAND_REVCOMP)  { delta = 0; }
    else { throw "ERROR..."; }

    extension_nucleotide_marked = (val>>(nt+delta))&1;

    return  extension_nucleotide_marked == 1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BranchingTerminator::mark (const Node& node)
{
    bool could_mark = false;

    // if it is a branching kmer, mark it directly (it may have no branching neighbor)
    if (is_indexed(node))
    {
        Value val;
        branching_kmers.get (node.kmer, val);
        branching_kmers.set (node.kmer, val|(1<<8));
        could_mark = true;
    }

    /** We try each direction (outcoming and incoming). */
    foreach_direction (dir)
    {
        /** We loop the neighbors edges of the current node. */
        Graph::Vector<Edge> neighbors = _graph.neighbors<Edge> (node, dir);

        /** We loop the branching neighbors. */
        for (size_t i=0; i<neighbors.size(); i++)
        {
            if (is_indexed(neighbors[i].to)==false)  { continue; }

            /** We mark this edge (reversed first, in order to have the branching as the 'from' node) */
            mark (neighbors[i].reverse());
            could_mark = true;
        }
    }

    if (could_mark)  {   assert(is_marked(actualNode) == true);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BranchingTerminator::is_marked (const Node& node) const
{
    Node actualNode = node; //  actualNode.strand = STRAND_FORWARD;

    // if it is a branching kmer, read marking directly (it may have no branching neighbor)
    if (is_indexed(actualNode))  {   return is_marked_branching(actualNode);  }

    /** We try each direction (outcoming and incoming). */
    foreach_direction (dir)
    {
        /** We loop the neighbors edges of the current node. */
        Graph::Vector<Edge> neighbors = _graph.neighbors<Edge> (actualNode, dir);

        for (size_t i=0; i<neighbors.size() && is_indexed(neighbors[i].to); i++)
        {
            if (is_marked (neighbors[i].reverse()) == true)  {  return true;  }
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
bool BranchingTerminator::is_marked_branching (const Node& node) const
{
    Value val;
    branching_kmers.get (node.kmer, val);
    return (val&(1<<8)) != 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BranchingTerminator::reset()
{
    branching_kmers.clear();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BranchingTerminator::is_indexed (const Node& node) const
{
    return branching_kmers.contains (node.kmer);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BranchingTerminator::dump ()
{
    Value checksum = 0;

    std::vector<Node::Type>& liste = branching_kmers.liste ;

    for (size_t i=0; i<liste.size(); i++)
    {
        Value v;  if (branching_kmers.get (liste[i],v))  { checksum += v; }
    }
    cout << "TERMINATOR:  nb=" << liste.size() << "  checksum=" << checksum << endl;
}

