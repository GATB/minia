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

#include <Terminator.hpp>

#include <cassert>

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
int BranchingTerminator::getDelta (const Edge& edge) const
{
         if (edge.direction == DIR_OUTCOMING && edge.from.strand == STRAND_FORWARD)  { return 0; }
    else if (edge.direction == DIR_OUTCOMING && edge.from.strand == STRAND_REVCOMP)  { return 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == STRAND_FORWARD)  { return 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == STRAND_REVCOMP)  { return 0; }
    else { return -1; }
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

    int delta = getDelta (edge);
    if (delta >= 0)
    {
        // set a 1 at the right NT & strand position
        val |= 1 << (edge.nt + delta);

        branching_kmers.set (edge.from.kmer,val); //was insert for Hash16
    }

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

    int delta = getDelta (edge);
    if (delta >= 0)
    {
        // set a 1 at the right NT & strand position
        extension_nucleotide_marked = (val>>(edge.nt+delta))&1;
    }

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

    /** We loop the neighbors edges of the current node. */
    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge> (node.kmer);

    /** We loop the branching neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& e = neighbors[i];

        if (is_indexed(e.to)==false)  { continue; }

        /** We mark this edge (reversed first, in order to have the branching as the 'from' node) */
        mark (_graph.reverse(e));

        could_mark = true;
    }

    if (could_mark)  {   assert(is_marked(node) == true);  }
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
    // if it is a branching kmer, read marking directly (it may have no branching neighbor)
    if (is_indexed(node))  {   return is_marked_branching(node);  }

    /** We loop the neighbors edges of the current node. */
    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge> (node.kmer);

    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& e = neighbors[i];

        if  (is_indexed(e.to)==false)  { continue; }

        if (is_marked (_graph.reverse(e)) == true)  {  return true;  }
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

    std::vector<Node::Value>& liste = branching_kmers.liste ;

    for (size_t i=0; i<liste.size(); i++)
    {
        Value v;  if (branching_kmers.get (liste[i],v))  { checksum += v; }
    }
    cout << "TERMINATOR:  nb=" << liste.size() << "  checksum=" << checksum << endl;
}

