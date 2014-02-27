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

#include <Traversal.hpp>
#include <NodeSelector.hpp>
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
Traversal* Traversal::create (
    const std::string&  type,
    const Graph&        graph,
    Terminator&         terminator,
    int                 max_len,
    int                 max_depth,
    int                 max_breadth
)
{
    Traversal* result = 0;

         if (type == "unitig")    { result = new SimplePathsTraversal (graph, terminator, max_len, max_depth, max_breadth); }
    else if (type == "monument")  { result = new MonumentTraversal    (graph, terminator, max_len, max_depth, max_breadth); }
    else if (type == "null")      { result = new NullTraversal        (graph, terminator, max_len, max_depth, max_breadth); }
    else                          { result = new MonumentTraversal    (graph, terminator, max_len, max_depth, max_breadth); }

    return result;
}

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
    int max_len,
    int max_depth,
    int max_breadth
)
    : graph(graph), terminator(terminator), stats(TraversalStats()), final_stats(TraversalStats()),
      maxlen      (max_len     == 0 ? Traversal::defaultMaxLen     : max_len),
      max_depth   (max_depth   == 0 ? Traversal::defaultMaxDepth   : max_depth),
      max_breadth (max_breadth == 0 ? Traversal::defaultMaxBreadth : max_breadth)
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
int Traversal::traverse (const Node& startingNode, Direction dir, Path& consensus)
{
    Node currentNode = startingNode;
    Node previousNode;

    int nnt = 0;

    bool looping = false;

    Path path;  path.resize (max_depth+1);

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
            consensus.push_back (path[i]);

            /** We update previous and current nodes. */
            previousNode = currentNode;

            /** We compute the neighbor node for the current nucleotide of the path.
             * WARNING: we 'build' here a node without checking that it belongs to the Debruijn graph, because
             * the transition nucleotide has been got through a call to Graph::neighbors, and therefore it should
             * be a trustable transition nucleotide. */
            currentNode  = graph.neighbor<Node> (currentNode, dir, path[i]);

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
            #     #  #     #  ###  #######  ###   #####
            #     #  ##    #   #      #      #   #     #
            #     #  # #   #   #      #      #   #
            #     #  #  #  #   #      #      #   #  ####
            #     #  #   # #   #      #      #   #     #
            #     #  #    ##   #      #      #   #     #
             #####   #     #  ###     #     ###   #####
*********************************************************************/

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
    Path& path,
    const Node& previousNode
)
{
    return  max (graph.simplePathAvance (node, dir, path[0]),  0);
}

/*********************************************************************
#     #  #######  #     #  #     #  #     #  #######  #     #  #######
##   ##  #     #  ##    #  #     #  ##   ##  #        ##    #     #
# # # #  #     #  # #   #  #     #  # # # #  #        # #   #     #
#  #  #  #     #  #  #  #  #     #  #  #  #  #####    #  #  #     #
#     #  #     #  #   # #  #     #  #     #  #        #   # #     #
#     #  #     #  #    ##  #     #  #     #  #        #    ##     #
#     #  #######  #     #   #####   #     #  #######  #     #     #
*********************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
MonumentTraversal::MonumentTraversal (
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
char MonumentTraversal::avance (
    const Node& node,
    Direction dir,
    bool first_extension,
    Path& consensus,
    const Node& previousNode
)
{
    // if we're on a simple path, just traverse it
    int is_simple_path = graph.simplePathAvance (node, dir, consensus[0]);
    if (is_simple_path > 0)  {  return 1;  }

    // the following function does:
    // * a bfs from the starting kmer, stopping when:
    //      - breadth > max_breadth
    //      or
    //      - depth > max_depth
    // * check if there a single end point
    // * computing all possible paths between start and end
    // * returns one flattened consensus sequence
    bool success = explore_branching (node, dir, consensus, previousNode);
    if (!success)
    { 
        stats.ended_traversals++;
        return 0;
    }

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
bool MonumentTraversal::explore_branching (
    const Node& node,
    Direction dir,
    Path& consensus,
    const Node& previousNode
)
{
    set<Node> all_involved_extensions;

    return explore_branching (node, dir, consensus, previousNode, all_involved_extensions);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool MonumentTraversal::explore_branching (
    const Node& startNode,
    Direction dir,
    Path& consensus,
    const Node& previousNode,
    std::set<Node>& all_involved_extensions
)
{
    Node endNode;

    // find end of branching, record all involved extensions (for future marking)
    // it returns false iff it's a complex bubble
    int traversal_depth = find_end_of_branching (dir, startNode, endNode, previousNode, all_involved_extensions);
    if (!traversal_depth)  
    {
        stats.couldnt_find_all_consensuses++;
        return false;
    }

    // find all consensuses between start node and end node
    bool success;
    set<Path> consensuses = all_consensuses_between (dir, startNode, endNode, traversal_depth+1, success);

    // if consensus phase failed, stop
    if (!success)  {  return false;  }

    consensus.resize (0);
    // validate paths, based on identity
    bool validated = validate_consensuses (consensuses, consensus);
    if (!validated)   
    {  
        stats.couldnt_validate_consensuses++;
        return false;  
    }

    // the consensuses agree, mark all the involved extensions
    // (corresponding to alternative paths we will never traverse again)
    mark_extensions (all_involved_extensions);

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
int MonumentTraversal::find_end_of_branching (
    Direction    dir,
    const Node&  startingNode,
    Node&        endNode,
    const Node&  previousNode,
    std::set<Node>& all_involved_extensions
)
{
    /** We need a branching frontline. */
    FrontlineBranching frontline (dir, graph, terminator, startingNode, previousNode, &all_involved_extensions);

    do  {
        bool should_continue = frontline.go_next_depth();
        if (!should_continue) 
        {
            if (frontline.stopped_reason == Frontline::MARKED)
                stats.couldnt_because_marked_kmer++;
            if (frontline.stopped_reason == Frontline::IN_BRANCHING_DEPTH)
                stats.couldnt_inbranching_depth++;
            if (frontline.stopped_reason == Frontline::IN_BRANCHING_BREADTH)
                stats.couldnt_inbranching_breadth++;
            if (frontline.stopped_reason == Frontline::IN_BRANCHING_OTHER)
                stats.couldnt_inbranching_other++;
            return 0;
        }

        // don't allow a depth too large
        if (frontline.depth() > max_depth)  
        {  
            stats.couldnt_traverse_bubble_depth++;
            return 0;  
        }

        // don't allow a breadth too large
        if (frontline.size()> max_breadth)  
        {  
            stats.couldnt_traverse_bubble_breadth++;
            return 0;
        }

        // stopping condition: frontline is either empty, or contains only 1 kmer
        // needs the kmer to be non-branching, in order to avoid a special case of bubble immediatly after a bubble
        // affects mismatch rate in ecoli greatly
        if (frontline.size() == 0)  
        {
            stats.couldnt_find_extension++;
            return 0;
        }

        // if (frontline.size() == 1) // longer contigs but for some reason, higher mismatch rate
        if (frontline.size() == 1 &&   (!terminator.isEnabled() || !terminator.is_branching(frontline.front().node)) )  {   break;  }
    }
    while (1);

    if (frontline.size()==1)
    {
        endNode = frontline.front().node;
        return frontline.depth();
    }

   return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void MonumentTraversal::mark_extensions (std::set<Node>& extensions_to_mark)
{
    if (terminator.isEnabled())
    {
        for(set<Node>::iterator it = extensions_to_mark.begin(); it != extensions_to_mark.end() ; ++it)
        {
            terminator.mark (*it);
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
set<Path> MonumentTraversal::all_consensuses_between (
    Direction    dir,
    const Node& startNode,
    const Node& endNode,
    int traversal_depth,
    set<Node::Value> usedNode,
    Path current_consensus,
    bool& success
)
{
    set<Path> consensuses;

    // find_end_of_branching and all_consensues_between do not always agree on clean bubbles ends
    // until I can fix the problem, here is a fix
    // to reproduce the problem: SRR001665.fasta 21 4
    if (traversal_depth < -1)
    {
        success = false;
        return consensuses;
    }

    if (startNode.kmer == endNode.kmer)// not testing for end_strand anymore because find_end_of_branching doesn't care about strands
    {
        consensuses.insert(current_consensus);
        return consensuses;
    }

    /** We retrieve the neighbors of the provided node. */
    Graph::Vector<Edge> neighbors = graph.neighbors<Edge> (startNode, dir);

    /** We loop these neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& edge = neighbors[i];

        // don't resolve bubbles containing loops
        // (tandem repeats make things more complicated)
        // that's a job for a gapfiller
        if (usedNode.find(edge.to.kmer) != usedNode.end())
        {
            success = false;
            return consensuses;
        }

        // generate extended consensus sequence
        Path extended_consensus(current_consensus);
        extended_consensus.push_back (edge.nt);

        // generate list of used kmers (to prevent loops)
        set<Node::Value> extended_kmers (usedNode);
        extended_kmers.insert (edge.to.kmer);

        // recursive call to all_consensuses_between
        set<Path> new_consensuses = all_consensuses_between (
            dir,
            edge.to,
            endNode,
            traversal_depth - 1,
            extended_kmers,
            extended_consensus,
            success
        );

        consensuses.insert (new_consensuses.begin(), new_consensuses.end());

        // mark to stop we end up with too many consensuses
        if (consensuses.size() > (unsigned int )max_breadth)  {    success = false;  }

        // propagate the stop if too many consensuses reached
        if (success == false)  {   return consensuses;  }
    }

    return consensuses;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
set<Path> MonumentTraversal::all_consensuses_between (
    Direction    dir,
    const Node& startNode,
    const Node& endNode,
    int traversal_depth,
    bool &success
)
{
    set<Node::Value> usedNode;
    usedNode.insert(startNode.kmer);
    Path current_consensus;
    success = true;

    return all_consensuses_between (dir, startNode, endNode, traversal_depth, usedNode, current_consensus, success);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool MonumentTraversal::validate_consensuses (set<Path>& consensuses, Path& result)
{
    bool debug = false;
    // compute mean and stdev of consensuses
    int mean = 0;
    int path_number = 0;
    for(set<Path>::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        //if (debug)  printf("bubble path %d: %s (len=%lu)\n",path_number,(*it).c_str(),(*it).length());
        mean+=(*it).size();
        path_number++;
    }
    mean/=consensuses.size();
    double stdev = 0;
    for(set<Path>::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        int consensus_length = (*it).size();
        stdev += pow(fabs(consensus_length-mean),2);
    }
    stdev = sqrt(stdev/consensuses.size());

    // don't traverse large bubbles
    if (mean > max_depth)
        return false;

    // don't traverse large deadends (here, having one consensus means the other paths were large deadends)
    if (consensuses.size() == 1 && mean > graph.getKmerSize()+1) // deadend length should be < k+1 (most have length 1, but have seen up to 10 in ecoli)
        return false;

    if (debug) printf("%lu-bubble mean %d, stdev %.1f\n",consensuses.size(),mean,stdev);

    // traverse bubbles if paths have roughly the same length
    if (stdev>mean/5)
        return false;

    // check that all consensuses are similar
    if (! all_consensuses_almost_identical(consensuses))
        return false;

    // if all good, an arbitrary consensus is chosen
    Path chosen_consensus = *consensuses.begin();
    int result_length = chosen_consensus.size();
    if  (result_length> max_depth) // it can happen that consensus is longer than max_depth, despite that we didn't explore that far (in a messy bubble with branchings inside)
        return false;

    /** We the the result consensus. */
    result = chosen_consensus;

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
bool MonumentTraversal::all_consensuses_almost_identical (set<Path>& consensuses)
{
    for (set<Path>::iterator it_a = consensuses.begin(); it_a != consensuses.end(); it_a++)
    {
        set<Path>::iterator it_b = it_a;
        advance(it_b,1);
        while (it_b != consensuses.end())
        {
            if (needleman_wunch(*it_a,*it_b) * 100 < consensuses_identity)
                return false;
            advance(it_b,1);
        }
    }
    return true;
}
