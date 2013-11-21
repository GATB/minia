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
#include <Frontline.hpp>

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
Traversal* Traversal::create (
    const std::string&  type,
    const Graph&        graph,
    Terminator&         terminator,
    INodeSelector*      selector
)
{
    Traversal* result = 0;

         if (type == "unitig")    { result = new SimplePathsTraversal (graph, terminator, selector); }
    else if (type == "monument")  { result = new MonumentTraversal    (graph, terminator, selector); }
    else                          { result = new MonumentTraversal    (graph, terminator, selector); }

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
    INodeSelector* selector,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : graph(graph), terminator(terminator), _selector(0),
      maxlen(1000000),max_depth(500),max_breadth(20)
{
    setSelector (selector);
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
    setSelector (0);
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
    return _selector->select (from, to);
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

    PATH path;  path.resize (max_depth+1);

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
    INodeSelector* selector,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : Traversal (graph, terminator, selector, maxlen, max_depth, max_breadth)
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
    PATH& path,
    const Node& previousNode
)
{
    return  max (GraphHelper(graph).simplePathAvance (node, dir, path[0]),  0);
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
    INodeSelector* selector,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : Traversal (graph, terminator, selector, maxlen, max_depth, max_breadth)
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
    PATH& consensus,
    const Node& previousNode
)
{
    // if we're on a simple path, just traverse it
    int is_simple_path = GraphHelper(graph).simplePathAvance (node, dir, consensus[0]);
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
    if (!success)  {  return 0; }

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
    PATH& consensus,
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
    PATH& consensus,
    const Node& previousNode,
    std::set<Node>& all_involved_extensions
)
{
    Node endNode;

    // find end of branching, record all involved extensions (for future marking)
    // it returns false iff it's a complex bubble
    int traversal_depth = find_end_of_branching (dir, startNode, endNode, previousNode, all_involved_extensions);
    if (!traversal_depth)  {   return false;  }

    // find all consensuses between start node and end node
    set<PATH> consensuses;
    bool success;
    consensuses = all_consensuses_between (dir, startNode, endNode, traversal_depth+1, success);

    // if consensus phase failed, stop
    if (!success)  {  return false;  }

    consensus.resize (0);
    // validate paths, based on identity
    bool validated = validate_consensuses (consensuses, consensus);
    if (!validated)   {  return false;  }

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
        if (!should_continue) {  return 0;  }

        // don't allow a depth too large
        if (frontline.depth() > max_depth)  {  return 0;  }

        // don't allow a breadth too large
        if (frontline.size()> max_breadth)  {  return 0;  }

        // stopping condition: frontline is either empty, or contains only 1 kmer
        // needs the kmer to be non-branching, in order to avoid a special case of bubble immediatly after a bubble
        // affects mismatch rate in ecoli greatly
        if (frontline.size() == 0)  {  return 0;  }

        // if (frontline.size() == 1) // longer contigs but for some reason, higher mismatch rate
        if (frontline.size() == 1 &&   !terminator.is_branching(frontline.front().node) )  {   break;  }
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
    for(set<Node>::iterator it = extensions_to_mark.begin(); it != extensions_to_mark.end() ; ++it)
    {
        terminator.mark (*it);
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
set<PATH> MonumentTraversal::all_consensuses_between (
    Direction    dir,
    const Node& startNode,
    const Node& endNode,
    int traversal_depth,
    set<Node> usedNode,
    PATH current_consensus,
    bool& success
)
{
    set<PATH> consensuses;

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
        if (usedNode.find(edge.to) != usedNode.end())
        {
            success = false;
            return consensuses;
        }

        Nucleotide& NT = edge.nt;

        if (hack)
        {
            if (edge.direction == DIR_INCOMING)
            {
                size_t span = graph.getKmerSize();

                if (edge.to.strand == STRAND_FORWARD)
                {
                    NT = ((Nucleotide) (edge.to.kmer[span-1]));
                    NT = (Nucleotide) incomingTable[NT];
                }
                else
                {
                    NT = reverse ((Nucleotide) (edge.to.kmer[0]));
                    NT = (Nucleotide) incomingTable[NT];
                }
            }
        }

        // generate extended consensus sequence
        PATH extended_consensus(current_consensus);
        extended_consensus.push_back (edge);

        // generate list of used kmers (to prevent loops)
        set<Node> extended_kmers (usedNode);
        extended_kmers.insert (edge.to);

        // recursive call to all_consensuses_between
        set<PATH> new_consensuses = all_consensuses_between (
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
set<PATH> MonumentTraversal::all_consensuses_between (
    Direction    dir,
    const Node& startNode,
    const Node& endNode,
    int traversal_depth,
    bool &success
)
{
    set<Node> usedNode;
    usedNode.insert(startNode);
    PATH current_consensus;
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
bool MonumentTraversal::validate_consensuses (set<PATH>& consensuses, PATH& result)
{
    bool debug = false;
    // compute mean and stdev of consensuses
    int mean = 0;
    int path_number = 0;
    for(set<PATH>::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        //if (debug)  printf("bubble path %d: %s (len=%d)\n",path_number,(*it).c_str(),(*it).length());
        mean+=(*it).size();
        path_number++;
    }
    mean/=consensuses.size();
    double stdev = 0;
    for(set<PATH>::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
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

    if (debug) printf("%d-bubble mean %d, stdev %.1f\n",consensuses.size(),mean,stdev);

    // traverse bubbles if paths have roughly the same length
    if (stdev>mean/5)
        return false;

    // check that all consensuses are similar
    if (! all_consensuses_almost_identical(consensuses))
        return false;

    // if all good, an arbitrary consensus is chosen
    PATH chosen_consensus = *consensuses.begin();
    int result_length = chosen_consensus.size();
    if  (result_length> max_depth) // it can happen that consensus is longer than max_depth, despite that we didn't explore that far (in a messy bubble with branchings inside)
        return false;

    result.resize(result_length);

    for (size_t i=0; i<chosen_consensus.size(); i++)
    {
        result[i] = chosen_consensus[i];
    }

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
bool MonumentTraversal::all_consensuses_almost_identical (set<PATH>& consensuses)
{
    for (set<PATH>::iterator it_a = consensuses.begin(); it_a != consensuses.end(); it_a++)
    {
        set<PATH>::iterator it_b = it_a;
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
