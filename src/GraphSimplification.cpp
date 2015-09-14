/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014-2015  INRIA
 *   Authors: R.Chikhi, G.Rizk, D.Lavenier, E.Drezen
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

#define DEBUG(a)   //a
#define TIME(a)   a

#include <vector>
#include <set>
#include <stack>
#include <chrono>
#include <GraphSimplification.hpp>

#define get_wtime() chrono::system_clock::now()
#define diff_wtime(x,y) (unsigned long)chrono::duration_cast<chrono::nanoseconds>(y - x).count()
#define DIR2STR(dir) ((dir==DIR_OUTCOMING) ? "outcoming" : "incoming")

using namespace std;

static const char* progressFormat0 = "removing tips,    pass %2d ";
static const char* progressFormat1 = "removing bubbles, pass %2d ";
static const char* progressFormat2 = "removing bulges,  pass %2d ";
static const char* progressFormat3 = "removing ec,      pass %2d ";

GraphSimplification::GraphSimplification(const Graph & graph, int nbCores)
        : _nbTipRemovalPasses(0), _nbBubbleRemovalPasses(0), _nbBulgeRemovalPasses(0), _nbECRemovalPasses(0), _graph(graph), 
        _nbCores(nbCores), _firstNodeIteration(true)
{
    // just a way to get number of nodes
    char buffer[128];
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);
    unsigned long nbNodes = itNode.size();

    interestingNodes.resize(nbNodes); // number of graph nodes // (!) this will alloc 1 bit per kmer.
    for (unsigned long i = 0; i < nbNodes; i++)
        interestingNodes[i] = false;
}


double GraphSimplification::getSimplePathCoverage(Node node, Direction dir, unsigned int *pathLenOut, unsigned int maxLength)
{
    Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (node, dir);
    unsigned long total_abundance = _graph.queryAbundance(node.kmer);
    unsigned int pathLen = 1;
    for (itNodes.first(); !itNodes.isDone(); itNodes.next())
    {
        total_abundance += _graph.queryAbundance((*itNodes).kmer);
        pathLen++;
        if (maxLength > 0 && pathLen >= maxLength)
            break;
    }
    *pathLenOut = pathLen;
    return ((double)total_abundance) / ((double)pathLen);
}

// gets the mean abundance of neighboring paths around a branching node (excluding the path that starts with nodeToExclude, e.g. the tip itself)
double GraphSimplification::getMeanAbundanceOfNeighbors(Node branchingNode, Node nodeToExclude)
{
    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(branchingNode);
    unsigned int nbNeighbors = 0;
    double meanNeighborsCoverage = 0;
    //DEBUG(cout << endl << "called getMeanAbudanceOfNeighbors for node " << _graph.toString(branchingNode) << " of degrees " << _graph.indegree(branchingNode) <<"/"<< _graph.outdegree(branchingNode)<< " excluding node  " <<  _graph.toString (nodeToExclude) << endl);
    for (size_t i = 0; i < neighbors.size(); i++)
    {
        Node neighbor = neighbors[i].to;
        if (neighbor == nodeToExclude) // (in gatb-core, Node == Node means Node.kmer == Node.kmer)
        {
            //DEBUG(cout << endl << "good, seen the node to exclude" << endl);
            continue; 
        }

        unsigned int pathLen;
        double simplePathCoverage = getSimplePathCoverage(neighbor, neighbors[i].direction, &pathLen, 100);
        meanNeighborsCoverage += simplePathCoverage;
        nbNeighbors++;

        DEBUG(cout << endl << "got simple path coverage for neighbor " << nbNeighbors  << " : " << " meancoverage: " <<simplePathCoverage << " over " << pathLen << " kmers" << endl);
    }
    meanNeighborsCoverage /= nbNeighbors;
    return meanNeighborsCoverage;
}

// this needs to be in Graph.cpp of gatb-core
string GraphSimplification::path2string(Direction dir, Path p, Node endNode)
{
    // naive conversion from path to string

    string p_str;
    if (dir == DIR_INCOMING)
    {
        // TODO: remove this code once gatb-core is fixed w.r.t DIR_INCOMING bug
        Node revstart = endNode;
        p_str = _graph.toString(revstart);
        for (size_t i = 0; i < p.size(); i++)
            p_str.push_back(p.ascii(p.size()-1-i));
    }
    else
    {
        p_str = _graph.toString(p.start);
        for (size_t i = 0; i < p.size(); i++)
            p_str.push_back(p.ascii(i));
    }
    return p_str;
}

string maybe_print(long value, string str)
{
    if (value == 0)
        return "";
    return std::to_string(value) + " " + str;
}


/* coverage of the simple path (stored in "nodes" vector)
      then compares it to coverage of other paths connected to the last node of it. */
bool GraphSimplification::satisfyRCTC(vector<Node> nodes, double RCTCcutoff)
{

    unsigned long mean_abundance = 0;
    for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
    {
        unsigned int abundance = _graph.queryAbundance((*itVecNodes).kmer);
        mean_abundance += abundance;
    }
    double meanTipAbundance = (double)mean_abundance / (double)(nodes.size());
    double stdevTipAbundance = 0;

    // get std dev, for debug only
    bool debugstdev = false;
    if (debugstdev)
    {
        for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
        {
            unsigned int abundance = _graph.queryAbundance((*itVecNodes).kmer);
            stdevTipAbundance += pow(fabs(abundance-meanTipAbundance),2);
        }
        stdevTipAbundance = sqrt(stdevTipAbundance/nodes.size());
    }

    // explore the other two or more simple paths connected to that path, to get an abundance estimate
    // but first, get the branching node(s) the tip is connected to 
    /* (it's weird when it's more than one branching node though, it's a situation like:
     * 
     *                ..o--o--o--o..
     *                    /
     *  nodes-> o--o--o--o
     *                    \
     *                ..o--o--o--o.. 
     *
     *  instead of the classical:
     *
     *  nodes-> o--o--o--o
     *                    \
     *                ..o--o--o--o..)*/

    Graph::Vector<Edge> connectedBranchingNodes = _graph.neighbors<Edge>(nodes.back());
    unsigned int nbBranchingNodes = 0;
    double meanNeighborsCoverage = 0;
    bool foundIt = false; // just a safety. can be removed later
    for (size_t j = 0; j < connectedBranchingNodes.size(); j++)
    {
        /* we should the second-to-last node from "nodes" as a neighbor, and skip it */
        if ((nodes.size() >= 2) && (connectedBranchingNodes[j].to == nodes[nodes.size() - 2]))
        {
            foundIt = true;
            continue;
        }
        meanNeighborsCoverage += getMeanAbundanceOfNeighbors(connectedBranchingNodes[j].to, nodes.back());
        nbBranchingNodes++;
    }
    if (foundIt == false && nodes.size() >= 2)
        cout << "WTF!!..!!!" << endl;

    if (nbBranchingNodes > 0)
        meanNeighborsCoverage /= nbBranchingNodes;

    bool isRCTC = (meanNeighborsCoverage > RCTCcutoff * meanTipAbundance);

    DEBUG(cout << endl << "RCTC test, over " << nbBranchingNodes << " connected nodes. Global mean neighbors coverage: " << meanNeighborsCoverage <<  " compared to mean tip abundance over "<< nodes.size() << " values : " << meanTipAbundance << (debugstdev? " stddev: " : "") << (debugstdev? to_string(stdevTipAbundance): "") << ", is RCTC satisfied? " << isRCTC << endl);

    return isRCTC;
}

/* okay let's analyze SPAdes 3.5 tip clipping conditions, just for fun: (following graph_simplifications.hpp and simplifications.info)
 *
 * * tc_lb is a coefficient, setting tip length to be max(read length, g*tc_lb) (see simplification_settings.hpp);
 *
 * * cb is a plain coverage upper bound (see basic_edge_conditions.hpp)
 * when it's auto, it's equal to ec_bound. (source: set_detected_coverage_bound(gp.ginfo.ec_bound() in graph_simplification.hpp)
 *
 * * * ec_bound is given by the error threshold found by the coverage model in kmer_coverage_model.cpp (see: genomic_info_filler.cpp)
 * it's the largest abundance value such that the probability of an erroneous kmer is less than 0.05. And, if it's smaller than the valley, it's the valley.
 * EC bound looks off in metagenomic assemblies (a value over 50.. seems high).
 *
 * * * in single cell mode, it's max(average edge coverage, graph theshold), and graph threshold is quite advanced! (omni_tools.hpp)
 * but already, in multicell mode, spades does a good job at tip clipping, so let's look at that first.
 *
 * * rctc is a "relative coverage tip condition" and the value given next to it is: max_relative_coverage
 * rctc examines whether the coverage of the tip is <= max_relative_coverage * (max_{over all neighbors}(coverage of neighbor) + 1)
 
 * actual tip clipping code for "tc", i.e. implementations of the conditions, is in tip_clipper
 * (omni/graph_processing_algorithm.hpp:EdgeRemovingAlgorithm() is just generically parsing conditions)
 *
 * in one condition, tc_lb = 3.5, cb is 1000000 (minia's coverages values don't go that far.. yet..), rctc 2
 * and in another, tc_lb is 10 and cb is auto
 *
 * some historical facts:
 * for CAMI we were more loose than SPAdes: topological tip had no cov criterion, wherehas it should have had rctc (like spades)
 * and long tc_lb10 tips should have the auto coverage bound, but instead they had rctc 2. 
 *
 * so TODO: make it more strict. but for now I'm focusing on EC.
 */
unsigned long GraphSimplification::removeTips()
{
    unsigned int k = _graph.getKmerSize();
    
    unsigned int maxTipLengthTopological = (unsigned int)((float)k * (3.5 - 1.0)); // aggressive with SPAdes length threshold, but no coverage criterion
    unsigned int maxTipLengthRCTC = (unsigned int)(k * 10); // experimental, SPAdes-like
    double RCTCcutoff = 2; // SPAdes-like

    unsigned long nbTipsRemoved = 0;

    /** We get an iterator over all nodes */
    char buffer[128];
    sprintf(buffer, progressFormat0, ++_nbTipRemovalPasses);
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // parallel stuff
    vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
    unsigned long nbNodes = itNode.size();
    nodesToDelete.resize(nbNodes); // number of graph nodes // (!) this will alloc 1 bit per kmer.
    for (unsigned long i = 0; i < nbNodes; i++)
        nodesToDelete[i] = false;

    bool haveInterestingNodesInfo = !_firstNodeIteration;
    _firstNodeIteration = false;

    dispatcher.iterate (itNode, [&] (Node& node)
    {
        unsigned long index = _graph.nodeMPHFIndex(node);

        if (haveInterestingNodesInfo)
            if (interestingNodes[index] == false)
                return; // no point in examining non-branching nodes, saves calls to in/out-degree, i.e. accesses to the minia datastructure

        if (_graph.isNodeDeleted(index)) { return; } // {continue;} // sequential and also parallel
        if (nodesToDelete[index]) { return; }  // parallel // actually not sure if really useful

        unsigned inDegree = _graph.indegree(node), outDegree = _graph.outdegree(node);

        if (!haveInterestingNodesInfo)
            interestingNodes[index] = interestingNodes[index] || (!(inDegree == 1 && outDegree == 1));

        /* tips have out/in degree of 0 on one side, and any non-zero degree on the other */
        if ((inDegree == 0 || outDegree == 0) && (inDegree != 0 || outDegree != 0))
        {
            bool isShortTopological = true;
            bool isShortRCTC = true;

            //DEBUG(cout << endl << "deadend node: " << _graph.toString (node) << endl);

            /** We follow the simple path to get its length */
            Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node.kmer); // so, it has one or more neighbors in a single direction
            
            /* it may appear that we're only going to follow its first neighbor, but in fact, neighbors[0].from is node.kmer */
            /* so, follow the simple path from this start tip node to the further node that has out-branching (out is w.r.t to the direction) */
            Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[0].from, neighbors[0].direction); //
            //DEBUG(cout << endl << "neighbors from: " << _graph.toString (neighbors[0].from) << " direction: " << neighbors[0].direction << endl);
            unsigned int pathLen = 1;
            vector<Node> nodes;
            nodes.push_back(node);
            /* get that putative tip length (stop at a max) */
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                nodes.push_back(*itNodes);
                if (k + pathLen >= maxTipLengthTopological) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                    isShortTopological = false;
                /* don't break here, tip might still be long enough for RCTC length*/

                if (k + pathLen >= maxTipLengthRCTC) 
                {
                    isShortRCTC= false;
                    break;
                }

                pathLen++;
            }

            if (node != neighbors[0].from)
                cout << "WTF!.!" << endl; // TODO remove that after a while

            // if it's not short, then no point in computing whether it's connected (this is a little optimization to a save the next neighbors call)
            if ( ! (isShortTopological || isShortRCTC) )
                return;
            
            // at this point, the last node in "nodes" is the last node of the tip.
            // check if it's connected to something. 
            // condition: degree > 1, because connected to the tip and to that "something"
            bool isConnected = (_graph.neighbors<Edge>(nodes.back()).size() > 1);
            if (pathLen == 1)
            {
                // special case: only a single tip node, check if it's not isolated
                isConnected |=  (_graph.indegree(node) != 0 || _graph.outdegree(node) != 0); 
            }

            bool isTopologicalShortTip = isShortTopological && isConnected; 
            bool isMaybeRCTCTip = isShortRCTC && isConnected;

            //DEBUG(cout << endl << "pathlen: " << pathLen << " last node " << _graph.toString(nodes.back()) << " neighbors in/out: " <<_graph.indegree(nodes.back()) << " " << _graph.outdegree(nodes.back()) << " istoposhorttip: " << isTopologicalShortTip << endl);

            bool isRCTCTip = false;
            if (!isTopologicalShortTip && isMaybeRCTCTip)
            {
                isRCTCTip = satisfyRCTC(nodes, RCTCcutoff);
            }

            bool isTip = isTopologicalShortTip || isRCTCTip; 


            if (isTip)
            {
                // delete it
                //

                //DEBUG(cout << endl << "TIP of length " << pathLen << " FOUND: " <<  _graph.toString (node) << endl);
                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                {
                    //DEBUG(cout << endl << "deleting tip node: " <<  _graph.toString (*itVecNodes) << endl);
                    //_graph.deleteNode(*itVecNodes); // sequential version

                    unsigned long index = _graph.nodeMPHFIndex(*itVecNodes); // parallel version
                    nodesToDelete[index] = true; // parallel version

                }

                // update interesting status of connected nodes
                Graph::Vector<Node> connectedBranchingNodes = _graph.neighbors<Node>(nodes.back());
                for (size_t j = 0; j < connectedBranchingNodes.size(); j++)
                {
                    unsigned long index = _graph.nodeMPHFIndex(connectedBranchingNodes[j]);
                    interestingNodes[index] = true;
                }


                __sync_fetch_and_add(&nbTipsRemoved, 1);


            }
        }
        // } // sequential
    }); // parallel

    // now delete all nodes, in sequential (shouldn't take long)
    for (unsigned long i = 0; i < nbNodes; i++)
    {
        if (nodesToDelete[i])
        {
           _graph.deleteNode(i);
        }
    }


    return nbTipsRemoved;
}

enum HMCP_Success { HMCP_DEADEND = 0, HMCP_FOUND_END = 1 , HMCP_MAX_DEPTH = -1, HMCP_LOOP = - 2};

/* note: the returned mean abundance does not include start and end nodes */
Path GraphSimplification::heuristic_most_covered_path(
        Direction dir, const Node startNode, const Node endNode, 
        int traversal_depth, int& success, double &abundance, bool most_covered, 
        unsigned int backtrackingLimit, Node *avoidFirstNode)
{
    set<Node::Value> usedNode;
    usedNode.insert(startNode.kmer);
    Path current_path;
    current_path.start = startNode;
    success = HMCP_DEADEND;
    vector<int> abundances; 
    unsigned long nbCalls = 0;

    Path res = heuristic_most_covered_path(dir, startNode, endNode, traversal_depth, current_path, usedNode, success, abundances, most_covered,
            backtrackingLimit, avoidFirstNode, nbCalls);

    abundance = 0;
    for (unsigned int i = 0; i < abundances.size(); i++){
        abundance += abundances[i];}
    abundance /= abundances.size();

    bool debug_abundances = false;
    if (debug_abundances)
    {
        cout << "abundance for path (is most covered path: " << most_covered << "): ";
        for (unsigned int i = 0; i < abundances.size(); i++)
            cout << abundances[i]<< " ";
        cout << ";"<<endl;
    }

    bool debug_nbcalls = false;
    if (debug_nbcalls)
        cout << "number of path-finding calls: " << nbCalls << endl;

    return res;
}
        
Path GraphSimplification::heuristic_most_covered_path(
        Direction dir, const Node startNode, const Node endNode, 
        int traversal_depth, Path current_path, set<Node::Value> usedNode, int& success, vector<int>& abundances, bool most_covered,
        unsigned int backtrackingLimit, Node *avoidFirstNode, unsigned long &nbCalls)
{
    // inspired by all_consensuses_between
    nbCalls++;
    
    if (traversal_depth < -1)
    {
        success = HMCP_MAX_DEPTH;
        return current_path;
    }

    if (startNode.kmer == endNode.kmer)
    {
        success = 1;
        return current_path;
    }

    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge> (startNode, dir);

    /** We loop these neighbors. */
    vector<std::pair<int, Edge> > abundance_node;
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& edge = neighbors[i];

        if (avoidFirstNode != NULL && edge.to.kmer == avoidFirstNode->kmer)
            continue;

        // don't resolve bubbles containing loops
        // (tandem repeats make things more complicated)
        // that's a job for a gapfiller
        if (usedNode.find(edge.to.kmer) != usedNode.end())
        {
            success = -2;
            return current_path;
        }

        unsigned int abundance = _graph.queryAbundance(neighbors[i].to.kmer);
        abundance_node.push_back(std::make_pair(abundance, edge));
    }

    std::sort(abundance_node.begin(), abundance_node.end()); // sort nodes by abundance
    if (most_covered) 
        std::reverse(abundance_node.begin(), abundance_node.end()); // decreasing abundances

    // traverse graph in abundance order, return most abundant path
    for (unsigned int i = 0; i < abundance_node.size(); i++)
    {
        Edge edge = abundance_node[i].second;

        // generate extended consensus sequence
        Path extended_path(current_path);
        extended_path.push_back (edge.nt);

        // generate list of used kmers (to prevent loops)
        set<Node::Value> extended_kmers (usedNode);
        extended_kmers.insert (edge.to.kmer);

        // extend abundances
        vector<int> extended_abundances (abundances);
        if (edge.to.kmer != endNode.kmer) // skip abundance of last node
            extended_abundances.push_back(abundance_node[i].first);

        // recursive call to all_consensuses_between
        Path new_path = heuristic_most_covered_path (
            dir,
            edge.to,
            endNode,
            traversal_depth - 1,
            extended_path,
            extended_kmers,
            success,
            extended_abundances,
            most_covered,
            backtrackingLimit,
            NULL, // no longer avoid nodes
            nbCalls
        );

        if ((success == 1)|| (backtrackingLimit > 0 && nbCalls >= backtrackingLimit)) // on success or if no more backtracking, return immediately
        {
            abundances = extended_abundances;
            return new_path; 
        }
    }

    return current_path;


}


/* bulge removal algorithm. mimics spades, which doesnt remove bubbles, but only bulges. looks as effective.
 * it's slow to do heuristic_find_most_covered path so i'm testing it with no backtracking
 *
 * see a-b-c here for an explanation of bulge removal: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791033/figure/f4/
 *
 * spades pops bulges based on something like the ratio between most examined simple path and a more covered path is (whether it is above 1.1).
 * so i'm actually doing just that. I recall checking spades source code to implement this. this was during CAMI.
 */ 
unsigned long GraphSimplification::removeBulges()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int coeff = 3;
    unsigned int additive_coeff = 100;
    unsigned int maxBulgeLength = std::max((unsigned int)((double)k * coeff), (unsigned int)(k + additive_coeff)); // SPAdes, exactly

    unsigned int backtrackingLimit = maxBulgeLength; // arbitrary

    // stats
    //
    unsigned long nbBulgesRemoved = 0;
    unsigned long nbSimplePaths = 0;
    unsigned long nbLongSimplePaths = 0;
    unsigned long nbShortSimplePaths = 0;
    unsigned long nbTopologicalBulges = 0;
    unsigned long nbFirstNodeDeleted = 0;
    unsigned long nbFirstNodeGraphDeleted = 0;
    unsigned long nbNoAltPathBulgesLoop = 0, nbNoAltPathBulgesDepth = 0,nbNoAltPathBulgesDeadend = 0;
    unsigned long nbBadCovBulges = 0;

    unsigned long timeAll = 0, timePathFinding = 0, timeFailedPathFinding = 0, timeLongestFailure = 0,
                  timeSimplePath = 0, timeDelete = 0, timePost = 0, timeVarious = 0;

    unsigned long longestFailureDepth = 0;

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, progressFormat2, ++_nbBulgeRemovalPasses);
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // parallel stuff
    vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
    unsigned long nbNodes = itNode.size();
    nodesToDelete.resize(nbNodes); // number of graph nodes // (!) this will alloc 1 bit per kmer.
    for (unsigned long i = 0; i < nbNodes; i++)
        nodesToDelete[i] = false;

    bool haveInterestingNodesInfo = !_firstNodeIteration;

    dispatcher.iterate (itNode, [&] (Node& node)
    {

      TIME(auto start_thread_t=get_wtime());

      unsigned long index = _graph.nodeMPHFIndex(node);

      // TODO think about cases where bulge suppression could make a node intresting
    /*  if (haveInterestingNodesInfo)
          if (interestingNodes[index] == false)
            return; // no pont in examining non-branching nodes, saves calls to in/out-degree, i.e. accesses to the minia datastructure
*/

      // need to search in both directions
      for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
      {
         if ((_graph.outdegree(node) >= 2 && dir == DIR_OUTCOMING) || (_graph.indegree(node) >= 2 && dir == DIR_INCOMING))
         {
            TIME(auto start_various_overhead_t=get_wtime());
            if (_graph.isNodeDeleted(index)) { return; } // {continue;} // sequential and also parallel
            if (nodesToDelete[index]) { return; }  // parallel // actually not sure if really useful

            DEBUG(cout << endl << "putative bulge node: " << _graph.toString (node) << endl);

            /** We follow the outgoing simple paths to get their length and last neighbor */
            Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node, dir);
            TIME(auto end_various_overhead_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

            // do everying for each possible short simple path that is neighbor of that node
            for (unsigned int i = 0; i < neighbors.size(); i++)
            {
                vector<Node> nodes;
                bool foundShortPath = false;
                unsigned int pathLen = 0;
            
                TIME(auto start_various_overhead_t=get_wtime());
                unsigned long index =_graph.nodeMPHFIndex(neighbors[i].to);
                if (_graph.isNodeDeleted(index)) { 
                     __sync_fetch_and_add(&nbFirstNodeGraphDeleted, 1);
                    continue;}
                if (nodesToDelete[index]) { 
                     __sync_fetch_and_add(&nbFirstNodeDeleted, 1);
                    continue;}

                TIME(auto end_various_overhead_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

                /* explore the simple path from that node */
                TIME(auto start_simplepath_t=get_wtime());
                Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[i].to, dir);
                DEBUG(cout << endl << "neighbors " << i+1 << "/" << neighbors.size() << " from: " << _graph.toString (neighbors[i].to) << " dir: " << DIR2STR(dir) << endl);
                bool isShort = true;
                pathLen = 0;
                nodes.push_back(neighbors[i].to);
                for (itNodes.first(); !itNodes.isDone(); itNodes.next())
                {
                    nodes.push_back(*itNodes);
                    if (k + pathLen++ >= maxBulgeLength) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                    {
                        __sync_fetch_and_add(&nbLongSimplePaths, 1);
                        isShort = false;
                        break;       
                    }
                }
                TIME(auto end_simplepath_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeSimplePath, diff_wtime(start_simplepath_t,end_simplepath_t)));

                __sync_fetch_and_add(&nbSimplePaths, 1);

                if (!isShort || pathLen == 0) // can't do much if it's pathLen=0, we don't support edge removal, only node removal
                    continue;
                
                __sync_fetch_and_add(&nbShortSimplePaths, 1);

                TIME(start_various_overhead_t=get_wtime());
                Graph::Vector<Edge> outneighbors = _graph.neighbors<Edge>(nodes.back(), dir);
                DEBUG(cout << "last node of simple path: "<< _graph.toString(nodes.back()) << " has indegree/outdegree: " <<_graph.indegree(nodes.back()) << "/" << _graph.outdegree(nodes.back()) << endl);

                if (outneighbors.size() == 0) // might still be a tip, unremoved for some reason
                    continue;

                Node endNode = outneighbors[0].to;
                DEBUG(cout << "endNode: " << _graph.toString(endNode) << endl);

                // at this point, the last node in "nodes" is the last node of a potential Bulge path, and endNode is hopefully a branching node right after.
                // check if it's connected to something that has in-branching. 
                bool isDoublyConnected = (_graph.indegree(endNode) > 1 && dir==DIR_OUTCOMING) || (dir==DIR_INCOMING && _graph.outdegree(endNode) > 1);

                bool isTopologicalBulge = isDoublyConnected;

                DEBUG(cout << "pathlen: " << pathLen << " istopobulge: " << isTopologicalBulge << endl);

                TIME(end_various_overhead_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

                if (!isTopologicalBulge)
                    continue;

                __sync_fetch_and_add(&nbTopologicalBulges, 1);

                unsigned int depth = std::max((unsigned int)(pathLen * 1.1),(unsigned int) 3); // following SPAdes
                double mean_abundance_most_covered;
                int success;
                Node startNode = node;

                TIME(auto start_pathfinding_t=get_wtime());

                Path heuristic_p_most = heuristic_most_covered_path(dir, startNode, endNode, depth+2, success, mean_abundance_most_covered,
                        true, // most covered
                        backtrackingLimit, // avoid backtracking
                        &(neighbors[i].to) // avoid that node
                        );

                TIME(auto end_pathfinding_t=get_wtime());
                TIME(__sync_fetch_and_add(&timePathFinding, diff_wtime(start_pathfinding_t,end_pathfinding_t)));
                TIME(auto start_post_t=get_wtime());

                if (success != 1)
                {
                    TIME(__sync_fetch_and_add(&timeFailedPathFinding, diff_wtime(start_pathfinding_t,end_pathfinding_t)));
                    TIME(if (diff_wtime(start_pathfinding_t,end_pathfinding_t) > timeLongestFailure) { timeLongestFailure = diff_wtime(start_pathfinding_t,end_pathfinding_t); longestFailureDepth = depth;});

                    if (success == -2)
                        __sync_fetch_and_add(&nbNoAltPathBulgesLoop, 1);
                    if (success == -1)
                         __sync_fetch_and_add(&nbNoAltPathBulgesDepth, 1);
                    if (success == 0)
                         __sync_fetch_and_add(&nbNoAltPathBulgesDeadend, 1);
                    continue;
                }

                DEBUG(cout << "alternative path is:  "<< path2string(dir, heuristic_p_most, endNode)<< " abundance: "<< mean_abundance_most_covered <<endl);

                bool debug = false;
                if (debug)
                {
                    double mean_abundance_least_covered;
                    Path heuristic_p_least = heuristic_most_covered_path(dir, startNode, endNode, depth+2, success, mean_abundance_least_covered,false);
                    DEBUG(cout << endl << "alternative least is: "<< path2string(dir, heuristic_p_least, endNode)<< " abundance: "<< mean_abundance_least_covered <<endl);
                }

                unsigned int dummyLen;
                double simplePathCoverage = getSimplePathCoverage(nodes[1], dir, &dummyLen);

                DEBUG(cout << "retraced bulge path over length: " << dummyLen << endl);

                bool isBulge =  simplePathCoverage * 1.1  <=  mean_abundance_most_covered;

                DEBUG(cout << "bulge coverages: " << simplePathCoverage<< "/" <<  mean_abundance_most_covered  << endl);

                if (!isBulge)
                {
                    __sync_fetch_and_add(&nbBadCovBulges, 1);
                    DEBUG(cout << "not a bulge due to coverage criterion" << endl);
                    continue;
                }

                // delete it
                //

                DEBUG(cout << endl << "BULGE of length " << pathLen << " FOUND: " <<  _graph.toString (node) << endl);
                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                {
                    //DEBUG(cout << endl << "deleting node " << _graph.toString(*itVecNodes) << endl);
                    unsigned long index = _graph.nodeMPHFIndex(*itVecNodes); // parallel version
                    nodesToDelete[index] = true; // parallel version
                }

                __sync_fetch_and_add(&nbBulgesRemoved, 1);
                TIME(auto end_post_t=get_wtime());
                TIME(__sync_fetch_and_add(&timePost, diff_wtime(start_post_t,end_post_t)));

            } // for neighbors
        } // if outdegree
      } // for direction
        TIME(auto end_thread_t=get_wtime());
        TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t)));
    }); // parallel
    
    // now delete all nodes, in sequential (shouldn't take long)
    TIME(auto start_nodedelete_t=get_wtime());
    for (unsigned long i = 0; i < nbNodes; i++)
    {
        if (nodesToDelete[i])
        {
           _graph.deleteNode(i);
        }
    }
    TIME(auto end_nodedelete_t=get_wtime());
    TIME(__sync_fetch_and_add(&timeDelete, diff_wtime(start_nodedelete_t,end_nodedelete_t)));

    cout << nbBulgesRemoved << " bulges removed. " << endl <<
        nbSimplePaths << "/" << nbLongSimplePaths << "+" <<nbShortSimplePaths << " any=long+short simple path examined across all threads, among them " <<
        nbTopologicalBulges << " topological bulges, " << nbFirstNodeDeleted << "+" << nbFirstNodeGraphDeleted << " were first-node duplicates." << endl;
    cout << nbNoAltPathBulgesDepth << "+" << nbNoAltPathBulgesLoop << "+" << nbNoAltPathBulgesDeadend << " without alt. path (complex+loop+deadend), " 
    << nbBadCovBulges << " didn't satisfy cov. criterion." << endl;

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);
    TIME(cout << "Timings: " << timeAll / unit << " CPUsecs total."<< endl);
    TIME(cout << "         " << timeSimplePath / unit << " CPUsecs simple path traversal." << endl);
    TIME(cout << "         " << timePathFinding / unit << "(/" << timePathFinding / unit << ") CPUsecs path-finding(/failed). Longest: " << timeLongestFailure / (unit/1000) << " CPUmillisecs (depth " << longestFailureDepth << ")." << endl);
    TIME(cout << "         " << timePost / unit << " CPUsecs topological bulge processing, " << timeDelete / unit << " CPUsecs nodes deletion." << endl);
    TIME(cout << "         " << timeVarious / unit << " CPUsecs various overhead." << endl);

    return nbBulgesRemoved;
}



/* Again, let's see what spades 3.5 does.
 *erroneous connection remover is:

 RemoveLowCoverageEdges
 calls config then calls:
  omnigraph::EdgeRemovingAlgorithm which is same as tc

  to_ec_lb in ../src/debruijn/simplification/simplification_settings.hpp
  is exactly like tip clipping but with a different length (2 * (max tip length with coeff 5) - 1, so that's exactly 10*min(k,readlen/2) - 1)

  icb is something that removes edges with coverage = (cov_bound * iter / nb_iters )
  where cov_bound is same as in cb

  and it's asking for AlternativesPresenceCondition:

    o   o                           o-->-O
   /     \                             /
   O-->---O    drawn differently:     /
                                    O-->-o

  so anyway, since we're not computing the coverage model like SPAdes does, I'm going to use RCTC 2

*/
unsigned long GraphSimplification::removeErroneousConnections()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int maxECLength = (unsigned int)((float)k * (10 - 1.0)) ;  // SPAdes mode 
    double RCTCcutoff = 2.0;

    unsigned long nbSimplePaths = 0;
    unsigned long nbLongSimplePaths = 0;
    unsigned long nbShortSimplePaths = 0;
    unsigned long nbTopologicalEC = 0;
    unsigned long nbECRemoved = 0;

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, progressFormat3, ++_nbECRemovalPasses);
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // parallel stuff
    vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
    unsigned long nbNodes = itNode.size();
    nodesToDelete.resize(nbNodes); // number of graph nodes // (!) this will alloc 1 bit per kmer.
    for (unsigned long i = 0; i < nbNodes; i++)
        nodesToDelete[i] = false;

    unsigned long timeAll = 0, timeSimplePath = 0;

    dispatcher.iterate (itNode, [&] (Node& node)
            {
            TIME(auto start_thread_t=get_wtime());

            unsigned long index = _graph.nodeMPHFIndex(node);


            /* TODO think about interestingnodes info, at the same time as we think for it for bulge removal. right now it's only implemented for tips. might speed bulges/EC removal up too */
            //if (interestingNodes[index] == false)
            //return; // no point in examining non-branching nodes, saves calls to in/out-degree, i.e. accesses to the minia datastructure

            if (_graph.isNodeDeleted(index)) { return; } // {continue;} // sequential and also parallel
            if (nodesToDelete[index]) { return; }  // parallel // actually not sure if really useful

            unsigned inDegree = _graph.indegree(node), outDegree = _graph.outdegree(node);

            // if (!haveInterestingNodesInfo)
            // interestingNodes[index] = interestingNodes[index] || (!(inDegree == 1 && outDegree == 1));

            /* ec nodes have out/in degree of 1 or more on one side, and of 2 or more on the other */
            if (!((inDegree >= 1 && outDegree > 1 ) || (inDegree > 1 && outDegree >=1 )))
                return;

            // need to search in both directions
            for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
            {

                if ((_graph.outdegree(node) >= 2 && dir == DIR_OUTCOMING) || (_graph.indegree(node) >= 2 && dir == DIR_INCOMING))
                {  
                    if (_graph.isNodeDeleted(node)) { return; } // {continue;} // sequential and also parallel
                    if (nodesToDelete[_graph.nodeMPHFIndex(node)]) { return; }  // parallel // actually not sure if really useful

                    DEBUG(cout << endl << "putative EC node: " << _graph.toString (node) << endl);

                    /** We follow the outcoming simple paths (so, if it's outdegree 2, we follow the outcoming simple paths.
                     * to get their length and last neighbor */
                    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node, dir);

                    // do everying for each possible short simple path that is neighbor of that node
                    for (unsigned int i = 0; i < neighbors.size(); i++)
                    {

                        unsigned long index =_graph.nodeMPHFIndex(neighbors[i].to);
                        if (_graph.isNodeDeleted(index)) { 
                            continue;}
                        if (nodesToDelete[index]) { 
                            continue;}

                        /* explore the simple path from that node */
                        vector<Node> nodes;
                        bool foundShortPath = false;
                        unsigned int pathLen = 0;
                        TIME(auto start_simplepath_t=get_wtime());
                        Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[i].to, dir);
                        DEBUG(cout << endl << "neighbors " << i+1 << "/" << neighbors.size() << " from: " << _graph.toString (neighbors[i].to) << " dir: " << DIR2STR(dir) << endl);
                        bool isShort = true;
                        pathLen = 0;
                        nodes.push_back(neighbors[i].to);
                        for (itNodes.first(); !itNodes.isDone(); itNodes.next())
                        {
                            nodes.push_back(*itNodes);
                            if (k + pathLen++ >= maxECLength) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                            {
                                __sync_fetch_and_add(&nbLongSimplePaths, 1);
                                isShort = false;
                                break;       
                            }
                        }
                        TIME(auto end_simplepath_t=get_wtime());
                        TIME(__sync_fetch_and_add(&timeSimplePath, diff_wtime(start_simplepath_t,end_simplepath_t)));

                        __sync_fetch_and_add(&nbSimplePaths, 1);

                        if (!isShort || pathLen == 0) // can't do much if it's pathLen=0, we don't support edge removal, only node removal
                        {
                            DEBUG(cout << "direction: " << DIR2STR(dir) << ", not an EC: foundShortPath: " << foundShortPath << " pathLen: " << pathLen << endl);
                            continue;
                        }

                        __sync_fetch_and_add(&nbShortSimplePaths, 1);

                        Graph::Vector<Edge> outneighbors = _graph.neighbors<Edge>(nodes.back(), dir);
                        DEBUG(cout << "last simple path node: "<< _graph.toString(nodes.back()) << " has " << outneighbors.size() << " outneighbors" << endl);

                        if (outneighbors.size() == 0) // might still be a tip, unremoved for some reason
                            continue;

                        Node endNode = outneighbors[0].to;
                        DEBUG(cout << "endNode: " << _graph.toString(endNode) << endl);

                        // at this point, the last node in "nodes" is the last node of a potential EC, and endNode is hopefully a branching node right after.
                        // check if it's connected to something that has in-branching and also an out neighbor. 
                        bool isDoublyConnected = (_graph.indegree(endNode) > 1 && _graph.outdegree(endNode) >= 1 && dir==DIR_OUTCOMING) || \
                                                 (dir==DIR_INCOMING && _graph.outdegree(endNode) > 1 && _graph.indegree(endNode) >= 1 );


                        bool isTopologicalEC = isDoublyConnected;

                        DEBUG(cout << "direction: " << DIR2STR(dir) << ", pathlen: " << pathLen << " last node neighbors size: " << _graph.neighbors<Edge>(nodes.back()).size() << " indegree outdegree: " <<_graph.indegree(node) << " " << _graph.outdegree(node) << " isDoublyConnected: " << isDoublyConnected << " isTopoEC: " << isTopologicalEC << endl);

                        if (isTopologicalEC)
                        {

                            bool isRCTC = satisfyRCTC(nodes, RCTCcutoff);

                            bool isEC = isRCTC;


                            if (isEC)
                            {
                                // delete it
                                //

                                DEBUG(cout << endl << "EC of length " << pathLen << " FOUND: " <<  _graph.toString (node) << endl);
                                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                                {
                                    //DEBUG(cout << endl << "deleting EC node: " <<  _graph.toString (*itVecNodes) << endl);
                                    //_graph.deleteNode(*itVecNodes); // sequential version

                                    unsigned long index = _graph.nodeMPHFIndex(*itVecNodes); // parallel version
                                    nodesToDelete[index] = true; // parallel version
                                }

                                __sync_fetch_and_add(&nbECRemoved, 1);

                            }
                        }
                    }
                }
            }
            TIME(auto end_thread_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t)));

            }); // parallel

    // now delete all nodes, in sequential (shouldn't take long)
    for (unsigned long i = 0; i < nbNodes; i++)
    {
        if (nodesToDelete[i])
        {
            _graph.deleteNode(i);
        }
    }

    bool debugThis = true;
    if (debugThis)
    {
        cout << nbECRemoved << " erroneous connections removed. " << endl;
        cout << nbSimplePaths << "/" << nbLongSimplePaths << "+" <<nbShortSimplePaths << " any=long+short simple path examined across all threads" << endl;
        double unit = 1000000000;
        cout.setf(ios_base::fixed);
        cout.precision(1);
        TIME(cout << "Timings: " << timeAll / unit << " CPUsecs total."<< endl);
    }

    return nbECRemoved;
}
