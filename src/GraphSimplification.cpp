/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
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

using namespace std;

static const char* progressFormat0 = "removing tips,    pass %2d ";
static const char* progressFormat1 = "removing bubbles, pass %2d ";
static const char* progressFormat2 = "removing bulges,  pass %2d ";

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

        //DEBUG(cout << endl << "tip neighbor " << nbNeighbors << " : " <<  _graph.toString (neighbor) << " direction: " << neighbors[i].direction << " meancoverage: " <<simplePathCoverage << " over " << pathLen << " kmers" << endl);
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
 *
 * actual tip clipping code for "tc", i.e. implementations of the conditions, is in tip_clipper
 * (omni/graph_processing_algorithm.hpp:EdgeRemovingAlgorithm() is just generically parsing conditions)
 *
 * in one condition, tc_lb = 3.5, cb is 1000000 (minia's coverages values don't go that far.. yet..), rctc 2
 * and in another, tc_lb is 10 and cb is auto
 */
unsigned long GraphSimplification::removeTips()
{
    unsigned int k = _graph.getKmerSize();
    //int maxTipLength = _graph.getKmerSize() + 2; // in line with legacyTraversal
    
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

    dispatcher.iterate (itNode, [&] (Node& node)
    {

        unsigned inDegree = _graph.indegree(node), outDegree =  _graph.outdegree(node);
        if ((inDegree == 0 || outDegree == 0) && (inDegree != 0 || outDegree != 0))
        {
            if (_graph.isNodeDeleted(node)) { return; } // {continue;} // sequential and also parallel
            if (nodesToDelete[_graph.nodeMPHFIndex(node)]) { return; }  // parallel // actually not sure if really useful

            //DEBUG(cout << endl << "deadend node: " << _graph.toString (node) << endl);

            /** We follow the simple path to get its length */
            Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node.kmer); // so, it has a single neighbor
            Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[0].from, neighbors[0].direction); // .from should be .to but doesn't matter here
            //DEBUG(cout << endl << "neighbors from: " << _graph.toString (neighbors[0].from) << " direction: " << neighbors[0].direction << endl);

            bool isShortTopological = true;
            bool isShortRCTC = true;
            unsigned int pathLen = 1;
            vector<Node> nodes;
            nodes.push_back(node);
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                nodes.push_back(*itNodes);
                if (k + pathLen >= maxTipLengthTopological) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                    isShortTopological = false;

                if (k + pathLen >= maxTipLengthRCTC) 
                {
                    isShortRCTC= false;
                    break;
                }

                pathLen++;
            }

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
                
                // coverage of the putative tip
                unsigned long mean_abundance = 0;
                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                {
                    unsigned int abundance = _graph.queryAbundance((*itVecNodes).kmer);
                    mean_abundance += abundance;
                }
                double meanTipAbundance = (double)mean_abundance / (double)(nodes.size());
                double stdevTipAbundance = 0;

                // for debug only, can be disabled to save time
                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                {
                    unsigned int abundance = _graph.queryAbundance((*itVecNodes).kmer);
                    stdevTipAbundance += pow(fabs(abundance-meanTipAbundance),2);
                }
                stdevTipAbundance = sqrt(stdevTipAbundance/nodes.size());


                // explore the other two or more simple paths connected to that deadend, to get an abundance estimate
                // but first, get the branching node(s) the tip is connected to (it's weird when it's more than one branching node though)
                Graph::Vector<Node> connectedBranchingNodes = _graph.neighbors<Node>(nodes.back(), neighbors[0].direction);
                unsigned int nbBranchingNodes = 0;
                double meanNeighborsCoverage = 0;
                for (size_t j = 0; j < connectedBranchingNodes.size(); j++)
                {
                    meanNeighborsCoverage += getMeanAbundanceOfNeighbors(connectedBranchingNodes[j], nodes.back());
                    nbBranchingNodes++;
                }
                if (nbBranchingNodes > 0)
                    meanNeighborsCoverage /= nbBranchingNodes;

                isRCTCTip = (meanNeighborsCoverage > RCTCcutoff * meanTipAbundance);
                
                //DEBUG(cout << endl << "RCTCTip test, over " << nbBranchingNodes << " connected nodes. Global mean neighbors coverage: " << meanNeighborsCoverage <<  " compared to mean tip abundance over "<< nodes.size() << " values : " << meanTipAbundance << " stddev: " << stdevTipAbundance <<", is RCTC tip? " << isRCTCTip << endl);
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

/* note: the returned mean abundance does not include start and end nodes */
Path GraphSimplification::heuristic_most_covered_path(
        Direction dir, const Node startNode, const Node endNode, 
        int traversal_depth, bool& success, double &abundance, bool most_covered, bool backtracking, Node *avoidFirstNode)
{
    set<Node::Value> usedNode;
    usedNode.insert(startNode.kmer);
    Path current_path;
    current_path.start = startNode;
    success = false;
    vector<int> abundances; 
    unsigned long nbCalls = 0;

    Path res = heuristic_most_covered_path(dir, startNode, endNode, traversal_depth, current_path, usedNode, success, abundances, most_covered,
            backtracking, avoidFirstNode, nbCalls);

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
        int traversal_depth, Path current_path, set<Node::Value> usedNode, bool& success, vector<int>& abundances, bool most_covered,
        bool backtracking, Node *avoidFirstNode, unsigned long &nbCalls)
{
    // inspired by all_consensuses_between
    nbCalls++;

    if (traversal_depth < -1)
    {
        success = false;
        return current_path;
    }

    if (startNode.kmer == endNode.kmer)
    {
        success = true;
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
            success = false;
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
            backtracking,
            NULL, // no longer avoid nodes
            nbCalls
        );

        if (success ||(!backtracking)) // on success or if backtracking is disallowed, return immediately, don't backtrack
        {
            abundances = extended_abundances;
            return new_path; 
        }
    }

    return current_path;


}

/* same treatment as tip clipping: let's take some information from spades. 
 * it's interesting to see that SPAdes 3.5 does not do identity-based bubble popping, at all! 
 * actually, it's not even doing bubble popping. it's bulge popping.
 * it pops bulges based on something like the ratio between most examined simple path and a more covered path is (whether it is above 1.1).
 * 
 * not sure if this ratio is better than identity.
 * anyhow it'd rather do less strict identity bubble popping. legacy minia has 90%, i lowered it to 80%. 
 * this present algo still failed to pop most bubbles in the buchnera population test (test/buchnera_test.sh), because of difficulty to enumerate all paths in a bubble
 * so i switched to a heuristic path finding algorithm (didn't look at spades')
 * still not ideal, as FrontlineReachable missed many bubbles. that led to a way-increased assembly size (1.3 Mbp instead of 600kbp), 
 * that can be decreased to 900k if FrontlineReachable check() is forced to always true (but also pops bad bubbles)
 * 
 * that being said, let's attempt to pop bulges instead. seems better! 
 */

unsigned long GraphSimplification::removeBubbles()
{
    unsigned long nbBubblesRemoved = 0;
    unsigned long nbCandidateBubbles = 0;
    unsigned long nbConsensusFailed = 0;
    unsigned long nbBubbleSimilarAbundance = 0;
    unsigned long nbValidationFailed = 0;
    unsigned long nbAlreadyPopped = 0;
    
    unsigned int k = _graph.getKmerSize();
    
    // constants are the same as legacyTraversal
    unsigned int max_depth = std::max(500, 3* (int) _graph.getKmerSize()); // legacytraversal with a small change for depth: the max with 3k-1 (a bit arbitrary..)
    unsigned int max_breadth = 20; 

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, progressFormat1, ++_nbBubbleRemovalPasses);
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);
    ISynchronizer* synchro = System::thread().newSynchronizer();

    Terminator& dummyTerminator = NullTerminator::singleton(); // Frontline wants one

    MonumentTraversal * traversal = new MonumentTraversal(
            _graph,
            dummyTerminator,
            10000000,
            max_depth,
            max_breadth
            );
    LOCAL(traversal);

    // sequential
    /** We loop over all nodes. */
    //for (itNode.first(); !itNode.isDone(); itNode.next())
    //{
      //  Node node = itNode.item();

    // parallel stuff
    vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
    unsigned long nbNodes = itNode.size();
    nodesToDelete.resize(nbNodes); // number of graph nodes
    for (unsigned long i = 0; i < nbNodes; i++)
        nodesToDelete[i] = false;

    set<Node::Value> bubblesAlreadyPopped; // endnodes of bubbles already popped (using kmers instead of Nodes because of possibly reversing them)

    dispatcher.iterate (itNode, [&] (Node& node)
    {

        if (_graph.outdegree(node) <= 1 && _graph.indegree(node) <= 1)
            return; // parallel
            // continue; // sequential

        // we're at the start of a branching
        Node startingNode = node;
        Node previousNode = node; // dummy, no consequence  
        set<Node> all_involved_extensions;

        // do not pop the same bubble twice (and possibly both paths..)
        // so this is a compromise if check-late-to-avoid-synchro-overhead
        {
            LocalSynchronizer local(synchro);
            if (bubblesAlreadyPopped.find(startingNode.kmer) == bubblesAlreadyPopped.end())
                bubblesAlreadyPopped.insert(startingNode.kmer);
            else
                return;
        }

        // pick a direction. for those complex nodes that might be forming two bubbles, let's just pop one of the two bubbles only, the other will follow on the other end
        Direction dir;
        if (_graph.outdegree(node) > 1)
            dir = DIR_OUTCOMING; // is the word "outcoming" really used in that context?
        else
        {
            //dir = DIR_INCOMING;
            // until traversal with dir_incoming is buggy in gatbcore, let's not use pop bubbles using DIR_INCOMING
            dir = DIR_OUTCOMING;
            startingNode = _graph.reverse(startingNode); 
            previousNode = startingNode;

            // that's a good fix, but maybe someday investigate why some bubbles are only popped in one direction and not the other. it's probably because of longer tips..
        }

        FrontlineReachable frontline (dir, _graph, dummyTerminator, startingNode, previousNode, &all_involved_extensions);
        // one would think using FrontlineBranching will pop more bubbles less conservatively. actually wasn't the case in my early tests. need to test more.

        bool cleanBubble = true;

        do  {
            bool should_continue = frontline.go_next_depth();
            if (!should_continue) 
            {
                cleanBubble = false;
                break;
            }

            // don't allow a depth too large
            if (frontline.depth() > max_depth)
            {  
                //    stats.couldnt_traverse_bubble_depth++;
                DEBUG(cout << endl << "Candidate bubble from node " <<  _graph.toString(startingNode) << " frontline exceeds depth" << endl);
                cleanBubble = false;
                break;
            }

            // don't allow a breadth too large
            if (frontline.size()> max_breadth)
            {  
                //    stats.couldnt_traverse_bubble_breadth++;
                DEBUG(cout << endl << "Candidate bubble from node " <<  _graph.toString(startingNode) << " frontline exceeds breadth" << endl);
                cleanBubble = false;
                break;
            }

            // stopping condition: frontline is either empty, or contains only 1 kmer
            // needs the kmer to be non-branching, in order to avoid a special case of bubble immediatly after a bubble
            // affects mismatch rate in ecoli greatly
            if (frontline.size() == 0)  
            {
                //    stats.couldnt_find_extension++;
                DEBUG(cout << endl << "Candidate bubble from node " <<  _graph.toString(startingNode) << " frontline empty" << endl);
                cleanBubble = false;
                break;
            }

            if (frontline.size() == 1) {break;}
        }
        while (1);

        if (frontline.size()!=1)
        {
            DEBUG(cout << endl << "Candidate bubble from node " <<  _graph.toString(startingNode) << " frontline ends of size " << frontline.size() << endl);
            cleanBubble = false;
        }

        cleanBubble &= frontline.isReachable();

        if (!cleanBubble)
            return; // parallel
        // continue; // sequential
        
        __sync_fetch_and_add(&nbCandidateBubbles, 1);

        Node startNode = startingNode; 
        Node endNode = frontline.front().node;
        int traversal_depth = frontline.depth();

        bool success;

        // not exploring all consensuses anymore, instead, we'll greedily search directly for the most abundant one
        // hope that it will be faster.
        // turns out it's also more accurate!
#if 0
        // code taken from Traversal
        //
        // find all consensuses between start node and end node
        set<Path> consensuses = traversal->all_consensuses_between (dir, startNode, endNode, traversal_depth+2, success);

        // if consensus phase failed, stop
        if (!success)
        {
            __sync_fetch_and_add(&nbConsensusFailed, 1);
            return;
        }

        Path consensus_exhaust;
        consensus.resize (0);
        // validate paths, based on identity
       bool validated = traversal->validate_consensuses (consensuses, consensus_exhaust);
        if (!validated){
            __sync_fetch_and_add(&nbValidationFailed, 1);
            return; 
        }

        // ready to pop the bubble and keep only the validated consensus
        DEBUG(cout << endl << "READY TO POP!\n");

        string p_str_exhaust = path2string(dir, consensus, endNode);
#endif

        // another possibly faster method, and possibly more accurate, than above
        double mean_abundance_most_covered;
        double mean_abundance_least_covered;
        Path heuristic_p_most = heuristic_most_covered_path(dir, startNode, endNode, traversal_depth+2, success, mean_abundance_most_covered);

        if (!success)
        {
            __sync_fetch_and_add(&nbConsensusFailed, 1);
            return;
        }

        // now get the least covered path
        Path heuristic_p_least = heuristic_most_covered_path(dir, startNode, endNode, traversal_depth+2, success, mean_abundance_least_covered, false);
 
        if (!success)
        {
            __sync_fetch_and_add(&nbConsensusFailed, 1);
            return;
        }

        //cout << "most/least covered path: " << mean_abundance_most_covered << "/" << mean_abundance_least_covered << endl;

        if (mean_abundance_most_covered < 1.1 * mean_abundance_least_covered)
        {
            __sync_fetch_and_add(&nbBubbleSimilarAbundance, 1);
            return;
        }

        string p_str_heur = path2string(dir, heuristic_p_most, endNode);

#if 0 
        // comparison with both bubble exporing codes (exhaustive but bounded (ie minia 1) vs greedy but unbounded)
        if (p_str != p_str_heur)
        {
            cout << "heuristic path doesn't agree with validated path: " << endl << p_str << endl << p_str_heur << endl;
        }
        else
            cout << "heuristic path AGREE!" << endl;
#endif

        string p_str = p_str_heur;
        Path consensus = heuristic_p_most;

        DEBUG(cout << "consensus string to keep: " << p_str << endl);

        for (size_t i = 0; i < consensus.size(); i++)
        {            
            Node node = _graph.buildNode((char *)(p_str.c_str()), i); 
            int nb_erased = all_involved_extensions.erase(node.kmer);
            //if (nb_erased != 1)
            //    cout << "error: wanted to keep kmer" << _graph.toString (node) << "but wasn't there" << endl;

            //DEBUG(cout << endl << "keeping bubble node: " <<  _graph.toString (node) << endl);
        }
        all_involved_extensions.erase(startNode.kmer);
        all_involved_extensions.erase(endNode.kmer);

        // do not pop the same bubble twice (and possibly both paths..)
        // same code as above but check end node too, important! 
        {
            LocalSynchronizer local(synchro);
            if (bubblesAlreadyPopped.find(endNode.kmer) == bubblesAlreadyPopped.end())
                bubblesAlreadyPopped.insert(endNode.kmer);
            else
            {
                nbAlreadyPopped++;
                return;
            }
        }

        for (set<Node>::iterator itVecNodes = all_involved_extensions.begin(); itVecNodes != all_involved_extensions.end(); itVecNodes++)
        {
            //DEBUG(cout << endl << "deleting bubble node: " <<  _graph.toString (*itVecNodes) << endl);
            // _graph.deleteNode(*itVecNodes); // sequential

            unsigned long index = _graph.nodeMPHFIndex(*itVecNodes); // parallel version
            nodesToDelete[index] = true; // parallel version
        }
        
        __sync_fetch_and_add(&nbBubblesRemoved, 1);

        // } // sequential
    }); // parallel


    // now delete all nodes, in sequential
    for (unsigned long i = 0; i < nbNodes; i++)
        if (nodesToDelete[i])
           _graph.deleteNode(i);

    bool debugThatPass = true;
    if (debugThatPass)
    {
        traversal->commit_stats();
        cout << nbCandidateBubbles << " candidate bubbles.\nAmong them, " << nbBubblesRemoved << " bubbles popped. " << endl;
        cout << "kept : " << maybe_print(nbConsensusFailed, "due to topology ; ") <<  
            // those belong to those kept due to topology
            maybe_print(traversal->final_stats.couldnt_consensus_negative_depth, "abnormal depth ; ")  <<
            maybe_print(traversal->final_stats.couldnt_consensus_loop, "loop in consensus generation ; ") << 
            maybe_print(traversal->final_stats.couldnt_consensus_amount, "too many paths ; ");

        cout << maybe_print(nbValidationFailed, "during path collapsing ; ") <<
            maybe_print(nbBubbleSimilarAbundance, "due to similar abundance ; ") << 
            maybe_print(nbAlreadyPopped, "found by another thread ") << 
            maybe_print(traversal->final_stats.couldnt_validate_bubble_mean_depth, "high mean length, ") << 
            maybe_print(traversal->final_stats.couldnt_validate_bubble_deadend, "just one long path, ") << 
            maybe_print(traversal->final_stats.couldnt_validate_bubble_stdev, "not roughly same length, ") << 
            maybe_print(traversal->final_stats.couldnt_validate_bubble_identity, "not identical enough, ") << 
            maybe_print(traversal->final_stats.couldnt_validate_bubble_long_chosen, "too long chosen consensus.") << endl;
    }

    return nbBubblesRemoved;
}


/* bulge removal algorithm. mimics spades, which doesnt remove bubbles, but only bulges. looks as effective.
 * it's slow to do heuristic_find_most_covered path so i'm testing it with no backtracking
 *
 * see a-b-c here for an explanation of bulge removal: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791033/figure/f4/
 */ 
unsigned long GraphSimplification::removeBulges()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int coeff = 3;
    unsigned int additive_coeff = 100;
    unsigned int maxBulgeLength = std::max((unsigned int)((double)k * coeff), (unsigned int)(k + additive_coeff)); // SPAdes, exactly

    // stats
    //
    unsigned long nbBulgesRemoved = 0;
    unsigned long nbSimplePaths = 0;
    unsigned long nbLongSimplePaths = 0;
    unsigned long nbShortSimplePaths = 0;
    unsigned long nbTopologicalBulges = 0;
    unsigned long nbFirstNodeDeleted = 0;
    unsigned long nbFirstNodeGraphDeleted = 0;
    unsigned long nbNoAltPathBulges = 0;
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

    dispatcher.iterate (itNode, [&] (Node& node)
    {

    auto start_thread_t=get_wtime();


      // need to search in both directions
      for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
      {
         if ((_graph.outdegree(node) >= 2 && dir == DIR_OUTCOMING) || (_graph.indegree(node) >= 2 && dir == DIR_INCOMING))
         {
            TIME(auto start_various_overhead_t=get_wtime());
            unsigned long index =_graph.nodeMPHFIndex(node);
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

                TIME(auto start_simplepath_t=get_wtime());
                Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[i].to, dir);
                DEBUG(cout << endl << "neighbors " << i << "/" << neighbors.size() << " from: " << _graph.toString (neighbors[i].to) << " dir: " << dir << endl);
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
                Node endNode = outneighbors[0].to;

                DEBUG(cout << endl << "last simple path node: "<< _graph.toString(nodes.back()) << " has " << outneighbors.size() << " outneighbors" << endl);

                // at this point, the last node in "nodes" is the last node of a potential Bulge path, and endNode is hopefully a branching node right after.
                // check if it's connected to something that has in-branching. 
                bool isDoublyConnected = (_graph.indegree(endNode) > 1 && dir==DIR_OUTCOMING) || (dir==DIR_INCOMING && _graph.outdegree(endNode) > 1);

                bool isTopologicalBulge = isDoublyConnected;

                DEBUG(cout << endl << "pathlen: " << pathLen << " last node neighbors size: " << _graph.neighbors<Edge>(nodes.back()).size() << " indegree/outdegree: " <<_graph.indegree(nodes.back()) << "/" << _graph.outdegree(nodes.back()) << " istopobulge: " << isTopologicalBulge << endl);

                TIME(end_various_overhead_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

                if (!isTopologicalBulge)
                    continue;

                __sync_fetch_and_add(&nbTopologicalBulges, 1);

                unsigned int depth = std::max((unsigned int)(pathLen * 1.1),(unsigned int) 3); // following SPAdes
                double mean_abundance_most_covered;
                bool success;
                Node startNode = node;

                TIME(auto start_pathfinding_t=get_wtime());

                Path heuristic_p_most = heuristic_most_covered_path(dir, startNode, endNode, depth+2, success, mean_abundance_most_covered,
                        true, // most covered
                        false, // avoid backtracking
                        &(neighbors[i].to) // avoid that node
                        );

                TIME(auto end_pathfinding_t=get_wtime());
                TIME(__sync_fetch_and_add(&timePathFinding, diff_wtime(start_pathfinding_t,end_pathfinding_t)));
                TIME(auto start_post_t=get_wtime());

                if (!success)
                {
                    TIME(__sync_fetch_and_add(&timeFailedPathFinding, diff_wtime(start_pathfinding_t,end_pathfinding_t)));
                    TIME(if (diff_wtime(start_pathfinding_t,end_pathfinding_t) > timeLongestFailure) { timeLongestFailure = diff_wtime(start_pathfinding_t,end_pathfinding_t); longestFailureDepth = depth;});
                    __sync_fetch_and_add(&nbNoAltPathBulges, 1);
                    continue;
                }

                DEBUG(cout << endl << "alternative path is:  "<< path2string(dir, heuristic_p_most, endNode)<< " abundance: "<< mean_abundance_most_covered <<endl);

                bool debug = false;
                if (debug)
                {
                    double mean_abundance_least_covered;
                    Path heuristic_p_least = heuristic_most_covered_path(dir, startNode, endNode, depth+2, success, mean_abundance_least_covered,false);
                    DEBUG(cout << endl << "alternative least is: "<< path2string(dir, heuristic_p_least, endNode)<< " abundance: "<< mean_abundance_least_covered <<endl);
                }

                unsigned int dummyLen;
                double simplePathCoverage = getSimplePathCoverage(nodes[1], dir, &dummyLen);

                DEBUG(cout << endl << "retraced bulge path over length: " << dummyLen << endl);

                bool isBulge =  simplePathCoverage * 1.1  <=  mean_abundance_most_covered;

                DEBUG(cout << endl << "bulge coverages: " << simplePathCoverage<< "/" <<  mean_abundance_most_covered  << endl);

                if (!isBulge)
                {
                    __sync_fetch_and_add(&nbBadCovBulges, 1);
                    DEBUG(cout << endl << "not a bulge due to similar coverage" << endl);
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
            // } // sequential
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
    cout << nbNoAltPathBulges << " without alt. path, " << nbBadCovBulges << " didn't satisfy cov. criterion." << endl;

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);
    TIME(cout << "Timings: " << timeAll / unit << " CPUsecs total."<< endl);
    TIME(cout << "         " << timeSimplePath / unit << " CPUsecs simple path traversal." << endl);
    TIME(cout << "         " << timePathFinding / unit << "(/" << timePathFinding / unit << ") CPUsecs path-finding(/failed). Longest: " << timeLongestFailure / unit << " CPUsecs (depth " << longestFailureDepth << ")." << endl);
    TIME(cout << "         " << timePost / unit << " CPUsecs topological bulge processing, " << timeDelete / unit << " CPUsecs nodes deletion." << endl);
    TIME(cout << "         " << timeVarious / unit << " CPUsecs various overhead." << endl);

    return nbBulgesRemoved;
}

