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

#include <vector>
#include <set>
#include <GraphSimplification.hpp>

using namespace std;

static const char* progressFormat0 = "Minia : removing tips,    pass %2d";
static const char* progressFormat1 = "Minia : removing bubbles, pass %2d";

double GraphSimplification::getSimplePathCoverage(Node node, Direction dir, unsigned int *pathLenOut, unsigned int maxLength)
{
    Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (node, dir);
    unsigned long mean_abundance = _graph.queryAbundance(node.kmer);
    unsigned int pathLen = 1;
    for (itNodes.first(); !itNodes.isDone(); itNodes.next())
    {
        Node node = *itNodes;
        unsigned int abundance = _graph.queryAbundance(node.kmer);
        mean_abundance += abundance;
        pathLen++;
        if (maxLength > 0 && pathLen >= maxLength)
            break;
    }
    *pathLenOut = pathLen;
    return (double)mean_abundance / (double)pathLen;
}

// gets the mean abundance of neighboring paths around a branching node (excluding the path that starts with nodeToExclude, e.g. the tip itself)
double GraphSimplification::getMeanAbundanceOfNeighbors(Node branchingNode, Node nodeToExclude)
{
    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(branchingNode);
    unsigned int nbNeighbors = 0;
    double meanNeighborsCoverage = 0;
    DEBUG(cout << endl << "called getMeanAbudanceOfNeighbors for node " << _graph.toString(branchingNode) << " of degrees " << _graph.indegree(branchingNode) <<"/"<< _graph.outdegree(branchingNode)<< " excluding node  " <<  _graph.toString (nodeToExclude) << endl);
    for (size_t i = 0; i < neighbors.size(); i++)
    {
        Node neighbor = neighbors[i].to;
        if (neighbor == nodeToExclude) // (in gatb-core, Node == Node means Node.kmer == Node.kmer)
        {
            DEBUG(cout << endl << "good, seen the node to exclude" << endl);
            continue; 
        }

        unsigned int pathLen;
        unsigned int simplePathCoverage = getSimplePathCoverage(neighbor, neighbors[i].direction, &pathLen, 100);
        meanNeighborsCoverage += simplePathCoverage;
        nbNeighbors++;

        DEBUG(cout << endl << "tip neighbor " << nbNeighbors << " : " <<  _graph.toString (neighbor) << " direction: " << neighbors[i].direction << " meancoverage: " <<simplePathCoverage << " over " << pathLen << " kmers" << endl);
    }
    meanNeighborsCoverage /= nbNeighbors;
    return meanNeighborsCoverage;
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

            DEBUG(cout << endl << "deadend node: " << _graph.toString (node) << endl);

            /** We follow the simple path to get its length */
            Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node.kmer); // so, it has a single neighbor
            Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[0].from, neighbors[0].direction);
            DEBUG(cout << endl << "neighbors from: " << _graph.toString (neighbors[0].from) << " direction: " << neighbors[0].direction << endl);

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

            bool isTopologicalTip = isShortTopological && isConnected; 
            bool isMaybeRCTCTip = isShortRCTC && isConnected;

            DEBUG(cout << endl << "pathlen: " << pathLen << " last node " << _graph.toString(nodes.back()) << " neighbors in/out: " <<_graph.indegree(nodes.back()) << " " << _graph.outdegree(nodes.back()) << " istopotip: " << isTopologicalTip << endl);

            bool isRCTCTip = false;
            if (!isTopologicalTip && isMaybeRCTCTip)
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
                    // maybe use that code later, or clean it up
                    //if (nodes.size() > 1 && connectedBranchingNodes[j] == nodes[nodes.size() - 2])
                    //meanNeighborsCoverage += getMeanAbundanceOfNeighbors(nodes.back(), ( (nodes.size() > 1) ? nodes[nodes.size() - 2] : nodes.back()));
                    meanNeighborsCoverage += getMeanAbundanceOfNeighbors(connectedBranchingNodes[j], nodes.back());
                    nbBranchingNodes++;
                }
                if (nbBranchingNodes > 0)
                    meanNeighborsCoverage /= nbBranchingNodes;

                isRCTCTip = (meanNeighborsCoverage > RCTCcutoff * meanTipAbundance);
                
                DEBUG(cout << endl << "RCTCTip test, over " << nbBranchingNodes << " connected nodes. Global mean neighbors coverage: " << meanNeighborsCoverage <<  " compared to mean tip abundance over "<< nodes.size() << " values : " << meanTipAbundance << " stddev: " << stdevTipAbundance <<", is RCTC tip? " << isRCTCTip << endl);
            }

            bool isTip = isTopologicalTip || isRCTCTip; 


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


unsigned long GraphSimplification::removeBubbles()
{
    unsigned long nbBubblesRemoved = 0;
    
    // constants are the same as legacyTraversal
    // small change for depth: the max with 3k-1 (a bit arbitrary..)
    unsigned int max_depth = std::max(500, 3* (int) _graph.getKmerSize());
    unsigned int max_breadth = 20; 

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, progressFormat1, ++_nbBubbleRemovalPasses);
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);
    ISynchronizer* synchro = System::thread().newSynchronizer();

    Terminator& dummyTerminator = NullTerminator::singleton(); // Frontline wants one

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
                break;

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

        if (!cleanBubble)
            return; // parallel
        // continue; // sequential

        Node startNode = startingNode; 
        Node endNode = frontline.front().node;
        int traversal_depth = frontline.depth();

        // code taken from Traversal
        MonumentTraversal * traversal = new MonumentTraversal(
                _graph,
                dummyTerminator,
                10000000,
                max_depth,
                max_breadth
                );
        LOCAL(traversal);

        // find all consensuses between start node and end node
        bool success;
        set<Path> consensuses = traversal->all_consensuses_between (dir, startNode, endNode, traversal_depth+1, success);

        // if consensus phase failed, stop
        if (!success)
            return; // parallel
        // continue; // sequential

        Path consensus;
        consensus.resize (0);
        // validate paths, based on identity
        /* small digression. it's interesting to see that SPAdes 3.5 does not do identity-based bubble popping, at all! 
         * it pops bubble based on something like the ratio between most covered and less covered path is (whether it is above 1.1).
         * not sure if it's betetr. anyhow it'd rather do less strict identity bubble popping. legacy minia has 90%, i'll lower to 85%. 
         * I still see 80% bubbles in coli (SRR001665_1.fa) with k=22 for instance.
         */
        bool validated = traversal->validate_consensuses (consensuses, consensus);
        if (!validated)   
            return; // parallel
        // continue; // sequential

        // ready to pop the bubble and keep only the validated consensus
        DEBUG(cout << endl << "READY TO POP!\n");

        // also taken from Traversal
        // naive conversion from path to string
        Path p = consensus;

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

        DEBUG(cout << "consensus string to keep: " << p_str << endl);

        for (size_t i = 0; i < p.size(); i++)
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
                return;
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


    return nbBubblesRemoved;
}


// maybe.. someday.. will finish implementing that. not persuaded it's fully important for CAMI
#if 0

/* Again, let's see what spades 3.5 does.
 *erroneous connection remover is:

 RemoveLowCoverageEdges
 calls config then calls:
  omnigraph::EdgeRemovingAlgorithm which is same as tc

  to_ec_lb in ../src/debruijn/simplification/simplification_settings.hpp
  is exactly like tip clipping but with a different length (2 * (max tip length with coeff 5) - 1, so that's exactly 10*min(k,readlen/2) - 1)

  icb is something that removes edges with coverage cov_bound / nb_iters * iter
  where cov_bound is same as in cb

  and it's asking for AlternativesPresenceCondition:

    o   o                           o-->-O
   /     \                             /
   O-->---O    drawn differently:     /
                                    O-->-o
*/
unsigned long GraphSimplification::removeErroneousConnections()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int maxTipLength = (unsigned int)((float)k * (10 - 1.0)) * 3   ;  // SPAdes mode

    unsigned long nbECRemoved = 0;

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, progressFormat0, ++_nbECRemovalPasses);
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

        if (_graph.outdegree(node) == 2)
        {
            if (_graph.isNodeDeleted(node)) { return; } // {continue;} // sequential and also parallel
            if (nodesToDelete[_graph.nodeMPHFIndex(node)]) { return; }  // parallel // actually not sure if really useful

            DEBUG(cout << endl << "putative EC node: " << _graph.toString (node) << endl);

            /** We follow the outgoing simple paths to get their length and last neighbor */
            Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node.kmer);

            vector<Node> nodes;
            bool foundShortPath = false;
            for (int i = 0; i < neighbors.size(); i++)
            {
                nodes.clear();
                Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[i].from, neighbors[i].direction);
                DEBUG(cout << endl << "neighbors from: " << _graph.toString (neighbors[0].from) << endl);
                bool isShort = true;
                unsigned int pathLen = 1;
                nodes.push_back(node);
                for (itNodes.first(); !itNodes.isDone(); itNodes.next())
                {
                    nodes.push_back(*itNodes);
                    if (k + pathLen++ >= maxECLength) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                    {
                        isShort = false;
                        break;       
                    }
                }
                if (isShort)
                {
                    foundShortPath = true;
                    break;
                }
            }

            if (!foundShortPath || pathLen == 1) // can't do much if it's pathLen=1, we don't support edge removal, only node removal
                return;

            // at this point, the last node in "nodes" is the last node of a potential EC.
            // check if it's connected to something that has in-branching and a out-neighbor. 
            bool isDoublyConnected = (_graph.outdegree(nodes.back) >= 1 && _graph.indegree(nodes.back) > 1);

            bool isTopologicalEC = isDoublyConnected;

            DEBUG(cout << endl << "pathlen: " << pathLen << " last node neighbors size: " << _graph.neighbors<Edge>(nodes.back()).size() << " indegree outdegree: " <<_graph.indegree(node) << " " << _graph.outdegree(node) << "istopotip: " << isTopologicalTip << endl);

            if (isTopologicalEC)
            {
                
#if 0
                // coverage criterion
                unsigned long mean_abundance = 0;
                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                {
                    unsigned int abundance = _graph.queryAbundance(node.kmer);
                    mean_abundance += abundance;
                }

                // explore the other two or more simple paths connected to that deadend to get a abundance estimate
                Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(nodes.back());
                for (size_t i = 0; i < neighbors.size(); i++)
                {
                    Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[i].to, neighbors[i].direction);
                    int pathLen = 0;
                    vector<Node> simplepath_nodes;
                    for (itNodes.first(); !itNodes.isDone(); itNodes.next())
                    {
                        simplepath_nodes.push_back(*itNodes);
                        if (pathLen++ >= 100)
                        {
                            isShort = false;
                            break;       
                        }
                    }
                }
#endif


                bool isTip = true; // TODO maybe: (mean_abundance < local_min_abundance); but in fact, we didn't use that in legacyTraversal. maybe it's a mistake.


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

#endif
