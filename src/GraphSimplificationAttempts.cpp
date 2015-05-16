
// those are unfinished attempt at coding some classical graph simplification algos:
// - erroneous edge removal
// - superbubbles


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

#if 0
/* superbubbles algorithm */
/* unfinished attempt */ /* made it in sequential only */
unsigned long GraphSimplification::removeSuperBubbles()
{
    unsigned long nbBubblesRemoved = 0;
    unsigned long nbCandidateBubbles = 0;
    unsigned long nbConsensusFailed = 0;
    unsigned long nbBubbleSimilarAbundance = 0;
    unsigned long nbValidationFailed = 0;
    unsigned long nbAlreadyPopped = 0;
    
    unsigned int k = _graph.getKmerSize();
    
    // constants are the same as legacyTraversal
    //unsigned int max_depth = std::max(500, 3* (int) _graph.getKmerSize()); // legacytraversal with a small change for depth: the max with 3k-1 (a bit arbitrary..)
    unsigned int coeff = 3;
    unsigned int additive_coeff = 100;
    unsigned int max_depth = std::max((unsigned int)((double)k * coeff), (unsigned int)(k + additive_coeff)); // SPAdes, exactly
    unsigned int max_breadth = 20; 

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, progressFormat1, ++_nbBubbleRemovalPasses);
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), buffer);

    // parallel stuff: create a dispatcher ; support atomic operations
/*    Dispatcher dispatcher (_nbCores);
    ISynchronizer* synchro = System::thread().newSynchronizer();*/

    MonumentTraversal * traversal = new MonumentTraversal(
            _graph,
            dummyTerminator,
            10000000,
            max_depth,
            max_breadth
            );
    LOCAL(traversal);


    vector<bool> visitedNodes;
    unsigned long nbNodes = itNode.size();
    visitedNodes.resize(nbNodes); // number of graph nodes
    for (unsigned long i = 0; i < nbNodes; i++)
        visitedNodes[i] = false;


    // sequential
    /** We loop over all nodes. */
    for (itNode.first(); !itNode.isDone(); itNode.next())
    {
        Node node = itNode.item();

    // parallel stuff
/*    vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
    unsigned long nbNodes = itNode.size();
    nodesToDelete.resize(nbNodes); // number of graph nodes
    for (unsigned long i = 0; i < nbNodes; i++)
        nodesToDelete[i] = false;

    set<Node::Value> bubblesAlreadyPopped; // endnodes of bubbles already popped (using kmers instead of Nodes because of possibly reversing them)

    dispatcher.iterate (itNode, [&] (Node& node)
    {
*/
        if (_graph.outdegree(node) <= 2)
            //return; // parallel
             continue; // sequential

        // we're at the start of a branching
        Node startingNode = node;
        Node previousNode = node; // dummy, no consequence  
        set<Node> all_involved_extensions;

        // do not pop the same bubble twice (and possibly both paths..)
        // so this is a compromise if check-late-to-avoid-synchro-overhead
/*        {
            LocalSynchronizer local(synchro);
            if (bubblesAlreadyPopped.find(startingNode.kmer) == bubblesAlreadyPopped.end())
                bubblesAlreadyPopped.insert(startingNode.kmer);
            else
                return;
        }
*/
        // pick a direction. for those complex nodes that might be forming two bubbles, let's just pop one of the two bubbles only, the other will follow on the other end
        Direction dir;
        if (_graph.outdegree(node) > 1)
            dir = DIR_OUTCOMING; // is the word "outcoming" really used in that context?
  /*      else
        {
            //dir = DIR_INCOMING;
            // until traversal with dir_incoming is buggy in gatbcore, let's not use pop bubbles using DIR_INCOMING
            dir = DIR_OUTCOMING;
            startingNode = _graph.reverse(startingNode); 
            previousNode = startingNode;

            // that's a good fix, but maybe someday investigate why some bubbles are only popped in one direction and not the other. it's probably because of longer tips..
        }
*/

        std::stack<Node> stack;

        bool cleanBubble = true;

        do  {

            Node node = stack.pop();
            
            unsigned long index = _graph.nodeMPHFIndex(node); 
            visitedNodes[index] = true;

            // TODO: continue here!

        }
        while (1);

        if (!cleanBubble)
            return; // parallel
        // continue; // sequential
        
        __sync_fetch_and_add(&nbCandidateBubbles, 1);

        /* now it's a clean bubble. find the most covered path */

        Node startNode = startingNode; 
        Node endNode = frontline.front().node;
        int traversal_depth = frontline.depth();

        bool success;

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
        /*{
            LocalSynchronizer local(synchro);
            if (bubblesAlreadyPopped.find(endNode.kmer) == bubblesAlreadyPopped.end())
                bubblesAlreadyPopped.insert(endNode.kmer);
            else
            {
                nbAlreadyPopped++;
                return;
            }
        }*/

        for (set<Node>::iterator itVecNodes = all_involved_extensions.begin(); itVecNodes != all_involved_extensions.end(); itVecNodes++)
        {
            //DEBUG(cout << endl << "deleting bubble node: " <<  _graph.toString (*itVecNodes) << endl);
            // _graph.deleteNode(*itVecNodes); // sequential

            unsigned long index = _graph.nodeMPHFIndex(*itVecNodes); // parallel version
            nodesToDelete[index] = true; // parallel version
        }
        
        __sync_fetch_and_add(&nbBubblesRemoved, 1);

         } // sequential
    //}); // parallel


    // now delete all nodes, in sequential
    /*for (unsigned long i = 0; i < nbNodes; i++)
        if (nodesToDelete[i])
           _graph.deleteNode(i);
*/

    bool debugThatPass = true;
    if (debugThatPass)
    {
        traversal->commit_stats();
        cout << nbCandidateBubbles << " candidate bubbles.\nAmong them, " << nbBubblesRemoved << " bubbles popped. " << endl;
        cout << "kept : " << maybe_print(nbConsensusFailed, "due to topology ; ") <<  

        cout << maybe_print(nbValidationFailed, "during path collapsing ; ") <<
            maybe_print(nbBubbleSimilarAbundance, "due to similar abundance ; ") << 
            maybe_print(nbAlreadyPopped, "found by another thread ") << 
    }

    return nbBubblesRemoved;
}

#endif
