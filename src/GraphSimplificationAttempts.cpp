
// those are unfinished attempt at coding some classical graph simplification algos:
// - bubbles
// - superbubbles


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

#if 0

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

        int success;

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
