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
#include <GraphSimplification.hpp>

using namespace std;

static const char* progressFormat0 = "removing tips          ";
static const char* progressFormat1 = "removing bubbles       ";

unsigned long GraphSimplification::removeTips()
{
    int maxTipLength = _graph.getKmerSize() + 2; // in line with legacyTraversal

    unsigned long nbTipsRemoved = 0;

    /** We get an iterator over all nodes . */
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), progressFormat0);

    /** We loop over all nodes. */
    for (itNode.first(); !itNode.isDone(); itNode.next())
    {
        Node node = itNode.item();

        if (_graph.indegree(node) + _graph.outdegree(node) == 1)
        {
            if (_graph.isNodeDeleted(node)) { continue; }
            //DEBUG(cout << endl << "deadend node: " << _graph.toString (node) << endl);

            /** We follow the simple path to get its length */
            Graph::Vector<Edge> neighbors = _graph.neighbors<Edge>(node.kmer); // so, it has a single neighbor
            Graph::Iterator <Node> itNodes = _graph.simplePath<Node> (neighbors[0].from, neighbors[0].direction);
            bool isShort = true;
            int pathLen = 2; // node and neighbors[0].to
            vector<Node> nodes;
            nodes.push_back(node);
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                nodes.push_back(*itNodes);
                if (pathLen++ >= maxTipLength)
                {
                    isShort = false;
                    break;       
                }
            }

            bool isTopologicalTip = isShort && (_graph.neighbors<Edge>(nodes.back()).size() > 1); // is not a deadend on the other side

            if (isTopologicalTip)
            {
                
                // coverage criterion
                unsigned long mean_abundance = 0;
                for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                {
                    unsigned int abundance = _graph.queryAbundance(node.kmer);
                    mean_abundance += abundance;
                }

#if 0
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

                    //DEBUG(cout << endl << "TIP of length " << pathLen << " FOUND: " <<  _graph.toString (node) << endl);
                    for (vector<Node>::iterator itVecNodes = nodes.begin(); itVecNodes != nodes.end(); itVecNodes++)
                    {
                        //DEBUG(cout << endl << "deleting tip node: " <<  _graph.toString (*itVecNodes) << endl);
                        _graph.deleteNode(*itVecNodes);
                    }
                    nbTipsRemoved++;
                }
            }
        }
    }
    return nbTipsRemoved;
}


unsigned long GraphSimplification::removeBubbles()
{
    unsigned long nbBubblesRemoved = 0;
    
    // constants are the same as legacyTraversal
    // small change for depth: the max with 3k-1 (a bit arbitrary..)
    int max_depth = std::max(500, _graph.getKmerSize() * 3);
    int max_breadth = 20; 

    /** We get an iterator over all nodes . */
    ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (_graph.iterator<Node>(), progressFormat1);

    Terminator& dummyTerminator = NullTerminator::singleton(); // Frontline wants one

    /** We loop over all nodes. */
    for (itNode.first(); !itNode.isDone(); itNode.next())
    {
        Node node = itNode.item();

        if (_graph.outdegree(node) <= 1 && _graph.indegree(node) <= 1)
            continue;

        // we're at the start of a branching
        Node startingNode = node;
        Node previousNode = node; // dummy, no consequence  
        set<Node> all_involved_extensions;

        // pick a direction. for those complex nodes that might be forming two bubbles, let's just pop one of the two bubbles only, the other will follow on the other end
        Direction dir;
        if (_graph.outdegree(node) > 1)
            dir = DIR_OUTCOMING; // is the word "outcoming" really used in that context?
        else
            dir = DIR_INCOMING;

        FrontlineReachable frontline (dir, _graph, dummyTerminator, startingNode, previousNode, &all_involved_extensions);
        // one would think using FrontlineBranching will pop more bubbles less conservatively. actually wasn't the case in my early tests. need to test more.

        bool cleanBubble = true;

        do  {
            bool should_continue = frontline.go_next_depth();
            if (!should_continue) 
                break;

            // don't allow a depth too large
            if ((int)frontline.depth() > max_depth)
            {  
                //    stats.couldnt_traverse_bubble_depth++;
                DEBUG(cout << endl << "Candidate bubble from node " <<  _graph.toString(startingNode) << " frontline exceeds depth" << endl);
                cleanBubble = false;
                break;
            }

            // don't allow a breadth too large
            if ((int)frontline.size()> max_breadth)
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

            if (frontline.size() == 1) // in legacyTraversal: longer contigs but for some reason, higher mismatch rate; TODO: verify it in the context of graph simplification
                {break;}
            //if (frontline.size() == 1 &&   (!terminator.isEnabled() || !terminator.is_branching(frontline.front().node)) )  {   break;  }
        }
        while (1);

        if (frontline.size()!=1)
        {
            DEBUG(cout << endl << "Candidate bubble from node " <<  _graph.toString(startingNode) << " frontline ends of size " << frontline.size() << endl);
            cleanBubble = false;
        }

        printf("node %s is a bubble? %d\n",_graph.toString(startingNode).c_str(),cleanBubble);

        if (!cleanBubble)
            continue;
       
        Node startNode = node; 
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
        if (!success)  {
            continue;
        }

        Path consensus;
        consensus.resize (0);
        // validate paths, based on identity
        bool validated = traversal->validate_consensuses (consensuses, consensus);
        if (!validated)   
        {  
            //stats.couldnt_validate_consensuses++;
            continue;
        }

        // ready to pop the bubble and keep only the validated consensus
        DEBUG(cout << endl << "READY TO POP!\n");

        // also taken from Traversal
        // naive conversion from path to string
        Path p = consensus;
        
        string p_str;

        if (dir == DIR_INCOMING)
        {
            // TODO: port that change into Traversal
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
            if (nb_erased != 1)
                cout << "error: wanted to keep kmer" << _graph.toString (node) << "but wasn't there" << endl;
            //DEBUG(cout << endl << "keeping bubble node: " <<  _graph.toString (node) << endl);
        }
        all_involved_extensions.erase(startNode.kmer);
        all_involved_extensions.erase(endNode.kmer);

        for (set<Node>::iterator itVecNodes = all_involved_extensions.begin(); itVecNodes != all_involved_extensions.end(); itVecNodes++)
        {
            //DEBUG(cout << endl << "deleting bubble node: " <<  _graph.toString (*itVecNodes) << endl);
            _graph.deleteNode(*itVecNodes);
        }

        nbBubblesRemoved++;
    }
    return nbBubblesRemoved;
}
