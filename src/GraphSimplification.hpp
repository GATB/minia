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

#ifndef _GATB_TOOLS_GRAPH_SIMPLIFICATION_HPP_
#define _GATB_TOOLS_GRAPH_SIMPLIFICATION_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/


/** */
class GraphSimplification : public SmartPointer
{
public:

    GraphSimplification (const Graph& graph, int nbCores);

    unsigned long removeTips();
    unsigned long removeBubbles();
    unsigned long removeBulges();

    double getSimplePathCoverage(Node node, Direction dir, unsigned int* pathLen = NULL, unsigned int maxLength = 0);
    double getMeanAbundanceOfNeighbors(Node branchingNode, Node nodeToExclude);

    int _nbTipRemovalPasses;
    int _nbBubbleRemovalPasses;
    int _nbBulgeRemovalPasses;

protected:
    const Graph&  _graph;
    int _nbCores;
    bool _firstNodeIteration;

    string path2string(Direction dir, Path p, Node endNode);

    Path heuristic_most_covered_path(Direction dir, const Node startingNode, const Node endingNode, 
                                    int traversal_depth, int& success, double& mean_abundance, bool most_covered = true,
                                    unsigned int backtrackingLimit = 0, Node *avoidFirstNode = NULL);
    Path heuristic_most_covered_path(Direction dir, const Node startingNode, const Node endingNode, 
                                    int traversal_depth, Path current_path, set<Node::Value> usedNode, int& success, vector<int>& abundances, bool most_covered,
                                    unsigned int backtrackingLimit, Node *avoidFirstNode, unsigned long &nbCalls);

    vector<bool> interestingNodes;
};

/********************************************************************************/

#endif 

