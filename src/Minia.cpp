/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

#include <Minia.hpp>
#include <NodeSelector.hpp>
#include <GraphSimplification.hpp>

#include <fstream>
#include <string>
using namespace std;

/********************************************************************************/

#define DEBUG(a)   //a

/********************************************************************************/

static const char* STR_TRAVERSAL_KIND  = "-traversal";
static const char* STR_STARTER_KIND    = "-starter";
static const char* STR_CONTIG_MAX_LEN  = "-contig-max-len";
static const char* STR_BFS_MAX_DEPTH   = "-bfs-max-depth";
static const char* STR_BFS_MAX_BREADTH = "-bfs-max-breadth";
static const char* STR_NO_LENGTH_CUTOFF = "-no-length-cutoff";
static const char* STR_FASTA_LINE_SIZE  = "-fasta-line";

static const char* progressFormat0 = "Assembly                               ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Minia::Minia () : Tool ("minia")
{
	// reinit the parser to get rid of options added by the Tool class, as we'll add them again in the Graph parser
	setParser (new OptionsParser ("minia")); 

    /** We add options specific to Minia (most important at the end). */
	OptionsParser* assemblyParser = new OptionsParser ("assembly");

	assemblyParser->push_front (new OptionOneParam (STR_FASTA_LINE_SIZE, "number of nucleotides per line in fasta output (0 means one line)",  false, "0"));
	assemblyParser->push_front (new OptionOneParam (STR_BFS_MAX_BREADTH, "maximum breadth for BFS",               false,  "0"         ));
	assemblyParser->push_front (new OptionOneParam (STR_BFS_MAX_DEPTH,   "maximum depth for BFS",                 false,  "0"         ));
	assemblyParser->push_front (new OptionOneParam (STR_CONTIG_MAX_LEN,  "maximum length for contigs",            false,  "0"         ));
	assemblyParser->push_front (new OptionOneParam (STR_STARTER_KIND,    "starting node ('best', 'simple')",      false,  "best"      ));
	assemblyParser->push_front (new OptionOneParam (STR_TRAVERSAL_KIND,  "traversal type ('contig', 'unitig')", false,  "contig"  ));
	assemblyParser->push_front (new OptionNoParam  (STR_NO_LENGTH_CUTOFF, "turn off length cutoff of 2*k in output sequences", false));
	assemblyParser->push_front (new OptionOneParam (STR_URI_INPUT,       "input reads (fasta/fastq/compressed)",   false));
	assemblyParser->push_front (new OptionOneParam (STR_URI_GRAPH,       "input graph file (hdf5)",                false));

    getParser()->push_back (assemblyParser);

    // when we input reads, dbgh5 is executed, so its options are needed here
    IOptionsParser* graphParser = Graph::getOptionsParser(false, true);

    // we hide the STR_URI_INPUT option, otherwise we would have it twice
    if (IOptionsParser* p = graphParser->getParser(STR_URI_INPUT))  {  p->setVisible(false); }

    getParser()->push_back(graphParser, 1);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Minia::execute ()
{
 	Graph graph;

    // graph to not construct branching nodes in the event of mphf != none
    if (getInput()->getStr(STR_MPHF_TYPE).compare("emphf") == 0)  { getInput()->setStr(STR_BRANCHING_TYPE,  "none");     }

	if (getInput()->get(STR_URI_GRAPH) != 0)
	{
		graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
	}
	else if (getInput()->get(STR_URI_INPUT) != 0)
    {
        graph = Graph::create (getInput());
    }
    else
    {
        throw OptionFailure (getParser(), "Specifiy -graph or -in");
	}

    /** We build the contigs. */
    assemble (graph);

    /** We gather some statistics. */
    getInfo()->add (1, getTimeInfo().getProperties("time"));
}

void Minia::assembleFrom(Node startingNode, Traversal *traversal, const Graph& graph, IBank *outputBank)
{

    Path consensusRight;
    Path consensusLeft;
    Sequence seq (Data::ASCII);

    /** We compute right and left extensions of the starting node. */
    unsigned int lenRight = traversal->traverse (startingNode,                DIR_OUTCOMING, consensusRight);
    unsigned int lenLeft  = traversal->traverse (graph.reverse(startingNode), DIR_OUTCOMING, consensusLeft);

    unsigned int lenTotal = graph.getKmerSize() + lenRight + lenLeft;

    /** We keep this contig if its size is long enough. */
    if ((unsigned long)lenTotal >= 2*graph.getKmerSize()+1 || isNoLengthCutoff)
    {
        /** We create the contig sequence. */
        buildSequence (graph, startingNode, lenTotal, nbContigs, consensusRight, consensusLeft, seq);

        /** We add the sequence into the output bank. */
        outputBank->insert (seq);

        nbContigs += 1;
        totalNt   += lenTotal;

        traversal->commit_stats();

        if (lenTotal > maxContigLen)      { maxContigLen      = lenTotal;   }
        if (lenLeft  > maxContigLenLeft)  { maxContigLenLeft  = lenLeft;    }
        if (lenRight > maxContigLenRight) { maxContigLenRight = lenRight;   }
    }
    else
    {
        traversal->revert_stats();
        nbSmallContigs++;
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
void Minia::assemble (const Graph& graph)
{
    TIME_INFO (getTimeInfo(), "assembly");

    string output = (getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) :
        System::file().getBaseName (
		  (getInput()->get(STR_URI_INPUT) ? getInput()->getStr(STR_URI_INPUT) : 
		   getInput()->getStr(STR_URI_GRAPH))
                                   ) 
                    )+ ".contigs.fa";

    /** We setup default values if needed. */
    if (getInput()->getInt (STR_CONTIG_MAX_LEN)  == 0)  { getInput()->setInt (STR_CONTIG_MAX_LEN,  Traversal::defaultMaxLen);     }
    if (getInput()->getInt (STR_BFS_MAX_DEPTH)   == 0)  { getInput()->setInt (STR_BFS_MAX_DEPTH,   Traversal::defaultMaxDepth);   }
    if (getInput()->getInt (STR_BFS_MAX_BREADTH) == 0)  { getInput()->setInt (STR_BFS_MAX_BREADTH, Traversal::defaultMaxBreadth); }

    /** We create the output bank. Note that we could make this a little bit prettier
     *  => possibility to save the contigs in specific output format (other than fasta).  */
    IBank* outputBank = new BankFasta (output);
    LOCAL (outputBank);

    /** We set the fasta line size. */
    BankFasta::setDataLineSize (getInput()->getInt (STR_FASTA_LINE_SIZE));

    bool hasMphf = (graph.getState() & Graph::STATE_MPHF_DONE);
    bool legacyTraversal = !hasMphf;

    Terminator *terminator;
    INodeSelector* starter = NULL;
    string traversalKind;
    bool simplifyGraph = false;

    // Legacy Minia traversal mode: index branching nodes to mark kmers, use branching nodes to start traversals
    if (legacyTraversal)
    { 
        /** We create the Terminator. */
        terminator = new BranchingTerminator(graph);

        /** We create the starting node selector according to the user choice. */
        starter = NodeSelectorFactory::singleton().create (getInput()->getStr(STR_STARTER_KIND), graph, *terminator);
        
        traversalKind = getInput()->getStr(STR_TRAVERSAL_KIND);
    }
    else
    {
        terminator = new MPHFTerminator(graph);
        traversalKind = "unitig"; // we output unitigs of the simplified graph or the original graph
        simplifyGraph = getInput()->getStr(STR_TRAVERSAL_KIND).compare("contig") == 0;
    }
    LOCAL (terminator);
    LOCAL (starter);

    /** We create the Traversal instance according to the user choice. */
    Traversal* traversal = Traversal::create (
        traversalKind,
        graph,
        *terminator,
        getInput()->getInt (STR_CONTIG_MAX_LEN),
        getInput()->getInt (STR_BFS_MAX_DEPTH),
        getInput()->getInt (STR_BFS_MAX_BREADTH)
    );
    LOCAL (traversal);

    nbContigs         = 0;
    nbSmallContigs    = 0;
    totalNt           = 0;
    maxContigLen      = 0;
    maxContigLenLeft  = 0;
    maxContigLenRight = 0;

    isNoLengthCutoff = getParser()->saw(STR_NO_LENGTH_CUTOFF);

    string tipRemoval = "", bubbleRemoval = "";

    if (legacyTraversal)
    {
        /** We get an iterator over the branching nodes. */
        ProgressGraphIterator<BranchingNode,ProgressTimerAndSystem> itBranching (graph.iterator<BranchingNode>(), progressFormat0);

        /** We loop over the branching nodes. */
        for (itBranching.first(); !itBranching.isDone(); itBranching.next())
        {
            DEBUG ((cout << endl << "-------------------------- " << graph.toString (itBranching.item()) << " -------------------------" << endl));

            Node startingNode;

            // keep looping while a starting kmer is available from this kmer
            // everything will be marked during the traversal()'s
            while (starter->select (itBranching.item(), startingNode) == true)
            {
                assembleFrom(startingNode, traversal, graph, outputBank);        

            } /* end of  while (starter->select() */

        } /* end of for (itBranching.first() */
    }
    else
    {
        /** We get an iterator over all nodes . */
        ProgressGraphIterator<Node,ProgressTimerAndSystem> itNode (graph.iterator<Node>(), progressFormat0);

        // if we want unitigs, then don't simplify the graph; else do it
        if (simplifyGraph)
        {
            GraphSimplification graphSimplification(graph);

            unsigned long nbTipsRemoved_1 = graphSimplification.removeTips();
            unsigned long nbTipsRemoved_2 = graphSimplification.removeTips();

            unsigned long nbBubblesRemoved_1 = graphSimplification.removeBubbles();
            unsigned long nbBubblesRemoved_2 = graphSimplification.removeBubbles();
           
            unsigned long nbTipsRemoved_3 = graphSimplification.removeTips();
            
            unsigned long nbBubblesRemoved_3 = graphSimplification.removeBubbles();
            
            tipRemoval = std::to_string(nbTipsRemoved_1) + " + " + std::to_string(nbTipsRemoved_2) + " + " + std::to_string(nbTipsRemoved_3) ;
            bubbleRemoval = std::to_string(nbBubblesRemoved_1) + " + " + std::to_string(nbBubblesRemoved_2) +  " + " + std::to_string(nbBubblesRemoved_3);
        }

        /** We loop over all nodes. */
        for (itNode.first(); !itNode.isDone(); itNode.next())
        {
            Node node = itNode.item();

            // in this setting it's very simple, we don't even need NodeSelector anymore. Just assemble from any non-deleted unmarked node
            if (terminator->is_marked (node))  {  continue;   }
            if (graph.isNodeDeleted(node)) { continue; }
            
            DEBUG ((cout << endl << "-------------------------- " << graph.toString (node) << " -------------------------" << endl));

            assembleFrom(node, traversal, graph, outputBank);
        }

    }

    /** We add the input parameters to the global properties. */
    getInfo()->add (1, getInput());

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "traversal",         "%s", getInput()->getStr(STR_TRAVERSAL_KIND).c_str());
    getInfo()->add (2, "using_mphf",    "%d", hasMphf);
    if (legacyTraversal)
        getInfo()->add (2, "start_selector",    "%s", starter->getName().c_str());
    getInfo()->add (2, "nb_contigs",         "%d", nbContigs);
    getInfo()->add (2, "nb_small_contigs_discarded","%d", nbSmallContigs);
    getInfo()->add (2, "nt_assembled",      "%ld", totalNt);
    getInfo()->add (2, "max_length",        "%d", maxContigLen);
    getInfo()->add (2, "max_length_left",   "%d", maxContigLenLeft);
    getInfo()->add (2, "max_length_right",  "%d", maxContigLenRight);

    if (legacyTraversal)
    {
        getInfo()->add (2, "debugging traversal stats");
        getInfo()->add (2, "large breadth",          "%d", traversal->final_stats.couldnt_traverse_bubble_breadth);
        getInfo()->add (2, "large depth",            "%d", traversal->final_stats.couldnt_traverse_bubble_depth);
        getInfo()->add (2, "marked kmer inside traversal",        "%d", traversal->final_stats.couldnt_because_marked_kmer);
        getInfo()->add (2, "traversal ends with dead-ends",             "%d", traversal->final_stats.couldnt_find_extension);
        getInfo()->add (2, "in-branching large depth",       "%d", traversal->final_stats.couldnt_inbranching_depth);
        getInfo()->add (2, "in-branching large breadth",    "%d", traversal->final_stats.couldnt_inbranching_breadth);
        getInfo()->add (2, "in-branching other",            "%d", traversal->final_stats.couldnt_inbranching_other);
        getInfo()->add (2, "couldn't validate consensuses", "%d", traversal->final_stats.couldnt_validate_consensuses);
    }
    else
    {
        getInfo()->add (2, "graph simpification stats");
        getInfo()->add (2, "tips removed",          "%s", tipRemoval.c_str());
        getInfo()->add (2, "bubbles removed",          "%s", bubbleRemoval.c_str());
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
void Minia::buildSequence (
    const Graph& graph,
    const Node& startingNode,
    size_t length,
    size_t nbContigs,
    const Path& consensusRight,
    const Path& consensusLeft,
    Sequence& seq
)
{
    /** Shortcuts. */
    Data&  data     = seq.getData();
    size_t lenRight = consensusRight.size();
    size_t lenLeft  = consensusLeft.size ();

    /** We set the sequence comment. */
    stringstream ss1;
    ss1 << nbContigs << "__len__" << length;
    seq._comment = ss1.str();

    /** We set the data length. */
    seq.getData().resize (length);

    size_t idx=0;

    /** We dump the left part. */
    for (size_t i=0; i<lenLeft;  i++)  {  data[idx++] = ascii (reverse(consensusLeft [lenLeft-i-1])); }

    /** We dump the starting node. */
    string node = graph.toString (startingNode);
    for (size_t i=0; i<node.size(); i++)  { data[idx++] = node[i]; }

    /** We dump the right part. */
    for (size_t i=0; i<lenRight; i++)  {  data[idx++] = ascii (consensusRight[i]); }
}
