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

#include <fstream>

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
static const char* STR_URI_GRAPH        = "-graph";
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

   // when we input reads, dbgh5 is executed, so its options are needed here
   getParser()->push_back(Graph::getOptionsParser(false));

    /** We add options specific to Minia (most important at the end). */
    getParser()->push_front (new OptionOneParam (STR_FASTA_LINE_SIZE, "number of nucleotides per fasta line (0 means one line)",  false, "0"));
    getParser()->push_front (new OptionOneParam (STR_BFS_MAX_BREADTH, "maximum breadth for BFS",               false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_BFS_MAX_DEPTH,   "maximum depth for BFS",                 false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_CONTIG_MAX_LEN,  "maximum length for contigs",            false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_STARTER_KIND,    "starting node ('best', 'simple')",      false,  "best"      ));
    getParser()->push_front (new OptionOneParam (STR_TRAVERSAL_KIND,  "traversal type ('monument', 'unitig')", false,  "monument"  ));
    getParser()->push_front (new OptionNoParam  (STR_NO_LENGTH_CUTOFF, "turn off length cutoff of 2*k in output sequences", false));
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT,       "input reads (fasta/fastq/compressed)",   false));
    getParser()->push_front (new OptionOneParam (STR_URI_GRAPH,       "input graph file (hdf5)",                false));
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
 if (getInput()->get(STR_VERSION) != 0){
        cout << "Minia from GATB version "<< STR_LIBRARY_VERSION << endl;
        return;
    }

	Graph graph;

	if (getInput()->get(STR_URI_GRAPH) != 0)
	{
		graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
	}
	else
	{

		if (getInput()->get(STR_URI_INPUT) != 0)
		{
			graph = Graph::create (getInput());

		}
		else
		{
			cout << "Specifiy -graph or -in \n";
			throw OptionFailure(*getParser());
		}
	}

    /** We build the contigs. */
    assemble (graph);

    /** We gather some statistics. */
    getInfo()->add (1, getTimeInfo().getProperties("time"));
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

    /** We get an iterator over the branching nodes. */
    ProgressGraphIterator<BranchingNode,ProgressTimerAndSystem> itBranching (graph.iterator<BranchingNode>(), progressFormat0);

    /** We create the Terminator. */
    BranchingTerminator terminator (graph);

    /** We create the starting node selector according to the user choice. */
    INodeSelector* starter = NodeSelectorFactory::singleton().create (getInput()->getStr(STR_STARTER_KIND), graph, terminator);
    LOCAL (starter);

    /** We create the Traversal instance according to the user choice. */
    Traversal* traversal = Traversal::create (
        getInput()->getStr(STR_TRAVERSAL_KIND),
        graph,
        terminator,
        getInput()->getInt (STR_CONTIG_MAX_LEN),
        getInput()->getInt (STR_BFS_MAX_DEPTH),
        getInput()->getInt (STR_BFS_MAX_BREADTH)
    );
    LOCAL (traversal);

    Path consensusRight;
    Path consensusLeft;

    u_int64_t nbContigs         = 0;
    u_int64_t nbSmallContigs    = 0;
    u_int64_t totalNt           = 0;
    u_int64_t maxContigLen      = 0;
    u_int64_t maxContigLenLeft  = 0;
    u_int64_t maxContigLenRight = 0;

    bool isNoLengthCutoff = getParser()->saw(STR_NO_LENGTH_CUTOFF);

    Sequence seq (Data::ASCII);

    /** We loop over the branching nodes. */
    for (itBranching.first(); !itBranching.isDone(); itBranching.next())
    {
        DEBUG ((cout << endl << "-------------------------- " << graph.toString (itBranching.item()) << " -------------------------" << endl));

        Node startingNode;

        // keep looping while a starting kmer is available from this kmer
        // everything will be marked during the traversal()'s
        while (starter->select (itBranching.item(), startingNode) == true)
        {
            /** We compute right and left extensions of the starting node. */
            int lenRight = traversal->traverse (startingNode,                DIR_OUTCOMING, consensusRight);
            int lenLeft  = traversal->traverse (graph.reverse(startingNode), DIR_OUTCOMING, consensusLeft);

            int lenTotal = graph.getKmerSize() + lenRight + lenLeft;

            /** We keep this contig if its size is long enough. */
            if (lenTotal >= 2*graph.getKmerSize()+1 || isNoLengthCutoff)
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

        } /* end of  while (starter->select() */

    } /* end of for (itBranching.first() */

    /** We add the input parameters to the global properties. */
    getInfo()->add (1, getInput());

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "traversal",         "%s", traversal->getName().c_str());
    getInfo()->add (2, "start_selector",    "%s", starter->getName().c_str());
    getInfo()->add (2, "nb_contigs",         "%d", nbContigs);
    getInfo()->add (2, "nb_small_contigs_discarded","%d", nbSmallContigs);
    getInfo()->add (2, "nt_assembled",      "%ld", totalNt);
    getInfo()->add (2, "max_length",        "%d", maxContigLen);
    getInfo()->add (2, "max_length_left",   "%d", maxContigLenLeft);
    getInfo()->add (2, "max_length_right",  "%d", maxContigLenRight);

    getInfo()->add (1, "debugging traversal stats");
    getInfo()->add (2, "couldn't validate consensuses", "%d", traversal->final_stats.couldnt_validate_consensuses);
    getInfo()->add (2, "large bubble breadth",          "%d", traversal->final_stats.couldnt_traverse_bubble_breadth);
    getInfo()->add (2, "large bubble depth",            "%d", traversal->final_stats.couldnt_traverse_bubble_depth);
    getInfo()->add (2, "stopped at marked kmer",        "%d", traversal->final_stats.couldnt_because_marked_kmer);
    getInfo()->add (2, "no kmer extension",             "%d", traversal->final_stats.couldnt_find_extension);
    getInfo()->add (2, "in-branchin large depth",       "%d", traversal->final_stats.couldnt_inbranching_depth);
    getInfo()->add (2, "in-branching large breadth",    "%d", traversal->final_stats.couldnt_inbranching_breadth);
    getInfo()->add (2, "in-branching other",            "%d", traversal->final_stats.couldnt_inbranching_other);
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
    ss1 << nbContigs << "__len__" << length << " ";
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
