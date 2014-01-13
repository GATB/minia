/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <Minia.hpp>
#include <Frontline.hpp>
#include <Terminator.hpp>
#include <Traversal.hpp>
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
    /** We add options specific to Minia (most important at the end). */
    getParser()->push_front (new OptionOneParam (STR_VERBOSE,         "verbosity level",                       false,  "1"));
    getParser()->push_front (new OptionOneParam (STR_BFS_MAX_BREADTH, "maximum breadth for BFS",               false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_BFS_MAX_DEPTH,   "maximum depth for BFS",                 false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_CONTIG_MAX_LEN,  "maximum length for contigs",            false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_STARTER_KIND,    "starting node ('best', 'simple')",      false,  "best"      ));
    getParser()->push_front (new OptionOneParam (STR_TRAVERSAL_KIND,  "traversal type ('monument', 'unitig')", false,  "monument"  ));
    getParser()->push_front (new OptionOneParam (STR_MAX_DISK,        "max disk space in MBytes",              false,  "0"         ));
    getParser()->push_front (new OptionOneParam (STR_MAX_MEMORY,      "max memory in MBytes",                  false,  "1000"      ));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT,      "output file",                           false));
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT,       "input file (likely a hdf5 file)",       true));
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
    /** We load the graph from the provided uri. */
    Graph graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

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

    string output = getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) :
        System::file().getBaseName (getInput()->getStr(STR_URI_INPUT)) + ".contigs";

    /** We setup default values if needed. */
    if (getInput()->getInt (STR_CONTIG_MAX_LEN)  == 0)  { getInput()->setInt (STR_CONTIG_MAX_LEN,  Traversal::defaultMaxLen);     }
    if (getInput()->getInt (STR_BFS_MAX_DEPTH)   == 0)  { getInput()->setInt (STR_BFS_MAX_DEPTH,   Traversal::defaultMaxDepth);   }
    if (getInput()->getInt (STR_BFS_MAX_BREADTH) == 0)  { getInput()->setInt (STR_BFS_MAX_BREADTH, Traversal::defaultMaxBreadth); }

    /** We create the output bank. Note that we could make this a little bit prettier
     *  => possibility to save the contigs in specific output format (other than fasta).  */
    IBank* outputBank = new BankFasta (output);
    LOCAL (outputBank);

    /** We get an iterator over the branching nodes. */
    ProgressIterator<BranchingNode,ProgressTimer> itBranching (graph.iterator<BranchingNode>(), "assembly");

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
    u_int64_t totalNt           = 0;
    u_int64_t maxContigLen      = 0;
    u_int64_t maxContigLenLeft  = 0;
    u_int64_t maxContigLenRight = 0;

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
            if (lenTotal >= 2*graph.getKmerSize()+1)
            {
                /** We create the contig sequence. */
                buildSequence (graph, startingNode, lenTotal, nbContigs, consensusRight, consensusLeft, seq);

                /** We add the sequence into the output bank. */
                outputBank->insert (seq);

                nbContigs += 1;
                totalNt   += lenTotal;

                if (lenTotal > maxContigLen)      { maxContigLen      = lenTotal;   }
                if (lenLeft  > maxContigLenLeft)  { maxContigLenLeft  = lenLeft;    }
                if (lenRight > maxContigLenRight) { maxContigLenRight = lenRight;   }
            }

        } /* end of  while (starter->select() */

    } /* end of for (itBranching.first() */

    /** We add the input parameters to the global properties. */
    getInfo()->add (1, getInput());

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "traversal",         "%s", traversal->getName().c_str());
    getInfo()->add (2, "start_selector",    "%s", starter->getName().c_str());
    getInfo()->add (2, "nb_contig",         "%d", nbContigs);
    getInfo()->add (2, "nt_assembled",      "%d", totalNt);
    getInfo()->add (2, "max_length",        "%d", maxContigLen);
    getInfo()->add (2, "max_length_left",   "%d", maxContigLenLeft);
    getInfo()->add (2, "max_length_right",  "%d", maxContigLenRight);
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
