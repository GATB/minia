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

static const char* STR_TRAVERSAL_KIND = "-traversal";
static const char* STR_STARTER_KIND   = "-starter";

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
    /** We add options specific to DSK. */
    getParser()->add (new OptionOneParam (STR_URI_DB,          "databank uri",                          true));
    getParser()->add (new OptionOneParam (STR_URI_OUTPUT,      "output file",                           false));
    getParser()->add (new OptionOneParam (STR_MAX_MEMORY,      "max memory in MBytes",                  false,  "1000"      ));
    getParser()->add (new OptionOneParam (STR_MAX_DISK,        "max disk space in MBytes",              false,  "0"         ));
    getParser()->add (new OptionOneParam (STR_TRAVERSAL_KIND,  "traversal type ('monument', 'unitig')", false,  "monument"  ));
    getParser()->add (new OptionOneParam (STR_STARTER_KIND,    "starting node ('simple')",              false,  "simple"    ));
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
    Graph graph = Graph::load (getInput()->getStr(STR_URI_DB));

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
        System::file().getBaseName (getInput()->getStr(STR_URI_DB)) + ".contigs";

    /** We create the output bank. Note that we could make this a little bit prettier
     *  => possibility to save the contigs in specific output format (other than fasta).  */
    IBank* outputBank = new Bank (output);
    LOCAL (outputBank);

    /** We get an iterator over the branching nodes. */
    ProgressIterator<BranchingNode,ProgressTimer> itBranching (graph.iterator<BranchingNode>(), "assembly");

    /** We create the Terminator. */
    BranchingTerminator terminator (graph);

    /** We create the Traversal instance according to the user choice. */
    Traversal* traversal = Traversal::create (getInput()->getStr(STR_TRAVERSAL_KIND), graph, terminator);
    LOCAL (traversal);

    std::vector<Nucleotide> consensusRight;
    std::vector<Nucleotide> consensusLeft;

    u_int64_t nbContigs         = 0;
    u_int64_t totalNt           = 0;
    u_int64_t maxContigLen      = 0;
    u_int64_t maxContigLenLeft  = 0;
    u_int64_t maxContigLenRight = 0;

    Sequence seq (Data::ASCII);

    /** We loop over the branching nodes.
     * IMPORTANT... by now, use a dispatcher with only 1 thread since the terminator is not thread safe. */
    Dispatcher(1).iterate (itBranching, [&] (const Node& node)
    {
        Node startingNode;

        DEBUG ((cout << endl << "-------------------------- " << graph.toString (node) << " -------------------------" << endl));

        // keep looping while a starting kmer is available from this kmer
        // everything will be marked during the traversal()'s
        while (traversal->findStartingNode (node, startingNode) == true)
        {
            /** We compute right and left extensions of the starting node. */
            int lenRight = traversal->traverse (startingNode, DIR_OUTCOMING, consensusRight);
            int lenLeft  = traversal->traverse (startingNode, DIR_INCOMING,  consensusLeft);

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
        }
    });

    //terminator.dump ();

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "uri",               "%s", getInput()->getStr(STR_URI_DB).c_str());
    getInfo()->add (2, "traversal",         "%s", traversal->getName().c_str());
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
    const std::vector<Nucleotide>& consensusRight,
    const std::vector<Nucleotide>& consensusLeft,
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
    string node = graph.toString (startingNode, startingNode.strand, 1);
    for (size_t i=0; i<node.size(); i++)  { data[idx++] = node[i]; }

    /** We dump the right part. */
    for (size_t i=0; i<lenRight; i++)  {  data[idx++] = ascii (consensusRight[i]); }
}
