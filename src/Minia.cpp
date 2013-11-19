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
    getParser()->add (new OptionOneParam (STR_URI_DB,          "databank uri",                         true));
    getParser()->add (new OptionOneParam (STR_URI_OUTPUT,      "output file",                          false));
    getParser()->add (new OptionOneParam (STR_MAX_MEMORY,      "max memory in MBytes",                 false,  "1000"  ));
    getParser()->add (new OptionOneParam (STR_MAX_DISK,        "max disk space in MBytes",             false,  "0"     ));
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

    ofstream file (output);

    /** We get an iterator over the branching nodes. */
    ProgressIterator<BranchingNode,ProgressTimer> itBranching (graph.iterator<BranchingNode>(), "assembly");

    /** We create the Terminator. */
    BranchingTerminator terminator (graph);

    SimplePathsTraversal traversal (graph, terminator);

    std::vector<Nucleotide> consensusRight;
    std::vector<Nucleotide> consensusLeft;

    u_int64_t nbContigs         = 0;
    u_int64_t totalNt           = 0;
    u_int64_t maxContigLen      = 0;
    u_int64_t maxContigLenLeft  = 0;
    u_int64_t maxContigLenRight = 0;

    /** We loop over the branching nodes.
     * IMPORTANT... by now, use a dispatcher with only 1 thread since the terminator is not thread safe. */
    Dispatcher(1).iterate (itBranching, [&] (const Node& node)
    {
        Node startingNode;

        DEBUG ((cout << endl << "-------------------------- " << graph.toString (node) << " -------------------------" << endl));

        while (traversal.findStartingNode (node, startingNode) == true)
        {
            int lenRight = traversal.traverse (startingNode, DIR_OUTCOMING, consensusRight);
            int lenLeft  = traversal.traverse (startingNode, DIR_INCOMING,  consensusLeft);

            int len = graph.getKmerSize() + lenRight + lenLeft;

            if (len >= 2*graph.getKmerSize()+1)
            {
                if (startingNode.strand==STRAND_FORWARD)
                {
                    dump (graph, startingNode, len, nbContigs, consensusRight, consensusLeft, file);
                }
                else
                {
                    dump (graph, startingNode, len, nbContigs, consensusLeft, consensusRight, file);
                }

                nbContigs += 1;
                totalNt   +=  len;

                if (len > maxContigLen)             { maxContigLen      = len;      }
                if (lenLeft    > maxContigLenLeft)  { maxContigLenLeft  = lenLeft;  }
                if (lenRight   > maxContigLenRight) { maxContigLenRight = lenRight; }
            }
        }
    });

    file.close ();

    //terminator.dump ();

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "uri",               "%s", getInput()->getStr(STR_URI_DB).c_str());
    getInfo()->add (2, "traversal",         "%s", traversal.getName().c_str());
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
void Minia::dump (
    const Graph& graph,
    const Node& startingNode,
    size_t length,
    size_t nbContigs,
    const std::vector<Nucleotide>& consensusRight,
    const std::vector<Nucleotide>& consensusLeft,
    std::ostream& file
)
{
    /** Shortcuts. */
    size_t lenRight = consensusRight.size();
    size_t lenLeft  = consensusLeft.size ();

    /** We dump the sequence comment. */
    file << ">" << nbContigs << "__len__" << length << " " << endl;

    /** We dump the left part. */
    for (size_t i=0; i<lenLeft;  i++)  {  file << ascii (reverse(consensusLeft [lenLeft-i-1])); }

    /** We dump the starting node. */
    file << graph.toString (startingNode, startingNode.strand, 1);

    /** We dump the right part. */
    for (size_t i=0; i<lenRight; i++)  {  file << ascii (consensusRight[i]);           }

    file << endl;
}
