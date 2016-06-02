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
#include <gatb/debruijn/impl/Simplifications.hpp>

#include <fstream>
#include <string>
#include <iomanip> // for setprecision
using namespace std;

/********************************************************************************/

#define DEBUG(a)   //a

/********************************************************************************/

static const char* STR_TRAVERSAL_KIND  = "-traversal";
static const char* STR_STARTER_KIND    = "-starter";
static const char* STR_CONTIG_MAX_LEN  = "-contig-max-len";
static const char* STR_BFS_MAX_DEPTH   = "-bfs-max-depth";
static const char* STR_BFS_MAX_BREADTH = "-bfs-max-breadth";
static const char* STR_FASTA_LINE_SIZE = "-fasta-line";
static const char* STR_KEEP_ISOLATED   = "-keep-isolated";

static const char* progressFormat0 = "Minia : assembly";

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
	assemblyParser->push_front (new OptionNoParam  (STR_KEEP_ISOLATED,   "keep short (<= max(2k, 150 bp)) isolated output sequences", false));
	assemblyParser->push_front (new OptionOneParam (STR_URI_INPUT,       "input reads (fasta/fastq/compressed) or hdf5 file",   false));

    getParser()->push_back (assemblyParser);

    // when we input reads, dbgh5 is executed, so its options are needed here
    IOptionsParser* graphParser = Graph::getOptionsParser(false);

    // we hide the STR_URI_INPUT option, otherwise we would have it twice
    if (IOptionsParser* p = graphParser->getParser(STR_URI_INPUT))  {  p->setVisible(false); }

    // we set the default value for the abundance min
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("2"); }

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
struct Parameter
{
    Parameter (Minia& minia) : minia(minia){}
    Minia&         minia;
};

// TODO refactor: can it be done that MiniaFunctor is a member of Minia class?
template<size_t span> 
struct MiniaFunctor  {  void operator ()  (Parameter parameter)
{
    Minia&       minia = parameter.minia;

    typedef GraphTemplate<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>> GraphFast;

    GraphFast graph;

    {
        TIME_INFO (minia.getTimeInfo(), "graph construction");

        //Warning if kmer size >128 cascading debloom does not work
        if(minia.getInput()->getInt(STR_KMER_SIZE)>128){
            minia.getInput()->get(STR_DEBLOOM_TYPE)->value="original";
        }
 
        // graph to not construct branching nodes
        minia.getInput()->setStr(STR_BRANCHING_TYPE,  "none");

        if (minia.getInput()->get(STR_URI_INPUT) != 0)
        {
            graph = GraphFast::create (minia.getInput());
        }
        else
        {
            throw OptionFailure (minia.getParser(), "Specifiy -in");
        }
    }

    // new    
    graph.precomputeAdjacency(minia.getInput()->getInt(STR_NB_CORES));

    /** We build the contigs. */
    minia.assemble<GraphFast, NodeFast<span>, EdgeFast<span>, GraphDataVariantFast<span> >(graph);

    /** We gather some statistics. */
    minia.getInfo()->add (1, minia.getTimeInfo().getProperties("time"));
}
};

void Minia::execute ()
{
    /** we get the kmer size chosen by the end user. */
    size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    /** We launch Minia with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<MiniaFunctor,Parameter> (kmerSize, Parameter (*this));
}


template <typename Graph_type, typename Node, typename Edge, typename GraphDataVariant>
void Minia::assembleFrom(Node startingNode, TraversalTemplate<Node,Edge,GraphDataVariant> *traversal, const Graph_type& graph, IBank *outputBank)
{

    Path_t<Node> consensusRight;
    Path_t<Node> consensusLeft;
    Sequence seq (Data::ASCII);

    /** We compute right and left extensions of the starting node. */
    unsigned int lenRight = traversal->traverse (startingNode,                DIR_OUTCOMING, consensusRight);
    bool isolatedLeft = traversal->deadend;
    Node rev_node = graph.reverse(startingNode);
    unsigned int lenLeft  = traversal->traverse (rev_node, DIR_OUTCOMING, consensusLeft);
    bool isolatedRight = traversal->deadend;

    unsigned int lenTotal = graph.getKmerSize() + lenRight + lenLeft;
    unsigned int isolatedCutoff = std::max(2*(unsigned int)graph.getKmerSize(), (unsigned int)150);

    /** We keep this contig if its not [shorter than isolatedCutoff and isolated (SPAdes-like criterion)], or if -keep-isolated passed */
    if (lenTotal > isolatedCutoff || (lenTotal <= isolatedCutoff && (!(isolatedLeft && isolatedRight))) || keepIsolatedTigs)
    {
        /** We create the contig sequence. */
        buildSequence<Graph_type,Node,Edge,GraphDataVariant> (graph, startingNode, lenTotal, nbContigs, consensusRight, consensusLeft, seq);

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
template <typename Graph_type, typename Node, typename Edge, typename GraphDataVariant>
void Minia::assemble (/*const, removed because Simplifications isn't const anymore*/ Graph_type& graph)
{
    TIME_INFO (getTimeInfo(), "assembly");

    string output = (getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) :
        System::file().getBaseName (getInput()->getStr(STR_URI_INPUT)) 
                    )+ ".contigs.fa";

    /** We setup default values if needed. */
    if (getInput()->getInt (STR_CONTIG_MAX_LEN)  == 0)  { getInput()->setInt (STR_CONTIG_MAX_LEN,  TraversalTemplate<Node,Edge,GraphDataVariant>::defaultMaxLen);     }
    if (getInput()->getInt (STR_BFS_MAX_DEPTH)   == 0)  { getInput()->setInt (STR_BFS_MAX_DEPTH,   TraversalTemplate<Node,Edge,GraphDataVariant>::defaultMaxDepth);   }
    if (getInput()->getInt (STR_BFS_MAX_BREADTH) == 0)  { getInput()->setInt (STR_BFS_MAX_BREADTH, TraversalTemplate<Node,Edge,GraphDataVariant>::defaultMaxBreadth); }

    /** We create the output bank. Note that we could make this a little bit prettier
     *  => possibility to save the contigs in specific output format (other than fasta).  */
    IBank* outputBank = new BankFasta (output);
    LOCAL (outputBank);

    /** We set the fasta line size. */
    BankFasta::setDataLineSize (getInput()->getInt (STR_FASTA_LINE_SIZE));

    TerminatorTemplate<Node,Edge,GraphDataVariant> *terminator;
    bool simplifyGraph = false;

    /* New Minia traversal mode: use MPHF to mark visited nodes. output all unitigs of simplified graph. doesn't care about -starter option */
    terminator = new MPHFTerminatorTemplate<Node,Edge,GraphDataVariant>(graph);
    string traversalKind = "unitig"; // we output unitigs of the simplified graph or the original graph
    simplifyGraph = getInput()->getStr(STR_TRAVERSAL_KIND).compare("contig") == 0;

    LOCAL (terminator);

    /** We create the Traversal instance according to the user choice. */
    TraversalTemplate<Node,Edge,GraphDataVariant>* traversal = TraversalTemplate<Node,Edge,GraphDataVariant>::create (
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

    hasMphf = true;
    keepIsolatedTigs = getParser()->saw(STR_KEEP_ISOLATED);

    string str_tipRemoval = "", str_bubbleRemoval = "", str_ECRemoval = "";


    /** We get an iterator over all nodes . */
    ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem,Node,Edge,GraphDataVariant> itNode (graph.Graph_type::iterator(), progressFormat0);

    // if we want unitigs, then don't simplify the graph; else do it
    if (simplifyGraph)
    {
        int nbCores = getInput()->getInt(STR_NB_CORES);
        bool verbose=true;
        Simplifications<Node,Edge,GraphDataVariant> graphSimplifications(graph, nbCores, verbose);

        graphSimplifications.simplify();

        str_tipRemoval = graphSimplifications.tipRemoval;
        str_bubbleRemoval = graphSimplifications.bubbleRemoval;
        str_ECRemoval = graphSimplifications.ECRemoval;
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

    /** We add the input parameters to the global properties. */
    getInfo()->add (1, getInput());

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "traversal",         "%s", getInput()->getStr(STR_TRAVERSAL_KIND).c_str());
    getInfo()->add (2, "using_mphf",    "%d", hasMphf); // should always be true
    getInfo()->add (2, "nb_contigs",         "%d", nbContigs);
    getInfo()->add (2, "nb_small_contigs_discarded","%d", nbSmallContigs);
    getInfo()->add (2, "nt_assembled",      "%ld", totalNt);
    getInfo()->add (2, "max_length",        "%d", maxContigLen);
    getInfo()->add (2, "max_length_left",   "%d", maxContigLenLeft);
    getInfo()->add (2, "max_length_right",  "%d", maxContigLenRight);

    getInfo()->add (2, "graph simpification stats");
    getInfo()->add (3, "tips removed",          "%s", str_tipRemoval.c_str());
    getInfo()->add (3, "bubbles removed",          "%s", str_bubbleRemoval.c_str());
    getInfo()->add (3, "EC removed",          "%s", str_ECRemoval.c_str());
    getInfo()->add (2, "assembly traversal stats");
    getInfo()->add (3, "no extension",             "%d", traversal->final_stats.couldnt_no_extension);
    getInfo()->add (3, "out-branching",       "%d", traversal->final_stats.couldnt_outbranching);
    getInfo()->add (3, "in-branching",       "%d", traversal->final_stats.couldnt_inbranching);

}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Graph_type, typename Node, typename Edge, typename GraphDataVariant>
void Minia::buildSequence (
    const Graph_type& graph,
    Node& startingNode,
    size_t length,
    size_t nbContigs,
    const Path_t<Node>& consensusRight,
    const Path_t<Node>& consensusLeft,
    Sequence& seq
)
{
    /** Shortcuts. */
    Data&  data     = seq.getData();
    size_t lenRight = consensusRight.size();
    size_t lenLeft  = consensusLeft.size ();

    /** We set the data length. */
    seq.getData().resize (length);

    string sequence = "";
    size_t idx=0;

    /** We dump the left part. */
    for (size_t i=0; i<lenLeft;  i++)  {  
        unsigned char c = ascii (reverse(consensusLeft [lenLeft-i-1]));
        data[idx++] = c; 
        sequence += c;
    }

    /** We dump the starting node. */
    string node = graph.toString (startingNode);
    for (size_t i=0; i<node.size(); i++)  { 
        data[idx++] = node[i]; 
        sequence += node[i];
    }

    /** We dump the right part. */
    for (size_t i=0; i<lenRight; i++)  {  
        unsigned char c = ascii (consensusRight[i]);
        data[idx++] = c;
        sequence += c;
    }

    // get coverage
    double coverage = 0;
    bool computeCoverage = true;
    if (computeCoverage)
    {
        for (unsigned int i = 0; i < length - graph.getKmerSize() + 1; i ++)
        {
            // taken from Traversal.cpp (might be good to factorize into something like getCoverage(string))
            Node node = graph.buildNode((char *)(sequence.c_str()), i); 
            /* I know that buildNode was supposed to be used for test purpose only,
             * but couldn't find anything else to transform my substring into a kmer */

            unsigned char abundance = graph.queryAbundance(node);
            coverage += (unsigned int)abundance;

        }
        coverage /= length - graph.getKmerSize() + 1;
    }

    /** We set the sequence comment. */
    stringstream ss1;
    // spades-like header (compatible with bandage) 
    ss1 << "NODE_"<< nbContigs + 1 << "_length_" << length << "_cov_" << fixed << std::setprecision(3) << coverage << "_ID_" << nbContigs;
    seq._comment = ss1.str();
}
