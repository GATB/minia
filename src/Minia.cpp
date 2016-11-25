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
#include <gatb/debruijn/impl/GraphUnitigs.hpp>

#include <fstream>
#include <string>
#include <iomanip> // for setprecision
using namespace std;

/********************************************************************************/

#define DEBUG(a)   //a

/********************************************************************************/

static const char* STR_TRAVERSAL_KIND  = "-traversal";
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

#ifdef GIT_SHA1
    std::cout << "Minia 3, git commit " << GIT_SHA1 << std::endl;
#endif

	// reinit the parser to get rid of options added by the Tool class, as we'll add them again in the Graph parser
	setParser (new OptionsParser ("minia")); 

    /** We add options specific to Minia (most important at the end). */
	OptionsParser* assemblyParser = new OptionsParser ("assembly");

	assemblyParser->push_front (new OptionOneParam (STR_FASTA_LINE_SIZE, "number of nucleotides per line in fasta output (0 means one line)",  false, "0"));
	assemblyParser->push_front (new OptionOneParam (STR_TRAVERSAL_KIND,  "traversal type ('contig', 'unitig')", false,  "contig"  ));
	assemblyParser->push_front (new OptionNoParam  (STR_KEEP_ISOLATED,   "keep short (<= max(2k, 150 bp)) isolated output sequences", false));
	assemblyParser->push_front (new OptionOneParam (STR_URI_INPUT,       "input reads (fasta/fastq/compressed) or hdf5 file",   false));

    getParser()->push_back (assemblyParser);

    // when we input reads, dbgh5 is executed, so its options are needed here
    IOptionsParser* graphParser = Graph::getOptionsParser(false);

    // we hide the STR_URI_INPUT option, otherwise we would have it twice
    if (IOptionsParser* p = graphParser->getParser(STR_URI_INPUT))  {  p->setVisible(false); }

    // we set the default value for the abundance min (2)
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("2"); }
    
    // through an undocumented environment variable, we can set the temporary path (useful to not have to specify -out-tmp every time)
    char *hidden_env_variable_out_tmp = std::getenv("MINIA_OUT_TMP");
    if (hidden_env_variable_out_tmp)
    {
        if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_URI_OUTPUT_TMP)))  {  p->setDefaultValue ((string)hidden_env_variable_out_tmp); }
    }

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

template<size_t span> 
struct MiniaFunctor  {  void operator ()  (Parameter parameter)
{
    Minia&       minia = parameter.minia;

    // selection of type of graph is done HERE
    //typedef GraphTemplate<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>> GraphType;
    typedef GraphUnitigsTemplate<span> GraphType;
    
    GraphType graph;
    
    {
        TIME_INFO (minia.getTimeInfo(), "graph construction");

        if (minia.getInput()->get(STR_URI_INPUT) != 0)
        {
            graph = GraphType::create (minia.getInput());
        }
        else
        {
            throw OptionFailure (minia.getParser(), "Specifiy -in");
        }
    }

    /** We build the contigs. */
    minia.assemble<GraphType, NodeGU, EdgeGU, span>(graph);

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


template <typename Graph_type, typename Node, typename Edge, size_t span>
void Minia::assembleFrom(Node startingNode, Graph_type& graph, IBank *outputBank)
{
    unsigned int isolatedCutoff = std::max(2*(unsigned int)graph.getKmerSize(), (unsigned int)150);

    bool isolatedLeft, isolatedRight;
    float coverage = 0;
    string sequence = graph.simplePathBothDirections(startingNode, isolatedLeft, isolatedRight, true, coverage);

    Sequence seq (Data::ASCII);
    seq.getData().setRef ((char*)sequence.c_str(), sequence.size());
    /** We set the sequence comment. */
    stringstream ss1;
    // spades-like header (compatible with bandage) 
    ss1 << "NODE_"<< nbContigs + 1 << "_length_" << sequence.size() << "_cov_" << fixed << std::setprecision(3) << coverage << "_ID_" << nbContigs;
    seq._comment = ss1.str();
    unsigned int lenTotal = sequence.size();
    if (lenTotal > isolatedCutoff || (lenTotal <= isolatedCutoff && (!(isolatedLeft && isolatedRight))) || keepIsolatedTigs)
    {
        outputBank->insert (seq);
        nbContigs += 1;
        totalNt   += lenTotal;
        if (lenTotal > maxContigLen)      { maxContigLen      = lenTotal;   }
    }
    else
        nbSmallContigs++;
    return;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Graph_type, typename Node, typename Edge, size_t span>
void Minia::assemble (/*const, removed because Simplifications isn't const anymore*/ Graph_type& graph)
{
    TIME_INFO (getTimeInfo(), "assembly");

    string output = (getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) :
        System::file().getBaseName (getInput()->getStr(STR_URI_INPUT)) 
                    )+ ".contigs.fa";

    /** We create the output bank. Note that we could make this a little bit prettier
     *  => possibility to save the contigs in specific output format (other than fasta).  */
    IBank* outputBank = new BankFasta (output);
    LOCAL (outputBank);

    /** We set the fasta line size. */
    BankFasta::setDataLineSize (getInput()->getInt (STR_FASTA_LINE_SIZE));

    bool simplifyGraph = false;

    string traversalKind = "unitig"; // we output unitigs of the simplified graph or the original graph
    simplifyGraph = getInput()->getStr(STR_TRAVERSAL_KIND).compare("contig") == 0;

    nbContigs         = 0;
    nbSmallContigs    = 0;
    totalNt           = 0;
    maxContigLen      = 0;

    keepIsolatedTigs = getParser()->saw(STR_KEEP_ISOLATED);

    string str_tipRemoval = "", str_bubbleRemoval = "", str_ECRemoval = "";

    /** We get an iterator over all nodes . */
    ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem> itNode (graph.Graph_type::iterator(), progressFormat0);

    // if we want unitigs, then don't simplify the graph; else do it
    if (simplifyGraph)
    {
        int nbCores = getInput()->getInt(STR_NB_CORES);
        bool verbose=true;
        Simplifications<Graph_type,Node,Edge> graphSimplifications(graph, nbCores, verbose);

        graphSimplifications.simplify();

        str_tipRemoval = graphSimplifications.tipRemoval;
        str_bubbleRemoval = graphSimplifications.bubbleRemoval;
        str_ECRemoval = graphSimplifications.ECRemoval;
    }

    //graph.debugPrintAllUnitigs(); // debugging

    /** We loop over all nodes. */
    for (itNode.first(); !itNode.isDone(); itNode.next())
    {
        Node node = itNode.item();

        if (graph.unitigIsMarked(node))  {  continue;   }
        if (graph.isNodeDeleted(node)) { continue; }

        DEBUG ((cout << endl << "-------------------------- " << graph.toString (node) << " -------------------------" << endl));

        assembleFrom<Graph_type, Node, Edge, span>(node, graph, outputBank);
    }

    /** We add the input parameters to the global properties. */
    getInfo()->add (1, getInput());

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "traversal",         "%s", getInput()->getStr(STR_TRAVERSAL_KIND).c_str());
    getInfo()->add (2, "nb_contigs",         "%d", nbContigs);
    getInfo()->add (2, "nb_small_contigs_discarded","%d", nbSmallContigs);
    getInfo()->add (2, "nt_assembled",      "%ld", totalNt);
    getInfo()->add (2, "max_length",        "%d", maxContigLen);

    getInfo()->add (2, "graph simpification stats");
    getInfo()->add (3, "tips removed",          "%s", str_tipRemoval.c_str());
    getInfo()->add (3, "bubbles removed",          "%s", str_bubbleRemoval.c_str());
    getInfo()->add (3, "EC removed",          "%s", str_ECRemoval.c_str());
    getInfo()->add (2, "assembly traversal stats");

}

