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

	OptionsParser* simplificationsParser = new OptionsParser ("graph simplifications");

    // this requires that Minia is compiled to support kmer values of 32 at least.
    Simplifications<GraphUnitigsTemplate<32>,NodeGU,EdgeGU> graphSimplifications(nullptr, 1, false); // get a graph simplifications object just to get default parameters

	simplificationsParser->push_back (new OptionNoParam  ("-no-bulge-removal", "ask to not perform bulge removal", false));
	simplificationsParser->push_back (new OptionNoParam  ("-no-tip-removal",   "ask to not perform tip removal", false));
	simplificationsParser->push_back (new OptionNoParam  ("-no-ec-removal",   "ask to not perform erroneous connection removal", false));

	simplificationsParser->push_back (new OptionOneParam ("-tip-len-topo-kmult", "remove all tips of length <= k * X bp",  false, to_string(graphSimplifications._tipLen_Topo_kMult)));
	simplificationsParser->push_back (new OptionOneParam ("-tip-len-rctc-kmult", "remove tips that pass coverage criteria, of length <= k * X bp",  false, to_string(graphSimplifications._tipLen_RCTC_kMult)));
	simplificationsParser->push_back (new OptionOneParam ("-tip-rctc-cutoff",    "tip relative coverage coefficient: mean coverage of neighbors >  X * tip coverage",  false, to_string(graphSimplifications._tipRCTCcutoff)));

	simplificationsParser->push_back (new OptionOneParam ("-bulge-len-kmult",    "bulges shorter than k*X bp are candidate to be removed",  false, to_string(graphSimplifications._bulgeLen_kMult)));
	simplificationsParser->push_back (new OptionOneParam ("-bulge-len-kadd",     "bulges shorter than k+X bp are candidate to be removed",  false, to_string(graphSimplifications._bulgeLen_kAdd)));
	simplificationsParser->push_back (new OptionOneParam ("-bulge-altpath-kadd", "explore up to k+X nodes to find alternative path",  false, to_string(graphSimplifications._bulgeAltPath_kAdd))); // TODO k should not appear in that equation
	simplificationsParser->push_back (new OptionOneParam ("-bulge-altpath-covmult", "bulges of coverage <= X*cov_altpath will be removed",  false, to_string(graphSimplifications._bulgeAltPath_covMult))); 

	simplificationsParser->push_back (new OptionOneParam ("-ec-len-kmult",       "EC shorter than k*X bp are candidates to be removed",  false, to_string(graphSimplifications._ecLen_kMult)));
	simplificationsParser->push_back (new OptionOneParam ("-ec-rctc-cutoff",     "EC relative coverage coefficient (similar in spirit as tip)",  false, to_string(graphSimplifications._ecRCTCcutoff)));

	getParser()->push_back (simplificationsParser);     

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
    string output = minia.assemble<GraphType, NodeGU, EdgeGU, span>(graph);

    // link contigs
    uint nb_threads = 1;  // doesn't matter because for now link_tigs is single-threaded
    bool verbose = true;
    link_tigs<span>(output, minia.k, nb_threads, minia.nbContigs, verbose);


    /** We gather some statistics. */
    minia.getInfo()->add (1, minia.getTimeInfo().getProperties("time"));
}
};

void Minia::execute ()
{
    /** we get the kmer size chosen by the end user. */
    k = getInput()->getInt (STR_KMER_SIZE);

    /** We launch Minia with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<MiniaFunctor,Parameter> (k, Parameter (*this));
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
    //ss1 << "NODE_"<< nbContigs + 1 << "_length_" << sequence.size() << "_cov_" << fixed << std::setprecision(3) << coverage << "_ID_" << nbContigs;
    // bcalm-like header (that can be converted to GFA)
    ss1 << nbContigs  << " LN:i:" << sequence.size() << " KC:i:" << (unsigned int)(coverage*(sequence.size()-k+1)) << " km:f:" << fixed << std::setprecision(3) << coverage ;
   
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
string Minia::assemble (/*const, removed because Simplifications isn't const anymore*/ Graph_type& graph)
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
        Simplifications<Graph_type,Node,Edge> graphSimplifications(&graph, nbCores, verbose);

        if (getParser()->saw("-no-tip-removal"))
            graphSimplifications._doTipRemoval = false;
        if (getParser()->saw("-no-bulge-removal"))
            graphSimplifications._doBulgeRemoval = false;
        if (getParser()->saw("-no-ec-removal"))
            graphSimplifications._doECRemoval = false;

		if (getParser()->saw("-tip-len-topo-kmult"))
			graphSimplifications._tipLen_Topo_kMult = getInput()->getDouble("-tip-len-topo-kmult");
		if (getParser()->saw("-tip-len-rctc-kmult"))
			graphSimplifications._tipLen_RCTC_kMult = getInput()->getDouble("-tip-len-rctc-kmult");
		if (getParser()->saw("-tip-rctc-cutoff"))
			graphSimplifications._tipRCTCcutoff = getInput()->getDouble("-tip-rctc-cutoff");

		if (getParser()->saw("-bulge-len-kmult"))
			graphSimplifications._bulgeLen_kMult = getInput()->getDouble("-bulge-len-kmult");
		if (getParser()->saw("-bulge-len-kadd"))
			graphSimplifications._bulgeLen_kAdd = getInput()->getDouble("-bulge-len-kadd");
		if (getParser()->saw("-bulge-altpath-kadd"))
			graphSimplifications._bulgeAltPath_kAdd = getInput()->getDouble("-bulge-altpath-kadd");
		if (getParser()->saw("-bulge-altpath-covmult"))
			graphSimplifications._bulgeAltPath_covMult = getInput()->getDouble("-bulge-altpath-covmult");

		if (getParser()->saw("-ec-len-kmult"))
			graphSimplifications._ecLen_kMult = getInput()->getDouble("-ec-len-kmult");
		if (getParser()->saw("-ec-rctc-cutoff"))
			graphSimplifications._ecRCTCcutoff = getInput()->getDouble("-ec-rctc-cutoff");
                                                                                                                                                                                                                  
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
    getInfo()->add (3, "bulges removed",          "%s", str_bubbleRemoval.c_str());
    getInfo()->add (3, "EC removed",          "%s", str_ECRemoval.c_str());
    getInfo()->add (2, "assembly traversal stats");

    return output;
}

