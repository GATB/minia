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

#include <gatb/gatb_core.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include <iomanip> // for setprecision

const int max_delta = 15; // important parameter

using namespace std;

static void file_copy(string from, string to) // http://stackoverflow.com/questions/10195343/copy-a-file-in-a-sane-safe-and-efficient-way
{
    std::ifstream  src(from, std::ios::binary);
    std::ofstream  dst(to,   std::ios::binary);
    dst << src.rdbuf();
}

static void file_append(string dest, string with_file)
{
    std::ofstream fout ( dest, std::ios_base::app);
    std::ifstream f2 ( with_file);
    fout << f2.rdbuf();
}

template <int span> using assembly_index_t = std::unordered_map<typename Kmer<span>::Type,uint64_t>;

typedef std::unordered_map<uint64_t,set<uint64_t>> connections_index_t ;
typedef std::unordered_map<uint64_t,set<string>> connections_t ;

template<int span>
void index_assembly(string assembly, int k, int nb_threads, bool verbose, assembly_index_t<span>& index)
{
    if (verbose)
        std::cout << "indexing assembly extremities" << std::endl;

    class IndexAssembly
    {
        typedef typename Kmer<span>::ModelCanonical ModelCanon;
        
        int k;
        ModelCanon modelCanon;
        uint64_t current_seq;
        assembly_index_t<span> &index;
        int _currentThreadIndex; // for later
        bool debug;

        public: 
        IndexAssembly(int k, assembly_index_t<span> &index) : 
                      k(k), modelCanon(k), current_seq(0),
                         index(index), _currentThreadIndex(-1), debug(false)
        {}
 
        void operator()     (const Sequence& sequence) {
            const string seq = sequence.toString(); 
            const string comment = sequence.getComment();

            const bool no_link_left = comment.find("L:-") == string::npos;
            const bool no_link_right = comment.find("L:+") == string::npos;

            if (no_link_left)
            {
                const string kmerBegin = seq.substr(0, k );
                const typename ModelCanon::Kmer kmmerBegin = modelCanon.codeSeed(kmerBegin.c_str(), Data::ASCII);
                ExtremityInfo e (current_seq, false, UNITIG_BEGIN); /* actually we dont know if its false/UNITIG_BEGIN or reverse/UNITIG_END but should think about it, for now, putting naive*/ 
                index[kmmerBegin.value()] = e.pack();
                if (debug)
                std::cout << "index left kmer " << modelCanon.toString(kmmerBegin.value()) << " of contig " << current_seq << std::endl; 
            }

            if (no_link_right)
            {
                const string kmerEnd = seq.substr(seq.size() - k , k );
                const typename ModelCanon::Kmer kmmerEnd = modelCanon.codeSeed(kmerEnd.c_str(), Data::ASCII);
                ExtremityInfo e (current_seq, false, UNITIG_END); 
                index[kmmerEnd.value()] = e.pack();
                if (debug)
                std::cout << "index right kmer " << modelCanon.toString(kmmerEnd.value()) << " of contig " << current_seq << std::endl; 
            }
            current_seq++;
        }
    };

    IndexAssembly indexAssembly(k, index);

    IBank *in = Bank::open (assembly);
    auto it = in->iterator();    
    for (it->first (); !it->isDone(); it->next())
        indexAssembly(it->item());
}

static void 
insert_connection(uint64_t previous_extremityinfo, uint64_t current_extremityinfo, connections_index_t &connections_index, connections_t &connections, int previous_pos, int current_pos, const string &seq, int k, bool verbose)
{
    /* there are two ways to write a connection and we need to write them both */
    ExtremityInfo pe(previous_extremityinfo); 
    ExtremityInfo ce(current_extremityinfo); 

    if (pe.unitig == ce.unitig) return; // one case where we actually don't insert
    
    connections_index[pe.pack()].insert(ce.pack());
    connections_index[ce.pack()].insert(pe.pack());

    // this the part of the read that is the connection sequence
    string connection = seq.substr(previous_pos,current_pos-previous_pos+k);
    connections[pe.pack()].insert(connection);
    connections[ce.pack()].insert(connection);

    //if (verbose)        std::cout << "found connection between contig " << pe.unitig << " and contig " << ce.unitig << " : " << seq.substr(previous_pos,current_pos-previous_pos+k) << std::endl; 
}


template<int span>
void search_reads(string reads, int k, int nb_threads, bool verbose, assembly_index_t<span>& assembly_index, connections_index_t &connections_index, connections_t &connections)
{
    if (verbose)
        std::cout << "searching reads for assembly extremities" << std::endl;

    std::mutex insertMutex;

    class SearchReads 
    {
        typedef typename Kmer<span>::ModelCanonical ModelCanon;
        
        int k;
        bool verbose;
        ModelCanon modelCanon;
        connections_index_t &connections_index;
        connections_t &connections;
        std::mutex &insertMutex;
        assembly_index_t<span> &assembly_index;
        int _currentThreadIndex; // for later

        public: 
        uint64_t nb_found_connections; 

        SearchReads(int k, bool verbose, assembly_index_t<span> &assembly_index, connections_index_t &connections_index, connections_t &connections, std::mutex &insertMutex) : 
                      k(k), verbose(verbose), modelCanon(k), connections_index(connections_index), connections(connections), insertMutex(insertMutex),
                         assembly_index(assembly_index), _currentThreadIndex(-1), nb_found_connections(0)
        {}
 
        void operator()     (/* not making it const because if itKmer */ Sequence& sequence) {
            const string seq = sequence.toString(); 
            const string comment = sequence.getComment();

            typename Kmer<span>::ModelCanonical::Iterator itKmer (modelCanon);
            itKmer.setData (sequence.getData());
            int previous_pos = -1;
            int current_pos = 0;
            uint64_t previous_extremityinfo = 0;
            for (itKmer.first(); !itKmer.isDone(); itKmer.next())
            {
                typedef typename Kmer<span>::Type Type;

                Type kmer =  itKmer.item().value();
                //std::cout << "testing kmer " << modelCanon.toString(kmer) << std::endl;
                auto hit = assembly_index.find(kmer) ;
                if (hit != assembly_index.end())
                {
                    uint64_t current_extremityinfo = hit->second;
                    if (previous_pos != -1 && current_pos - previous_pos < max_delta)
                    {
                        insertMutex.lock();
                        insert_connection(previous_extremityinfo, current_extremityinfo, connections_index, connections, previous_pos, current_pos, seq, k, verbose);
                        nb_found_connections++;
                        insertMutex.unlock();
                    }
                    previous_pos    = current_pos;
                    previous_extremityinfo = current_extremityinfo;
                }
                current_pos++;
            }
        }
    };

    SearchReads searchReads(k, verbose, assembly_index, connections_index, connections, insertMutex);

    IBank *in = Bank::open (reads);
    Dispatcher dispatcher (nb_threads);
    dispatcher.iterate (in->iterator(), searchReads);
    /*auto it = in->iterator();     // single threaded version
    for (it->first (); !it->isDone(); it->next())
        searchReads(it->item());
    */
    if (verbose) std::cout << "found " << searchReads.nb_found_connections << " connections" << std::endl;
}

/* examines connections_index to see how contigs are connected to each other
 * returns true for "compactable" nodes, those which have single out-neighbor and single-inneighbor in the connection graph
 * oh indeed, this is a compaction problem all over again
 * well no wonder why i'm using bglue after then
 *
 * unpure function: updates seqs_to_glue
 */
static bool check_connections(uint64_t packed, connections_index_t &connections_index, connections_t &connections, set<uint64_t> &seqs_to_glue, bool given_verbose)
{
    bool verbose = false; // = given_verbose

    if (connections_index.find(packed) == connections_index.end())
        return false;

    ExtremityInfo current(packed);

    set<uint64_t> &links = connections_index[packed];
    if (links.size() != 1)
    {
        if (verbose) std::cout << "bad conn, contig: " << current.unitig << " connections: " << links.size() << std::endl;
        return false;
    }

    ExtremityInfo other(*connections_index[packed].begin());

    // check the incoming links of that other
    if (connections_index[other.pack()].size() != 1)
    {
        if (verbose) std::cout << "bad conn, other contig: " << current.unitig << " connections: " << links.size() << std::endl;
        return false;
    }

    // all good, we keep this one
    // and to avoid redundancy, we delete the other
    if (verbose) std::cout << "clearing connections of contig " << other.unitig << std::endl;
    connections[other.pack()].clear();
    seqs_to_glue.insert(other.pack());
    
    if (verbose) std::cout << "good conn: " << current.unitig << " connections: " << links.size() << std::endl;
    
    return connections[packed].size() == 1; // may return false if it was previously cleared
}


// there are not so many places where i need a sequence revcomp. GraphUnitigs and here.
static
char revcomp (char s) {
    if (s == 'A') return 'T';
    else if (s == 'C') return 'G';
    else if (s == 'G') return 'C';
    else if (s == 'T') return 'A';
    else if (s == 'a') return 't';
    else if (s == 'c') return 'g';
    else if (s == 'g') return 'c';
    else if (s == 't') return 'a';
    return 'N';
}

static string revcomp (const string &s) {
    string rc;
    for (signed int i = s.length() - 1; i >= 0; i--) {rc += revcomp(s[i]);}
    return rc;
}

/* naive, cause I'm lazy */
static void 
somehow_merge(string &seq, const string &extension, unsigned int k, bool verbose)
{
    const string rev_extension = revcomp(extension);

    unsigned last_pos = seq.size() - k;
    const vector<unsigned int> positions {0, last_pos};
    for (auto pos: positions)
    {
        const string extrem = seq.substr(pos,k);
        const vector<string> extensions {extension, rev_extension};
        for (auto e: extensions)
        {
            if (pos == 0)
            {
                /*     [     ]---seq----
                 *  -e-[     ]
                 */
                if (e.substr(e.size()-k) == extrem)
                {
                    seq = e.substr(0,e.size()-k) + seq;
                    return;
                }
            }
            else
            {
                /*    ----seq-----[      ]
                 *                [      ]--e--
                 */
                if (e.substr(0,k) == extrem)
                {
                    seq = seq + e.substr(k);
                    return;
                }
            }
        }
    }
    std::cout << "oops, couldn't somehow merge " << seq << " with " << extension << std::endl;
}
static void 
extend_assembly_with_connections(const string assembly, int k, int nb_threads, bool verbose, connections_index_t &connections_index, connections_t &connections, BankFasta &out, BankFasta &glue)
{
    if (verbose)
        std::cout << "extending assembly with unambiguous connections" << std::endl;
    IBank *in = Bank::open (assembly);
    auto it = in->iterator();    
    uint64_t tig_index = 0;
    set<uint64_t> seqs_to_glue;

    for (it->first (); !it->isDone(); it->next())
    {
        const Sequence &sequence = it->item();
        string seq = sequence.toString(); 
        const string comment = sequence.getComment();

        ExtremityInfo beg(tig_index, false, UNITIG_BEGIN);
        ExtremityInfo end(tig_index, false, UNITIG_END);
        bool lmark, rmark;

        // see if the contig was already previously marked as one to glue
        lmark = (seqs_to_glue.find(beg.pack()) != seqs_to_glue.end());
        rmark = (seqs_to_glue.find(end.pack()) != seqs_to_glue.end());

        //if (verbose) if (lmark || rmark) std::cout << "not extending " << tig_index << " bc previously set to glue with another extremity" << std::endl;
        
        if ((!lmark) && check_connections(beg.pack(), connections_index, connections, seqs_to_glue, verbose)) // order of evaluation matters, as check_connection is unpure
        {
            const string extension = *(connections[beg.pack()].begin());
            somehow_merge(seq,extension,k, verbose);
            lmark = true;
        }

        if ((!rmark) && check_connections(end.pack(), connections_index, connections, seqs_to_glue, verbose))
        {
            const string extension = *connections[end.pack()].begin();
            somehow_merge(seq,extension, k, verbose);
            rmark = true;
        }

        Sequence s (Data::ASCII);
        s.getData().setRef ((char*)seq.c_str(), seq.size());
        s._comment = string(lmark?"1":"0")+string(rmark?"1":"0"); //We set the sequence comment.
        s._comment += " ";
        
        if (lmark || rmark)
            glue.insert(s); 
        else 
        {
            s._comment = comment;
            out.insert(s);
        }

        tig_index++;
    }
}

template<int span>
void merci(int k, string reads, string assembly, int nb_threads, bool verbose)
{
    if (verbose)
        cout << "k: " << k << " threads: " << nb_threads << " reads: " << reads << " assembly: " << assembly << endl;


    // make a copy of the assembly, and find links
    string linked_assembly = assembly + ".linked";
    file_copy(assembly, linked_assembly);
    uint64_t nb_tigs = 0;
    link_tigs<span>( linked_assembly, k, nb_threads, nb_tigs, verbose);

    // index extremities of assembly that do not have links
    assembly_index_t<span> index;
    index_assembly<span>(linked_assembly, k, nb_threads, verbose, index);

    // now pass the reads and find putative connections
    connections_index_t connections_index; // associates a contig to a connection
    connections_t connections; // a connection is basically a chunk of read linking two indexed kmers close to each other
    search_reads<span>(reads,k,nb_threads,verbose, index, connections_index, connections);

    // pass the assembly again and slightly extend some contigs with unambiguous connections
    System::file().remove (assembly+".merci");
    BankFasta out (assembly+".merci");
    BankFasta glue(assembly+".glue.glue"); // bglue will glue the .glue.glue to .glue
    extend_assembly_with_connections(assembly, k, nb_threads, verbose, connections_index, connections, out, glue);
    glue.flush();

    // glue what needs to be glued. magic, we're re-using bcalm code
    bglue<span> (nullptr /*no storage*/, assembly+".glue", k, nb_threads, verbose);
    
    // append glued to merci
    out.flush();
    file_append(assembly+".merci", assembly+".glue");
}

class Merci : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    Merci () : Tool("merci")
    {

#ifdef GIT_SHA1
        std::cout << "Merci for Minia 3, git commit " << GIT_SHA1 << std::endl;
#endif

        // reinit the parser to get rid of options added by the Tool class
        setParser (new OptionsParser ("merci")); 
        OptionsParser* assemblyParser = new OptionsParser ("assembly");
        assemblyParser->push_front (new OptionOneParam ("-reads",       "input reads (fasta/fastq/compressed)",   true));
        assemblyParser->push_front (new OptionOneParam ("-assembly",    "assembly to improve",   true));
        assemblyParser->push_front (new OptionOneParam (STR_KMER_SIZE,    "kmer size",   true));
        assemblyParser->push_front (new OptionOneParam (STR_NB_CORES,    "number of threads",   false, "0")); // apparently its needed when there is no parser
        assemblyParser->push_front (new OptionOneParam (STR_VERBOSE,  "verbosity level", false, "0"));
        getParser()->push_back (assemblyParser);
    }

    struct Parameter
    {
        string reads, assembly;
        int k,t;
        bool verbose;
    };
    
    template<size_t span>
    struct MerciFunctor
    {
        void operator() (Parameter parameter)
        {
            merci<span>(parameter.k,parameter.reads,parameter.assembly,parameter.t,parameter.verbose);
        }
    };

    /** \copydoc Tool::execute. */
    void  execute ()
    {
        /** we get the kmer size chosen by the end user. */
        int kmerSize = getInput()->getInt (STR_KMER_SIZE);
        string reads = getInput()->getStr ("-reads");
        string assembly = getInput()->getStr ("-assembly");
        int nb_threads = getInput()->getInt(STR_NB_CORES);
        bool verbose =  getInput()->getInt("-verbose") == 1; 
        /** We launch Merci with the correct Integer implementation according to the choosen kmer size. */
        Integer::apply<MerciFunctor,Parameter> (kmerSize, Parameter{reads,assembly,kmerSize, nb_threads, verbose});
        
        getInfo()->add (1, getTimeInfo().getProperties("time"));
    }


private:


};


/********************************************************************************/

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We execute dsk. */
        Merci().run (argc, argv);
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }

    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

