// comments added a posteriori, 1 year later:
// script that keeps reads if they match at least one kmer from a contig file
// probably this was done for multi-k optimization? 

// to compile it is a bit tricky: go to dsk repository, put that script in dsk/utils folder and use that CMakeLists file: https://github.com/GATB/dsk/blob/12de3e1f49e04f372b80773e7c4bc829c28e5657/utils/CMakeLists.txt

#include <gatb/gatb_core.hpp>
#include <fstream>
#include <set>

// We use the required packages
using namespace std;


static /* probably important that it's static here too */
char revcomp (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	else if (s == 'a') return 't';
	else if (s == 'c') return 'g';
	else if (s == 'g') return 'c';
	else if (s == 't') return 'a';
	return 'X';
}

static string revcomp (const string &s) {
	string rc;
	for (signed int i = s.length() - 1; i >= 0; i--) {rc += revcomp(s[i]);}
	return rc;
}

/********************************************************************************/
class FilterByPreviousContigs : public Tool
{
    public:

        /** */
        FilterByPreviousContigs () : Tool ("filter_by_previous_contigs")
    {
        getParser()->push_front (new OptionOneParam ("-previous-contigs",   "previous contigs"  , true));
        getParser()->push_front (new OptionOneParam ("-reads",              "reads"  , true));
        getParser()->push_front (new OptionOneParam (STR_KMER_SIZE,  "k-mer size used to generate previous contigs"  , true));
        getParser()->push_front (new OptionOneParam ("-out",                "output file that will contain filtered reads",           true));
    }

        /** */
        void execute ()
        {
            size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

            Integer::apply<Functor,Parameter> (kmerSize, Parameter(*this, kmerSize));
        }

        struct Parameter
        {
            Parameter (FilterByPreviousContigs& tool, size_t kmerSize) : tool(tool), kmerSize(kmerSize)  {}
            FilterByPreviousContigs& tool;
            size_t     kmerSize;
        };

        template<size_t span> struct Functor  {  void operator ()  (Parameter parameter)
            {
                FilterByPreviousContigs& tool     = parameter.tool;
                size_t                   kmerSize = parameter.kmerSize;
                unsigned int k = kmerSize;

                typedef typename Kmer<span>::Count Count;

                typename Kmer<span>::ModelCanonical model (kmerSize);

                BankFasta bank(tool.getInput()->getStr("-previous-contigs"));
                IBank* bankReads = Bank::open(tool.getInput()->getStr("-reads"));
                BankFasta bankOut(tool.getInput()->getStr("-out"));

                BankFasta::Iterator itCtg (bank);
                set<string> kmers;

                std::cout << "indexing previous contigs" << std::endl;

                for (itCtg.first(); !itCtg.isDone(); itCtg.next())
                {
                    string seq = itCtg->toString();
                    string kmerBegin = seq.substr(0, k );
                    string kmerEnd = seq.substr(seq.size() - k , k );

                    kmers.insert(kmerBegin);
                    kmers.insert(kmerEnd);
                    kmers.insert(revcomp(kmerBegin));
                    kmers.insert(revcomp(kmerEnd));
                }
                
                std::cout << "filtering reads" << std::endl;

                Iterator<Sequence>* itReads = bankReads->iterator();
                uint64_t nb_reads = 0, nb_kept = 0;

                for (itReads->first(); !itReads->isDone(); itReads->next())
                {
                    bool contains = false;
                    string seq = (*itReads)->toString();
                    nb_reads++;
                    int found_pos = 0;
                    for (int i = 0; i < seq.size() - k; i++)
                    {
                        string kmer = seq.substr(i, k );
                        contains = (kmers.find(kmer) != kmers.end());
                        if (contains) 
                        {
                            found_pos = i;
                            break;
                        }
                    }
                    if (contains)
                    {
                        Sequence s (Data::ASCII);
                        s.getData().setRef ((char*)seq.c_str(), seq.size());
                        s._comment = (*itReads)->getComment() + " contains previous-k kmer at pos " + to_string(found_pos);
                        bankOut.insert(s);
                        nb_kept++;
                    }
                }

                std::cout << "done filtering reads, kept " << nb_kept << "/" << nb_reads << std::endl;

                /** We gather some statistics. */
                tool.getInfo()->add (1, "stats");
                tool.getInfo()->add (2, "kmer_size", "%ld", kmerSize);
            }};
};

/********************************************************************************/
/*                       Dump solid kmers in ASCII format                       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /*    try
          {*/
    FilterByPreviousContigs().run (argc, argv);
    /*    }
          catch (Exception& e)
          {
          std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
          return EXIT_FAILURE;
          }*/ 
    // gdb likes having clean exceptions
}

