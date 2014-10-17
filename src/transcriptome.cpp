/*    
Copyright (C) 2014 Gavin Wilson, Lincoln Stein

This file is part of RNASequel.

RNASequel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RNASequel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RNASequel.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "transcriptome.hpp"
#include "build_transcriptome.hpp"
#include "junction.hpp"
#include "fasta_index.hpp"
#include "tokenizer.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <unordered_set>
#include "gtf.hpp"
#include "seed_iterator.hpp"
#include "reader.hpp"
#include "timer.hpp"

namespace po = boost::program_options;

using namespace std;
using namespace rnasequel;

class SpliceSiteStrand {
    public:
        struct splice_strand {
            splice_strand() : tid(-1), pos(0) {

            }

            splice_strand(int tid, unsigned int pos) : tid(tid), pos(pos) {

            }

            int          tid;
            unsigned int pos;

            bool operator<(const splice_strand & s) const {
                return tid < s.tid || (tid == s.tid && pos < s.pos);
            }
        };

        void add_junction(const Junction & j){ 
            auto ret_lft = lfts_.insert(make_pair(splice_strand(j.tid, j.lft), j.strand));
            if(!ret_lft.second){
                if(ret_lft.first->second != j.strand){
                    ret_lft.first->second = BOTH;
                }
            }
            auto ret_rgt = rgts_.insert(make_pair(splice_strand(j.tid, j.rgt), j.strand));
            if(!ret_rgt.second){
                if(ret_rgt.first->second != j.strand){
                    ret_rgt.first->second = BOTH;
                }
            }
        }

        Strand find_lft(int tid, unsigned int p) const {
            auto it = lfts_.find(splice_strand(tid, p));
            if(it == lfts_.end()) return BOTH;
            return it->second;
        }

        Strand find_rgt(int tid, unsigned int p) const {
            auto it = rgts_.find(splice_strand(tid, p));
            if(it == rgts_.end()) return BOTH;
            return it->second;
        }

    private:
        std::map<splice_strand, Strand> lfts_;
        std::map<splice_strand, Strand> rgts_;
};

union SpliceSiteType {
    struct {
	uint8_t donor:1;
	uint8_t acceptor:1;
	uint8_t plus:1;
	uint8_t minus:1;
	uint8_t pad_:4;
    };

    uint8_t site_val;

    SpliceSiteType() : site_val(0) {

    }

    Strand infer(){
        if((donor && acceptor) || (plus && minus) || (!plus && !minus)) return UNKNOWN;
        return plus ? PLUS : MINUS;
    }
};

inline std::ostream & operator<<(std::ostream & os, SpliceSiteType ss){
    std::string s;
    if(ss.donor)    s.append("D");
    if(ss.acceptor) s.append("A");
    if(ss.plus)     s.append("+");
    if(ss.minus)    s.append("-");
    os << s;
    return os;
}

class SpliceSites {
    public:
	void add_lft(unsigned int pos, Strand strand){
	    if(strand == PLUS){
		sites_[pos].donor = true;
		sites_[pos].plus  = true;
	    }else{
		sites_[pos].acceptor = true;
		sites_[pos].minus    = true;
	    }
	}

	void add_rgt(unsigned int pos, Strand strand){
	    if(strand == PLUS){
		sites_[pos].acceptor = true;
		sites_[pos].plus     = true;
	    }else{
		sites_[pos].donor = true;
		sites_[pos].minus = true;
	    }
	}

	size_t size() const {
	    return sites_.size();
	}

	SpliceSiteType find(unsigned int pos) const {
	   std::map<unsigned int, SpliceSiteType>::const_iterator it = sites_.find(pos);
	   if(it == sites_.end()) return SpliceSiteType();
	   return it->second;
	}

    private:
	std::map<unsigned int, SpliceSiteType> sites_;
};

class SpliceSiteMap {
    public:
	typedef std::map<std::string, SpliceSites> ChromMap;

	void build_from_model(const Model & m){
	    for(Model::const_iterator it = m.begin(); it != m.end(); it++){
		SpliceSites & ptr = sites_[it->first];
		const Model::gene_list & g = it->second;
		for(size_t i = 0; i < g.size(); i++){
		    for(size_t j = 0; j < g[i].transcripts().size(); j++){
			const Transcript & tx = g[i].transcripts()[j];
			for(size_t k = 1; k < tx.exons().size(); k++){
			    ptr.add_lft(tx.exons()[k - 1].rgt(), tx.strand());
			    ptr.add_rgt(tx.exons()[k].lft(), tx.strand());
			}
		    }
		}
	    }
	}

	SpliceSiteType find(const std::string & chrom, unsigned int pos) const {
	    ChromMap::const_iterator it = sites_.find(chrom);
	    if(it == sites_.end()){
		return SpliceSiteType();
	    }
	    return it->second.find(pos);
	}

    private:
	ChromMap sites_;
};

class SpliceStrand {
    public:
	SpliceStrand(const std::string & sites) {
	    sites_.resize(256, UNKNOWN);
            Tokenizer::token_t toks;
            std::string cpy = sites;
            Tokenizer::get(cpy, ',', toks);

	    for(size_t i = 0; i < toks.size(); i++){
		PackedSequence seq(toks[i]);
		uint8_t plus  = code(seq[0], seq[1], seq[2], seq[3]);
		uint8_t minus = code(seq[3].cmpl(), seq[2].cmpl(), seq[1].cmpl(), seq[0].cmpl());
                //std::cout << "    Seq: " << seq << " + code = " << static_cast<int>(plus) << " - code = " << static_cast<int>(minus) << "\n";
		assert(sites_[plus] == UNKNOWN && sites_[minus] == UNKNOWN);
                if(sites_[plus] != UNKNOWN || sites_[minus] != UNKNOWN){
                    std::cerr << "Error ambiguous splice site seq = " << seq
                              << " plus = " << packed2string<uint8_t, 2>(plus, 4)
                              << " minus = " << packed2string<uint8_t, 2>(minus, 4)
                              << "\n";
                    exit(0);
                }
		sites_[plus]  = PLUS;
		sites_[minus] = MINUS;
	    }
	}

	uint8_t code(PackedBase b1, PackedBase b2, PackedBase b3, PackedBase b4) const {
	    return (b1.nt4() << 6) | (b2.nt4() << 4) | (b3.nt4() << 2) | b4.nt4();
	}

	Strand strand(uint8_t site) const {
	    return sites_[site];
	}

	Strand strand(PackedBase b1, PackedBase b2, PackedBase b3, PackedBase b4) const {
	    return strand(code(b1, b2, b3, b4));
	}

	bool canonical(uint8_t site) const {
	    return strand(site) != UNKNOWN;
	}

	bool canonical(PackedBase b1, PackedBase b2, PackedBase b3, PackedBase b4) const {
	    return canonical(code(b1, b2, b3, b4));
	}

    private:
	std::vector<Strand>         sites_;

};

void tx_init_options(int argc, char *argv[], po::variables_map & vm) {
    po::options_description generic("Arguments");
    generic.add_options()
    ("ref,r", po::value< string >(), "Reference sequence fasta file")
    ("out,o", po::value< string >(), "Output prefix")
    ("gtf,g", po::value< string >(), "Gene annotation GTF file (optional)")
    ("skip,s", po::value< string >()->default_value("MT,chrM"), "Comma separated list of chromosomes to skip")
    ("bam,b", po::value< string >(), "The bam file for de novo junctions (optional)")
    ("read-size,n", po::value< unsigned int >(), "Read Size")
    ("max-iter", po::value< unsigned int >()->default_value(1000), "Maximum number of graph iterations before giving up on a locus")
    ("debug,d", "Whether the reads are stranded or not")
    ("help,h", "help message")
    ;

    po::options_description filter("Junction Filtering");
    filter.add_options()
    ("min-intron", po::value< unsigned int >()->default_value(21), "Mininum intron size")
    ("max-intron", po::value< unsigned int >()->default_value(500000), "Maxinum intron size")
    ("min-unique", po::value< unsigned int >()->default_value(2), "Minimum number of unique reads mapping across the junction")
    ("min-count", po::value< unsigned int >()->default_value(2), "Minimum number of reads mapping across the junction (includes repeats if enabled)")
    ("min-nc-unique", po::value< unsigned int >()->default_value(5), "Minimum number of unique reads mapping across the junction for non-canonical junctions")
    ("min-nc-count", po::value< unsigned int >()->default_value(5), "Minimum number of reads mapping across the junction (includes repeats if enabled) for non-canonical junctions")
    ("canonical", po::value< string>()->default_value("GTAG,GCAG,GCTG,GCAA,GCGG,GTTG,GTAA,ATAC,ATAA,ATAG"), "Canonical splice sites for strand inference")
    ("use-repeats", "Include junctions contained in the primary alignment of multi-mapped pairs")
    ("min-length", po::value<unsigned int>()->default_value(8), "Minimum exonic alignment length")
    ("end-cutoff", po::value<unsigned int>()->default_value(10), "End cutoff ie. junctions are counted as being int the end of the read if they are in the last end-cutoff bases")
    ("end-count", po::value< unsigned int >()->default_value(1), "Require at least this many reads outside of the ends of the junction")
    ;


    po::options_description cmdline_options;
    cmdline_options.add(generic).add(filter);

    po::options_description visible;
    visible.add(generic).add(filter);

    //po::positional_options_description pd;
    //pd.add("juncs",1);

    po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << "Usage: " << endl;
        cout << "rnasequel " << string(argv[0]) << " [options]\n" << endl;
        cout << "Extracts junctions from an optional bam file and merges them with\nannotated junctions in a GTF file (optional)\n";
        cout << "A transcriptome index is then generated using the specified read size\n";
        cout << "**Note only proper read pairs are used for to identify de novo splice junctions\n";

        cout << visible << "\n";
        exit(0);
    }

    bool error = false;

    if(vm.count("ref") == 0) {
        cout << "The reference index must be specified\n";
        error = true;
    }

    if(vm.count("read-size") == 0) {
        cout << "Read size must be specified\n";
        error = true;
    }

    if(vm.count("out") == 0) {
        cout << "Output file must be specified\n";
        error = true;
    }

    if(vm.count("bam") == 0 && vm.count("gtf") == 0){
        cout << "There must at least be a GTF or bam file specified\n";
        error = true;
    }

    if(error) exit(1);
}

struct JunctionCount {
    JunctionCount() : count(0), ends(0) {

    }

    std::unordered_set<uint64_t> positions;
    unsigned int                 count;
    unsigned int                 ends;
};

bool is_skip(CigarElement c) { return c.op == REF_SKIP; }

int rnasequel::build_transcriptome(int argc, char * argv[]) {
    Timer ti("Total time building the transcriptome");
    po::variables_map vm;
    tx_init_options(argc,argv,vm);

    FastaIndex fi(vm["ref"].as<string>());
    fi.load_all(false);
    unsigned int read_size  = vm["read-size"].as<unsigned int>();

    std::vector<bool>   tskips(fi.size(), false);


    {
	string skip = vm["skip"].as<string>();
	Tokenizer::token_t skips;
	Tokenizer::get(skip, ',', skips);
	for(size_t i = 0; i < skips.size(); i++){
	    FastaIndex::const_iterator it = fi.get_entry(skips[i]);
	    if(it != fi.end()) tskips[it->tid] = true;
	}
    }

    JunctionSet pjuncs, mjuncs, ajuncs;
    std::map<Junction, Strand> strand_map;
    SpliceSiteStrand strander;
    SpliceSiteMap splice_sites;

    if(vm.count("gtf") > 0) {
        int min_intron = vm["min-intron"].as<unsigned int>();
	Model m;
	GTF::parse(vm["gtf"].as<string>(), m);
        splice_sites.build_from_model(m);
        size_t kept = 0;
	for(Model::iterator it = m.begin(); it != m.end(); it++){
	    unsigned int tid = fi.tid(it->first);
	    if(tskips[tid]) continue;
	    for(size_t i = 0; i < it->second.size(); i++){
		for(Gene::iterator it2 = it->second[i].begin(); it2 != it->second[i].end(); it2++){
		    Transcript::junc_iterator jit = it2->junc_it();
		    while(jit.next()){
                        int isize = jit->rgt() - jit->lft() - 1;
                        if(isize < min_intron) continue;
                        Junction j(tid, jit->lft(), jit->rgt(), jit->strand());
			if(it2->strand() == PLUS){
                            strander.add_junction(j);
			    if(pjuncs.insert(j).second) kept++;
			}else{
                            strander.add_junction(j);
			    if(mjuncs.insert(j).second) kept++;
			}
		    }
		}
	    }
	}

        cout << "  Kept " << kept << " annotated junctions from the gtf file\n";
    }

    if(vm.count("bam") > 0){
        bool use_repeats = vm.count("use-repeats") > 0;
        SpliceStrand infer_strand(vm["canonical"].as<string>());
        std::map<Junction, JunctionCount> counts;
        BamReader bin(vm["bam"].as<string>());
        BamRead r;

        size_t min_length = vm["min-length"].as<unsigned int>();
        unsigned int end_cutoff = vm["end-cutoff"].as<unsigned int>();

        std::vector<int> tidmap(bin.header().size());
        for(size_t i = 0; i < bin.header().size(); i++){
            tidmap[i] = fi[bin.header()[i]].tid;
        }

        cout << "Extracting novel splice junctions" << endl;
        Timer ti2("Time extracting junctions");
        size_t total = 0;
        while(bin.get_read(r)){
            total++;
            if(total % 0x40000 == 0){
                time_t e = ti2.elapsed();
                if(e > 0){
                    std::cerr << "\33[2K\rTotal processed: " << total 
                        << ",  reads / second: " << (total / e);
                    std::cerr.flush();
                }
            }

            if(!r.aligned() || !r.flag.paired || !r.flag.proper_pair || r.flag.secondary){
                continue;
            }
            
            size_t num_juncs = std::count_if(r.cigar.begin(), r.cigar.end(), is_skip);
            if(num_juncs == 0 || (!use_repeats && r.repeat())) continue;

            uint64_t block = 0;
            if(r.flag.paired && r.flag.proper_pair){
                int tlen = r.tlen();
                if(tlen < 0){
                    block = (static_cast<uint64_t>(r.rgt() + r.cigar.back_clipped() - (tlen - 1)) << 32) | static_cast<uint64_t>(r.lft() - r.cigar.front_clipped() + tlen - 1);
                }else{
                    block = (static_cast<uint64_t>(r.lft()) << 32) | static_cast<uint64_t>(r.lft() + tlen - 1);
                }
            }

            size_t junc_num = 1;
            SeedIterator<Cigar::const_iterator> bit(r.cigar.begin(), r.cigar.end(), r.lft(), 0);
            Seed l = bit();

            int tid = tidmap[r.tid()];
            const PackedSequence & ref = fi[tid].seq;
            unsigned int qrgt = r.qrgt();

            while(bit.next()){
                Junction j(tid, l.rrgt(), bit().rlft(), UNKNOWN);
                Strand s1 = splice_sites.find(r.tname(), j.lft).infer();
                Strand s2 = splice_sites.find(r.tname(), j.rgt).infer();
                Strand strand;

                if((s1 == UNKNOWN && s2 == UNKNOWN) || (s1 != UNKNOWN && s2 != UNKNOWN && s1 != s2)){
                    strand = infer_strand.strand(ref[j.lft + 1], ref[j.lft + 2], ref[j.rgt - 2], ref[j.rgt - 1]);
                }else{
                    strand = s1 == UNKNOWN ? s2 : s1;
                }

                bool end = l.qrgt() < end_cutoff || r.qlft() >= (qrgt - end_cutoff);

                j.strand = strand;

                if((junc_num != 1 || l.score() >= static_cast<int>(min_length)) 
                    && (junc_num != num_juncs || bit().score() >= static_cast<int>(min_length))){

                    JunctionCount & jc = counts[j];
                    jc.count++;
                    jc.positions.insert(block);
                    if(end) jc.ends++;
                }
                l = bit();
                junc_num++;
            }
        }


        unsigned int min_intron    = vm["min-intron"].as<unsigned int>();
        unsigned int max_intron    = vm["max-intron"].as<unsigned int>();
        unsigned int min_unique    = vm["min-unique"].as<unsigned int>();
        unsigned int min_count     = vm["min-count"].as<unsigned int>();
        unsigned int end_count     = vm["end-count"].as<unsigned int>();
        unsigned int min_nc_unique    = vm["min-nc-unique"].as<unsigned int>();
        unsigned int min_nc_count     = vm["min-nc-count"].as<unsigned int>();

        size_t kept = 0, total_novel = 0;
        size_t nc_kept = 0, nc_total_novel = 0;

        for(auto const & j : counts){
            const Junction & junc    = j.first;
            const JunctionCount & jc = j.second;
            unsigned int isize = junc.isize();
            const PackedSequence & ref = fi[junc.tid].seq;
            bool canonical = infer_strand.canonical(ref[junc.lft + 1], ref[junc.lft + 2], ref[junc.rgt - 2], ref[junc.rgt - 1]);

            bool found = false;
            if(junc.strand == PLUS){
                found = pjuncs.find(junc) != pjuncs.end();
            }else if(junc.strand == MINUS){
                found = mjuncs.find(junc) != mjuncs.end();
            }else{
                found = ajuncs.find(junc) != ajuncs.end();
            }
            
            if(!found) {
                if(canonical) total_novel++;
                else          nc_total_novel++;
            }

            if((canonical && (jc.count < min_count || jc.positions.size() < min_unique)) ||
               (!canonical && (jc.count < min_nc_count || jc.positions.size() < min_nc_unique)) || 
                isize < min_intron || isize > max_intron || (jc.count - jc.ends) < end_count){
                continue;
            }

            if(!found) {
                if(junc.strand == PLUS){
                    pjuncs.insert(junc);
                }else if(junc.strand == MINUS){
                    mjuncs.insert(junc);
                }else{
                    ajuncs.insert(junc);
                }
                if(canonical) kept++;
                else          nc_kept++;
            }
        }
        cout << "  Kept " << kept << " out of " << total_novel << " novel junctions\n";
        cout << "  Kept " << nc_kept << " out of " << nc_total_novel << " non-canonical novel junctions\n";
    }

    for(auto const & j : ajuncs){
        Strand lft = strander.find_lft(j.tid, j.lft);
        Strand rgt = strander.find_rgt(j.tid, j.rgt);
        Strand strand = BOTH;

        if(lft != BOTH && rgt != BOTH && lft == rgt){
            strand = lft;
        }else if(lft == BOTH && rgt != BOTH){
            strand = rgt;
        }else if(lft != BOTH && rgt == BOTH){
            strand = lft;
        }

        if(strand == PLUS){
            pjuncs.insert(j);
        }else if(strand == MINUS){
            mjuncs.insert(j);
        }else{
            pjuncs.insert(j);
            mjuncs.insert(j);
        }
    }

    Transcriptome builder(vm["out"].as<string>(), fi, read_size, vm["max-iter"].as<unsigned int>(), false);
    builder.process_junctions(pjuncs, PLUS);
    builder.process_junctions(mjuncs, MINUS);

    return 0;
}
