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

#include "merge.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <boost/thread/thread.hpp>
#include <boost/utility.hpp>
#include "read_pairs.hpp"
#include "pair_builder.hpp"
#include "fragment_size.hpp"
#include "pair_junctions.hpp"
#include "read_pair.hpp"
#include "gtf.hpp"
#include "intervals.hpp"
#include "stranded.hpp"
#include "report.hpp"
#include "resolve_fragments.hpp"
#include "splice_trim.hpp"
#include "pair_output.hpp"

namespace po = boost::program_options;

using namespace std;
using namespace rnasequel;

void merge_init_options(int argc, char *argv[], po::variables_map & vm) {
    po::options_description generic("Arguments");
    generic.add_options()
    ("ref,r", po::value< string >(), "The indexed reference prefix")
    ("gtf,g", po::value< string >(), "GTF file (optional)")
    ("fragments,f", po::value<string>(), "Transcriptome index prefix")
    ("output,o", po::value<string>(), "Output Prefix")
    ("threads,t", po::value< unsigned int >()->default_value(4), "Number of threads to use for processing")
    ("help,h", "help message")
    ;

    po::options_description fragments("Fragment Size Distribution Estimation");
    fragments.add_options()
    ("min-exon", po::value< unsigned int >()->default_value(250), "The minimum exon size to use when estimating the fragment size distribution")
    ("confidence,c", po::value< double >()->default_value(0.99), "The confidence interval for the fragment size distribution (from the left tail to this)")
    ("max-fragment,m", po::value<unsigned int>()->default_value(1500), "The maximum fragment size")
    ("max-dist", po::value<int>()->default_value(1000000), "The maximum distance between two pairs for examination")
    ("first-strand", "The library is first stranded")
    ("second-strand", "The library is first stranded")
    ("obs", po::value<int>()->default_value(1000000), "Number of fragment size observations to use, -1 to check use all pairs")
    ("min-obs", po::value<unsigned int>()->default_value(100000), "Minimum number of fragment size observations required")
    ("max-repeat", po::value<unsigned int>()->default_value(10), "Maximum number of repeat pairs")
    ("frag-sizes,F", po::value< string >(), "Previously estimated fragment sizes")
    //("min-unique,u", po::value< unsigned int >()->default_value(3), "Minimum number of unique reads mapping across the junction")
    //("end-prop,e", po::value< double >()->default_value(0.85), "Proportion of reads that fall within the end region of a read (defined when junctions are extracted)")
    ;
    po::options_description alignment("Alignment Filtering");
    alignment.add_options()
    ("canonical-motifs", po::value< string>()->default_value("GTAG,GCAG,GCTG,GCAA,GCGG,GTTG,GTAA,ATAC,ATAA,ATAG"), "Canonical Splice Sites (Plus strand orientation)")
    ("match,m", po::value<int>()->default_value(3), "Match score")
    ("mismatch,i", po::value<int>()->default_value(-3), "Mismatch score")
    ("gap-open,G", po::value<int>()->default_value(-8), "Gap open penalty")
    ("gap-ext,E", po::value<int>()->default_value(-1), "Gap extend")
    ("score-GTAG", po::value<int>()->default_value(-3), "Score penalty for GT-AG splice junctions")
    ("score-canonical", po::value<int>()->default_value(-6), "Score penalty for canonical junctions (see canonical-motifs option)")
    ("score-non-canonical", po::value<int>()->default_value(-9))
    ("score-diff", po::value<unsigned int>()->default_value(6), "Maximum score difference for singleton alignment (paired alignments can differ by as much as 2 x score-diff")
    ("min-length", po::value<unsigned int>()->default_value(8), "Minimum exonic alignment length")
    ("intron-trim", po::value<unsigned int>()->default_value(6), "Number of base pairs to trim off of alignments that overlap splice sites")
    ("big-intron-size", po::value<unsigned int>()->default_value(64000), "Minimum intron size to include an intron penalty in the alignment score")
    ("big-intron-penalty", po::value<int>()->default_value(-12), "Scaling factor, intron penalty is -1 * (log2(isize) + penalty)")
    ("min-score", po::value< double >()->default_value(2.0), "Cutoff factor for pre-filtering")
    ("max-discordant-dist", po::value<int>()->default_value(10000), "For uniquely mapped discordant pairs output them if they are less than this distance apart")
    ("max-gene-dist", po::value<int>()->default_value(0), "If the program fails to estimate the fragment size and the pair falls within a gene consider it concordant if their distance is less than this")
    ("max-fallback-dist", po::value<int>()->default_value(0), "If a genome GTF file is not specified a pair is considered proper if their distance is less than this apart")
    ("max-edit-dist", po::value<int>()->default_value(15), "Maximum edit distance for a single read")
    ("score-bonus", po::value<int>()->default_value(10), "Bonus score scaled on the maximum value of the normalized empircal fragment size distribution (ie. max density gets this score everything else lower)")
    ;

    po::options_description hidden("Hidden Options");
    hidden.add_options()
    ("refs1", po::value< string >(), "Read 1 Reference Alignment Bam File")
    ("juncs1", po::value< string >(), "Read 1 Junction Alignment Bam File")
    ("refs2", po::value< string >(), "Read 2 Reference Alignment Bam File")
    ("juncs2", po::value< string >(), "Read 2 Junction Alignment Bam File")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(fragments).add(alignment).add(hidden);

    po::options_description visible;
    visible.add(generic).add(alignment).add(fragments);

    po::positional_options_description pd;
    pd.add("refs1",1);
    pd.add("juncs1",1);
    pd.add("refs2",1);
    pd.add("juncs2",1);

    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pd).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << "Usage: " << endl;
        cout << "rnasequel " << string(argv[0]) << " [options] <ref1.bam> <ref2.bam> <juncs1.bam> <juncs2.bam>\n" << endl;
        cout << visible << "\n";
        exit(0);
    }

    bool error = false;

    if(vm.count("output") == 0) {
        cout << "An output prefix must be specified\n";
        error = true;
    }

    if(vm.count("refs1") == 0) {
        cout << "A reference alignment bam file for read 1 must be specified\n";
        error = true;
    }

    if(vm.count("juncs1") == 0){
        cout << "A transcriptome alignment bam file for read 1 be specified\n";
        error = true;
    }

    if(vm.count("refs2") == 0) {
        cout << "A reference alignment bam file for read 2 must be specified\n";
        error = true;
    }

    if(vm.count("juncs2") == 0){
        cout << "A transcriptome alignment bam file for read 2 be specified\n";
        error = true;
    }

    if(vm.count("ref") == 0) {
        cout << "The indexed reference genome prefix must be specified\n";
        error = true;
    }

    if(vm.count("fragments") == 0) {
        cout << "The fragment index must be specified\n";
        error = true;
    }

    if(error) exit(1);
}

int rnasequel::merge_alignments(int argc, char *argv[]) {
    po::variables_map vm;
    merge_init_options(argc,argv,vm);

    bool debug = false;

    FastaIndex fi(vm["ref"].as<string>());
    fi.load_all();

    ResolveFragments rf;
    SpliceTrimmer strimmer(vm["intron-trim"].as<unsigned int>());
    Model model;
    int fb_dist = vm["max-fallback-dist"].as<int>();
    if(vm.count("gtf") > 0){
        GTF::parse(vm["gtf"].as<string>(), model);
        fb_dist = 0;
        model.sort_genes();
    }
    Stranded stranded(vm.count("first-strand") ? Stranded::FIRST_STRAND : (vm.count("second-strand") ? Stranded::SECOND_STRAND : Stranded::UNSTRANDED));
    GeneIntervals gene_intervals;
    gene_intervals.build(model);
    SizeDist      size_dist(vm["confidence"].as<double>(), vm["max-fragment"].as<unsigned int>());
    PairJunctions pjuncs(size_dist, vm["min-length"].as<unsigned int>());
    EstimateDist estimate_dist(model, vm["min-exon"].as<unsigned int>());

    {

        const size_t STEP = 5000;
        PairedReader reader(vm["refs1"].as<string>(), vm["juncs1"].as<string>(), 
                            vm["refs2"].as<string>(), vm["juncs2"].as<string>(), STEP);

        rf.open(reader.tx1_header(), reader.ref1_header(), vm["fragments"].as<string>());
        {
            Timer ti("Building the splice junction maps");
            pjuncs.set_model(model);

            for(auto const & m : rf.fragment_map()){
                for(auto const & s : m.second){
                    for(size_t i = 1; i < s.size(); i++){
                        pjuncs.add_junction(m.first, PosBlock(s[i - 1].rgt, s[i].lft, s.strand()));
                        strimmer.add_junction(m.first, s[i - 1].rgt, s[i].lft);
                    }
                }
            }

            pjuncs.prepare();
            strimmer.merge();
        }
        //FragmentSize fragment_size(pjuncs, estimate_dist, size_dist, stranded, gene_intervals, 
        //        vm["max-gene-dist"].as<int>(), vm["max-dist"].as<int>(), fb_dist, vm["score-bonus"].as<int>());



        //("score-GTAG", po::value<int>()->default_value(-3), "Score penalty for GT-AG splice junctions")
        //("score-canonical", po::value<int>()->default_value(-6), "Score penalty for canonical junctions (see canonical-motifs option)")
        //("score-non-canonical", po::value<int>()->default_value(-12))
        if(vm.count("frag-sizes") == 0){
            std::vector<PairEstimator*> threads(vm["threads"].as<unsigned int>());
            for(size_t i = 0; i < threads.size(); i++){
                threads[i] = new PairEstimator;
                threads[i]->set_debug(debug);
                threads[i]->set_params(vm["min-score"].as<double>(), vm["min-length"].as<unsigned int>());
                threads[i]->score_filter.set_scores(vm["match"].as<int>(), vm["mismatch"].as<int>(),
                                           vm["gap-open"].as<int>(), vm["gap-ext"].as<int>(), vm["canonical-motifs"].as<string>(),
                                           vm["score-GTAG"].as<int>(), vm["score-canonical"].as<int>(), 
                                           vm["score-non-canonical"].as<int>());
                threads[i]->score_filter.set_seqs(fi, reader.ref1_header());
                threads[i]->trimmer = &strimmer;
                threads[i]->pair_factory.set_stranded(stranded);
                threads[i]->rf = &rf;
                threads[i]->score_filter.set_intron_penalties(vm["big-intron-size"].as<unsigned int>(), vm["big-intron-penalty"].as<int>());
                threads[i]->score_diff = vm["score-diff"].as<unsigned int>();
                threads[i]->score_filter.set_filtering_params(vm["min-score"].as<double>(), vm["score-diff"].as<unsigned int>() * 2, vm["max-edit-dist"].as<int>());
                threads[i]->dist.init(vm["confidence"].as<double>(), vm["max-fragment"].as<unsigned int>());
                threads[i]->fsize.init(pjuncs, estimate_dist, threads[i]->dist, stranded, gene_intervals, 
                                        vm["max-gene-dist"].as<int>(), vm["max-dist"].as<int>(), fb_dist, vm["score-bonus"].as<int>());
            }
            size_t total = 0;
            Timer ti("Total Time Estimating Fragment Sizes");
            int min_obs = vm["obs"].as<int>();
            while(reader.load_input() > 0){
                size_t start   = 0;
                size_t num_per = reader.count() / threads.size();
                size_t extra   = reader.count() - num_per * threads.size();
                for(size_t i = 0; i < threads.size(); i++){
                    size_t end = start + num_per;
                    if(extra > 0){
                        end++;
                        extra--;
                    }
                    //cout << "i: " << i << " Start: " << start << " end: " << end << " sz = " << reader.count() << "\n";
                    threads[i]->init(reader.input.begin() + start, reader.input.begin() + end);
                    start = end;
                }

                for(size_t i = 1; i < threads.size(); i++){
                    threads[i]->start();
                }
                threads[0]->operator()();
                //cout << "0: " << threads[0]->num_passed << " " << threads[0]->num_used << "\n";

                int obs = threads[0]->num_passed;
                size_t unq = threads[0]->num_used;
                for(size_t i = 1; i < threads.size(); i++){
                    threads[i]->join();
                    unq += threads[i]->num_used;
                    obs += threads[i]->num_passed;
                    //cout << i << ": " << threads[i]->num_passed << " " << threads[i]->num_used << "\n";
                }

                if(min_obs > 0 && obs >= min_obs){
                    break;
                }

                total = reader.total();
                if(total % 1000000 == 0){
                    time_t e = ti.elapsed();
                    if(e > 0){
                        cout  << "\33[2K\rObservations: " << obs << " [" << std::setprecision(2) << std::fixed << (100.0 * obs / min_obs) << "% complete] Total processed: " << total << " " << (total / e) << " pair groups / second";
                        cout.flush();
                    }
                }
            }
            cout << "\n\nTotal: " << reader.total() << "\n";
            size_t obs = 0;
            for(size_t i = 0; i < threads.size(); i++){
                obs += threads[i]->num_passed;
                size_dist += threads[i]->dist;
                delete threads[i];
            }
            if(obs < vm["min-obs"].as<unsigned int>()){
                cout << "Error only observed " << obs << " observations\n";
                return 1;
            }
            size_dist.normalize();
            {
                string fname = vm["output"].as<string>() + "-dist.txt";
                ofstream out(fname.c_str());
                size_dist.save_normed(out);
            }
            cout << "Used: " << size_dist.count() << " observations to calculate the fragment size distribution\n\n";
        }else{
            string fdist = vm["frag-sizes"].as<string>();
            ifstream in(fdist.c_str());
            size_dist.load(in);
            in.close();
            cout << "Loaded: " << size_dist.count() << " observations to calculate the fragment size distribution\n\n";
        }
        cout << "\n\n";
    }

    {
        const size_t STEP = 10000;
        PairedReader reader(vm["refs1"].as<string>(), vm["juncs1"].as<string>(), 
                            vm["refs2"].as<string>(), vm["juncs2"].as<string>(), STEP);

        size_t N = max(vm["threads"].as<unsigned int>() - 1, 1U);
        PairOutput output_worker(vm["output"].as<string>() + ".bam", reader.ref1_header(), N);
        std::vector<PairResolver*> threads(N);
        for(size_t i = 0; i < threads.size(); i++){
            threads[i] = new PairResolver;
            threads[i]->set_debug(debug);
            threads[i]->set_params(vm["min-score"].as<double>(), vm["min-length"].as<unsigned int>());
            threads[i]->score_filter.set_scores(vm["match"].as<int>(), vm["mismatch"].as<int>(),
                                       vm["gap-open"].as<int>(), vm["gap-ext"].as<int>(), vm["canonical-motifs"].as<string>(),
                                       vm["score-GTAG"].as<int>(), vm["score-canonical"].as<int>(), 
                                       vm["score-non-canonical"].as<int>());
            threads[i]->score_filter.set_intron_penalties(vm["big-intron-size"].as<unsigned int>(), vm["big-intron-penalty"].as<int>());
            threads[i]->score_filter.set_seqs(fi, reader.ref1_header());
            threads[i]->score_filter.set_filtering_params(vm["min-score"].as<double>(), vm["score-diff"].as<unsigned int>() * 2, vm["max-edit-dist"].as<int>());
            threads[i]->trimmer = &strimmer;
            threads[i]->pair_factory.set_stranded(stranded);
            threads[i]->rf = &rf;
            threads[i]->score_diff = vm["score-diff"].as<unsigned int>();
            threads[i]->max_repeat = vm["max-repeat"].as<unsigned int>();
            threads[i]->max_dist   = vm["max-discordant-dist"].as<int>();
            //threads[i]->dist.init(vm["confidence"].as<double>(), vm["max-fragment"].as<unsigned int>());
            threads[i]->fsize.init(pjuncs, estimate_dist, size_dist, stranded, gene_intervals, 
                                    vm["max-gene-dist"].as<int>(), vm["max-dist"].as<int>(), fb_dist, vm["score-bonus"].as<int>());
        }
        if(debug) threads[0]->set_debug(true);
        size_t total = 0;
        Timer ti("Total Time Merging Pairs");
        while(reader.load_input() > 0){
            size_t start   = 0;
            size_t num_per = reader.count() / threads.size();
            size_t extra   = reader.count() - num_per * threads.size();
            for(size_t i = 0; i < threads.size(); i++){
                size_t end = start + num_per;
                if(extra > 0){
                    end++;
                    extra--;
                }
                //cout << "i: " << i << " Start: " << start << " end: " << end << " sz = " << reader.count() << "\n";
                threads[i]->init(reader.input.begin() + start, reader.input.begin() + end);
                start = end;
            }

            for(size_t i = 1; i < threads.size(); i++){
                threads[i]->start();
            }
            threads[0]->operator()();
            size_t unique = threads[0]->counts.unique_pair;
            size_t totalp = threads[0]->counts.total;
            //cout << "0: " << threads[0]->num_passed << " " << threads[0]->num_used << "\n";

            for(size_t i = 1; i < threads.size(); i++){
                threads[i]->join();
                unique += threads[i]->counts.unique_pair;
                totalp += threads[i]->counts.total;
                //cout << i << ": " << threads[i]->num_passed << " " << threads[i]->num_used << "\n";
            }

            // Make sure we are done writing this round of output
            output_worker.join();
            for(size_t i = 0; i < threads.size(); i++){
                output_worker.pools[i].clear();
                output_worker.pools[i].swap(threads[i]->output);
            }

            output_worker.start();
            total = reader.total();
            if(total % 1000000 == 0){
                time_t e = ti.elapsed();
                if(e > 0){
                    cout  << "\33[2K\r Unique = " << unique << " [" << std::fixed << std::setprecision(2) << (100.0 * unique / totalp) << "] Total Test = " << totalp 
                          << " Total processed: " << total << " " << (total / e) << " pair groups / second";
                    cout.flush();
                }
            }
        }
        // Write out remaining reads
        output_worker.start();
        output_worker.join();
        output_worker.close();
        PairResolver::OutputCounts counts;
        for(size_t i = 0; i < threads.size(); i++){
            counts += threads[i]->counts;
            delete threads[i];
        }
        std::cout << "\n\n";
        Report<unsigned int> rep(cout);
        rep
            ("Unique Pairs", counts.unique_pair)
            ("Repeat Pairs", counts.repeat_pair)
            ("Filtered Repeat Pairs", counts.xrepeat_pair)
            ("R1 Singleton", counts.r1_single)
            ("R1 Repeat", counts.r1_repeat)
            ("R1 Filtered Repeat", counts.r1_xrepeat)
            ("R2 Singleton", counts.r2_single)
            ("R2 Repeat", counts.r2_repeat)
            ("R2 Filtered Repeat", counts.r2_xrepeat)
            ("Discordant Pairs", counts.discordant_pair)
            ("Discordant Unique", counts.discordant_single)
            ("Unaligned", counts.unaligned)
            .report();
    }

    return 0;
}
