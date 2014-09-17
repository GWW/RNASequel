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

#include "gtf.hpp"
#include <fstream>
#include "tokenizer.hpp"
#include "timer.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

using namespace std;
using namespace rnasequel;

std::string skip_sources[] = {};

typedef boost::unordered_map< std::string, size_t >      GeneMap;
typedef boost::unordered_map< std::string, size_t >      TranscriptMap;
typedef GeneMap::iterator g_iterator;
typedef TranscriptMap::iterator tx_iterator;

void GTF::parse(const std::string & file, Model & m) {
    boost::iostreams::filtering_istream in;
    std::ifstream                       ifs;
    bool gz = file.rfind(".gz") == (file.length() - 3);
    ifs.open(file.c_str(), std::ios_base::in | std::ios_base::binary);

    if(!ifs) {
        std::cout << "Error opening the gtf file `" << file << "` for reading\n";
        exit(1);
    }

    if(gz) {
        in.push(boost::iostreams::gzip_decompressor());
    }
    in.push(ifs);

    Timer ti("Total GTF parsing time");

    {

    Timer ti2("Time parsing the GTF file");

    GeneMap       gmap;
    TranscriptMap tmap;
    size_t gene_ptr = 0, tx_ptr = 0;

    std::string line, gene_id, tx_id, chrom_id, last_gene, last_tx;
    std::vector<char *> tokens;
    std::vector< std::pair<char *, char *> > tag_tokens;
    Model::iterator chrom_it;
    unsigned int line_no = 0;


    while(getline(in, line)) {
	if(line[0] == '#' || line[0] == 'G' || line[0] == 'H'){
	    line_no++;
	    continue;
	}

	Tokenizer::get(line, '\t', tokens);

	if(tokens.size() < 9 || tokens.size() > 10){
	    std::cout << "Error malformed GTF file: " << file << " at line: " << line_no << " number of tokens = " << tokens.size() << "\n";
	    exit(1);
	}

	if(strcmp(tokens[2], "exon") != 0 
		&& strcmp(tokens[2], "stop_codon") != 0
		&& strcmp(tokens[2], "start_codon") != 0) continue;

	char * tags = tokens[8];

	char * gid = NULL;
	char * tid = NULL;
	char * pid = NULL;
	char * gname = NULL;

        //std::cout << "  tags = '" << tags << "'\n";
	//char * tname = NULL;

	while(*tags != '\0'){
            if(*tags == '\0') break;
	    while(*tags == ' ') tags++;
            if(*tags == '\0') break;
	    char * keyp = tags;
	    while(*(++tags) != ' ');
	    *tags = '\0';
	    tags += 2;
	    char * valp = tags;
	    while(*(++tags) != '"');
	    *tags = '\0';
	    tags++;
	    if(*tags != ';'){
		std::cout << "Error malformed GTF file: " << file << " at line: " << line_no << " keyp = " << keyp << " valp = " << valp << "\n";
		exit(1);
	    }

	    if(strcmp("gene_name", keyp) == 0){
		gname = valp;
	    }else if(strcmp("gene_id", keyp) == 0){
		gid = valp;
	    }else if(strcmp("transcript_id", keyp) == 0){
		tid = valp;
	    //}else if(strcmp("transcript_name", keyp) == 0){
		//tname = valp;
	    }else if(strcmp("protein_id", keyp) == 0){
		pid = valp;
	    }

	    tags++;
	}

	if(gid == NULL || tid == NULL){
	    std::cout << "Error malformed GTF file: " << file << " at line: " << line_no << "\n";
	    exit(1);
	}

	chrom_id = tokens[0];
	chrom_it = m.chrom(chrom_id);
	gene_id  = gid;
	tx_id    = tid;

	if(last_gene != gene_id){
	    pair<g_iterator, bool> res = gmap.insert(GeneMap::value_type(gene_id, 0));
	    if(res.second){
		res.first->second = chrom_it->second.size();
		chrom_it->second.push_back(Gene());
		Gene & gene = chrom_it->second.back();
		gene.id()        = gid;
		gene.name()      = gname == NULL ? "" : gname;
		gene.strand()    = char2strand[static_cast<size_t>(*tokens[6])];
                gene.ref()       = chrom_id;
	    }
	    last_gene = gene_id;
	    gene_ptr = res.first->second;
	}

	Gene & gene = chrom_it->second[gene_ptr];

	if(last_tx != tx_id){
	    pair<tx_iterator, bool> res = tmap.insert(TranscriptMap::value_type(tx_id, 0));
	    if(res.second){
		res.first->second = gene.transcripts().size();
		gene.transcripts().push_back(Transcript());
		Transcript & tx = gene.transcripts().back();
		tx.source() = tokens[1];
		tx.strand() = char2strand[static_cast<size_t>(*tokens[6])];
		tx.id()     = tid;
		if(pid != NULL) tx.pid() = pid;
	    }
	    last_tx = tx_id; 
	    tx_ptr = res.first->second;
	}

	Transcript & tx = gene.transcripts()[tx_ptr];

	if(strcmp(tokens[2], "exon") == 0){
	    tx.add_exon(atoi(tokens[3]) - 1, atoi(tokens[4]) - 1);
	}else if(strcmp(tokens[2], "start_codon") == 0){
	    tx.set_coding();
	    if(tx.strand() == PLUS){
		tx.cds_start() = atoi(tokens[3]) - 1;
	    }else{
		tx.cds_end()   = atoi(tokens[4]) - 1;
	    }
	}else if(strcmp(tokens[2], "stop_codon") == 0){
	    tx.set_coding();
	    if(tx.strand() == PLUS){
		tx.cds_end()   = atoi(tokens[4]) - 1;
	    }else{
		tx.cds_start() = atoi(tokens[3]) - 1;
	    }
	}

	line_no++;

    }

    }

    {

    Timer ti2("Time sorting the exons / transcripts / genes");
    for(Model::iterator it = m.begin(); it != m.end(); ++it){
	for(size_t i = 0; i < it->second.size(); ++i){
	    Gene & g = it->second[i];
	    //cout << g.id() << " size: " << g.transcripts().size() << "\n";
	    for(Gene::iterator tx_it = g.begin(); tx_it != g.end(); ++tx_it){
		//cout << "    " << tx_it->id() << " size: " << tx_it->exons().size() << "\n";
		g.lft() = std::min(g.lft(), tx_it->lft());    
		g.rgt() = std::max(g.rgt(), tx_it->rgt());    
		if(tx_it->cds_start() == MAX_POS) tx_it->cds_start() = 0;
		if(tx_it->cds_end() == MIN_POS)   tx_it->cds_end() = 0;
		tx_it->sort_exons();
	    }
	    g.sort_transcripts();
	}
    }

    m.sort_genes();

    }

    /*
    for(Model::iterator it = m.begin(); it != m.end(); ++it){
	for(size_t i = 0; i < it->second.size(); ++i){
	    Gene & g = it->second[i];
            cout << g << "\n";
        }
    }
    */
}
