#include <iostream>
#include <cmath>
#include <fstream>
#include <set>
#include <map>
#include <algorithm>
#include <cstring>

using namespace std;

#include <tktools.hxx>
#include <seq_gene.hxx>

using namespace tktools;
using namespace tktools::io;
using namespace tktools::util;
using namespace tktools::bio;
using namespace tkbio;

//double twocolor_probe::_inv_log_2 = 1.0 / log(2.0);

//const float twocolor_probe::NO_DATA = - FLT_MAX;

int seq_gene::maximum_size = 0;

seq_gene::seq_gene(const string& name, const string& symbol, int geneid) {
    _id = name;
    _symbol = symbol;
    _geneid = geneid;
    _chrm = -1;
    _pseudo = false;
}

void seq_gene::output(ostream& ost) const {
    ost << "GeneID:" << _geneid << "\t" << _symbol << "\t" << chromosome() << ":" << _orientation << ":" << _start << "-" << _end;
}

seq_gene::~seq_gene() {
    // for (vector<exon*>::iterator it = _exons.begin(); it != _exons.end(); it++) {
    //     delete *it;
    // }
}

bool seq_gene::contains(int chromosome_code, int start, int end) const {
    int left = 0;
    int right = _exons.size();
    if (chromosome_code != _chrm || start > _end || end < _start) {
        return false;
    }
    //cout << start << "-" << end << " // " << _start << "-" << _end << endl;
    for (;;) {
        if (left == right) {
            return false;
        }
        int index = (left + right) / 2;
        const exon& e = _exons[index];
        //cout << index << "/" << _exons.size() << " : " << e._start << "-" << e._end << endl;
        if (e._start > end) {
            right = index;
        } else if (e._end < start) {
            left = index + 1;
        } else {
            //if (e._start <= end && e._end >= start) {
            return true;
        }
    }
}

void seq_gene::add_exon(exon::Feature feature, int start, int end) {
    _exons.push_back(exon(feature, start, end));
}

void seq_gene::set_location(int chrm, char orientation, int start, int end) {
    _chrm = chrm;
    _orientation = orientation;
    _start = start;
    _end = end;
//cout << _chrm << " : " << chromosome() << endl;
}

void seq_gene::set_location(const char* chrm, char orientation, int start, int end) {
if (strncmp(chrm, "chr", 3) == 0) {
        chrm += 3;
    }
//    _chrm = convert_chromosome_to_code(chrm);
set_location(convert_chromosome_to_code(chrm), orientation, start, end);
//    _orientation = orientation;
//    _start = start;
//    _end = end;
}


int seq_gene::get_distance(int chromosome, int position) const {
    if (chromosome != _chrm) {
        return 0x8ffffff;
    }
    if (_orientation == '+') {
        return position - _start;
    } else {
        return _end - position;
    }
}

void seq_gene::set_symbol(const string& symbol) {
    _symbol = symbol;
}

void seq_gene::serialize(ostream& ost) const {
    ost << genbank() << "\t" << symbol() << "\t" << geneid() << "\t";
    ost << chromosome() << "\t" << orientation() << "\t" << start() << "\t" << end();
    ost <<  "\t" << (is_pseudo() ? "pseudo" : "coding");
    if (_exons.size() > 0) {
        ost << "\t";
        for (int i = 0; i < (int)_exons.size(); i++) {
            if (i > 0) ost << ",";
            char feature;
            switch (_exons[i]._feature) {
            case exon::UTR:
                feature = 'U'; break;
            case exon::CDS:
                feature = 'C'; break;
            default:
                feature = 'x'; break;
            }
            ost << feature << ":" << _exons[i]._start << ":" << _exons[i]._end;
        }
    }
    ost << "\n";

}

bool seq_gene::compare_position(const seq_gene* lhs, const seq_gene* rhs) {
    if (lhs->_chrm != rhs->_chrm) {
        return lhs->_chrm < rhs->_chrm;
    } else {
        return lhs->_start < rhs->_start;
//        return (lhs->_start + lhs->_end) < (rhs->_start + rhs->_end);
        // int pl = lhs->_orientation == '+' ? lhs->_start : lhs->_end;
        // int pr = rhs->_orientation == '+' ? rhs->_start : rhs->_end;
        // return pl < pr;
    }
}

bool seq_gene::compare_tss_position(const seq_gene* lhs, const seq_gene* rhs) {
    if (lhs->_chrm != rhs->_chrm) {
        return lhs->_chrm < rhs->_chrm;
    } else {
        int pl = lhs->_orientation == '+' ? lhs->_start : lhs->_end;
        int pr = rhs->_orientation == '+' ? rhs->_start : rhs->_end;
        return pl < pr;
    }
}

namespace {
    const char* cache_file = ".seq_gene";
    template<typename T> T MAX(T a, T b) { return a > b ? a : b; }

    string get_cache_filename(const char* original_file, const char* assembly) {
        string fn = string(cache_file) + ".";
        for (int i = 0, len = strlen(assembly); i < len; i++) {
            char c = assembly[i];
            if (c >= 'A' && c <= 'Z') {
                c += (char)('a' - 'A');
            } else if ((c >= 'a' && c <= 'z') ||(c >= '0' && c <= '9')) {
            } else {
                continue;
            }
            fn += c;
        }
        return fn;
    }

    vector<seq_gene*> load_cache(const char* filename, const char* original_file=NULL, const char* assembly=NULL) throw (logic_error) {
        vector<seq_gene*> genes;
        ifstream file_in(filename);
        //cerr << filename << endl;
        if (file_in.is_open() == false) {
            //throw invalid_argument(string("cannot open ") + string(filename));
            return genes;
        }
        while (!file_in.eof()) {
            string line;
            getline(file_in, line);
            //cout << line << endl;
            if (line.c_str()[0] == '#') {
                vector<string> prop = split_items(line, '=');
                if (prop[0] == "#original" && original_file != NULL) {
                    if (prop[1] != original_file) {
                        throw invalid_argument(string("filename mismatches : ") + prop[1] + string(" != ") + original_file);
                    }
                }
                if (prop[0] == "#assembly" && assembly != NULL) {
                    if (prop[1] != assembly) {
                        throw invalid_argument(string("assembly mismatches : ") + prop[1] + string(" != ") + assembly);
                    }
                }
            } else {
                vector<string> items = split_items(line, '\t');
                if (items.size() > 6) {
                    seq_gene* g = new seq_gene(items[0], items[1], atoi(items[2].c_str()));
                    g->set_location(items[3].c_str(), items[4].c_str()[0], atoi(items[5].c_str()), atoi(items[6].c_str()));
                    g->set_pseudo(items.size() > 7 && items[7] == "pseudo");
                    if (items.size() > 8) {
                        vector<string> exons = split_items(items[8], ',');
                        for (int i = 0; i < (int)exons.size(); i++) {
                            vector<string> exon = split_items(exons[i], ':');
                            if (exon.size() == 3) {
                                exon::Feature feature;
                                if (exon[0] == "U") {
                                    feature = exon::UTR;
                                } else if (exon[0] == "C") {
                                    feature = exon::CDS;
                                } else if (exon[0] == "x") {
                                    feature = exon::OTHERS;
                                } else {
                                    throw invalid_argument(string("unkown exon data ") + exons[i]);
                                }
                                int start = atoi(exon[1].c_str());
                                int end = atoi(exon[2].c_str());
                                g->add_exon(feature, start, end);
                            }
                        }
                    }
                    genes.push_back(g);
                }
            }
        }
        file_in.close();
        return genes;
    }

    void save_cache(const char* filename, const vector<seq_gene*>& genes, const char* original_file=NULL, const char* assembly=NULL) throw (logic_error) {
        ofstream file_out(filename);
        if (file_out.is_open() == false) {
            throw invalid_argument(string("cannot open ") + string(filename));
        }
        if (original_file != NULL) {
            file_out << "#original=" << original_file << "\n";
        }
        if (assembly != NULL) file_out << "#assembly=" << assembly << "\n";
        for (vector<seq_gene*>::const_iterator it = genes.begin(); it != genes.end(); it++) {
            const seq_gene* g = *it;
            g->serialize(file_out);
            //g->serialize(cerr);//file_out);
            
        }
        file_out.close();
    }


    ////////////////////////////
    /////////// spliced_gene
    vector<spliced_gene*> load_cache_sp(const char* filename_cache, const char* original_file, const char* assembly=NULL) throw (logic_error) {
        vector<spliced_gene*> genes;
        ifstream file_in(filename_cache);
        if (file_in.is_open() == false) {
            return genes;
        }
        while (!file_in.eof()) {
            string line;
            getline(file_in, line);
            if (line.c_str()[0] == '#') {
                vector<string> prop = split_items(line, '=');
                if (prop[0] == "#original" && original_file != NULL) {
                    if (prop[1] != original_file) {
                        throw invalid_argument(string("filename mismatches : ") + prop[1] + string(" != ") + original_file);
                    }
                }
                if (prop[0] == "#assembly" && assembly != NULL) {
                    if (prop[1] != assembly) {
                        throw invalid_argument(string("assembly mismatches : ") + prop[1] + string(" != ") + assembly);
                    }
                }
            } else {
                vector<string> items = split_items(line, '\t');
                if (items.size() > 7) {
                    spliced_gene* g = new spliced_gene(items[0], items[1], atoi(items[2].c_str()));
                    g->set_location(items[3].c_str(), items[4].c_str()[0], atoi(items[5].c_str()), atoi(items[6].c_str()));
                    vector<string> sec = split_items(items[7], ',');
                    for (int i = 0; i < (int)sec.size(); i++) {
                        //cout << sec[i] << endl;
                        vector<string> elem = split_items(sec[i], ':');
                        if (elem.size() >= 3) {
                            g->set_section(atoi(elem[1].c_str()), atoi(elem[2].c_str()), elem[0].c_str());
                        }
                    }
                    g->set_pseudo(items.size() > 8 && items[8] == "pseudo");
                    genes.push_back(g);
                }
            }
        }
        
        return genes;
    }
    

    // void save_cache_sp(const char* filename, const vector<spliced_gene*>& genes, const char* original_file=NULL, const char* assembly=NULL) throw (logic_error) {
    //     ofstream file_out(filename);
    //     if (file_out.is_open() == false) {
    //         throw invalid_argument(string("cannot open ") + string(filename));
    //     }
    //     if (original_file != NULL) {
    //         file_out << "#original=" << original_file << "\n";
    //     }
    //     if (assembly != NULL) file_out << "#assembly=" << assembly << "\n";
    //     for (vector<spliced_gene*>::const_iterator it = genes.begin(); it != genes.end(); it++) {
    //         const spliced_gene* g = *it;
    //         file_out << g->genbank() << "\t" << g->symbol() << "\t" << g->geneid() << "\t";
    //         file_out << g->chromosome() << "\t" << g->orientation() << "\t" << g->start() << "\t" << g->end() << "\t";
    //         for (int i = 0; i < g->count(); i++) {
    //             const splice_section& ss = g->get_feature(i);
    //             if (i > 0) file_out << ",";
    //             file_out << ss.feature() << ":" << ss.start() << ":" << ss.end();
    //         }
    //         file_out <<  "\t" << (g->is_pseudo() ? "pseudo" : "coding") << "\n";
    //     }
    //     file_out.close();
    // }

}

string seq_gene::chromosome() const {
    return tktools::bio::convert_code_to_chromosome(_chrm);
}



// void seq_gene::insert_probe(vector<seq_gene*>& genes, const target_probe* probe, int distance_up, int distance_down) {
//     int left = 0;
//     int right = (int)genes.size();
//     int center = 0;
//     int pchr = probe->chromosome_code();
//     int ppos = probe->position();
//     for (;;) {
//         center = (left + right) >> 1;
//         if (center == left) break;
//         seq_gene* gene = genes[center];
//         if (gene->_chrm != pchr) {
//             if (gene->_chrm < pchr) {
//                 left = center + 1;
//             } else {
//                 right = center;
//             }
//         } else {
//             int tss_left, tss_right;
//             if (gene->_orientation == '+') {
//                 tss_left = gene->_start + distance_up;
//                 tss_right = gene->_start + distance_down;
//             } else {
//                 tss_left = gene->_end + distance_down;
//                 tss_right = gene->_end + distance_up;
//             }
//             if (tss_left <= ppos && ppos <= tss_right) {
//                 break;
//             } else if (tss_right < ppos) {
//                 left = center + 1;
//             } else {
//                 right = center;
//             }
//         }
//     }
//     int margin = distance_down - distance_up;
//     int index = center;
//     while (index >= 0) {
//         seq_gene* gene = genes[index];
//         if (pchr == 1) {
//             cout << index << "  " << pchr << ":" << ppos << " " << gene->chromosome() << ":" << gene->_orientation << ":" << gene->_start << ":" << gene->_end << " GeneID:" << gene->_geneid << " " << gene->_symbol << endl;
//         }
//         if (gene->_chrm != pchr) break;
//         int tss_left, tss_right;
//         if (gene->_orientation == '+') {
//             tss_left = gene->_start + distance_up;
//             tss_right = gene->_start + distance_down;
//         } else {
//             tss_left = gene->_end + distance_down;
//             tss_right = gene->_end + distance_up;
//         }
//         if (tss_right + margin <= ppos) break;
//         if (tss_left <= ppos && tss_right > ppos) {
//             gene->probes.push_back(probe);
//             if (pchr == 1) {
//                 gene->serialize(cout);
//                 cout << endl;
//             }
//         }
//         index--;
//     }
//     index = center + 1;
//     while (index <= (int)genes.size()) {
//         seq_gene* gene = genes[index];
//         if (gene->_chrm != pchr) break;
//         int tss_left, tss_right;
//         if (gene->_orientation == '+') {
//             tss_left = gene->_start + distance_up;
//             tss_right = gene->_start + distance_down;
//         } else {
//             tss_left = gene->_end + distance_down;
//             tss_right = gene->_end + distance_up;
//         }
//         if (tss_left - margin > ppos) break;
//         if (tss_left <= ppos && tss_right > ppos) {
//             gene->probes.push_back(probe);
            
//             if (pchr == 1) {
//                 gene->serialize(cout);
//                 cout << endl;
//             }
// //            cout << index << "  " << pchr << ":" << ppos << " " << gene->chromosome() << ":" << gene->_orientation << ":" << gene->_start << ":" << gene->_end << " GeneID:" << gene->_geneid << " " << gene->_symbol << endl;

//         }
//         index++;
//     }
// }

vector<seq_gene*> seq_gene::load(const char* filename, const char* assembly, bool verbose) throw (logic_error) {
    ifstream file_in(filename);
    if (file_in.is_open() == false) {
        throw invalid_argument(string("cannot open ") + string(filename));
    }
    if (verbose) cerr << "loading " << filename << endl;
    int col_chr = 1;
    int col_start = 2;
    int col_end = 3;
    int col_ori = 4;
    int col_feature = 9;
    int col_id = 10;
    int col_type = 11;
    int col_group = 12;
    int col_max = 12;
    int col_tr = 13;
    string filename_cache = get_cache_filename(filename, assembly);
    if (verbose) cerr << "trying to load " << filename_cache << endl;
    vector<seq_gene*> genes = load_cache(filename_cache.c_str(), filename, assembly);
    for (vector<seq_gene*>::const_iterator it = genes.begin(); it != genes.end(); it++) {
        maximum_size = MAX(maximum_size, (*it)->get_size());
    }
//    cerr << genes.size() << endl;
    if (genes.size() > 100) {
        if (verbose) cerr << "load from cache\n";
        return genes;
    }
//    map<string,pair<int,string> > rna2gene;
    map<int,string> geneid2symbol;
    set<int> pseudo;
    int num_lines = 0;
    map<string,vector<exon> > exons;
    while (!file_in.eof()) {
        string line;
//        cerr << line << endl;
        getline(file_in, line);
        num_lines++;
        vector<string> items = split_items(line, '\t');
        if (line.c_str()[0] == '#') {
            for (int i = 0; i < (int)items.size(); i++) {
//                cerr << i << " : " << items[i] << endl;
                if (items[i] == "chromosome") {
                    col_chr = i;
                } else if (items[i] == "chr_start") {
                    col_start = i;
                } else if (items[i] == "chr_end") {
                    col_end = i;
                } else if (items[i] == "feature_name") {
                    col_feature = i;
                } else if (items[i] == "feature_id") {
                    col_id = i;
                } else if (items[i] == "feature_type") {
                    col_type = i;
                } else if (items[i] == "group_label") {
                    col_group = i;
                } else if (items[i] == "transcript") {
                    col_tr = i;
                }
            }
            col_max = MAX(MAX(MAX(MAX(MAX(MAX(col_chr, col_start), col_end), col_feature), col_id), col_type), col_group);
        } else if ((int)items.size() > col_max) {
            int geneid = items[col_id].size() > 7 ? atoi(items[col_id].c_str() + 7) : -1;
            //cerr << items[col_group] << endl;
            if (items[col_type] == "GENE") {
                if (items[col_id].size() > 7) {
                    //int geneid = atoi(items[col_id].c_str() + 7);
                    if (geneid > 0 && geneid2symbol.find(geneid) == geneid2symbol.end()) {
                        geneid2symbol[geneid] = items[col_feature];
                    }
                }
                //genes.push_back(new seq_gene(
            } else if (items[col_type] == "PSEUDO") {
                if (geneid > 0) {//items[col_id].size() > 7) {
//                    int geneid = atoi(items[col_id].c_str() + 7);
//                    if (pseudo.find(geneid) == pseudo.end()) {
                    pseudo.insert(geneid);
                    // if (geneid == 12578) {
                    //     cerr << line << endl;
                    //     throw invalid_argument(line);
                    // }
//                    }
                }
            } else if (items[col_group] == assembly) {
                if (items[col_type] == "UTR" || items[col_type] == "CDS") {
                    const string& refseq = items[col_tr];
                    int es = atoi(items[col_start].c_str());
                    int ee = atoi(items[col_end].c_str());

//                    cout << refseq << " // " << items[col_type] << ":" << es << "-" << ee << endl;

                    map<string,vector<exon> >::iterator eit = exons.find(refseq);
                    exon::Feature feature = items[col_type] == "UTR" ? exon::UTR : exon::CDS;
                    if (eit == exons.end()) {
                        vector<exon> exons_for_rna;
                        exons_for_rna.push_back(exon(feature, es, ee));
                        exons[refseq] = exons_for_rna;
                    } else {
                        eit->second.push_back(exon(feature, es, ee));
                    }
                } else if (items[col_type] == "RNA") {
                    map<int,string>::const_iterator it = geneid2symbol.find(geneid);
                    seq_gene* gene = NULL;
                    if (it == geneid2symbol.end()) {
                        gene = new seq_gene(items[col_feature], "", geneid);
                    } else {
                        gene = new seq_gene(items[col_feature], it->second, geneid);
                    }
                    gene->set_location(items[col_chr].c_str(), items[col_ori].c_str()[0], atoi(items[col_start].c_str()), atoi(items[col_end].c_str()));
                    maximum_size = MAX(maximum_size, gene->get_size());
                    genes.push_back(gene);
                    if (verbose) {
                        //gene->output(cerr);
                        cerr << num_lines << "   " << genes.size() << "   " << geneid2symbol.size() << "  " << items[col_chr] << "   " << gene->symbol() << "                   \r";
                    }
                    //cout << gene->symbol() << "\t" << gene->geneid() << "\t" << gene->genbank() << endl;
                }
            }
        }
    }
    file_in.close();

    if (verbose) cerr << "rearrange genes symbol and set exons\n";
    for (vector<seq_gene*>::iterator it = genes.begin(); it != genes.end(); it++) {
        seq_gene* g = (*it);
        if ((*it)->symbol() == "") {
            int geneid = g->geneid();
            if (geneid > 0 && geneid2symbol.find(geneid) != geneid2symbol.end()) {
                g->set_symbol(geneid2symbol[geneid]);
                if (verbose) cerr << "rename GeneID:" << geneid << " as " << geneid2symbol[geneid] << endl;
            }
        }
        if (pseudo.find(g->geneid()) != pseudo.end()) {
            g->set_pseudo();
        }
        map<string,vector<exon> >::const_iterator sit = exons.find(g->_id);
        if (sit != exons.end()) {
            for (vector<exon>::const_iterator eit = sit->second.begin(); eit != sit->second.end(); eit++) {
                g->add_exon(eit->_feature, eit->_start, eit->_end);
            }
        }
    }
    
    if (verbose) cerr << "saving to cache file\n";
    sort(genes.begin(), genes.end(), seq_gene::compare_tss_position);
    save_cache(filename_cache.c_str(), genes, filename, assembly);
    return genes;
}

int seq_gene::count_exons() const {
    int index = 0;
    if (_orientation == '+') {
        for (int i = 0; i < (int)_exons.size(); i++) {
            if (i == 0 || _exons[i]._start - _exons[i-1]._end > 1) {
                index++;
            }
        }
    } else {
        for (int i = (int)_exons.size() - 1; i >= 0; i--) {
            if (i == (int)_exons.size() - 1 || _exons[i]._start - _exons[i-1]._end > 1) {
                index++;
            }
        }
    }
    return index;
}

int seq_gene::count_bases() const {
    int bases = 0;
    for (int i = (int)_exons.size() - 1; i >= 0; i--) {
        bases += _exons[i].end() - _exons[i].start() + 1;
    }
    return bases;
}

int seq_gene::get_exon_index(int position) const {//, bool utr_sensitive) const {
    int index = 0;
    if (position < _start || position > _end) return -1;
    if (_orientation == '+') {
        for (int i = 0; i < (int)_exons.size(); i++) {
            if (i == 0 || _exons[i]._start - _exons[i-1]._end > 1) {
                index++;
            }
            int es = _exons[i]._start;
            int ee = _exons[i]._end;
//            cout << es << "-" << ee << "  <= " << position << " : " << ((es <= position && position <= ee) ? "inside" : "outside") << endl;
            if (es <= position && position <= ee) {
                return index;
            }
        }
    } else {
        for (int i = (int)_exons.size() - 1; i >= 0; i--) {
            if (i == (int)_exons.size() - 1 || _exons[i]._start - _exons[i-1]._end > 1) {
                index++;
            }
            int es = _exons[i]._start;
            int ee = _exons[i]._end;
//            cout << es << "-" << ee << "  <= " << position << " : " << ((es <= position && position <= ee) ? "inside" : "outside") << endl;
            if ( es <= position && position <= ee) {
                return index;
            }
        }
    }
    return -1;
}

vector<const seq_gene*> seq_gene::find_genes(const vector<seq_gene*>& genes, const string& chromosome, int start, int end) {
    return find_genes(genes, tktools::bio::convert_chromosome_to_code(chromosome.c_str()), start, end);
}

vector<const seq_gene*> seq_gene::find_genes(const vector<seq_gene*>& genes, int chromosome_code, int start, int end) {
    int left = 0;
    int right = (int)genes.size();
    int index = 0;
    vector<const seq_gene*> selected;
    if (right == 0) {
        return selected;
    }
    for (;;) {
        index = (left + right) / 2;
        if (left == right) {
            break;
        }
        const seq_gene* g = genes[index];
        //cout << chromosome_code << ":" << start << " " << g->symbol() << "\t" << g->chromosome() << ":" << g->start() << "-" << g->end() << endl;
        if (g->_chrm < chromosome_code) {
            left = index + 1;
        } else if (g->_chrm > chromosome_code) {
            right = index;
        } else if (g->_start > end) {
            right = index;
        } else if (g->_start + maximum_size < start) {
            left = index + 1;
        } else {
            break;
        }
    }
    for (int i = index; i >= 0; i--) {
        const seq_gene* g = genes[i];
        if (g->_chrm != chromosome_code || g->_end + maximum_size < start) {
            break;
        } else if (g->contains(chromosome_code, start, end)) {
            //} else if (g->_end >= start && g->_start <= end) {
            //cout << *g << endl;
            selected.push_back(g);
        }
    }
    for (int i = index + 1; i < (int)genes.size(); i++) {
        const seq_gene* g = genes[i];
        //cout << *g << endl;
        if (g->_chrm != chromosome_code || g->_start > end) {
            break;
        } else if (g->contains(chromosome_code, start, end)) {//g->_end >= start && g->_start <= end) {
//        } else if (g->_end >= start && g->_start <= end) {
            //cout << *g << endl;
            selected.push_back(g);
        }
    }
    return selected;
}

splice_section splice_section::out_of_range;

splice_section::splice_section(int start, int end, splice_section::section_feature feature) throw (invalid_argument) {
    _start = start;
    _end = end;
    _feature = feature;
    _value = 0;
    if (feature < 1 || feature > 4) {
        throw invalid_argument("given feature is not supported");
    }
}

int spliced_gene::maximum_size = 0;

vector<spliced_gene*> spliced_gene::load(const char* filename, const char* assembly, bool verbose) throw (logic_error) {
    ifstream file_in(filename);
    if (file_in.is_open() == false) {
        throw invalid_argument(string("cannot open ") + string(filename));
    }
    if (verbose) cerr << "loading " << filename << endl;
    int col_chr = 1;
    int col_start = 2;
    int col_end = 3;
    int col_ori = 4;
    int col_feature = 9;
    int col_id = 10;
    int col_type = 11;
    int col_group = 12;
    int col_max = 13;
    int col_tr = 13;
//    int col_tr = 13;
    string filename_cache = get_cache_filename(filename, assembly) + ".sp";
    if (verbose) cerr << "trying to load " << filename_cache << endl;
    vector<spliced_gene*> genes = load_cache_sp(filename_cache.c_str(), filename, assembly);
    for (vector<spliced_gene*>::const_iterator it = genes.begin(); it != genes.end(); it++) {
        maximum_size = MAX(maximum_size, (*it)->get_size());
    }
//    cerr << genes.size() << endl;
    if (genes.size() > 100) {
        if (verbose) cerr << "load from cache\n";
        return genes;
    }
//    map<string,pair<int,string> > rna2gene;
    map<int,string> geneid2symbol;
    map<string,spliced_gene*> parents;
    set<int> pseudo;
    int num_lines = 0;
    while (!file_in.eof()) {
        string line;
//        cerr << line << endl;
        getline(file_in, line);
        num_lines++;
        vector<string> items = split_items(line, '\t');
        if (line.c_str()[0] == '#') {
            for (int i = 0; i < (int)items.size(); i++) {
//                cerr << i << " : " << items[i] << endl;
                if (items[i] == "chromosome") {
                    col_chr = i;
                } else if (items[i] == "chr_start") {
                    col_start = i;
                } else if (items[i] == "chr_end") {
                    col_end = i;
                } else if (items[i] == "feature_name") {
                    col_feature = i;
                } else if (items[i] == "feature_id") {
                    col_id = i;
                } else if (items[i] == "feature_type") {
                    col_type = i;
                } else if (items[i] == "group_label") {
                    col_group = i;
                } else if (items[i] == "transcript") {
                    col_tr = i;
                }
            }
            col_max = MAX(MAX(MAX(MAX(MAX(MAX(MAX(col_chr, col_start), col_end), col_feature), col_id), col_type), col_group), col_tr);
        } else if ((int)items.size() > col_max) {
            int geneid = items[col_id].size() > 7 ? atoi(items[col_id].c_str() + 7) : -1;
            if (items[col_type] == "GENE") {
                if (items[col_id].size() > 7) {
                    if (geneid > 0 && geneid2symbol.find(geneid) == geneid2symbol.end()) {
                        geneid2symbol[geneid] = items[col_feature];
                    }
                }
            } else if (items[col_type] == "PSEUDO") {
                if (geneid > 0) {
                    pseudo.insert(geneid);
                    if (geneid == 12578) {
                        cerr << line << endl;
                        throw invalid_argument(line);
                    }
                }
            } else if (items[col_group] == assembly) {
                int start = atoi(items[col_start].c_str());
                int end = atoi(items[col_end].c_str());
                if (items[col_type] == "RNA") {
                    map<int,string>::const_iterator it = geneid2symbol.find(geneid);
                    spliced_gene* gene = NULL;
                    if (it == geneid2symbol.end()) {
                        gene = new spliced_gene(items[col_feature], "", geneid);
                    } else {
                        gene = new spliced_gene(items[col_feature], it->second, geneid);
                    }
                    gene->set_location(items[col_chr].c_str(), items[col_ori].c_str()[0], start, end);
                    maximum_size = MAX(maximum_size, gene->get_size());
                    genes.push_back(gene);
                    parents[items[col_feature]] = gene;
                    if (verbose) cerr << num_lines << "   " << genes.size() << "   " << geneid2symbol.size() << "  " << pseudo.size() << "   " << gene->symbol() << "                   \r";
                } else if (items[col_type] == "CDS" || items[col_type] == "UTR") {
                    map<string,spliced_gene*>::iterator it = parents.find(items[col_tr]);
                    if (it == parents.end()) {
//                        cerr << "cannot find parent gene for " << items[col_feature] << endl;
                    } else {
                        it->second->set_section(start, end, items[col_type].c_str());
                    }
//                } else if (items[col_type] == "UTR") {
                }
            } 
        }
    }
    if (verbose) cerr << "rearrange genes symbol\n";
//        for (vector<seq_gene*>::iterator it = gen
//        for (map<int,string>::const_iterator it = geneid2symbol
    
//        cerr << "set pseudogene\n";
    for (vector<spliced_gene*>::iterator it = genes.begin(); it != genes.end(); it++) {
        spliced_gene* g = (*it);
        if ((*it)->symbol() == "") {
            int geneid = g->geneid();
            if (geneid > 0 && geneid2symbol.find(geneid) != geneid2symbol.end()) {
                g->set_symbol(geneid2symbol[geneid]);
                if (verbose) cerr << "rename GeneID:" << geneid << " as " << geneid2symbol[geneid] << endl;
            }
        }
        if (pseudo.find(g->geneid()) != pseudo.end()) {
            g->set_pseudo();
        }
        g->set_introns();
    }
    file_in.close();
    if (verbose) cerr << "saving to cache file\n";
    sort(genes.begin(), genes.end(), seq_gene::compare_tss_position);
    vector<seq_gene*> g_;
    for (vector<spliced_gene*>::iterator it = genes.begin(); it != genes.end(); it++) {
        g_.push_back(*it);
    }
    save_cache(filename_cache.c_str(), g_, filename, assembly);
    return genes;
}

splice_section::splice_section(int start, int end, const char* feature) throw (invalid_argument) {
    _start = start;
    _end = end;
    _feature = get_feature_code(feature);
    _value = 0;
    if (_feature != splice_section::exon 
        || _feature != splice_section::intron 
        || _feature != splice_section::utr) {
        throw invalid_argument(string(feature) + " is not supported");
    }
}

const char* splice_section::feature() const {
    return get_feature_label(_feature);
}

splice_section::section_feature splice_section::get_feature_code(const char* feature) {
    if (strcmp(feature, "CDS") == 0 || strcmp(feature, "exon") == 0) {
        return exon;
    } else if (strcmp(feature, "intron") == 0) {
        return intron;
    } else if (strcmp(feature, "UTR") == 0) {
        return utr;
    } else {
        return unknown;
    }
}

const char* splice_section::get_feature_label(section_feature code) {
    switch (code) {
    case exon:
        return "exon";
    case intron:
        return "intron";
    case utr:
        return "UTR";
    default:
        return "unknown";
    }
}

spliced_gene::spliced_gene(const string& name, const string& symbol, int geneid) 
    : seq_gene(name, symbol, geneid) {
}

namespace {
    bool _section_comparator(const splice_section& lhs, const splice_section& rhs) {
        return lhs.start() < rhs.start();
    }
}

void spliced_gene::set_introns() {
    sort(_sections.begin(), _sections.end(), _section_comparator);
    int end = _sections.size();
    for (int i = 1; i < end; i++) {
        const splice_section& left = _sections[i-1];
        const splice_section& right = _sections[i];
        if (left.end() + 1 < right.start()) {
            _sections.push_back(splice_section(left.end() + 1, right.start() - 1, splice_section::intron));
        }
    }
    sort(_sections.begin(), _sections.end(), _section_comparator);
}

void spliced_gene::set_section(int start, int end, const char* feature) {
    try {
        splice_section::section_feature f = splice_section::get_feature_code(feature);
        if (f > 0) {
            set_section(start, end, f);
        }
    } catch (exception& e) {
        throw;
    }
}
void spliced_gene::set_section(int start, int end, splice_section::section_feature feature) {
    _sections.push_back(splice_section(start, end, feature));
}

bool tkbio::operator == (const splice_section& lhs, const splice_section& rhs) {
    return lhs._start == rhs._start && lhs._end == rhs._end;
}

splice_section& spliced_gene::get_section_at_index(int index) {
    if (index < 0 || index >= (int)_sections.size()) {
        return splice_section::out_of_range;
    } else {
        return _sections[index];
    }
}

splice_section& spliced_gene::get_section_at_position(int position) {
    for (int i = 0; i < (int)_sections.size(); i++) {
        splice_section& s = _sections[i];
        if (s._start <= position && position <= s._end) {
            return s;
        }
    }
    return splice_section::out_of_range;
}

void spliced_gene::serialize(ostream& file_out) const {
    const spliced_gene* g = this;
    file_out << g->genbank() << "\t" << g->symbol() << "\t" << g->geneid() << "\t";
    file_out << g->chromosome() << "\t" << g->orientation() << "\t" << g->start() << "\t" << g->end() << "\t";
    for (int i = 0; i < g->count(); i++) {
        const splice_section& ss = g->get_feature(i);
        if (i > 0) file_out << ",";
        file_out << ss.feature() << ":" << ss.start() << ":" << ss.end();
    }
    file_out <<  "\t" << (g->is_pseudo() ? "pseudo" : "coding") << "\n";
}

ostream& operator << (ostream& ost, const seq_gene& gene) {
    gene.serialize(ost);
    return ost;
}
