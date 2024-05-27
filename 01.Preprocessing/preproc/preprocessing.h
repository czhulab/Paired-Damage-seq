//
// Chenxu, 6.22.2023
// Module for preprocessing
//

#ifndef REACHTOOLS2_PREPROCESSING_H
#define REACHTOOLS2_PREPROCESSING_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

class preprocessing {
private:

public:
  static void combine_doubleCut(string input);
  static void combine_384plex(string input);
  static void combine_384plex_pe(string input);
	static void convert(string input);
	static void convert_pe(string input);
};

class read2_dc { // #8bp barcodes
private:
	int align_score(string str1, string str2){
		int score = 0;
    for(int i = 0; i < str1.length(); ++i){
      if(str2[i] == 'N')continue;
      str2[i] == str1[i] ? (score += 2) : (score -= 1);
    }
    return score;
	}
	
public:
	int bc1, bc2, bc4;
	string bc, rawline, sbc1, sbc2, sbc4, bsbc4, type, umi;
	int dock;
	bool valid;

	void init(string line){
		rawline = line;
		dock = -1;
		valid = false;
		return;
	}

	void trim(){
		// NNNNNNNNNNNNNNNNNNnnnnGTGGCCGATGTTTCGGTGCGAACTCAGACCNNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCGXXNNNNN
		//"NNNNNNNNNNNNNNNNNNGTGGCCGATGTTTCGGTGCGAACTCAGACCNNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCGXNNNNYY";
		// 1 determine additional Ns

		if(rawline.length() < 96)return;
		dock = 0;
		int t = 0;
		int cur_s = 10;
		string bait = "GTGGCCGATGTTTCG";
		for(int i = -2; i < 7; ++i){
			string qu = rawline.substr(18+i, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			t = i;
		}
		//if(cur_s<15)return;
		if(cur_s<13)return;
		dock = 1;
		sbc1 = rawline.substr(10, 8);
		umi = rawline.substr(0, 10);
		if(t==-1){
			umi = "N" + rawline.substr(0, 9);
		}
		else if(t==-2){
			umi = "NN" + rawline.substr(0, 8);
		}

		//2nd bc
		bait = "ATCCACGTGCTTGAG";
		cur_s = 10;
		int tt = t;
		for(int i = -2; i < 3; ++i){
			string qu = rawline.substr(56+i+tt, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			tt = t + i;
		}
		//if(cur_s<13)return;
		if(cur_s<12)return;
		t = tt;
		sbc2 = rawline.substr(48+t, 8);
		dock = 2;

		//4th BC
		bait = "AGGCCAGAGCATTCG";
		cur_s = 10;
		for(int i = -2; i < 3; ++i ){
			string qu = rawline.substr(71+i+tt, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			tt = t + i;
		}
		if(cur_s<11)return;
		t = tt;
		if(rawline.length() < 89+t)return;
		sbc4 = rawline.substr(87+t, 4);
		dock = 4;		

		cur_s = 0;
		string b1 = rawline.substr(86+t, 1);
		string b2 = rawline.substr(91+t, 2);

		type = "n";
		if(b1 == "A"){
			if(b2 == "GC"){
				type = "h";
			}
			else if(b2 == "AT" || b2 == "TC" || b2 == "CA" || b2 == "TA" || b2 == "AC" || b2 == "TT"){
				type = "a";
			}
		}
		else if(b1 == "T"){
			type = "r";
		}
		return;
	}


	bool is_valid(){
		return valid;
	}
	int where_dock(){
		return dock;
	}
};


class read2_dna { // #8bp barcodes
private:
	int align_score(string str1, string str2){
		int score = 0;
    for(int i = 0; i < str1.length(); ++i){
      if(str2[i] == 'N')continue;
      str2[i] == str1[i] ? (score += 2) : (score -= 1);
    }
    return score;
	}
	
public:
	int bc1, bc2, bc4;
	string bc, rawline, sbc1, sbc2, sbc4, bsbc4, type, umi;
	int dock;
	bool valid;

	void init(string line){
		rawline = line;
		dock = -1;
		valid = false;
		return;
	}

	void trim(){
		// NNNNNNNNNNNNNNNNNNnnnnGTGGCCGATGTTTCGGTGCGAACTCAGACCNNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCGXXNNNNN
		//"NNNNNNNNNNNNNNNNNNGTGGCCGATGTTTCGGTGCGAACTCAGACCNNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCGXNNNNYY";
		// 1 determine additional Ns

		if(rawline.length() < 95)return;
		dock = 0;
		int t = 0;
		int cur_s = 10;
		string bait = "GTGGCCGATGTTTCG";
		for(int i = -2; i < 7; ++i){
			string qu = rawline.substr(18+i, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			t = i;
		}
		//if(cur_s<15)return;
		if(cur_s<13)return;
		dock = 1;
		sbc1 = rawline.substr(10, 8);
		umi = rawline.substr(0, 10);
		if(t==-1){
			umi = "N" + rawline.substr(0, 9);
		}
		else if(t==-2){
			umi = "NN" + rawline.substr(0, 8);
		}

		//2nd bc
		bait = "ATCCACGTGCTTGAG";
		cur_s = 10;
		int tt = t;
		for(int i = -2; i < 3; ++i){
			string qu = rawline.substr(56+i+tt, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			tt = t + i;
		}
		//if(cur_s<13)return;
		if(cur_s<12)return;
		t = tt;
		sbc2 = rawline.substr(48+t, 8);
		dock = 2;

		//4th BC
		bait = "AGGCCAGAGCATTCG";
		cur_s = 10;
		for(int i = -2; i < 3; ++i ){
			string qu = rawline.substr(71+i+tt, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			tt = t + i;
		}
		if(cur_s<11)return;
		t = tt;
		if(rawline.length() < 89+t)return;
		sbc4 = rawline.substr(87+t, 4);
		dock = 4;		

		cur_s = 0;
		string b1 = rawline.substr(86+t, 1);
		string b2 = rawline.substr(91+t, 1);

		type = "n";

		if(b1=="A" || b2 == "A"){
			type = "d";
			if(b1=="T"||b2=="C")type="n";
		}
		if(b1=="T" || b2 == "C"){
			type = "r";
			if(b1=="A" || b2 == "A")type="n";
		}


		return;
	}


	bool is_valid(){
		return valid;
	}
	int where_dock(){
		return dock;
	}
};





#endif //REACHTOOLS2_PREPROCESSING_H