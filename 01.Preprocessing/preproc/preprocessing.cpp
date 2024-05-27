// 06/26/2023
// Chenxu Zhu

#include "cxstring.h"
#include "preprocessing.h"



void preprocessing:: combine_384plex(string r2){
  int total = 0;
  int pass = 0;
  int atac = 0;
  int rna = 0;
  // int histone = 0;
  int und = 0;

  string s1;
  string s2;
  string s3;

  // open log file to save demultiplexing results
  s1 = "cat - > ";
  s2 = r2 + "_type_report.txt";
  s3 = s1 + s2;
  FILE * outType;
  outType = popen(s3.c_str(), "w");


  // open Read1 and Read2 files, must in gz format.
  s1 = "zcat ";
  s2 = r2 + "_R1.fq.gz";
  s3 = s1 + s2;
  FILE *readfile1;
  readfile1 = popen(s3.c_str(), "r");

  s1 = "zcat ";
  s2 = r2 + "_R2.fq.gz";
  s3 = s1 + s2;
  FILE *readfile2;
  readfile2 = popen(s3.c_str(), "r");

  // open output file three modalities
  // PAT, pA-Tn5 derived reads
  // TN5, Tn5 derived reads
  // RNA, RT derived reads
  // For DNA modality, file sizes for _RNA.fq.gz should be small; for RNA modality, file sizes for _PAT.fq.gz & _Tn5.fq.gz should be small.
  // s1 = "gzip - > ";
  // s2 = r2 + "_combined_PAT.fq.gz";
  // s3 = s1 + s2;
  // FILE *outPAT;
  // outPAT = popen(s3.c_str(), "w");

  s1 = "gzip - > ";
  s2 = r2 + "_combined_TN5.fq.gz";
  s3 = s1 + s2;
  FILE *outTN5;
  outTN5 = popen(s3.c_str(), "w");

  s1 = "gzip - > ";
  s2 = r2 + "_combined_RNA.fq.gz";
  s3 = s1 + s2;
  FILE *outRNA;
  outRNA = popen(s3.c_str(), "w");

  // init read lines
  char buffer[2000];
  fqline in_line1;
  fqline in_line2;

  while(fgets(buffer, sizeof(buffer), readfile1)){
    ++total;
    string line1(buffer);
    fgets(buffer, sizeof(buffer), readfile2);
    string line2(buffer);

    line1 = cxstring::chomp(line1);
    line2 = cxstring::chomp(line2);

    in_line1.read_part_record(readfile1, line1);
    in_line2.read_part_record(readfile2, line2);

    read2_dna read_2; // init a read2_dna object - name defined in preprocessing.h - not only for DNA.
    read_2.init(in_line2.seq);
    read_2.trim(); // get barcodes and type

    string new_seq = read_2.sbc1 + read_2.sbc2 + read_2.sbc4;
    string umi = read_2.umi;

    in_line2.seq = new_seq;
    in_line2.qual = in_line2.qual.substr(0, in_line2.seq.length()); // not a valid qual value, just to fill the fastq format

    if(in_line2.seq.length()!=20)continue; // 8 + 8 + 4 . but not necessary to do here. 
    pass++;

    string type = read_2.type;
    if(type == "d")atac++;
    if(type == "r")rna++;
    // if(type == "h")histone++;
    if(type == "n")und++;

    // proc new fastq line format
    string a;
    stringstream as;

    as << in_line2.readname;
    as >> a;
    in_line2.readname = a + "_" + umi + "_" + in_line1.seq + "_" + in_line1.qual; // replace ":" with "_" in reachtools2 to avoid incompatible platforms

    if(type == "d")in_line2.write_record(outTN5);

    if(type == "r")in_line2.write_record(outRNA);
  
  }
  pclose(readfile1);
  pclose(readfile2);

  pclose(outTN5);
  pclose(outRNA);

  // full barcode ratio
  float pass_f = (float)pass;
  float total_f = (float)total;
  float aaa = pass_f * 10000.f / total_f;
  int bbb = (int)aaa;
  float r_fbc = (float)bbb / 100.f;

  // atac ratio
  float atac_f = (float)atac;
  aaa = atac_f * 10000.f / total_f;
  bbb = (int)aaa;
  float r_atac = (float)bbb / 100.f;

  // rna ratio
  float rna_f = (float)rna;
  aaa = rna_f * 10000.f / total_f;
  bbb = (int)aaa;
  float r_rna = (float)bbb / 100.f;

  // und ratio
  float und_f = (float)und;
  aaa = und_f * 10000.f / total_f;
  bbb = (int)aaa;
  float r_und = (float)bbb / 100.f;

  fputs(("Total #:\t" + cxstring::int2str(total) + "\n").c_str(), outType);
  fputs(("FullBC #:\t" + cxstring::int2str(pass) + "\n").c_str(), outType);
  fputs(("DNA assigned#:\t" + cxstring::int2str(atac) + "\n").c_str(), outType);
  fputs(("RNA assigned #:\t" + cxstring::int2str(rna) + "\n").c_str(), outType);
  fputs(("Unknown #:\t" + cxstring::int2str(und) + "\n").c_str(), outType);
  fputs(("FullBC %:\t" + cxstring::float2str(r_fbc) + "\n").c_str(), outType);
  fputs(("DNA assigned %:\t" + cxstring::float2str(r_atac) + "\n").c_str(), outType);
  fputs(("RNA assigned %:\t" + cxstring::float2str(r_rna) + "\n").c_str(), outType);
  fputs(("Unknown %:\t" + cxstring::float2str(r_und) + "\n").c_str(), outType);


  pclose(outType);


}



void preprocessing:: convert(string prefix){
 	// processing Paired-seq2 with 2-round ligation
	int total = 0;
	int pass = 0;
	string s1 = "cat ";
	string s2 = prefix;
	string s3 = s1 + s2;
	FILE * inbam;
	inbam = popen(s3.c_str(), "r");
	s1 = "gzip - > ";
	s2 = prefix.substr(0, prefix.length()-4) + "_cov.fq.gz";
	s3 = s1 + s2;
	FILE * fout;
	fout = popen(s3.c_str(), "w");
	samline align_line;
	fqline fastq_line;
	char buffer[2000];
	while(fgets(buffer, sizeof(buffer), inbam)){
		string line(buffer);
		line=cxstring::chomp(line);
		if(line.substr(0, 1) == "@"){
			continue;
		}
		++total;
		align_line.init(line);
		if(align_line.chr == "*")continue;
		vector<string> tmp = cxstring::split(align_line.readname, "_");
		// fastq_line.readname = "@" + tmp[0] + ":" + tmp[1] + ":" + tmp[2] + ":" + tmp[3] + ":" + tmp[4] + ":" + tmp[5] + ":" + tmp[6] + ":" +  align_line.chr + ":" + tmp[7];
    fastq_line.readname = "@" + tmp[0]+":"+align_line.chr + ":" + tmp[1];
		fastq_line.seq = tmp[2];
		fastq_line.qual = tmp[3];
		fastq_line.mark = "+";
		fastq_line.write_record(fout);
		++pass;
	}
	pclose(inbam);
	pclose(fout);
	cout << total << " reads processed." << endl;
	cout << pass << " mapped reads." << endl;
	return;
}


