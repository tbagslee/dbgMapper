/*
 * main.c
 *
 *  Created on: Sep 12, 2019
 *      Author: pluto
 */

#include "analyse_dBG.h"
#include "basic.h"
#include "Binary_Search.h"
#include "read.h"
#include "load_DBG_full.h"
#include "seeding.h"
#include "prealignment.h"
#include "method.h"
#include "print.h"
#include "verification.h"


int main(int argc, char** argv)
{

//	const char *p_dbg_path, *p_read_path;
//	char *kmertest;
//	uint32_t kmerlen;
//	uint8_t tau;
//	uint32_t threshold;
//
//	for(int32_t i=1;i<argc;i=i+2)
//	{
//		if(argv[i][0] == '-' && argv[i][1] == 'R')//
//		{
//			p_dbg_path = argv[i+1];
//		}
//		if(argv[i][0] == '-' && argv[i][1] == 'r')//
//		{
//			p_read_path = argv[i+1];
//		}
//		if(argv[i][0] == '-' && argv[i][1] == 'L')//
//		{
//			kmerlen = atoi(argv[i+1]);
//		}
//		if(argv[i][0] =='-' && argv[i][1] == 't')//
//		{
//			tau = atoi(argv[i+1]);
//		}
//
//		if(argv[i][0] =='-' && argv[i][1] == 'H')//
//		{
//			threshold = atoi(argv[i+1]);
//		}
//
//	}
////	gen_ed_readfile(p_read_path, 0, 10000);
////	cout << "generate fininhed" << endl;
////	getchar();
//
//	struct timeval tvs,tve;
//
////	1)read FMindex
////	char nindex[] = "./nindex";
//	char rindex[] = "./rindex";
////	read_bfile2index(nindex,gnFMidx,0);
//	read_bfile2index(rindex,grFMidx,0);
//	gPair[0] = make_pair('O',&grFMidx);
////	gPair[1] = make_pair('I',&gnFMidx);
//
//	//2)generate non-branching kmers and sort them,load dbg index
//	get_para(&gbit_para, kmerlen);
//	char *ref;
//	gen_dBG_index(gbit_para, gdBGindex, p_dbg_path, 1, &ref);
//	gettimeofday(&tvs,NULL);
//	test_extime(p_read_path, threshold, tau);
////	get_result(p_read_path, ref, gnFMidx, tau, threshold);


//	3) test myers 和 pre_alignment
	char *pattern = "ATACTTACGTACGTCAAGATCAAGATCAAGATCAAGATCAAGAAAGCGCTT";
//	char *pattern = "GATTCTTCCTGGTATTTTTCAAGACACTCAAATTAAAATTTAAACTTCAGAGGTTATTATGGCTTTGGAGATGGTACAGCCATTCTGAATAAGAATTTGG";
	char *text = "AATCAAGATCAAGATCAAGATCAAGA";
	int mapping_end_position;
	int a = banded_edit_distance(pattern, text, strlen(text), 8,&mapping_end_position);
	cout << a <<endl;
	cout<<mapping_end_position<<endl;
	cout<<strlen(pattern)<<endl;
    int b = pre_alignment(text, pattern, 25, 25, 3);
    cout << b << endl;
//
//	gettimeofday(&tve,NULL);
//	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
//	cout << "total time is: "<< span << endl;
////	free_dBGindex(gdBGindex);
	return 0;
}
