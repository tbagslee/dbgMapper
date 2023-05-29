/*
 * verification.cpp
 *
 *  Created on: May 25, 2020
 *      Author: pluto
 */

#include "verification.h"

void pathCan2posCan(list<pth_candidate> &lcand, const sFMindex &fmindx, list<pos_candidate> &lposCan)
{

	for(auto ite = lcand.begin(); ite != lcand.end(); ++ite)
	{
		struct pos_candidate tmp(ite->cur_seed_id);

		//1)calculate the suffix interval of the candidate path
		uint32_t *saarry = calc_SArangeSeq(fmindx, ite->p_path);

		//2)calculate the value of sa[i]
		uint32_t sa_i;
		for(uint32_t j = saarry[0]; j <= saarry[1]; ++j)
		{
			sa_i = calc_SA(fmindx,j);

			//3)sa[i] += pathCan.cur_start_pos
			sa_i += (uint32_t)(ite->cur_start_pos);

			//4)tmp.ref_pos=sa[i];
			tmp.ref_pos = sa_i;

			//5)posCan.push_back(tmp);
			lposCan.push_back(tmp);
		}
		if(saarry != nullptr)
		{
			free(saarry);
			saarry = nullptr;
		}
	}
}

void verification(const vector<ms_seed> &vseed, const list<pos_candidate> &lposCan, list<ms_result> &lrslt,
				const char *read, const char *ref, uint32_t tau)
{
	vector<ms_seed>::const_iterator ites = --vseed.end();
	list<pos_candidate>::const_iterator itec = lposCan.begin();
	uint32_t read_len = ites->end_pos + 1;
	struct ms_result rsltmp;
	char *readtmp = const_cast<char *>(read);
	char *reftmp = const_cast<char *>(ref);
	char *readcur;
	uint32_t seed_len;
	uint32_t linel,liner;
	for(;itec != lposCan.end(); ++itec)
	{
		uint32_t tautmp = tau;
		readcur = readtmp + vseed[itec->seed_id].start_pos;
		seed_len = strlen(vseed[itec->seed_id].seed_segment);
		uint32_t remain_len = read_len - seed_len;
		if(!itec->seed_id)	//first seed of read and expand to right
		{
			readcur += seed_len;
			liner = min_lastlineEd(readcur, reftmp+itec->ref_pos+seed_len, remain_len, remain_len+tautmp, tautmp);
			if(liner != -1)
			{
				rsltmp.ref_start = itec->ref_pos;
				rsltmp.ref_end = rsltmp.ref_start + seed_len - 1 + liner;
				rsltmp.ed = tau - tautmp;
				lrslt.push_back(rsltmp);
			}
		}
		else if(itec->seed_id == vseed.size()-1)	//last seed of read and expand to left
		{
			readcur = readtmp;
			char *read_r = new char[remain_len+1]();
			strncpy(read_r,readcur,remain_len);
			reverseq(read_r);
			char *reftmp_r = new char[remain_len+tautmp+1]();
			strncpy(reftmp_r,reftmp+itec->ref_pos-remain_len-tautmp,remain_len+tautmp);
			reverseq(reftmp_r);
			linel = min_lastlineEd(read_r, reftmp_r, remain_len, remain_len+tautmp, tautmp);
			if(linel != -1)
			{
				rsltmp.ref_start = itec->ref_pos - linel;   //bug
				rsltmp.ref_end = itec->ref_pos + seed_len - 1;
				rsltmp.ed = tau - tautmp;
				lrslt.push_back(rsltmp);
			}
			delete [] read_r;
			delete [] reftmp_r;
		}
		else	//middle seed of read and expand to bidirection
		{
			remain_len -= vseed[itec->seed_id].start_pos;
			readcur += seed_len;
			liner = min_lastlineEd(readcur, reftmp+itec->ref_pos+seed_len, remain_len, remain_len+tautmp, tautmp);
			if(liner != -1)
			{
				readcur = readtmp;
//				rextlen -= liner;
				remain_len = vseed[itec->seed_id].start_pos;
				char *read_r = new char[remain_len+1]();
				strncpy(read_r,readcur,remain_len);
				reverseq(read_r);
				char *reftmp_r = new char[remain_len+tautmp+1]();
				strncpy(reftmp_r,reftmp+itec->ref_pos-remain_len-tautmp,remain_len+tautmp);
				reverseq(reftmp_r);
				linel = min_lastlineEd(read_r, reftmp_r, remain_len, remain_len+tautmp, tautmp);
				if(linel != -1)
				{
					rsltmp.ref_start = itec->ref_pos - linel; //bug
					rsltmp.ref_end = itec->ref_pos + seed_len-1 + liner;
					rsltmp.ed = tau - tautmp;
					lrslt.push_back(rsltmp);
				}
				delete[] read_r;
				delete[] reftmp_r;
			}
		}
	}

}

void get_result(const char *read_path, const char *ref, const sFMindex &fmindx, uint8_t tau, uint32_t threhld)
{
	//1)get readline num and set ptr to NULL
	uint32_t readline = get_file_linenums(read_path);
	if(readline == 0)
	{
		printf("%s contains no reads, program return...\n",read_path);
		return;
	}
	char **p2_read = (char**)malloc(sizeof(char*) * readline);

	//2)put every line in read file to p2_read[i]
	FILE *fp = fopen(read_path, "r");
	for(int32_t i = 0; i < readline; ++i)
	{
		size_t len = 0;
		p2_read[i] = NULL;
		getline(&p2_read[i], &len, fp);
	}
	fclose(fp);

	//3)generate seed from every read
	vector<PH_Node> vPHN;
	generate_PHNArray(vPHN, tau);

	//4)generate exact matching candidates
	vector<ms_seed> vseed;
	list<ms_result> lrslt;
//	ofstream outfile;
//	outfile.open("FreqMore1wReadRecord", ios::out);
//	vector<uint32_t> vcandfreq;
//	vcandfreq.resize(readline);

	struct timeval tvs,tve;


	double canditime = 0;
	double postime = 0;
	double vertime = 0;

	uint32_t lowerHreadnum = 0;
	for(int32_t i = 0; i < readline; ++i)
	{

		generate_seeds(p2_read[i], tau, vseed);

//		for(uint32_t j = 0; j < vseed.size(); ++j)
//		{
//			uint32_t sa_i;
//			uint32_t * sainterval = calc_SArangeSeq(fmindx,vseed[j].seed_segment);
//			for(uint32_t k = sainterval[0]; k <= sainterval[1]; ++k)
//			{
//				sa_i = calc_SA(fmindx,k);
//			}
//			if(sainterval != nullptr)
//			{
//				free(sainterval);
//				sainterval = nullptr;
//			}
//		}

//
		gettimeofday(&tvs,NULL);
		list<pth_candidate> lcand;
		uint32_t candinum = init_candidate(vseed, fmindx, lcand);
		gettimeofday(&tve,NULL);
		canditime += tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		if(candinum < threhld)
		{
			lowerHreadnum++;
			gettimeofday(&tvs,NULL);
			list<pos_candidate> lposCan;
			pathCan2posCan(lcand, fmindx, lposCan);
			gettimeofday(&tve,NULL);
			postime += tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;

			gettimeofday(&tvs,NULL);
			verification(vseed, lposCan, lrslt, p2_read[i], ref, tau);
			gettimeofday(&tve,NULL);
			vertime += tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		}
//		if(candinum > 10000)
//		{
//			outfile << p2_read[i];
//		}
//		iter_candi2PHTroot(vseed, vPHN, lcand);
		free_pthcandi(lcand);
		free_seeds(vseed);
		vseed.clear();
	}
	cout << "candidate number lower " << threhld << " is:" << lowerHreadnum << endl;
	cout << "calc interval time:" << canditime << " calc sa[i] time:" << postime << " calc verification time:" << vertime << endl;
//	sort(vcandfreq.begin(), vcandfreq.end(),greater<uint32_t>());
//	for(int32_t i = 0; i < readline; ++i)
//	{
//		outfile << vcandfreq[i] << endl;
//	}
//	outfile.close();

	for(int32_t i = 0; i < readline; ++i)
	{
		free(p2_read[i]);
		p2_read[i] = nullptr;
	}
	free(p2_read);

}
