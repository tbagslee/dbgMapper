/*
 * seeding.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: bio
 */

//不同的seeding的方法
//seeding的操作是将原read序列拆成包含tau+1个seeds的划分
//目标是拆分出的seed的目标
//1）不能精确匹配的seed数量越多越好
//2）精确匹配的seed的长度越长越好

#include "seeding.h"
#include "Hash.h"
#include "Binary_Search.h"
#include "method.h"
#include "print.h"
int startid, endid;
std::pair<char, const struct sFMindex *> gPair[2];

void generate_seeds(const char *read, uint8_t tau, vector<ms_seed> &vseed)
{
	uint32_t readlen = strlen(read);
	if(!readlen)
	{
		return;
	}
	if(read[readlen-1] == '\n')
	{
		//deal '\n' character
		readlen--;
	}
	uint8_t segnum = tau + 1;
	vseed.resize(segnum);
	struct ms_seed seedtmp;
	uint32_t start = 0;
	uint32_t cnt = 0;
	while(segnum)
	{
		assert(segnum != 0);
		uint32_t lentmp = (double)readlen / segnum + 0.5;
		seedtmp.seed_segment = new char[lentmp+1]();
		strncpy(seedtmp.seed_segment,read+start,lentmp);
		seedtmp.start_pos = start;
		seedtmp.end_pos = start + lentmp - 1;
		vseed[cnt++] = seedtmp;
		start += lentmp;
		readlen -= lentmp;
		--segnum;
	}
//	uint32_t longcnt = readlen % segnum;
//	uint32_t len = readlen / segnum;
//	for(uint32_t i = 0; i < segnum; ++i)
//	{
//		uint32_t lentmp = len + (longcnt != 0 ? 1 : 0);
//		if(longcnt != 0)
//		{
//			longcnt--;
//		}
//		seedtmp.seed_segment = new char[lentmp + 1]();
//		strncpy(seedtmp.seed_segment,read+start,lentmp);
//		seedtmp.start_pos = start;
//		seedtmp.end_pos = start + lentmp - 1;
//		vseed.push_back(seedtmp);
//		start += lentmp;
//	}

}

void free_seeds(vector<ms_seed> & vseed)
{
	for(uint32_t i = 0; i != vseed.size(); ++i)
	{
		if(vseed[i].seed_segment != nullptr)
		{
			delete[] vseed[i].seed_segment;
			vseed[i].seed_segment = nullptr;
		}
	}
}

void generate_PHNArray(vector<PH_Node> &PHArray, uint8_t tau)
{
	uint32_t level = ceil(log(tau+1) / log(2));
	uint32_t nodenum = pow(2,level+1);
	PHArray.resize(nodenum);
	for(uint32_t i = 0; i < nodenum; ++i)
	{
		PHArray[i].start_seed_id = -1;
		PHArray[i].end_seed_id = -1;
		PHArray[i].tau = -1;
	}
	PHArray[1].start_seed_id = 0;
	PHArray[1].end_seed_id = tau;
	PHArray[1].tau = tau;
	uint8_t cnt = 0;
	for(uint32_t i = 2; i < nodenum; ++i)
	{
		if(PHArray[i/2].tau)
		{
			if(i % 2 == 0)
			{
				if(PHArray[i/2].tau == 1 || PHArray[i/2].tau == 3 || PHArray[i/2].tau == 7)
				{
					PHArray[i].tau = PHArray[i/2].tau / 2;
				}
				else
				{
					PHArray[i].tau = (PHArray[i/2].tau + 1) / 2;
				}
				PHArray[i].start_seed_id = PHArray[i/2].start_seed_id;
				PHArray[i].end_seed_id = PHArray[i].start_seed_id + PHArray[i].tau;
			}
			else
			{
				PHArray[i].tau = PHArray[i/2].tau - PHArray[i-1].tau - 1;
				PHArray[i].start_seed_id = PHArray[i-1].end_seed_id + 1;
				PHArray[i].end_seed_id = PHArray[i].start_seed_id + PHArray[i].tau;
			}
			if(PHArray[i].tau == 0)
			{
				++cnt;
				if(cnt == tau + 1)
				{
					break;
				}
			}
		}

	}
//	int p = 1;
//	for(uint32_t i = 1; i < nodenum; ++i)
//	{
//		cout << + static_cast<uint8_t>(PHArray[i].tau) << " " << + static_cast<uint8_t>(PHArray[i].start_seed_id) << "-" << \
//				+ static_cast<uint8_t>(PHArray[i].end_seed_id) << "\t";
//		if(i == pow(2,p) - 1)
//		{
//			cout << endl;
//			p++;
//		}
//	}
}

uint32_t find_PHNAindex(const vector<PH_Node> &PHArray, uint8_t start_id, uint8_t end_id)
{
	for(uint32_t i = 1; i < PHArray.size(); ++i) {	//if tau <= 5
		if(PHArray[i].start_seed_id == start_id && PHArray[i].end_seed_id == end_id) {
			return i;
		}
		if(PHArray[i].start_seed_id == -1 && PHArray[i].end_seed_id == -1 && PHArray[i].tau == -1) {
			break;
		}
	}
	return 0; //not found
}

void calc_edarray(struct TPTnode *pnode, char *seq, uint8_t tau)
{
	pnode->edarry = (uint8_t *)malloc(sizeof(uint8_t)*(2*tau+1));
	memset(pnode->edarry, tau+1, 2*tau+1);
	uint8_t flag;
	int upper = pnode->level - tau;
	int seqlen = strlen(seq);
	int start = 0;
	if(upper > 1)
	{
		start = upper - 1;
	}
	uint8_t left, top ,leftop;
	for(int i = 0; i < 2*tau+1; i++)
	{
		flag = 1;
		if(upper + i <= 0)
		{
			if(upper + i == 0)
			{
				pnode->edarry[i] = pnode->level;
			}
			continue;
		}

		if(upper + i == seqlen + 1)
		{
			break;
		}
		if(pnode->c == seq[start++])
		{
			flag = 0;
		}
		leftop = pnode->p_parent->edarry[i] + flag;
		if(!i)
		{
			left = pnode->p_parent->edarry[i+1] + 1;
			pnode->edarry[i] = min(leftop,left);
			continue;
		}
		top = pnode->edarry[i-1] + 1;
		if(i == 2*tau)
		{
			pnode->edarry[i] = min(leftop,top);
			break;
		}
		left = pnode->p_parent->edarry[i+1] + 1;
		pnode->edarry[i] = min(leftop,top);
		pnode->edarry[i] = min(left,pnode->edarry[i]);
	}
}

void calc_edarray_iter(char *seq, char c, uint8_t **p2_parent, uint8_t lev, uint8_t tau)
{
	uint8_t *p_cur = (uint8_t*) malloc(sizeof(uint8_t) * (2 * tau + 1));
	memset(p_cur, tau+1, 2*tau + 1);
	uint8_t flag;
	int upper = lev - tau;
	int seqlen = strlen(seq);
	int start = 0;		//start position of alignseq
	if (upper > 1) {
		start = upper - 1;
	}
	uint8_t left, top, leftop;
	for (int32_t i = 0; i < 2 * tau + 1; i++) {
		flag = 1;
		if (upper + i <= 0) {
			if(upper + i == 0)
			{
				p_cur[i] = lev;
			}
			continue;
		}
		if (upper + i == seqlen + 1) {
			break;
		}
		if (c == seq[start++]) {
			flag = 0;
		}
		leftop = (*p2_parent)[i] + flag;
		if (!i) {
			left = (*p2_parent)[i + 1] + 1;
			p_cur[i] = min(leftop, left);
			continue;
		}
		top = p_cur[i - 1] + 1;
		if (i == 2 * tau) {
			p_cur[i] = min(leftop, top);
			break;
		}
		left = (*p2_parent)[i + 1] + 1;
		p_cur[i] = min(leftop, top);
		p_cur[i] = min(left, p_cur[i]);
	}
	free(*p2_parent);
	*p2_parent = p_cur;
}

void calc_kmer_edVec(struct dBG_extpara &ext_para, struct TPTnode *pnode, const struct pth_candidate *pcandi)
{
	uint8_t lev = 1;	//this might be wrong
	if(ext_para.dir == 'I')
	{
		for(int i = pcandi->cur_start_pos-1; i >= 0; --i)
		{
			calc_edarray_iter(ext_para.alignseq, *(pcandi->p_path+i), &(pnode->edarry), lev++, ext_para.tau);
			pnode->level = pcandi->cur_start_pos;
		}
	}
	else
	{
		uint8_t path_len = strlen(pcandi->p_path);
		for(int i = pcandi->cur_end_pos+1; i < path_len; ++i)
		{
			calc_edarray_iter(ext_para.alignseq, *(pcandi->p_path+i), &(pnode->edarry), lev++, ext_para.tau);
			pnode->level = path_len - pcandi->cur_end_pos - 1;
		}
	}
}


void init_rootnode(struct TPTnode *pnode, struct dBG_extpara &ext_para, struct pth_candidate *pcandi)   //mod == 0 edvec 0~tau   mod != 0 edvec calc the cur_seed_id + 1 ~ parent's end_seed_id
{
	pnode->c = 'R';
	pnode->offset = 0;
	pnode->edarry = (uint8_t*)malloc(sizeof(uint8_t)*(2*ext_para.tau+1));
	memset(pnode->edarry, 0, 2*ext_para.tau+1);
	pnode->p_parent = nullptr;
	for(int i = 0 ; i < 4; i++)
	{
		pnode->p_child[i] = nullptr;
	}
	if(pcandi == nullptr)
	{
		pnode->level = 0;
		int tmp = 0;
		for(int i = ext_para.tau+1; i < 2*ext_para.tau+1; i++)
		{
			pnode->edarry[i] = ++tmp;
		}
	}
	else
	{
		calc_kmer_edVec(ext_para, pnode, pcandi);
	}
	char *tmpseq = (char *)malloc(sizeof(char)*(strlen(ext_para.orignseq)+1));
	memset(tmpseq,'\0',strlen(ext_para.orignseq)+1);
	strcpy(tmpseq,ext_para.orignseq);
	if(ext_para.dir == 'O')
	{
		reverseq(tmpseq);
	}
	pnode->saarry = calc_SArangeSeq(*ext_para.pFMidx,tmpseq);
	if(tmpseq != nullptr)
	{
		free(tmpseq);
		tmpseq = nullptr;
	}
}

void push_vseq_upward(struct TPTnode *pnode, list<char *>& lseq);

bool init_childnode(const struct dBG_extpara &ext_para, struct TPTnode *pnodec, struct TPTnode *pnodep)  //返回真表示继续扩展
{
	pnodec->p_parent = pnodep;
	pnodec->level = pnodep->level + 1;
	calc_edarray(pnodec, ext_para.alignseq, ext_para.tau);		//
	pnodec->saarry = calc_SArangeChar(*ext_para.pFMidx,pnodep->saarry,pnodec->c);
//	if(pnodec->saarry[1] >= pnodec->saarry[0])
//	{
		for(int i = 0; i < 2*ext_para.tau+1; i++)
		{
			if(pnodec->edarry[i] <= ext_para.tau)
			{
				for(uint32_t ii = 0 ; ii < 4; ii++)
				{
					pnodec->p_child[ii] = nullptr;
				}
				return true;
			}
		}
//	}
	return false;
}

char *push_vseq_upward(struct TPTnode *pnode)
{
	int len = 0;
	//add char to string
	struct TPTnode *save = pnode;
	do{
		len++;
		pnode = pnode->p_parent;
	}while(pnode != nullptr);

	char *str = new char[len]();
	int i = 0;
	pnode = save;
	do{
		str[i++] = pnode->c;
		pnode = pnode->p_parent;
	}while(pnode->c != 'R');
	//reverse string
//	reverseq(str);
	return str;
}

void show_nearby_alilen(struct dBG_extpara &ext_para,struct TPTnode *pnode)
{
	//print upward from level == extlen - tau  to level == extlen + tau
	int pos = pnode->level - ext_para.alignlen;
	if(pos >= -ext_para.tau && pos <= ext_para.tau)
	{
		if(pnode->edarry[ext_para.tau-pos] <= ext_para.tau)
		{
			onlyprint_upward(pnode);
		}
	}
}

void find_child_mined(struct dBG_extpara &ext_para,struct TPTnode *pnode, uint8_t &mint)
{
	int pos = pnode->level - ext_para.alignlen;
	if(pnode->edarry[ext_para.tau-pos] < mint)
	{
		mint = pnode->edarry[ext_para.tau-pos];
	}
	for(int i = 0; i < 4; ++i)
	{
		if(pnode->p_child[i] != nullptr)
		{
			find_child_mined(ext_para, pnode->p_child[i], mint);
		}
	}
}


void find_saveupward(struct dBG_extpara &ext_para,struct TPTnode *pnode, list<pair<char*, uint8_t> >& lpair)
{
	uint8_t min_ed = 255;
	find_child_mined(ext_para, pnode, min_ed);
	if(min_ed <= ext_para.tau)
	{
		pair<char*, uint8_t> p = make_pair(push_vseq_upward(pnode), min_ed);
		lpair.push_back(p);
	}
}

void ext_treenode(struct dBG_extpara &ext_para, struct TPTnode *pnode, uint32_t extlen, uint32_t &cnt)//, list<pair<char*, uint8_t> > &lpair)
{
	if(extlen != 0)
	{
		char *original_save = ext_para.orignseq;
		bool extflag = false;
		if(pnode->offset > 0 && pnode->offset < strlen(ext_para.orignseq) - gbit_para.kmer1Len/2) //offset > 1?
		{
			if(ext_para.dir == 'I')
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->offset = pnode->offset - 1;
				pnode->p_child[0]->c = ext_para.orignseq[pnode->p_child[0]->offset];  //
			}
			else
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->offset = pnode->offset + 1;
				pnode->p_child[0]->c = ext_para.orignseq[pnode->offset+gbit_para.kmer1Len/2];
			}
			extflag = init_childnode(ext_para, pnode->p_child[0], pnode);
			if(extflag == true)
			{
				ext_para.onunipath = true;
				//store results
//				show_nearby_alilen(ext_para, pnode->p_child[0]);
				ext_treenode(ext_para, pnode->p_child[0], extlen-1, cnt); //, lpair);
			}
			else
			{
				free(pnode->p_child[0]);
				pnode->p_child[0] = nullptr;
				cnt++;
			}
//			if(pnode->level == ext_para.alignlen - ext_para.tau)
//			{
//				find_saveupward(ext_para, pnode, lpair);
//				destory_extree(pnode);
//			}

		}
		else
		{
			//initial the current k-mer
			uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * gbit_para.kmer64Len);
			uint32_t seqlen = gbit_para.kmer1Len/2;
			char *seqe = (char *)malloc(sizeof(char) * (seqlen+1));
			memset(seqe,0,seqlen+1);
			if(pnode->c == 'R')
			{
				strcpy(seqe,ext_para.orignseq);
			}
			else
			{
				if(true == ext_para.onunipath)
				{
					strncpy(seqe,ext_para.orignseq+pnode->offset,seqlen);
				}
				else
				{
					if(ext_para.dir == 'I')
					{
						seqe[0] = pnode->c;
						strncpy(seqe+1,ext_para.orignseq,seqlen-1);
					}
					else
					{
						strncpy(seqe,ext_para.orignseq+1,seqlen-1);
						seqe[seqlen-1] = pnode->c;
					}
				}
			}
			//calculate the k-mer hash value
			cal_hash_value_directly_256bit(seqe,hashvalue_tmp,gbit_para);
			uint64_t bkmerfindret = ULLONG_MAX,ukmerfindret = ULLONG_MAX;
			//dealing with the case of branching node
			bkmerfindret = Tfind_arrindexN(gdBGindex.p2_bkmer, gdBGindex.p_branchedkmer, \
					hashvalue_tmp,gbit_para.kmer64Len);
			if(bkmerfindret != ULLONG_MAX)
			{
				uint8_t adinfo = gdBGindex.p_branchedkmerad[bkmerfindret];
				if(ext_para.dir == 'I')
				{
					adinfo >>= 4;
				}
				for(uint32_t i = 0; i < 4; i++)
				{
					if(adinfo & 0x01)
					{
						pnode->p_child[i] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
						pnode->p_child[i]->offset = 0;
						switch(i)
						{
							case 0 :
								pnode->p_child[i]->c = 'A';break;
							case 1 :
								pnode->p_child[i]->c = 'C';break;
							case 2 :
								pnode->p_child[i]->c = 'G';break;
							case 3 :
								pnode->p_child[i]->c = 'T';break;
						}
						extflag = init_childnode(ext_para, pnode->p_child[i], pnode);
						if(extflag == true)
						{
							ext_para.onunipath = false;
							ext_para.orignseq = seqe;
							//store results
							ext_treenode(ext_para, pnode->p_child[i], extlen-1,cnt); //, lpair);
						}
						else
						{
							free(pnode->p_child[i]);
							pnode->p_child[i] = nullptr;
						}
					}
					adinfo >>= 1;
				}
				if(pnode->p_child[0] == nullptr && pnode->p_child[1] == nullptr &&\
						pnode->p_child[2] == nullptr && pnode->p_child[3] == nullptr)
				{
					cnt++;
				}
//				if(pnode->level == ext_para.alignlen - ext_para.tau)
//				{
//					find_saveupward(ext_para, pnode, lpair);
//					destory_extree(pnode);
//				}
			}
			else
			{
				// add else case
				//dealing with the case of uni-path k-mer
				ukmerfindret = Tfind_arrindexN(gdBGindex.p2_ukmer, gdBGindex.p_unbranchedkmer, \
						hashvalue_tmp,gbit_para.kmer64Len);
				if(ukmerfindret != ULLONG_MAX)//如果是unipath上的kmer  判断在unipath上的offset 根据是向出度方向还是入度方向扩展来找
				{
					uint32_t ukmerid = gdBGindex.p_unbranchedkmerid[ukmerfindret];
					char *pos = nullptr;
					pos = strstr(gdBGindex.upath_arr[ukmerid],seqe);
					pnode->offset = pos - gdBGindex.upath_arr[ukmerid];
					char ch;
					if(ext_para.dir == 'I')
					{
//						if(!pos)
//						{
//							cout << "unipath reached reference beginning..." << endl;
//							return 0;
//						}
						ch = *(pos-1);
						pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
						pnode->p_child[0]->offset = pnode->offset - 1;
						pnode->p_child[0]->c = ch;
					}
					else
					{
						ch = *(pos + strlen(seqe));
						pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
						pnode->p_child[0]->offset = pnode->offset + 1;
						pnode->p_child[0]->c = ch;
					}
					extflag = init_childnode(ext_para, pnode->p_child[0], pnode);
					if(extflag == true)
					{
						//add a condition
						ext_para.onunipath = true;
						ext_para.orignseq = gdBGindex.upath_arr[ukmerid];
						//store results
//						show_nearby_alilen(ext_para, pnode->p_child[0]);
						ext_treenode(ext_para, pnode->p_child[0], extlen-1, cnt); //, lpair);

					}
					else
					{
						free(pnode->p_child[0]);
						pnode->p_child[0] = nullptr;
						cnt++;
					}
//					if(pnode->level == ext_para.alignlen - ext_para.tau)
//					{
//						find_saveupward(ext_para, pnode, lpair);
//						destory_extree(pnode);
//					}
				}
			}
			free(seqe);
			free(hashvalue_tmp);
		}
		ext_para.orignseq = original_save;
	}
	else {cnt++;}
}

void ext_fromroot(struct dBG_extpara &ext_para, struct TPTnode *proot, list<pair<char*, uint8_t> > &lpair)
{
//	ext_treenode(ext_para, proot, ext_para.alignlen + ext_para.tau, lpair);
}

uint32_t init_candidate(const vector<ms_seed> & vseed, const sFMindex &fmindx, list<pth_candidate> &lcand)
{
	uint32_t candinum = 0;
	for(uint32_t i = 0; i != vseed.size(); ++i)
	{
		pth_candidate candtmp(i);
		uint8_t seedlen = strlen(vseed[i].seed_segment);
		candtmp.p_path = new char[seedlen + 1]();
		strncpy(candtmp.p_path, vseed[i].seed_segment, seedlen);
		candtmp.cur_end_pos = seedlen - 1;
		candtmp.sa_ary[0] = calc_SArangeSeq(fmindx,vseed[i].seed_segment);
		candtmp.sa_ary[1] = (uint32_t *)malloc(sizeof(uint32_t) * 2);
		candtmp.sa_ary[1][0] = candtmp.sa_ary[0][0];
		candtmp.sa_ary[1][1] = candtmp.sa_ary[0][1];
		if(candtmp.sa_ary[0][0] <= candtmp.sa_ary[0][1])
		{
			candinum += candtmp.sa_ary[0][1] - candtmp.sa_ary[0][0] + 1;
			lcand.push_back(candtmp);
		}
		else
		{
			delete[] candtmp.p_path;
			candtmp.p_path = nullptr;
			free(candtmp.sa_ary[0]);
			candtmp.sa_ary[0] = nullptr;
			free(candtmp.sa_ary[1]);
			candtmp.sa_ary[1] = nullptr;
		}
	}
	return candinum;
}

void free_pthcandi(list<pth_candidate> &lcand)
{
	for(auto it = lcand.begin(); it != lcand.end();)
	{
		if(it->p_path != nullptr)
		{
			delete[] it->p_path;
			it->p_path = nullptr;
 		}
		if(it->sa_ary[0] != nullptr)
		{
			free(it->sa_ary[0]);
			it->sa_ary[0] = nullptr;
		}
		if(it->sa_ary[1] != nullptr)
		{
			free(it->sa_ary[1]);
			it->sa_ary[1] = nullptr;
		}
		it = lcand.erase(it);
	}
	lcand.clear();
}

int candi2candi(const vector<ms_seed> & vseed, const vector<PH_Node> &PHArray, list<pth_candidate> &lcand, uint32_t i)
{
	//erase candidate in vector
	list<pth_candidate>::iterator itel = lcand.begin();
	for(uint32_t ii = 0; ii < i; ii++)
	{
		itel++;
	}
	struct pth_candidate candidate = *itel;

	//find candidate in PHTree
	int PHTidx = find_PHNAindex(PHArray, candidate.start_seed_id, candidate.end_seed_id);
	if( /*(!PHTidx)|| */PHTidx == 1)
	{
		return 1;
	}
	if(!PHTidx)
	{
		if(itel->p_path != nullptr)
		{
			delete[] itel->p_path;
			itel->p_path = nullptr;
		}
		if(itel->sa_ary[0] != nullptr)
		{
			free(itel->sa_ary[0]);		//double free ???
			itel->sa_ary[0] = nullptr;
		}
		if(itel->sa_ary[1] != nullptr)
		{
			free(itel->sa_ary[1]);
			itel->sa_ary[1] = nullptr;
		}
		lcand.erase(itel);
		return 0;
	}

	struct dBG_extpara ext_set(gPair[PHTidx % 2]);
	struct PH_Node PHParent = PHArray[PHTidx/2];

	ext_set.tau = PHParent.tau - candidate.tau[PHTidx % 2];		//how to deal PHParent.tau < candidate.tau  exit directly
	if(ext_set.tau > 128)
	{
		if(itel->p_path != nullptr)
		{
			delete[] itel->p_path;
			itel->p_path = nullptr;
		}
		if(itel->sa_ary[0] != nullptr)
		{
			free(itel->sa_ary[0]);		//double free ???
			itel->sa_ary[0] = nullptr;
		}
		if(itel->sa_ary[1] != nullptr)
		{
			free(itel->sa_ary[1]);
			itel->sa_ary[1] = nullptr;
		}
		lcand.erase(itel);
		return 0;
	}

	uint8_t fullseqlen = 0;
	uint8_t node_extlen = 0;
	uint8_t candpathlen = strlen(candidate.p_path);
	ext_set.orignseq = new char[gbit_para.kmer1Len/2+1]();
	if(ext_set.dir == 'I')
	{
		strncpy(ext_set.orignseq, candidate.p_path, gbit_para.kmer1Len/2);
		for(int i = PHParent.start_seed_id; i < candidate.cur_seed_id; ++i)
		{
			uint8_t seedseglen = strlen(vseed[i].seed_segment);
			fullseqlen += seedseglen;
			if(i < candidate.start_seed_id)
			{
				node_extlen += seedseglen;
			}
		}
		ext_set.alignseq = new char[fullseqlen+1]();
		for(int i = PHParent.start_seed_id; i < candidate.cur_seed_id; ++i)
		{
			strcat(ext_set.alignseq,vseed[i].seed_segment);
		}
		reverseq(ext_set.alignseq);
	}
	else
	{
		strncpy(ext_set.orignseq, candidate.p_path+candpathlen-gbit_para.kmer1Len/2, gbit_para.kmer1Len/2);
		for(int i = candidate.cur_seed_id + 1; i <= PHParent.end_seed_id; ++i)
		{
			uint8_t seedseglen = strlen(vseed[i].seed_segment);
			fullseqlen += seedseglen;
			if(i > candidate.end_seed_id)
			{
				node_extlen += seedseglen;
			}
		}
		ext_set.alignseq = new char[fullseqlen+1]();
		for(int i = candidate.cur_seed_id + 1; i <= PHParent.end_seed_id; ++i)
		{
			strcat(ext_set.alignseq,vseed[i].seed_segment);
		}
	}
	ext_set.alignlen = fullseqlen;
	struct TPTnode rootnode;
	if( (ext_set.dir == 'I' && candidate.cur_seed_id == candidate.start_seed_id) || (ext_set.dir == 'O' && candidate.cur_seed_id == candidate.end_seed_id))
	{
		init_rootnode(&rootnode, ext_set, nullptr);
	}
	else
	{
		init_rootnode(&rootnode, ext_set, &candidate);
	}
//		if(candidate.start_seed_id == startid && candidate.end_seed_id == endid)
//		{
//			cout << "===========================================================================" << endl;
//			cout << candidate.p_path << " " << (int)candidate.cur_end_pos << endl;
//			cout << "tau0, tau1:" << (int)candidate.tau[0] << "," << (int)candidate.tau[1] << " ext.tau"
//					<< (int)ext_set.tau << endl;
//			cout << "rootnode edarray:";
//			for(int i = 0; i < 2*ext_set.tau + 1; ++i)
//			{
//				cout << (int)rootnode.edarry[i] << " ";
//			}
//			cout << endl;
//			cout << "ext_set.alignlen:" << (int)ext_set.alignlen << " " << ext_set.alignseq << endl;
//			cout << "===========================================================================" << endl;
//		}


	//call ext_fromroot
	list<pair<char*, uint8_t> > lext_pair;
	ext_fromroot(ext_set, &rootnode, lext_pair);
	destory_extree(&rootnode);
	//store results
	if(lext_pair.size())
	{
		int ptroffset = ext_set.alignlen - node_extlen;
		ext_set.alignseq += ptroffset;
		ext_set.alignlen = node_extlen;
		for(auto ite = lext_pair.begin(); ite != lext_pair.end(); ++ite)
		{
			char *ext_seq = ite->first;
			uint8_t ext_len = strlen(ext_seq);
			uint8_t n_pathlen = ext_len + candpathlen;
			char *new_path = new char[n_pathlen + 1]();
			if(ext_set.dir == 'I')
			{
				strcat(new_path, ext_seq);
				strcat(new_path, itel->p_path);
				candidate.tau[0] = ite->second;
				candidate.sa_ary[1] = (uint32_t*)malloc(sizeof(uint32_t) * 2);
				candidate.sa_ary[1][0] = itel->sa_ary[1][0];
				candidate.sa_ary[1][1] = itel->sa_ary[1][1];
				candidate.p_path = new_path;
				candidate.start_seed_id = PHParent.start_seed_id;
				candidate.cur_start_pos = itel->cur_start_pos + ext_len;
				candidate.cur_end_pos = itel->cur_end_pos + ext_len;
				uint8_t leftlen = candidate.cur_end_pos+1;
				char *leftseq = new char[leftlen+1]();
				strncpy(leftseq, new_path, leftlen);
				candidate.sa_ary[0] = calc_SArangeSeq(*ext_set.pFMidx, leftseq);
				delete[] leftseq;
			}
			else
			{
				reverseq(ext_seq);
				candidate.sa_ary[0] = (uint32_t*)malloc(sizeof(uint32_t) * 2);
				candidate.sa_ary[0][0] = itel->sa_ary[0][0];
				candidate.sa_ary[0][1] = itel->sa_ary[0][1];
				candidate.tau[1] = ite->second;
				strcat(new_path, itel->p_path);
				strcat(new_path, ext_seq);
				candidate.sa_ary[1] = calc_SArangeSeq(*ext_set.pFMidx, new_path+itel->cur_start_pos);
				candidate.p_path = new_path;
				candidate.end_seed_id = PHParent.end_seed_id;
			}
			ext_set.alignseq -= ptroffset;
			lcand.push_back(candidate);
//			if(candidate.start_seed_id == startid && candidate.end_seed_id == endid)
//			{
//				cout << "===========================================================================" << endl;
//				cout << candidate.p_path << " " << (int)candidate.cur_end_pos << endl;
//				cout << "tau0, tau1:" << (int)candidate.tau[0] << "," << (int)candidate.tau[1] << " ext.tau"
//						<< (int)ext_set.tau << endl;
//				cout << "ext_set.alignlen:" << (int)ext_set.alignlen << " " << ext_set.alignseq << endl;
//				cout << "===========================================================================" << endl;
//			}
//			cout << "IN:start~cur~end:" << (int)candidate.start_seed_id << "~"
//					<< (int)candidate.cur_seed_id << "~" << (int)candidate.end_seed_id << endl;
		}
		//free memory
		for(auto ite = lext_pair.begin(); ite != lext_pair.end(); ++ite)
		{
			if(ite->first != nullptr)
			{
				delete[] ite->first;
				ite->first = nullptr;
			}
		}

	}
//	else
//	{
//		cout << "OUT:start~cur~end:" << (int)candidate.start_seed_id << "~"
//							<< (int)candidate.cur_seed_id << "~" << (int)candidate.end_seed_id << endl;
//	}
	if(itel->p_path != nullptr)
	{
		delete[] itel->p_path;
		itel->p_path = nullptr;
	}
	if(itel->sa_ary[0] != nullptr)
	{
		free(itel->sa_ary[0]);		//double free ???
		itel->sa_ary[0] = nullptr;
	}
	if(itel->sa_ary[1] != nullptr)
	{
		free(itel->sa_ary[1]);
		itel->sa_ary[1] = nullptr;
	}
	lcand.erase(itel);
	return 2;

}

void iter_candi2PHTroot(const vector<ms_seed> &vseed, const vector<PH_Node> &PHArray, list<pth_candidate> &lcand)
{
	int ret;
	for(uint32_t i = 0; i < lcand.size(); ++i)
	{
		ret = candi2candi(vseed, PHArray, lcand, i);
		if(ret == 2)
		{
			i = -1;
		}
	}
	if(!lcand.size())
	{
		cout << __FUNCTION__ << "~lcand.size() == 0\n";
	}
}

void destory_extree(struct TPTnode *pnode)
{
	if(pnode != nullptr)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(pnode->p_child[i] != nullptr)
			{
				destory_extree(pnode->p_child[i]);
				free(pnode->p_child[i]);
				pnode->p_child[i] = nullptr;
			}
		}
		if(pnode->saarry != nullptr)
		{
			free(pnode->saarry);
			pnode->saarry = nullptr;
		}
		if(pnode->edarry != nullptr)
		{
			free(pnode->edarry);
			pnode->edarry = nullptr;
		}
	}
}

//--------------------------------2021---------------------------------------//
void test_extime(const char *read_path, uint32_t ext2length, uint8_t tau)
{
	//1)get readline num and set ptr to NULL
	uint32_t readline = get_file_linenums(read_path);
	if(readline == 0)
	{
		printf("%s contains no reads, program return...\n",read_path);
		return;
	}
	char **p2_read = (char**)malloc(sizeof(char*) * readline);
	char **p2_origin = (char**)malloc(sizeof(char*) * readline);
	char **p2_align = (char**)malloc(sizeof(char*) * readline);

	uint32_t extlen = ext2length - gbit_para.kmer1Len / 2;
	//2)put every line in read file to p2_read[i]
	FILE *fp = fopen(read_path, "r");
	for(int32_t i = 0; i < readline; ++i)
	{
		size_t len = 0;
		p2_read[i] = NULL;
		getline(&p2_read[i], &len, fp);

		p2_origin[i] = (char*)malloc(sizeof(char)*(gbit_para.kmer1Len / 2 + 1));
		memset(p2_origin[i], 0, gbit_para.kmer1Len / 2 + 1);

		p2_align[i] = (char*)malloc(sizeof(char)*(extlen + 1));
		memset(p2_align[i], 0, ext2length - gbit_para.kmer1Len / 2 + 1);

		strncpy(p2_origin[i], p2_read[i], gbit_para.kmer1Len / 2);
		strncpy(p2_align[i], p2_read[i] + gbit_para.kmer1Len / 2, extlen + 1);
	}
	fclose(fp);

	struct timeval tvs,tve;
	double testime = 0;
	struct dBG_extpara ext_set(gPair[0]);
	ext_set.tau = tau;
	ext_set.alignlen = ext2length - gbit_para.kmer1Len / 2;
	struct TPTnode rootnode;
	uint32_t totalcnt = 0;
	for(int32_t i = 0; i < readline; ++i)
	{
		uint32_t cnt = 0;
		ext_set.orignseq = p2_origin[i];
		ext_set.alignseq = p2_align[i];
		gettimeofday(&tvs,NULL);
		init_rootnode(&rootnode, ext_set, nullptr);
		ext_treenode(ext_set, &rootnode, extlen, cnt);//, lpair);
		totalcnt += cnt;
		gettimeofday(&tve,NULL);
		testime += tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	}
	cout << __FUNCTION__ << ":total time : " << testime	\
		 << "\textend to " << ext2length << " total count:" << totalcnt << endl;
	for(int32_t i = 0; i < readline; ++i)
	{
		free(p2_read[i]);
		p2_read[i] = nullptr;

		free(p2_origin[i]);
		p2_origin[i] = nullptr;

		free(p2_align[i]);
		p2_align[i] = nullptr;
	}
	free(p2_read);
	free(p2_origin);
	free(p2_align);

}
