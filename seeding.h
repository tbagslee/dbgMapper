/*
 * seeding.h
 *
 *  Created on: Sep 18, 2019
 *      Author: bio
 */

#ifndef SEEDING_H_
#define SEEDING_H_

#include "basic.h"
#include "FMindex_ExactMatch.h"
#include "load_DBG_full.h"
#include <utility>
#include <list>
#define nullptr NULL

extern int startid, endid;
extern std::pair<char, const struct sFMindex *> gPair[2];

struct ms_seed{
	uint8_t start_pos;
	uint8_t end_pos;
	char *seed_segment;
};

struct candiPath{
	char *p_candidatepath;
	uint8_t start_seed_id;
	uint8_t end_seed_id;
	uint8_t *sa_ary;

};

struct PH_Node{
	uint8_t start_seed_id;
	uint8_t end_seed_id;
	uint8_t tau;
};

struct pth_candidate
{
	uint8_t start_seed_id;
	uint8_t cur_seed_id;
	uint8_t end_seed_id;
	uint8_t cur_start_pos;
	uint8_t cur_end_pos;
	uint8_t tau[2];
	char *p_path;
	uint32_t *sa_ary[2];
	pth_candidate() : start_seed_id(0), cur_seed_id(0), end_seed_id(0), cur_start_pos(0), cur_end_pos(0), tau{0,0}, p_path(nullptr), sa_ary{nullptr, nullptr}{}
	pth_candidate(uint32_t p) : start_seed_id(p), cur_seed_id(p), end_seed_id(p), cur_start_pos(0), cur_end_pos(0), tau{0,0}, p_path(nullptr), sa_ary{nullptr, nullptr}{}
};

struct pos_candidate
{
	uint8_t seed_id;
	uint32_t ref_pos;
	pos_candidate() : seed_id(0), ref_pos(0) {};
	pos_candidate(uint8_t s) : seed_id(s), ref_pos(0) {};
};


struct ms_result
{
	uint32_t ref_start;
	uint32_t ref_end;
	uint32_t ed;
};

struct dBG_extpara
{	//dBG extension parameter
	bool onunipath;		//whether extending on unipath
	char dir;			//direction of extending
	uint8_t tau;		//ed threshold of extending
	uint8_t alignlen;	//length of alignseq
	char *orignseq;		//kmer to be located in dBG
	char *alignseq;		//aligning sequence during extending
	const struct sFMindex *pFMidx;	//according to dir pFMidx select normal FMindex or reverse FMindex
	dBG_extpara(std::pair<char, const struct sFMindex *> &qpr) : onunipath(false), dir(qpr.first), pFMidx(qpr.second) {}
};

struct TPTnode
{	//TPTnode is node of the tree which generated when kmer extending in dBG
	char c;							//character in sequence including {A,C,G,T}
	uint8_t level;					//level of node i.e. depth of tree
	uint16_t offset;				//offset of the kmer in unipath
	struct TPTnode * p_parent;		//parent of node
	struct TPTnode * p_child[4];	//child of node branched node have more than one children while node in unipath only one
	uint8_t * edarry;				//edarray of node length : 2tau + 1  including ed of specify node nearby
	uint32_t * saarry;				//appearance frequency of sequence from root to current node
};

struct seed_segment
{
	char * seedseq;
	uint32_t * saarry;
//	uint32_t *
};


void generate_PHNArray(vector<PH_Node> &PHArray, uint8_t tau);
uint32_t find_PHNAindex(const vector<PH_Node> &PHArray, uint8_t start_id, uint8_t end_id);
void generate_seeds(const char *read, uint8_t tau, vector<ms_seed> & vseed);
void free_seeds(vector<ms_seed> &vseed);
void free_pthcandi(list<pth_candidate> &lcand);
uint32_t init_candidate(const vector<ms_seed> &vseed, const sFMindex &fmindx, list<pth_candidate> &lcand);
int candi2candi(const vector<ms_seed> &vseed, const vector<PH_Node> &PHArray, list<pth_candidate> &vcand, uint32_t i);
void iter_candi2PHTroot(const vector<ms_seed> &vseed, const vector<PH_Node> &PHArray, list<pth_candidate> &lcand);
void init_rootnode(struct TPTnode *pnode,const struct dBG_extpara &ext_para, struct pth_candidate *pcandi = nullptr);
void ext_treenode(struct dBG_extpara &ext_para, struct TPTnode *pnode, uint32_t extlen, uint32_t &cnt);
void ext_fromroot(struct dBG_extpara &ext_para, struct TPTnode *proot, list<pair<char*, uint8_t> > &lpair);
void destory_extree(struct TPTnode *pnode);
void test_extime(const char *read_path, uint32_t ext2length, uint8_t tau);
#endif /* SEEDING_H_ */
