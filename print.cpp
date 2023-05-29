/*
 * mapperStruct.cpp
 *
 *  Created on: May 9, 2020
 *      Author: pluto
 */

#include "print.h"

void printbytevector(uint8_t *data, int length) {
	int i;
	for (i = 0; i < length; i++) {
		int m;
		for (m = 0; m < 8; m++) {
			if (data[i] & (1ULL << m))
				printf("1");
			else
				printf("0");
		}
	}
}

void print128_bit(__m128i var) {
	uint8_t *val = (uint8_t*) &var;
	printbytevector(val, 16);

	printf("\n");

}

void onlyprint_upward(struct TPTnode *pnode)
{
	string str;
	//add char to string
	do{
		str += pnode->c;
		pnode = pnode->p_parent;
	}while(pnode != nullptr);

	//reverse string
	reverse(str.begin(), str.end());
	cout << str << endl;
}

void print_extree(const struct TPTnode &node,char *seq)
{
	seq[strlen(seq)] = node.c;
	if(node.p_child[0] == nullptr && node.p_child[1] == nullptr && node.p_child[2] == nullptr && node.p_child[3] == nullptr)
	{
		printf("%s\n",seq);
	}
	else
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(node.p_child[i] != nullptr)
			{
				print_extree(*node.p_child[i],seq);
			}
		}
	}
	seq[strlen(seq)-1] = '\0';
}


void print_seed(const vector<ms_seed> &vseed)
{
	for(uint32_t i = 0; i != vseed.size(); ++i)
	{
		cout << +static_cast<uint8_t>(vseed[i].start_pos) << ":" << +static_cast<uint8_t>(vseed[i].end_pos) << ":" << vseed[i].seed_segment << endl;
	}
}

void print_candidate(const list<pth_candidate> & lcand)
{
	for(auto it = lcand.begin(); it != lcand.end(); ++it)
	{
		cout << "cur_seed_id:" << (uint32_t)it->cur_seed_id << ":" << it->p_path << endl;
	}
}

void print_poscand(const list<pos_candidate> & lcand)
{
	for(auto it = lcand.begin(); it != lcand.end(); ++it)
	{
		cout << "seed_id:" << (uint32_t)it->seed_id << " refpos:" << it->ref_pos << endl;
	}
}

void print_result(const list<ms_result> & lrslt)
{
	list<ms_result>::const_iterator ite = lrslt.begin();
	for(; ite != lrslt.end(); ++ite)
	{
		cout << ite->ref_start << "~" << ite->ref_end << " ed:"<< ite->ed << endl;
	}
}

