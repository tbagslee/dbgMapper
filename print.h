/*
 * mapperStruct.h
 *
 *  Created on: May 9, 2020
 *      Author: pluto
 */

#ifndef PRINT_H_
#define PRINT_H_

#include "basic.h"
#include "seeding.h"
#include "FMindex_ExactMatch.h"
#include <tmmintrin.h>
#include <list>

void print128_bit(__m128i var);
void print_seed(const vector<ms_seed> & vseed);
void print_candidate(const list<pth_candidate> & lcand);
void print_poscand(const list<pos_candidate> & lcand);
void print_result(const list<ms_result> & lrslt);
void onlyprint_upward(struct TPTnode *pnode);
void print_extree(const struct TPTnode &node,char *seq);

#endif /* PRINT_H_ */
