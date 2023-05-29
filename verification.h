/*
 * verification.h
 *
 *  Created on: May 25, 2020
 *      Author: pluto
 */

#ifndef VERIFICATION_H_
#define VERIFICATION_H_


#include "basic.h"
#include "method.h"
#include "print.h"
#include <list>

void pathCan2posCan(list<pth_candidate> &lcand, const sFMindex &fmindx, list<pos_candidate> &lposCan);
void verification(const vector<ms_seed> &vseed, const list<pos_candidate> &lposCan, list<ms_result> &lrslt,
				const char *read, const char *ref, uint32_t tau);
void get_result(const char *read_path, const char *ref, const sFMindex &fmindx, uint8_t tau, uint32_t threhld);

#endif /* VERIFICATION_H_ */
