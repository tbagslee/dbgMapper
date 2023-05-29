/*
 * prealignment.h
 *
 *  Created on: May 25, 2020
 *      Author: pluto
 */

#ifndef PREALIGNMENT_H_
#define PREALIGNMENT_H_

#include "basic.h"
#include <xmmintrin.h>
#include <tmmintrin.h>
#include <emmintrin.h>
#include <nmmintrin.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration.hpp>
#include <boost/preprocessor/arithmetic.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

#ifndef __aligned__
	#define __aligned__ __attribute__((aligned(16)))
#endif

#define SSE_BIT_LENGTH		128
#define SSE_BYTE_NUM		BOOST_PP_DIV(SSE_BIT_LENGTH, 8)
#define BASE_SIZE1			1
#define SSE_BASE_NUM1		BOOST_PP_DIV(SSE_BIT_LENGTH, BASE_SIZE1)

int banded_edit_distance(const char *pattern, const char *text, int text_length, const int error_threshold,int *mapping_end_position);
void sse3_convert2bit1(char *str, uint8_t *bits0, uint8_t *bits1);
int pre_alignment(char* read, char* ref, int start, int end, int max_error);
extern uint8_t read_vec0_t[SSE_BYTE_NUM] __aligned__;
extern uint8_t read_vec1_t[SSE_BYTE_NUM] __aligned__;


#endif /* PREALIGNMENT_H_ */
