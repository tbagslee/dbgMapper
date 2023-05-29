/*
 * method.cpp
 *
 *  Created on: Feb 16, 2019
 *      Author: bio
 */
#include "method.h"

uint32_t min_lastlineEd(char *a,char * b,int32_t x,int32_t y, uint32_t &tau)  //b is ref   1。read 长度不同时间   2.逻辑是否有问题  3.candidate数量
{
	uint32_t r;
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
	}

	for(int32_t i=0;i<x+1;++i)
	{
		m[i][0]=i;
	}
	for(int32_t i=0;i<y+1;++i)
	{
		m[0][i]=i;
	}

	for(int32_t i=1;i<x+1;++i)
	{
		for(int32_t j=1;j<y+1;++j)
		{
			int32_t tmp=0;
			if(a[i-1]!=b[j-1])
			{
				tmp=1;
			}
			m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
			m[i][j]=min(m[i][j],m[i-1][j]+1);
		}
	}
	r = m[x][0];
	uint32_t col = 0;
	for(int32_t i=1;i<y+1;++i)
	{
		if(m[x][i] <= r)
		{
			r = m[x][i];
			col = i;
		}
	}
	for(int32_t i=0;i<x+1;++i)
	{
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
	if(r <= tau)
	{
		tau -= r;
		return col;
	}
	else
	{
		return -1;
	}
}


uint32_t **originalEd(char *a,char * b,int32_t x,int32_t y)
{
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
	}

	for(int32_t i=0;i<x+1;i++)
	{
		m[i][0]=i;
	}
	for(int32_t i=0;i<y+1;i++)
	{
		m[0][i]=i;
	}

	for(int32_t i=1;i<x+1;i++)
	{
		for(int32_t j=1;j<y+1;j++)
		{
			int32_t tmp=0;
			if(a[i-1]!=b[j-1])
			{
				tmp=1;
			}
			m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
			m[i][j]=min(m[i][j],m[i-1][j]+1);
		}
	}
	return m;
}

void free_edmatrix(uint32_t ***m, int32_t x)
{
	if(*m != NULL)
	{
		for(int32_t i=0;i<x+1;i++)
		{
			free((*m)[i]);
			(*m)[i] = NULL;
		}
		free(*m);
		*m = NULL;
	}
}

uint32_t rangeEd(char *a,char * b,int32_t x,int32_t y,int32_t ed)
{
	if((x-y > ed) || (x-y < -ed))
	{
		return (uint32_t)ed+1;
	}
	uint32_t r;
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
	}

	for(int32_t i=0;i<x+1;i++)
	{
		m[i][0]=i;
	}
	for(int32_t i=0;i<y+1;i++)
	{
		m[0][i]=i;
	}

	int32_t y_start;
	int32_t y_end;
	int32_t mid_ed=ed/2;
	int32_t l_bound;
	int32_t r_bound;
	int32_t y1=1;
	if(x<=y)
	{
		l_bound=mid_ed;
		r_bound=(y-x)+mid_ed;
	}
	else
	{
		l_bound=(x-y)+mid_ed;
		r_bound=mid_ed;
	}
	for(int32_t i=1;i<x+1;i++)
	{
		y_start=max(y1,i-l_bound);
		y_end=min(i+r_bound,y);
		int32_t line_check = ed+1;
		for(int32_t j=y_start;j<=y_end;j++)
		{
			uint32_t tmp=0;
			if(j==i-l_bound && j==i+r_bound)
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=m[i-1][j-1]+tmp;
			}
			else if(j==i-l_bound)
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=min(m[i-1][j]+1,m[i-1][j-1]+tmp);
			}
			else if(j==i+r_bound)
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
			}
			else
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
				m[i][j]=min(m[i][j],m[i-1][j]+1);
			}
			if(m[i][j] < line_check)
			{
				line_check = m[i][j];
			}
		}
		if(line_check == (ed+1))
		{
			for(int32_t i=0;i<x+1;i++)
			{
				free(m[i]);
			}
			free(m);
			return ed+1;
		}
	}
	r = m[x][y];
	for(int32_t i=0;i<x+1;i++)
	{
		free(m[i]);
	}
	free(m);
	return r;
}
