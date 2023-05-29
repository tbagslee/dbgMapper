#include "basic.h"
#include "method.h"
#define nullptr NULL

struct bit256KmerPara gbit_para;

void reverseq(char *seq)
{
	char *front = seq;
	char *rear = seq + strlen(seq) - 1;
	char tmp;
	while(front < rear)
	{
		tmp = *front;
		*front = *rear;
		*rear = tmp;
		front++;
		rear--;
	}
}

void kmercpy(uint64_t *dest, const uint64_t * from, uint32_t unitper)
{
	if(dest == NULL || from == NULL)
	{
		printf("kmercpy error\n");
		exit(1);
	}
	for(uint32_t i = 0; i < unitper; i++)
	{
		dest[i] = from[i];
	}
}

void get_para(struct bit256KmerPara *para1,uint32_t kmer_length)
{
	para1->kmer1Len = kmer_length*2;
	para1->remainer1to64 = para1->kmer1Len%64;
	para1->kmer64Len = para1->kmer1Len/64+(para1->remainer1to64?1:0);
	if(para1->remainer1to64==0)
	{
		para1->codefor1 = 0;
	}
	else
	{
		para1->codefor1 = 0;
		for(uint32_t i=0;i<para1->remainer1to64-1;i++)
		{
			para1->codefor1 = para1->codefor1|1;
			para1->codefor1 = para1->codefor1<<1;
		}
		para1->codefor1 = para1->codefor1|1;
	}
}


void fill_char_with_four_char(uint8_t* current,char*p)
{
	for(uint32_t i=0;i<3;i++)
	{
		switch(p[i])
		{
			case 'A':
				*current=*current<<2;
				break;
			case 'C':
				*current=*current|1;
				*current=*current<<2;
				break;
			case 'G':
				*current=*current|2;
				*current=*current<<2;
				break;
			case 'T':
				*current=*current|3;
				*current=*current<<2;
				break;
			default:
				*current=*current<<2;
				break;
		}
	}
	switch(p[3])
	{
		case 'A':
			break;
		case 'C':
			*current=*current|1;
			break;
		case 'G':
			*current=*current|2;
			break;
		case 'T':
			*current=*current|3;
			break;
		default:
			break;
	}
}
void ReadSeq_bit(uint8_t **seq1,uint64_t *seq_length,char* p_ref)
{
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;

	maxsize=pow(2,20);
	addsize=pow(2,20);

	uint8_t *seq;//=seq1;
	seq=(uint8_t*) malloc (sizeof(uint8_t)*maxsize);
	cursize=maxsize;

	uint32_t buffer_size=_BufferSize_;
	char buffer_line[_BufferSize_];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint64_t len=0;
	uint32_t cycle=0;
	uint32_t number_left=0;
	uint32_t number_buffer_char;
	char array_left[4];
	while (1)
	{
		if(fgets(buffer_line,buffer_size-number_left,fp)==NULL)
		{
			break;
		}
		if(buffer_line[0]=='>')
		{
			continue;
		}
		else
		{
			char buffer_line_tmp[_BufferSize_];
			if(number_left!=0)
			{
				strcpy(buffer_line_tmp,array_left);
				strcpy(buffer_line_tmp+number_left,buffer_line);
			}
			else
			{
				strcpy(buffer_line_tmp,buffer_line);
			}

			uint32_t buffer_line_char_number=strlen(buffer_line_tmp);
			if(buffer_line_tmp[buffer_line_char_number-1]=='\n')
			{
				buffer_line_char_number--;
			}
			number_left=buffer_line_char_number%4;
			number_buffer_char=buffer_line_char_number/4;
			if(number_left!=0)
			{
				for(uint32_t i=0;i<number_left;i++)
				{
					array_left[i]=buffer_line_tmp[buffer_line_char_number-(number_left-i)];
				}
				array_left[number_left]='\0';
			}

			if(len+number_buffer_char<cursize)
			{
				for(uint32_t i=0;i<number_buffer_char;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
			}
			else
			{
				seq=(uint8_t*) realloc (seq,sizeof(uint8_t)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_line_char_number/4;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
				cout <<"add 1024*1024 byte for seq: " << cycle++ <<endl;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len*4+number_left;
	if(number_left!=0)
	{
		for(uint32_t i=number_left;i<4;i++)
		array_left[i]='A';
		uint8_t tmp;
		fill_char_with_four_char(&tmp,array_left);
		seq[len]=tmp;
		len++;
	}

	seq1[0]=seq;

	uint8_t* seq_shift1=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift1[i]=seq[i]<<2;
		seq_shift1[i]=seq_shift1[i]|(seq[i+1]>>6);
	}
	seq_shift1[len-1]=seq[len-1]<<2;
	seq1[1]=seq_shift1;

	uint8_t* seq_shift2=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift2[i]=seq_shift1[i]<<2;
		seq_shift2[i]=seq_shift2[i]|(seq_shift1[i+1]>>6);
	}
	seq_shift2[len-1]=seq_shift1[len-1]<<2;
	seq1[2]=seq_shift2;

	uint8_t* seq_shift3=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift3[i]=seq_shift2[i]<<2;
		seq_shift3[i]=seq_shift3[i]|(seq_shift2[i+1]>>6);
	}
	seq_shift3[len-1]=seq_shift2[len-1]<<2;
	seq1[3]=seq_shift3;

	cout << "the length of seq is: " << *seq_length << endl;
}

void ReadSeq(char **seq1,uint32_t *seq_length,const char* p_ref)
{
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;

	maxsize=pow(2,20);
	addsize=pow(2,20);

	char *seq;
	seq=(char*) malloc (sizeof(char)*maxsize);
	cursize=maxsize;

	uint32_t buffer_size=256;
	char buffer_line[256] = {0};

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint32_t len=0;
	while (fgets(buffer_line,buffer_size,fp)!=NULL)
	{
		if(buffer_line[0]=='>')
			continue;
		else
		{
			if(len+buffer_size<cursize)
			{
				for(uint32_t i=0;i<buffer_size;i++)
				{
					if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
					{
						break;
					}
					if(buffer_line[i]>='a')
					{
						buffer_line[i] -= 32;
					}
					if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						continue;
					}
					seq[len]=buffer_line[i];
					len++;
				}
			}
			else
			{
				seq=(char*) realloc (seq,sizeof(char)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_size;i++)
				{
					if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
					{
						break;
					}
					if(buffer_line[i]>='a')
					{
						buffer_line[i] -= 32;
					}
					if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						continue;
					}
					seq[len]=buffer_line[i];
					len++;
				}
			}
		}
		memset(buffer_line,0,buffer_size*sizeof(char));
	}
	*seq_length=len;
	*seq1=seq;
//	cout << "the length of seq is: " << len << endl;
}


char *RefSaveAsSeq(const char *p_ref)
{
	char *ref = nullptr;
	uint32_t ref_len = 0;
	ReadSeq(&ref, &ref_len, p_ref);
	cout << "ReadSeq done..." << endl;

	char ref_name[16] = {0};
	const char *p_name = p_ref;
	while(p_name != p_ref + strlen(p_ref))
	{
		if(*p_name == '.')
		{
			break;
		}
		strncat(ref_name,p_name,1);
		p_name++;
	}

	char *RefOutName = new char[strlen("RefSave")+strlen(ref_name)+2]();
	sprintf(RefOutName,"%s_%s","RefSave",ref_name);

	cout << "Ref fwriting..." << endl;
	FILE* RefFile;
	RefFile = fopen(RefOutName,"w");
	fwrite(ref, sizeof(char)*ref_len, 1, RefFile);
	fclose(RefFile);

	cout << "Ref fwrite finished!" << endl;
	return RefOutName;
}

void SubstitudeOtherChar(const char* p_ref, uint32_t line_len)
{
	//open a .fa file and substitude the character who is not included by {A,C,G,T},jump across the line start with '>'
	string line;
	string strsave;
	ifstream filein(p_ref);
	ofstream outfile;
	char outname[32] = {0};
	sprintf(outname, "sub_%s", p_ref);
	outfile.open(outname, ios::out);

	if(filein)
	{
		while(getline(filein, line))
		{
			if(line[0] == '>')
			{
				if(strsave.size())
				{
					uint32_t looptime = strsave.size() / line_len;
					if(strsave.size() % line_len)
					{
						looptime++;
					}
					uint32_t strinx = 0;
					for(uint32_t i = 0; i < looptime; ++i)
					{
						string segment(strsave, strinx, 50);
						outfile << segment << endl;
						strinx += 50;
					}
					strsave.clear();
				}
				outfile << line << endl;
			}
			else
			{
				for(size_t i = 0; i < line.size(); ++i)
				{
					if(line[i] >= 'a')
					{
						line[i] -= 32;
					}
					if(line[i] != 'A' && line[i] != 'C' && line[i] != 'G' && line[i] != 'T')
					{
						continue;
					}
					strsave += line[i];
				}
			}
		}
		outfile.close();
	}
	else
	{
		cout << __FUNCTION__ << " open " << p_ref << " failed..." << endl;
	}

}

uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len)
{
	uint32_t r=2;
	for(uint32_t i=0;i<len;i++)
	{
		if(a[i]<b[i])
		{
			r=0;
			break;
		}
		else
		{
			if(a[i]>b[i])
			{
				r=1;
				break;
			}
		}
	}
	return r;
}
void cal_hash_value_directly_256bit(char *seq,uint64_t * current,\
		struct bit256KmerPara para)
{
	uint64_t tmp;
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len;i++)
	{
		char* loop_tmp=k_mer_temp+32*i;
		if(i==para.kmer64Len-1&&para.remainer1to64!=0)
		{
			tmp=0;
			for(uint32_t j=0;j<para.remainer1to64/2-1;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[para.remainer1to64/2-1])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
		else
		{
			tmp=0;
			for(uint32_t j=0;j<31;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[31])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
	}
}
void cal_hash_value_indirectly_256bit(char *seq,uint64_t* original,uint64_t* current,\
		struct bit256KmerPara para)
{
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len-1;i++)
	{
		if(i==para.kmer64Len-2&&para.remainer1to64!=0)
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>(para.remainer1to64-2));
		}
		else
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>62);
		}
	}
	if(para.remainer1to64==0)
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
	else
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		current[para.kmer64Len-1]=current[para.kmer64Len-1]&para.codefor1;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
}

uint32_t get_file_linenums(const char *filepath)
{
	char c,lc;
	int32_t num = 0;
	FILE *fp = fopen(filepath, "r");
	if(fp == NULL)
	{
		return 0;
	}
    while((c = fgetc(fp)) != EOF) //read each character
    {
        if(c == '\n') //accumulate line num
        {
        	num++;
        }
        lc = c; 	//save last character
    }
    if(lc != '\n')	//deal last line
    {
    	num++;
    }
	fclose(fp);
	return num;
}


char replace_rule(char a)
{
	switch(a)
	{
		case 'A':
			return 'C';
			break;
		case 'C':
			return 'G';
			break;
		case 'G':
			return 'T';
			break;
		case 'T':
			return 'A';
			break;
		default:
			return '\0';
			break;
	}
}

char rand_generate_char()
{
	switch(rand()%4)
	{
		case 0:
			return 'C';
			break;
		case 1:
			return 'G';
			break;
		case 2:
			return 'T';
			break;
		case 3:
			return 'A';
			break;
		default:
			return '\0';
	}
}

void ed_operate(char * line,int len)   //计划将函数修改为   void ed_operate(char * line,int len, int &pos, int &type)  其中pos为编辑操作位置  type为标记操作类型
{
	int pos_ed=rand()%len;
	int type_ed=rand()%3; //0:sub;1:del;2:insert

	switch(type_ed)
	{
		case 0:
			line[pos_ed]=replace_rule(line[pos_ed]);
			break;
		case 1:
			for(int k=pos_ed;k<len-1;k++)
			{
				line[k]=line[k+1];
			}
			line[len-1]='\0';
			len--;
			break;
		case 2:
			for(int k=len-1;k>=pos_ed;k--)
			{
				line[k+1]=line[k];
			}
			line[pos_ed]=rand_generate_char();
			len++;
			line[len]='\0';
			break;
		default:
			break;
	}
}

void gen_ed_readfile(const char *filepath,uint8_t edstart, uint32_t line)
{
	/*
	 * filepath is path of input file
	 * edstart indicate that generate ed from edstart to edstart + 5
	 * line indicate the read number of each ed
	 * */
	const int READLENMAX = 128;
	const int READLEN = 40;

	char *ref_seq = nullptr;
	uint32_t reflength = 0;

	char filetype[8] = {0};
	const char *p_ref = filepath;
	while(p_ref != filepath + strlen(filepath))
	{
		if(*p_ref == '.')
		{
			strcpy(filetype, p_ref);
		}
		p_ref++;
	}

	if(strcmp(filetype, ".fa") == 0 || strcmp(filetype, ".fasta") == 0)
	{
		ReadSeq(&ref_seq, &reflength, filepath);
		cout << "generate the read seq from primitive file!" << endl;
	}
	else
	{
		FILE* ref_file;
		ref_file = fopen(filepath, "r");
		if(ref_file == NULL)
		{
			cout << __FUNCTION__ << ":" << filepath << "open failed!" << endl;
			exit(-1);
		}
		fseek(ref_file,0,2);

		// get the length of reference
		reflength = ftell(ref_file)/(sizeof(char));
		fseek(ref_file,0,0);

		// get reference seq
		ref_seq = (char*)malloc(sizeof(char)*(reflength + 1));
		memset(ref_seq, 0, reflength + 1);
		size_t readcnt = fread(ref_seq,sizeof(char),reflength,ref_file);

		if(readcnt != reflength)
		{
			cout << __FUNCTION__ <<":fread block is not adequate!" << endl;
		}
	}

    ofstream outfile;
    char outname[64] = {0};
    strcat(outname,"readof_");
	const char *p_name = filepath;
	while(p_name != filepath + strlen(filepath))
	{
		if(*p_name == '.')
		{
			break;
		}
		strncat(outname,p_name,1);
		p_name++;
	}
    outfile.open(outname,ios::out);

    //two char array to save read
    char linea[128] = {0};
    char lineb[128] = {0};

    uint32_t pos;
//    for(uint8_t ed = edstart; ed < edstart+5; ++ed)
//    {
		for(uint32_t i = 0; i < line; ++i)
		{
			pos = rand()%(reflength - READLENMAX);
			strncpy(linea, pos + ref_seq, READLEN);
			strcpy(lineb, linea);
//			while(rangeEd(linea,lineb,strlen(linea),strlen(lineb),ed+1) != ed)
//			{
//				ed_operate(lineb,strlen(lineb));
//			}
			outfile << lineb << endl;
		}
//    }
    outfile.close();
    free(ref_seq);
}
