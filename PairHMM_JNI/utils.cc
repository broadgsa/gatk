#include "headers.h"
#include "template.h"

uint8_t ConvertChar::conversionTable[255];
float (*g_compute_full_prob_float)(testcase *tc, float* before_last_log) = 0;
double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log) = 0;

using namespace std;

bool is_avx_supported()
{
  int ecx = 0, edx = 0, ebx = 0;
  __asm__("cpuid"
      : "=b" (ebx),
      "=c" (ecx),
      "=d" (edx)
      : "a" (1)
      );
  return ((ecx >> 28)&1) == 1;
}

bool is_sse42_supported()
{
  int ecx = 0, edx = 0, ebx = 0;
  __asm__("cpuid"
      : "=b" (ebx),
      "=c" (ecx),
      "=d" (edx)
      : "a" (1)
      );
  return ((ecx >> 20)&1) == 1;
}

int normalize(char c)
{
	return ((int) (c - 33));
}

int read_testcase(testcase *tc, FILE* ifp)
{
	char *q, *i, *d, *c, *line = NULL;
	int _q, _i, _d, _c;
	int x, size = 0;
	ssize_t read;

	read = getline(&line, (size_t *) &size, ifp == 0 ? stdin : ifp);
	if (read == -1)
		return -1;


	tc->hap = (char *) malloc(size);
	tc->rs = (char *) malloc(size);
	q = (char *) malloc(size);
	i = (char *) malloc(size);
	d = (char *) malloc(size);
	c = (char *) malloc(size);

	if (sscanf(line, "%s %s %s %s %s %s\n", tc->hap, tc->rs, q, i, d, c) != 6)
		return -1;

	tc->haplen = strlen(tc->hap);
	tc->rslen = strlen(tc->rs);
	assert(tc->rslen < MROWS);
	tc->ihap = (int *) malloc(tc->haplen*sizeof(int));
	tc->irs = (int *) malloc(tc->rslen*sizeof(int));

	//tc->q = (int *) malloc(sizeof(int) * tc->rslen);
	//tc->i = (int *) malloc(sizeof(int) * tc->rslen);
	//tc->d = (int *) malloc(sizeof(int) * tc->rslen);
	//tc->c = (int *) malloc(sizeof(int) * tc->rslen);

	for (x = 0; x < tc->rslen; x++)
	{
		_q = normalize(q[x]);
		_i = normalize(i[x]);
		_d = normalize(d[x]);
		_c = normalize(c[x]);
		tc->q[x] = (_q < 6) ? 6 : _q;
		tc->i[x] = _i;
		tc->d[x] = _d;
		tc->c[x] = _c;
		tc->irs[x] = tc->rs[x];
	}
	for (x = 0; x < tc->haplen; x++)
	  tc->ihap[x] = tc->hap[x];


	free(q);
	free(i);
	free(d);
	free(c);
	free(line);

	return 0;
}

unsigned MAX_LINE_LENGTH = 65536;
int convToInt(std::string s)
{
  int i;
  std::istringstream strin(s);
  strin >> i;
  return i;
}

void tokenize(std::ifstream& fptr, std::vector<std::string>& tokens)
{
  int i = 0;
  std::string tmp;
  std::vector<std::string> myVec;
  vector<char> line;
  line.clear();
  line.resize(MAX_LINE_LENGTH);
  vector<char> tmpline;
  tmpline.clear();
  tmpline.resize(MAX_LINE_LENGTH);
  myVec.clear();

  while(!fptr.eof())
  {
    i = 0;
    bool still_read_line = true;
    unsigned line_position = 0;
    while(still_read_line)
    {
      fptr.getline(&(tmpline[0]), MAX_LINE_LENGTH);
      if(line_position + MAX_LINE_LENGTH > line.size())
	line.resize(2*line.size());
      for(unsigned i=0;i<MAX_LINE_LENGTH && tmpline[i] != '\0';++i,++line_position)
	line[line_position] = tmpline[i];
      if(fptr.eof() || !fptr.fail()) 
      {
	still_read_line = false;
	line[line_position++] = '\0';
      }
    }
    std::istringstream kap(&(line[0]));

    while(!kap.eof())
    {
      kap >> std::skipws >> tmp;
      if(tmp != "")
      {
	myVec.push_back(tmp);
	++i;
	//std::cout <<tmp <<"#";
      }
      tmp = "";
    }
    //std::cout << "\n";
    if(myVec.size() > 0)
      break;
  }
  tokens.clear();
  //std::cout << "Why "<<myVec.size()<<"\n";
  tokens.resize(myVec.size());
  for(i=0;i<(int)myVec.size();++i)
    tokens[i] = myVec[i];
  line.clear();
  tmpline.clear();
}

int read_mod_testcase(ifstream& fptr, testcase* tc, bool reformat)
{
  static bool first_call = true;
  vector<string> tokens;
  tokens.clear();
  tokenize(fptr, tokens);
  if(tokens.size() == 0)
    return -1;
  tc->hap = new char[tokens[0].size()+2];
  tc->haplen = tokens[0].size();
  memcpy(tc->hap, tokens[0].c_str(), tokens[0].size());
  tc->rs = new char[tokens[1].size()+2];
  tc->rslen = tokens[1].size();
  //cout << "Lengths "<<tc->haplen <<" "<<tc->rslen<<"\n";
  memcpy(tc->rs, tokens[1].c_str(),tokens[1].size());
  assert(tokens.size() == 2 + 4*(tc->rslen));
  assert(tc->rslen < MROWS);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->q[j] = convToInt(tokens[2+0*tc->rslen+j]);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->i[j] = convToInt(tokens[2+1*tc->rslen+j]);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->d[j] = convToInt(tokens[2+2*tc->rslen+j]);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->c[j] = convToInt(tokens[2+3*tc->rslen+j]);
 
  if(reformat)
  {
    ofstream ofptr;
    ofptr.open("reformat/debug_dump.txt",first_call ? ios::out : ios::app);
    assert(ofptr.is_open());
    ofptr << tokens[0] << " ";
    ofptr << tokens[1] << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->q[j]+33));
    ofptr << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->i[j]+33));
    ofptr << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->d[j]+33));
    ofptr << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->c[j]+33));
    ofptr << " 0 false\n";

    ofptr.close();
    first_call = false;
  }


  return tokens.size();
}

void debug_dump(string filename, string s, bool to_append, bool add_newline)
{
  ofstream fptr;
  fptr.open(filename.c_str(), to_append ? ofstream::app : ofstream::out);
  fptr << s;
  if(add_newline)
    fptr << "\n";
  fptr.close();
}
