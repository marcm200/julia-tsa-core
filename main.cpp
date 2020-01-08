/*

	My (Marc Meidlinger, July-October 2019)	implementation of 
	the ingenious trustworthy Julia set algorithm from the article:
	
	"Images of Julia sets that you can trust" by Luiz Henrique de Figueiredo, Diego Nehab,
	Jorge Stolfi, Joao Batista Oliveira from 2013
	
	"Rigorous bounds for polynomial Julia sets" by Luiz Henrique de Figueiredo, Diego Nehab,
	Jorge Stolfi, Joao Batista Oliveira 

*/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include "time.h"

// used floating type
// comment out or in what is needed

//#define _DOUBLE
//#define _LONGDOUBLE
//#define _QUADMATH
#define _FPA
//#define _F107


typedef signed long long VLONG;
typedef uint8_t BYTE;
typedef unsigned long long UVLONG;
typedef uint32_t DDBYTE;
typedef DDBYTE *PDDBYTE;
typedef uint16_t DBYTE;
typedef uint16_t *PDBYTE;

// chunk size depe4ndend on operating system
// win32 => 512 MB
//#define _CHUNK512

#ifdef _CHUNK512
const uint64_t CHUNKSIZE=( (uint64_t)1 << 27 );
#else
// win64 => 1 GB
const uint64_t CHUNKSIZE=( (uint64_t)1 << 30 );
#endif

int32_t SAVE128KIMAGE=0;

#ifdef _QUADMATH
#include "quadmath.h"
typedef __float128 NTYP;
typedef __float128 *Pfloat128;
const char NNTYPSTR[]="f128_";
const char NTS[]="QD";
#endif

#ifdef  _LONGDOUBLE
typedef long double NTYP;
const char NNTYPSTR[]="ld_";
const char NTS[]="LD";
#endif

#ifdef _DOUBLE
typedef double NTYP;
const char NNTYPSTR[]="";
const char NTS[]="D";
#endif

#ifdef _F107
// https://mrob.com/pub/math/f161.html#source
#include "f107_o.cpp"
typedef f107_o NTYP;
typedef f107_o *Pf107;
inline void	minimaxF107AB(f107_o&,f107_o&,const f107_o,const f107_o);
inline void	minimaxF107ABCD(f107_o&,f107_o&,const f107_o,const f107_o,const f107_o,const f107_o);
inline void	minimaxF107AB(Pf107&,Pf107&,f107_o*,f107_o*);
inline void	minimaxF107ABCD(Pf107&,Pf107&,f107_o*,f107_o*,f107_o*,f107_o*);
const char NNTYPSTR[]="f107_";
const char NTS[]="F1";
f107_o operator*(int,const f107_o);
#endif

const int MAXZAEHLER=16;
const VLONG ZWEIHOCH32=( (UVLONG)1 << 32 );
const UVLONG MAXDDBYTE=ZWEIHOCH32-1;
const VLONG INT32MAX=(((int64_t)1 << 32) - 1);
const UVLONG MAXREFPOINTS=( (uint64_t)1 << 31 );

#ifdef _FPA
struct FPA {
	// UVLONG or DDBYTE not performance-relevant
	int8_t vorz; // -1,0,1
	DDBYTE a,b,c,d; // a + b*2^-32 + c*2^-64 + d*2^-96
	
	FPA();
	FPA(const double);
	FPA(const VLONG);
	
	void set_vlong(const VLONG);
	void setNull(void);
	void squareTo(FPA&);
	void checkNull(void);
	void copyFrom(const FPA);
	void set_double32(const double);
	char* str(char*);
	double convert_to_double(void);

	void shiftLeft(const int32_t); // Multiplication
	FPA& operator=(const FPA);
};

typedef FPA NTYP;
typedef FPA *PFPA;
const char NNTYPSTR[]="fpa_";
const char NTS[]="FP";

bool operator<(const FPA,const FPA);
bool operator>(const FPA,const FPA);
bool operator!=(const FPA,const FPA);
bool operator==(const FPA,const FPA);
FPA operator-(const FPA,const FPA);
FPA operator+(const FPA,const FPA);
FPA operator*(const FPA,const FPA);
FPA operator*(const int,const FPA);
inline void FPA_sub_abs_ZAB(FPA&,const FPA,const FPA);
inline void FPA_add_abs_ZAB(FPA&,const FPA,const FPA);
// const FPA beibehalten, dann kann a=a+b dort
// ausgeführt werden
inline void FPA_add_ZAB(FPA&,const FPA,const FPA);
inline void FPA_sub_ZAB(FPA&,const FPA,const FPA);
inline void FPA_sub2tpos_ZAB(FPA&,const FPA,const FPA);
inline void FPA_sub_abs_ovgl_ZAB(FPA&,const FPA,const FPA);
inline int FPA_vgl_abs(const FPA,const FPA);
inline int FPA_vgl(const FPA,const FPA);
inline int FPA_vgl_abs(PFPA,PFPA);
inline int FPA_vgl(PFPA,PFPA);
inline void FPA_mul_ZAuvlong(FPA&,const FPA,const UVLONG);
inline void FPA_mul_ZAB(FPA&,const FPA,const FPA);
inline void minimaxFPAMIMAAB(FPA&,FPA&,const FPA,const FPA);
inline void minimaxFPAMIMAAB(PFPA&,PFPA&,PFPA,PFPA);
inline void minimaxFPAMIMAABCD(FPA&,FPA&,const FPA,const FPA,const FPA,const FPA);
inline void minimaxFPAMIMAABCD(PFPA&,PFPA&,PFPA,PFPA,PFPA,PFPA);
inline VLONG floorFPA(const FPA);
#endif

typedef uint16_t DBYTE;
typedef DBYTE *PDBYTE;


// const definitions

// pixel/square colors in 2bit format
const uint32_t SQUARE_GRAY=0b00;
const uint32_t SQUARE_WHITE=0b01;
const uint32_t SQUARE_BLACK=0b10;
const uint32_t SQUARE_GRAY_POTENTIALLY_WHITE=0b11;
const DBYTE COLORRED=4;

#define CF(AR,FF) \
const uint32_t AR[] = {\
	(uint32_t)(FF),\
	(uint32_t)(FF) << 2,\
	(uint32_t)(FF) << 4,\
	(uint32_t)(FF) << 6,\
	(uint32_t)(FF) << 8,\
	(uint32_t)(FF) << 10,\
	(uint32_t)(FF) << 12,\
	(uint32_t)(FF) << 14,\
	(uint32_t)(FF) << 16,\
	(uint32_t)(FF) << 18,\
	(uint32_t)(FF) << 20,\
	(uint32_t)(FF) << 22,\
	(uint32_t)(FF) << 24,\
	(uint32_t)(FF) << 26,\
	(uint32_t)(FF) << 28,\
	(uint32_t)(FF) << 30\
};

// constructing an array indexed by relative pixel number to access the bits 
// in the 32bit-integer
CF(ARRAY_SQUARE_GRAY,SQUARE_GRAY)
CF(ARRAY_SQUARE_WHITE,SQUARE_WHITE)
CF(ARRAY_SQUARE_BLACK,SQUARE_BLACK)
CF(ARRAY_SQUARE_GRAYPOTW,SQUARE_GRAY_POTENTIALLY_WHITE)

const uint32_t uint32max=0b11111111111111111111111111111111;
const int32_t MAXPTR=2048;

enum { 
	CMD_CALC=1,
	CMD_PERIOD
};

enum { 
	FUNC_Z2C=0,FUNC_Z2AZC=1,FUNC_Z3AZC=2,
	FUNC_Z4AZC=3,FUNC_Z5AZC=4,FUNC_Z6AZC=5,
	FUNCANZ
};

const char funcname[][32] = {
	"Z2C","Z2AZC","Z3AZC","Z4AZC","Z5AZC","Z6AZC"
};

// clearing a specific pixel's color via bitwise AND-mask
const uint32_t COLOR_CLEARMASK[] = {
	(uint32max - (0b11)),
	(uint32max - (0b11 << 2)),
	(uint32max - (0b11 << 4)),
	(uint32max - (0b11 << 6)),
	(uint32max - (0b11 << 8)),
	(uint32max - (0b11 << 10)),
	(uint32max - (0b11 << 12)),
	(uint32max - (0b11 << 14)),
	(uint32max - (0b11 << 16)),
	(uint32max - (0b11 << 18)),
	(uint32max - (0b11 << 20)),
	(uint32max - (0b11 << 22)),
	(uint32max - (0b11 << 24)),
	(uint32max - (0b11 << 26)),
	(uint32max - (0b11 << 28)),
	(uint32max - (0b11 << 30))
};

#define CFALL(FF) \
	( \
		((FF) << 30) | ((FF) << 28) | ((FF) << 26) | ((FF) << 24) |\
		((FF) << 22) | ((FF) << 20) | ((FF) << 18) | ((FF) << 16) |\
		((FF) << 14) | ((FF) << 12) | ((FF) << 10) | ((FF) << 8) |\
		((FF) << 6)  | ((FF) << 4)  | ((FF) << 2)  | (FF)  \
	)
	
// if all 16 consecutive pixels of a 32bit integer have the same color
const uint32_t SQUARE_GRAY_16_CONSECUTIVE=CFALL(SQUARE_GRAY);
const uint32_t SQUARE_GRAYPOTW_16_CONSECUTIVE=CFALL(SQUARE_GRAY_POTENTIALLY_WHITE);
const uint32_t SQUARE_WHITE_16_CONSECUTIVE=CFALL(SQUARE_WHITE);
const uint32_t SQUARE_BLACK_16_CONSECUTIVE=CFALL(SQUARE_BLACK);
const int32_t BASEDENOMINATOR=25;
const int64_t DENOM225=( (int64_t)1 << BASEDENOMINATOR );
#ifdef _CHUNK512
const int MAXFATOUCOMPONENTS=8192;
#else
const int MAXFATOUCOMPONENTS=65500;
#endif
const int MAXCYCLES=110;
const int MAXPERIODICPOINTS=1024;
const int FATOUCOMPONENTCOLOROFFSET=24;


// structs

// palette-entry for bitmaps
struct RGB4 {
	uint8_t B,G,R,alpha; // in the order of the bitmap file structure
};

struct Palette4 {
	RGB4 rgbs[256]; // 8-bit bitmap palette
	
	void setPaletteRGB(const int,const int,const int,const int);
};

// a rectangle in the complex plane - used for a square and its bounding box
struct PlaneRect {
	NTYP x0,x1,y0,y1;
};

// for the reverse cell graph: position of one tile
struct Parent {
	uint16_t BX,BY; 
};

typedef Parent *PParent;

struct Int2 {
	int x,y;
};

typedef Int2 *PInt2;

struct DFSPunkt {
	int x,y;
	DBYTE tiefe;
};

struct ListeFIFO {
	int next_readpos,next_writepos;
	int allokiereAmStueck;
	Int2* werte;
	
	ListeFIFO();
	virtual ~ListeFIFO();
	void start(void);
	char read(int&,int&);
	char write(const int,const int);
};

struct ListeDFS {
	int anz;
	int allokiereAmStueck;
	DFSPunkt* werte;
	
	ListeDFS();
	virtual ~ListeDFS();
	void start(void);
	char read(int&,int&,DBYTE&);
	char write(const int,const int,const DBYTE);
};

struct Cycle {
	int len;
	DBYTE immediateBasinColorIdx;
	DBYTE attractionBasinColorIdx;
	DBYTE fatouidx0,fatouidx1; // Index auf ibfcomponents = immedaite basin fatou components
};

// pixel-coordinate rectangle
struct ScreenRect {
	int32_t x0,x1,y0,y1;
};

struct FatouComponent {
	// scrc can overlap with a second or more other
	// fatou components

	ScreenRect scrc;
	DBYTE currentOrbitColorIdxTemp; // < 256 => set, >= 256: temporary
	int32_t inCycleNbr; // < 0 => none
	int8_t isimmediate; // flag
};

// one reverse cell graph tile and its parents (preimages before the iteration)
struct RevCGBlock {
	int32_t howmany; 
	// flag, whether this tile has to be checked (again) for gray pixels and the bounding
	// boxes' hits after one iteration
	int8_t tovisit;
	// how many parents are needed here? First pass result
	int32_t memused;
	// parents as variable size array
	Parent* parent;
		
	RevCGBlock();
	void addParent(const int32_t,const int32_t);
};

struct Gray_in_row {
	int32_t g0,g1;
	int32_t mem0,mem1;
};

struct ColorPalette {
	int32_t anz;
	RGB4* rgbs;

	ColorPalette();
	virtual ~ColorPalette();
	void setlen(const int32_t);
	void setInterval(const double,const double,const int32_t,const int32_t,const int32_t,const int32_t,const int32_t,const int32_t);
	void getColor(const double,int32_t&,int32_t&,int32_t&);
};

struct DByteMgr {
	DBYTE* current;
	int32_t allocatedIdx,freeFromIdx,allocatePerBlockIdx;
	PDBYTE ptr[MAXPTR];
	int32_t anzptr;
	
	DByteMgr();
	virtual ~DByteMgr();
	DBYTE* getMemory(const int32_t);
};

struct ArrayDDByteManager {
	DDBYTE* current;
	int allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PDDBYTE ptr[MAXPTR];
	int anzptr;
	
	ArrayDDByteManager ();
	virtual ~ArrayDDByteManager ();
	PDDBYTE getMemory(const int);
};

// main object
struct Data5 {
	uint32_t** zeilen;
	Gray_in_row* memgrau;
	RevCGBlock* revcgYX;
	ArrayDDByteManager* datamgr;
		
	Data5();
	virtual ~Data5();
		
	void saveBitmap4(const char*);
	void saveBitmap4_twd(const char*,const int32_t);
	
	void saveRaw(const char*);
	int32_t readRawBlowUp(void);
	int32_t readtovisit(const char*);
	void savetovisit(const char*);
};

struct ParentManager {
	Parent* lastallocated;
	int32_t memused;
	int32_t freefrom;
	int32_t anzptr;
	PParent ptr[MAXPTR];
	
	ParentManager();
	virtual ~ParentManager();
	
	Parent* getParentSpace(const int32_t);
};


// globals

VLONG ctrbbxfa=0;
int8_t interiorpresent=0;
int64_t checkclockatbbxcount0=10000000;
int64_t checkclockatbbxadd=(1 << 26);
int32_t CLOCK1HOUR=CLOCKS_PER_SEC*3600;
int8_t _PERIODICPOINTS=0;
FILE *flog=NULL;
Cycle* cycles=NULL;
FatouComponent* ibfcomponents=NULL;
int32_t anzibf=0;
int32_t anzcycles=0;
ColorPalette basinpal;
void (*getBoundingBoxfA)(PlaneRect&,PlaneRect&) = NULL;
int32_t _FUNC;
ParentManager* parentmgr=NULL;
int32_t REFINEMENTLEVEL=0;
NTYP seedC0re,seedC1re,seedC0im,seedC1im; 
NTYP FAKTORAre,FAKTORAim;
Data5 *data5;
NTYP scaleRangePerPixel,scalePixelPerRange;
int32_t scalePixelPerRangeExponent=0;
int64_t countsquares_white,countsquares_gray;
int64_t countsquares_black,countsquares_graypotw;
// region in the complex plane where all the currently still gray squares reside
NTYP planegrayx0,planegrayx1;
NTYP planegrayy0,planegrayy1;
// region in screen coordinates
int32_t encgrayx0,encgrayx1;
int32_t encgrayy0,encgrayy1;
// width of screen in pixel, must be a power of 2
int32_t SCREENWIDTH;
// number of 32bit-integers a screen row holds, equal to (SCREENWIDTH >> 4)
int32_t RANGE0=-2,RANGE1=2;
// complete holds the same values as range just in a different numerical type
NTYP COMPLETE0,COMPLETE1;
// for the low-resolution reverse cell graph working on (usually) 64x64 pixel squares or bigger
int32_t REVCGBITS,REVCGBLOCKWIDTH;
int32_t REVCGmaxnumber,REVCGmaxnumberQ;
// variables for bit precision checks


// forward declarations

// main function to compute the Julia set
void compute(void);
void freeRevCGMem(void);
// constructing the low-level reverse cell-graph
void construct_static_reverse_cellgraph(void);
// squares that directly hit the special exterior outside COMPLETE
void find_special_exterior_hitting_squares(void);
// white and potentially-white cells wiul
void propagate_definite(void);
void propagate_potw(void);
int color_changeS32(const DDBYTE,const DDBYTE,const DDBYTE,const DDBYTE);
void copy_pixel_to_2x2grid(const uint32_t,uint32_t*);

static inline int32_t scrcoord_as_lowerleft(const NTYP&);
void write2(FILE*,const uint8_t,const uint8_t);
void write4(FILE*,const uint8_t,const uint8_t,const uint8_t,const uint8_t);
static inline int minimumI(const int,const int);
static inline int maximumI(const int,const int);
static inline int maximumI(const int,const int,const int,const int);

inline NTYP minimumD(const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP);
inline NTYP minimumD(const NTYP,const NTYP,const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP,const NTYP,const NTYP);

// defines used as small expressions

#define SETDATA5BYMEM_MY(MM,YY,WW32) \
{\
	if ( \
		((MM) >= data5->memgrau[YY].mem0) &&\
		((MM) <= data5->memgrau[YY].mem1)\
	) {\
		data5->zeilen[YY][MM - data5->memgrau[YY].mem0]=WW32;\
	} else {\
		if ( (WW32) != SQUARE_WHITE_16_CONSECUTIVE ) { \
			LOGMSG4("Implementation Error. SET MM=%i YY=%i WW=%i\n",MM,YY,WW32);\
		}\
	}\
}

#define GETDATA5BYMEM_MY(MM,YY,ERG) \
{\
	if ( \
		((MM) >= data5->memgrau[YY].mem0) &&\
		((MM) <= data5->memgrau[YY].mem1)\
	) {\
		ERG=data5->zeilen[YY][MM - data5->memgrau[YY].mem0];\
	} else {\
		ERG=SQUARE_WHITE_16_CONSECUTIVE;\
	}\
}

#define LOGMSG(TT) \
{\
	fprintf(flog,TT); fflush(flog);\
	printf(TT);\
}

#define LOGMSG2(TT,AA) \
{\
	fprintf(flog,TT,AA); fflush(flog);\
	printf(TT,AA);\
}

#define LOGMSG3(TT,AA,BB) \
{\
	fprintf(flog,TT,AA,BB); fflush(flog);\
	printf(TT,AA,BB);\
}

#define LOGMSG4(TT,AA,BB,CC) \
{\
	fprintf(flog,TT,AA,BB,CC); fflush(flog);\
	printf(TT,AA,BB,CC);\
}

#define LOGMSG5(TT,AA,BB,CC,DD) \
{\
	fprintf(flog,TT,AA,BB,CC,DD); fflush(flog);\
	printf(TT,AA,BB,CC,DD);\
}

#define LOGMSG6(TT,AA,BB,CC,DD,EE) \
{\
	fprintf(flog,TT,AA,BB,CC,DD,EE); fflush(flog);\
	printf(TT,AA,BB,CC,DD,EE);\
}

#define SQUARE_LIES_ENTIRELY_IN_SPECEXT(BBX) \
	(\
		(BBX.x1 < COMPLETE0) ||\
		(BBX.x0 > COMPLETE1) ||\
		(BBX.y1 < COMPLETE0) ||\
		(BBX.y0 > COMPLETE1)\
	)
	
#define SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(BBX) \
	(\
		(BBX.x1 < planegrayx0) ||\
		(BBX.x0 > planegrayx1) ||\
		(BBX.y1 < planegrayy0) ||\
		(BBX.y0 > planegrayy1)\
	)

#define SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(BBX) \
	(\
		(BBX.x1 < planegrayx1) &&\
		(BBX.x0 > planegrayx0) &&\
		(BBX.y1 < planegrayy1) &&\
		(BBX.y0 > planegrayy0)\
	)

#define GET_SINGLE_CELLCOLOR_XY(XX,YY,FF) \
{\
	int mem=(XX) >> 4;\
	int bpos=(XX) & 0b1111;\
	DDBYTE w;\
	GETDATA5BYMEM_MY(mem,YY,w);\
	FF=(w >> (bpos << 1)) &0b11;\
}

#define SET_SINGLE_CELLCOLOR_XYFF(XX,YY,FF) \
{\
	int mem=(XX) >> 4;\
	int bpos=(XX) & 0b1111;\
	DDBYTE w;\
	GETDATA5BYMEM_MY(mem,YY,w);\
	w &= COLOR_CLEARMASK[bpos];\
	w |= ( (FF) << (bpos << 1) );\
	SETDATA5BYMEM_MY(mem,YY,w);\
}

#define GET_SINGLE_PIXELCOLOR_FROM_4BYTEINTEGER(WW,BPOS) \
	( ( (WW) >> (BPOS) ) & (uint32_t)3 )
	
#define SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(WW,CLEARMASKE,SETF) \
	((WW & CLEARMASKE) | SETF)


// routines

char* upper(char* s) {
	if (!s) return NULL;
	
	for(int32_t i=0;i<strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

inline int32_t maximumI(const int32_t a,const int32_t b) {
	if (a > b) return a;
	return b;
}

inline int32_t minimumI(const int32_t a,const int32_t b) {
	if (a < b) return a;
	return b;
}

inline int32_t maximumI(const int32_t a,const int32_t b,const int32_t c,const int32_t d) {
	int32_t m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}

// RevCGBlock
// static reverse cell graph at the given SCREENWIDTH and REVCGBLOCKWIDTH

RevCGBlock::RevCGBlock() {
	howmany=0;
	memused=0;
}

void RevCGBlock::addParent(const int ax,const int ay) {
	if (
		(ax < 0) ||
		(ay < 0) ||
		(ax >= REVCGmaxnumber) ||
		(ay >= REVCGmaxnumber) 
	) {
		LOGMSG3("Implementation error. Parent at %i,%i not valid.\n",ax,ay);
		return; 
	}
	
	if (parent == NULL) {
		parent=parentmgr->getParentSpace(memused);
		howmany=0;
		if (!parent) {
			LOGMSG("Memory failure for parent.\n");
			exit(99);
		}
	}
	
	parent[howmany].BX=ax;
	parent[howmany].BY=ay;
	howmany++;
}

// struct data5
// main struct: handles the pixels/squares and the storing and reading

// saves the currently computed square data (color) as a raw file for future continuation
// of computing or blowing-up for the next refinement level
void Data5::savetovisit(const char* afn) {
	FILE *f=fopen(afn,"wb");
	fwrite(&REVCGmaxnumberQ,1,sizeof(REVCGmaxnumberQ),f);
	for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
		fwrite(&revcgYX[i].tovisit,1,sizeof(revcgYX[i].tovisit),f);
	}
	fclose(f);
}

int32_t Data5::readtovisit(const char* afn) {
	FILE *f=fopen(afn,"rb");
	if (!f) return 0;
	
	int32_t a;
	fwrite(&a,1,sizeof(a),f);
	if (a != REVCGmaxnumberQ) {
		// maybe wrong source ?
		fclose(f);
		return 0;
	}
	
	for(uint32_t i=0;i<REVCGmaxnumberQ;i++) {
		fread(&revcgYX[i].tovisit,1,sizeof(revcgYX[i].tovisit),f);
	}
	fclose(f);
	
	return 1;
}

void Data5::saveRaw(const char* afn) {
	// raw data - new structure
	char fn[1024];
	sprintf(fn,"%s.raw",afn);
	FILE *f=fopen(fn,"wb");
	DDBYTE w=SCREENWIDTH;
	fwrite(&w,sizeof(w),1,f);
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if (memgrau[y].g0 <= memgrau[y].g1) {
			// memg0 <= g0 <= g1 <= memg1
			int32_t m0=SCREENWIDTH >> 4,m1=0;
			for(int32_t mem=memgrau[y].mem0;mem<=memgrau[y].mem1;mem++) {
				DDBYTE w;
				GETDATA5BYMEM_MY(mem,y,w);
				if (w != SQUARE_WHITE_16_CONSECUTIVE) {
					if (mem < m0) m0=mem;
					if (mem > m1) m1=mem;
				}
			}

			int32_t laenge=m1-m0+1;
			if (laenge>0) {
				fwrite(&m0,1,sizeof(m0),f);
				fwrite(&laenge,1,sizeof(laenge),f);
				fwrite(
					&zeilen[y][m0 - memgrau[y].mem0],
					laenge,
					sizeof(DDBYTE),
				f);
			} else {
				// empty row
				int32_t start=0,laenge=0;
				fwrite(&start,1,sizeof(start),f);
				fwrite(&laenge,1,sizeof(laenge),f);
			}
		} else {
			// empty row
			int32_t start=0,laenge=0;
			fwrite(&start,1,sizeof(start),f);
			fwrite(&laenge,1,sizeof(laenge),f);
		}
	} // y
	
	fclose(f);
}

// does a row contain interior cells
int32_t interiorinrow(DDBYTE* aw,const int32_t alen) {
	for(int32_t i=0;i<alen;i++) {
		DDBYTE w=aw[i];
		for(int bit=0;bit<16;bit++) {
			if ((w & 0b11) == SQUARE_BLACK) return 1;
			w >>= 2;
		}
	}
	return 0;
}

int32_t Data5::readRawBlowUp(void) {
	FILE *f=fopen("_in.raw","rb");
	if (!f) {
		printf("No stored data found. Computation starts anew.\n");
		return 0;
	}
	
	DDBYTE savedlen;
	fread(&savedlen,sizeof(savedlen),1,f);
	// zeilen-Ptr are not initialised here
	encgrayx0=encgrayy0=SCREENWIDTH-1;
	encgrayx1=encgrayy1=0;
	
	DDBYTE *eine=new DDBYTE[SCREENWIDTH];
	
	printf("reading stored data ");
	
	VLONG memused=0;

	if ((int32_t)savedlen == SCREENWIDTH) {
		// 1 auf 1
		
		int32_t start,laenge;
		for(int32_t y=0;y<SCREENWIDTH;y++) {
			fread(&start,1,sizeof(start),f);
			fread(&laenge,1,sizeof(laenge),f);
			if (laenge<=0) {
				zeilen[y]=NULL;
				memgrau[y].g0=SCREENWIDTH;
				memgrau[y].g1=0;
				memgrau[y].mem0=(SCREENWIDTH >> 4);
				memgrau[y].mem1=0;
			} else {
				zeilen[y]=data5->datamgr->getMemory(laenge);
				memused += (laenge*sizeof(DDBYTE));
				if (!zeilen[y]) {
					LOGMSG("Speicherfehler. ReadRaw\n");
					exit(99);
				}
				fread(zeilen[y],laenge,sizeof(DDBYTE),f);
				// is there at least 1 black cell ?
				if (interiorpresent<=0) {
					interiorpresent=interiorinrow(zeilen[y],laenge);
				}
				memgrau[y].mem0=start;
				memgrau[y].mem1=start+laenge-1;
				memgrau[y].g0=memgrau[y].mem0 << 4;
				memgrau[y].g1=((memgrau[y].mem1+1) << 4)-1;
				if (y < encgrayy0) encgrayy0=y;
				if (y > encgrayy1) encgrayy1=y;
				if (memgrau[y].g0 < encgrayx0) encgrayx0=memgrau[y].g0;
				if (memgrau[y].g1 > encgrayx1) encgrayx1=memgrau[y].g1;
			}
		}
	} else if ( (int32_t)savedlen == (SCREENWIDTH >> 1) ) {
		int32_t start,laenge;
		for(int32_t yread=0;yread<(SCREENWIDTH-1);yread+=2) {
			fread(&start,1,sizeof(start),f);
			fread(&laenge,1,sizeof(laenge),f);
			if (laenge<=0) {
				zeilen[yread]=zeilen[yread+1]=NULL;
				memgrau[yread].g0=memgrau[yread+1].g0=SCREENWIDTH;
				memgrau[yread].g1=memgrau[yread+1].g1=0;
				memgrau[yread].mem0=memgrau[yread+1].mem0=SCREENWIDTH >> 4;
				memgrau[yread].mem1=memgrau[yread+1].mem1=0;
			} else {
				// one cell into a 2x2 grid - refinement process
				DDBYTE readlaenge=laenge;
				laenge <<= 1;
				start <<= 1;
				zeilen[yread]=data5->datamgr->getMemory(laenge);
				zeilen[yread+1]=data5->datamgr->getMemory(laenge);
				memused += (2*laenge*sizeof(DDBYTE));
				if (!zeilen[yread+1]) {
					LOGMSG("Memory error. ReadRaw\n");
					exit(99);
				}
				fread(eine,readlaenge,sizeof(DDBYTE),f);
				if (interiorpresent<=0) {
					interiorpresent=interiorinrow(eine,readlaenge);
				}
				
				DDBYTE mem=0;
				for(DDBYTE k=0;k<readlaenge;k++) {
					// substitute POTWGRAU with GRAY
					// as potentially-white-information cannot be blowed-up
					DDBYTE ziel[2];
					copy_pixel_to_2x2grid(eine[k],ziel);
					data5->zeilen[yread][mem]=ziel[0];
					data5->zeilen[yread][mem+1]=ziel[1];
					data5->zeilen[yread+1][mem]=ziel[0];
					data5->zeilen[yread+1][mem+1]=ziel[1];
					mem += 2;
				} 
				
				memgrau[yread].mem0=memgrau[yread+1].mem0=start;
				memgrau[yread].mem1=memgrau[yread+1].mem1=start+laenge-1;
				memgrau[yread].g0=memgrau[yread+1].g0=(memgrau[yread].mem0 << 4);
				memgrau[yread].g1=memgrau[yread+1].g1=((memgrau[yread].mem1+1) << 4)-1;

				if (yread < encgrayy0) encgrayy0=yread;
				if ( (yread+1) > encgrayy1) encgrayy1=(yread+1);
				if (memgrau[yread].g0 < encgrayx0) encgrayx0=memgrau[yread].g0;
				if (memgrau[yread].g1 > encgrayx1) encgrayx1=memgrau[yread].g1;
				
			} // 2x2 gird
		} // yread 
	} else {
		LOGMSG("ReadBlowup. Wrong resolution. File ignored.\n");
		fclose(f);
		
		return 0;
	}
	
	fclose(f);
	delete[] eine;
	
	printf("\n  %I64d GB pixel memory allocated\n",1+(memused >> 30));
	
	// encgray 
	planegrayx0=encgrayx0*scaleRangePerPixel + COMPLETE0;
	planegrayy0=encgrayy0*scaleRangePerPixel + COMPLETE0;
	planegrayx1=(encgrayx1+16)*scaleRangePerPixel + COMPLETE0;
	planegrayy1=(encgrayy1+16)*scaleRangePerPixel + COMPLETE0;
	
	return 1;
}

void Data5::saveBitmap4_twd(const char* afn,const int atwdexp) {
	// saves a trustworthily downsized version of the image: 16-fold. 
	// image format is: 8 bit Bitmap
	// 2^TWDEXP x 2^TWDEXP gives one final pixel
	
	int32_t _TWDEXPONENT=atwdexp;
	if (atwdexp<0) {
		// adjust exponentn so final image is at most 2^16 x 2^16
		while ( (SCREENWIDTH >> _TWDEXPONENT) > 65536) _TWDEXPONENT++;
	}

	int32_t bytes_per_row = SCREENWIDTH >> _TWDEXPONENT;
	uint8_t* rgbz=new uint8_t[bytes_per_row];
	uint32_t off
	=		14 // size of file header
		+	40 // size of bitmap header
		+	256*4 // palette entries
	;
	uint32_t filelen
		=	off
		+	(bytes_per_row*bytes_per_row);
	;
	
	char tmp[1024];

	// setting up the 256-entry palette
	RGB4 pal[256];
	for(int32_t i=0;i<256;i++) pal[i].R=pal[i].G=pal[i].B=pal[i].alpha=0;
	pal[SQUARE_GRAY].R=127;
	pal[SQUARE_GRAY].G=127;
	pal[SQUARE_GRAY].B=127;
	pal[SQUARE_BLACK].R=0;
	pal[SQUARE_BLACK].G=0;
	pal[SQUARE_BLACK].B=0;
	pal[SQUARE_WHITE].R=255;
	pal[SQUARE_WHITE].G=255;
	pal[SQUARE_WHITE].B=255;

	sprintf(tmp,"%s_2_%i-fold.bmp",afn,_TWDEXPONENT);
	FILE *fbmp=fopen(tmp,"wb");
	if (!fbmp) {
		fprintf(stderr,"Error. twd16\n");
		exit(99);
	}
	write2(fbmp,66,77); 
	fwrite(&filelen,1,sizeof(filelen),fbmp);
	write4(fbmp,0,0,0,0); 
	fwrite(&off,1,sizeof(off),fbmp); 
	write4(fbmp,40,0,0,0); 
	fwrite(&bytes_per_row,sizeof(bytes_per_row),1,fbmp);
	fwrite(&bytes_per_row,sizeof(bytes_per_row),1,fbmp);
	write2(fbmp,1,0);
	write2(fbmp,8,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,19,11,0,0);
	write4(fbmp,19,11,0,0);
	write4(fbmp,0,1,0,0);
	write4(fbmp,0,0,0,0);

	for(int32_t i=0;i<256;i++) write4(fbmp,pal[i].B,pal[i].G,pal[i].R,pal[i].alpha);
	
	int32_t TWDSTEP=(1 << _TWDEXPONENT);
			
	// 16x16 pixels into one: 16 rows and each with one 32bit-integer
	for(int32_t y=0;y<SCREENWIDTH;y+=TWDSTEP) {
		int32_t xzeile=-1;
		for(int32_t x=0;x<SCREENWIDTH;x+=TWDSTEP) {
			xzeile++;
			
			int32_t finalf=-1;
			for(int dy=0;dy<TWDSTEP;dy++) {
				for(int dx=0;dx<TWDSTEP;dx++) {
					int f;
					GET_SINGLE_CELLCOLOR_XY(x+dx,y+dy,f);
					if (finalf<0) finalf=f;
					else if (finalf != f) {
						finalf=SQUARE_GRAY;
						break;
					} else if (f == SQUARE_GRAY_POTENTIALLY_WHITE) {
						finalf=SQUARE_GRAY;
						break;
					}
					
					if (finalf == SQUARE_GRAY) break;
				}

				if (finalf == SQUARE_GRAY) break;
			} // dy

			rgbz[xzeile]=finalf;
		} // x
		
		fwrite(rgbz,bytes_per_row,sizeof(uint8_t),fbmp);
	} // y
	
	fclose(fbmp);
				
	delete[] rgbz;
}

void Data5::saveBitmap4(const char* afn) {
	// save the current pixels in several 4-bit-Bitmaps (each maximal 2 GB in size)

	int32_t width_of_one_image=SCREENWIDTH;
	const int64_t max_filesize_per_bitmap=(int64_t)(1) << 31;
	while (width_of_one_image > 16) {
		int64_t memory_used=width_of_one_image;
		// since storing 2 pixels per bitmap entry
		memory_used *= (width_of_one_image >> 1); // 2 Pixel pro uint8_t
		if (memory_used > max_filesize_per_bitmap) {
			// image width too high => halven in
			width_of_one_image >>= 1;
		} else break;
	} // while
	
	countsquares_gray=countsquares_black=0;
	countsquares_white=countsquares_graypotw=0;
	
	// 2 resulting pixels in bitmap per byte
	int32_t bytes_per_row=(width_of_one_image >> 1);
	if ((bytes_per_row % 4) != 0) {
		bytes_per_row=(((bytes_per_row >> 2) + 1) << 2);
	}
	uint8_t* rgbz=new uint8_t[bytes_per_row];
	uint32_t off
	=		14 // file header
		+	40 // bitmap header
		+	16*4 // palette entries, 16 here
	;
	uint32_t filelen
		=	off
		+	(bytes_per_row*width_of_one_image);
	;
	
	char tmp[1024];
	RGB4 pal[16];
	for(int32_t i=0;i<16;i++) pal[i].R=pal[i].G=pal[i].B=pal[i].alpha=0;
	pal[SQUARE_GRAY].R=127;
	pal[SQUARE_GRAY].G=127;
	pal[SQUARE_GRAY].B=127;
	pal[SQUARE_BLACK].R=0;
	pal[SQUARE_BLACK].G=0;
	pal[SQUARE_BLACK].B=0;
	pal[SQUARE_WHITE].R=255;
	pal[SQUARE_WHITE].G=255;
	pal[SQUARE_WHITE].B=255;
	
	int32_t cy=0;
	for(int32_t fy=0;fy<SCREENWIDTH;fy+=width_of_one_image) {
		int32_t cx=0;
		for(int32_t fx=0;fx<SCREENWIDTH;fx+=width_of_one_image) {
			sprintf(tmp,"%s_Y%02iX%02i.bmp",afn,cy,cx);
			FILE *fbmp=fopen(tmp,"wb");
			if (!fbmp) {
				fprintf(stderr,"Error. saveBitmap4\n");
				exit(99);
			}

			write2(fbmp,66,77);
			fwrite(&filelen,1,sizeof(filelen),fbmp);
			write4(fbmp,0,0,0,0);
			fwrite(&off,1,sizeof(off),fbmp);
			write4(fbmp,40,0,0,0);
	
			uint32_t w = width_of_one_image;
			fwrite(&w,sizeof(w),1,fbmp);
			fwrite(&w,sizeof(w),1,fbmp);
			write2(fbmp,1,0);
			write2(fbmp,4,0);
			write4(fbmp,0,0,0,0);
			write4(fbmp,0,0,0,0);
			write4(fbmp,19,11,0,0);
			write4(fbmp,19,11,0,0);
			write4(fbmp,16,0,0,0);
			write4(fbmp,0,0,0,0);

			for(int32_t i=0;i<16;i++) write4(fbmp,pal[i].B,pal[i].G,pal[i].R,pal[i].alpha);
			
			const int32_t yfende=fy+width_of_one_image;
			for(int32_t y=fy;y<yfende;++y) {
				int32_t mem=fx >> 4;
				int32_t xpos=0;

				for(int32_t x=fx;x<(fx+width_of_one_image);x+=16) {
					//uint32_t w=data5->pixelrow[y][mem];
					uint32_t w;
					GETDATA5BYMEM_MY(mem,y,w);
					
					int32_t bit=4;
					uint8_t erg=0;
					for(int32_t bp=0;bp<32;bp+=2) {
						switch (GET_SINGLE_PIXELCOLOR_FROM_4BYTEINTEGER(w,bp)) {
							case SQUARE_WHITE:
								erg |= (SQUARE_WHITE << bit);
								countsquares_white++;
								break;
							case SQUARE_BLACK:
								erg |= (SQUARE_BLACK << bit);
								countsquares_black++;
								break;
							case SQUARE_GRAY_POTENTIALLY_WHITE:
								erg |= (SQUARE_GRAY << bit);
								countsquares_graypotw++;
								break;
							default:
								erg |= (SQUARE_GRAY << bit);
								countsquares_gray++;
								break;
						} // switch
						bit -= 4;
						if (bit < 0) {
							bit=4;
							rgbz[xpos]=erg;
							erg=0;
							xpos++;
						}
					} // bp
					
					mem++;
				} // x
			
				fwrite(rgbz,bytes_per_row,1,fbmp);
			} // y
		
			fclose(fbmp);
			cx++;
		} // fx
		cy++;
	} // fy
	
	delete[] rgbz;

	LOGMSG2("%I64d interior",countsquares_black);
	LOGMSG2(" %I64d exterior",countsquares_white);
	LOGMSG3(" %I64d (%.0lf%%) gray squares\n",
		countsquares_gray+countsquares_graypotw,
		100.0*(double)(countsquares_gray+countsquares_graypotw)/(countsquares_black+countsquares_graypotw+countsquares_gray+countsquares_white));

	if (countsquares_graypotw>0) {
		LOGMSG2("(of which %I64d are gray but potentially white)\n",countsquares_graypotw);
	}
}

Data5::Data5() {
	printf("initialising main object ...\n");
	
	memgrau=new Gray_in_row[SCREENWIDTH];
	zeilen=new uint32_t*[SCREENWIDTH];
	revcgYX=new RevCGBlock[REVCGmaxnumber*REVCGmaxnumber];
	datamgr=new ArrayDDByteManager;
}

Data5::~Data5() {
	delete datamgr;
	if (revcgYX) delete[] revcgYX;
	delete[] memgrau;
	delete[] zeilen;
}

// one pixel transforming into a 2x2 grid, gray or gray-potentially-white will both
// be set to gray. works on a 32bit-integer, hence 16 consecutive bits directly
void copy_pixel_to_2x2grid(const uint32_t q,uint32_t* erg) {
	int32_t eidx=0;
	int32_t zbit=0;
	erg[0]=erg[1]=0;
	
	for(int32_t qbit=0;qbit<32;qbit+=2) {
		uint32_t w=((q >> qbit) & 0b11);
		if (w == SQUARE_GRAY_POTENTIALLY_WHITE) w=SQUARE_GRAY;
		erg[eidx] |= (w << zbit);
		zbit += 2;
		erg[eidx] |= (w << zbit);
		zbit += 2;
		if (zbit >= 32) {
			zbit=0;
			eidx=1;
		}
	}
}

static inline int32_t scrcoord_as_lowerleft(const NTYP& a) {
	// calculating the screen coordinte of the pixel that contains the coordinate
	// if the coordinate lies on an edge/corner (and belongs to more than one pixel)
	// the pixel where it lies on the left,bottom edge/corner is returned
	// check for out-of-COMPLETE before flooring values
	// to remain in the representable range

	int32_t w;
	
	#ifdef _QUADMATH
	if (a <= COMPLETE0) return 0; // left
	if (a >= COMPLETE1) return (SCREENWIDTH-1);
	w=(int)floorq( (a - COMPLETE0) * scalePixelPerRange );
	if (w >= SCREENWIDTH) return (SCREENWIDTH-1);
	return w;
	#endif

	#ifdef _F107
	if (a <= COMPLETE0) return 0; // left
	if (a >= COMPLETE1) return (SCREENWIDTH-1);
	w= (int)floor( (a - COMPLETE0) * scalePixelPerRange );
	if (w >= SCREENWIDTH) return (SCREENWIDTH-1);
	return w;
	#endif
	
	#ifdef _FPA
	if (a.vorz<0) {
		if (a.a >= (-RANGE0)) return 0;
	} else if (a.vorz==0) return (SCREENWIDTH >> 1);
	else if (a.a >= RANGE1) return (SCREENWIDTH-1);

	FPA b=(a-COMPLETE0);
	// scalePixelPerRange is a power of 2
	// so bit-shifting can be used as multiplication
	//b=b*scalePixelPerRange;
	b.shiftLeft(scalePixelPerRangeExponent);
	// Floor
	VLONG wfl=floorFPA(b);

	if (
		(wfl < -INT32MAX) ||
		(wfl > INT32MAX)
	) {
		printf("Implementation error. FPA. out-of-range floor scrcoord %I64d\n",wfl);
		exit(99);
	}
	if (wfl >= SCREENWIDTH) return (SCREENWIDTH-1);
	
	return (int32_t)wfl;
	#endif
	
	#ifndef _FPA
	if (a <= COMPLETE0) return 0; // left
	if (a >= COMPLETE1) return (SCREENWIDTH-1);
	w=(int)floor( (a - COMPLETE0) * scalePixelPerRange );
	if (w >= SCREENWIDTH) return (SCREENWIDTH-1);
	return w;
	#endif
}

#ifndef _FPA
inline NTYP minimumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	return m;
}

inline NTYP minimumD(const NTYP a,const NTYP b) {
	if (a < b) return a;
	return b;
}

inline NTYP maximumD(const NTYP a,const NTYP b) {
	if (a > b) return a;
	return b;
}

// maximum of 4 doubles. often used by the boundingbox function
inline NTYP maximumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}
#endif

// IMPORTANT FUNCTION
// - function that performs the iteration of the quadratic case: z := z*z + seedC
// - function can be replaced by any other boundiong box generating routine to accommodate
// for other types of iteration formulas like z^3+c, z^6+c etc.
// - CAVE: Care must be taken that the C++ double floating point type used here can handle
// the produced numbers to absolute accuracy
// - expression is semi-automatically generated and left un-optimized for ease of substitution

// iterting function z^2+c
void getBoundingBoxfA_z2c(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;

	fA.x0=minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+seedC0re;
	fA.x1=maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+seedC1re;
	fA.y0=2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedC0im;
	fA.y1=2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedC1im;
}

// z^2+A*z+c
void getBoundingBoxfA_z2azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-maximumD(A.y0*A.y0,A.y1*A.y1);
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-minimumD(A.y0*A.y0,A.y1*A.y1);
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
}

// z^3+A*z+c
void getBoundingBoxfA_z3azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;

	#ifdef _F107
	// optimized version - path of calculation
	NTYP x02=A.x0*A.x0;
	NTYP y02=A.y0*A.y0;
	NTYP x12=A.x1*A.x1;
	NTYP y12=A.y1*A.y1;
	NTYP mi1,ma1;
	minimaxF107AB(mi1,ma1,x02,x12);
	NTYP mi2,ma2;
	minimaxF107AB(mi2,ma2,y02,y12);
	NTYP mi3,ma3;
	minimaxF107AB(mi3,ma3,FAKTORAre*A.x0,FAKTORAre*A.x1);
	NTYP mi4,ma4;
	minimaxF107AB(mi4,ma4,FAKTORAim*A.y0,FAKTORAim*A.y1);
	NTYP mi5,ma5;
	minimaxF107AB(mi5,ma5,FAKTORAre*A.y0,FAKTORAre*A.y1);
	NTYP mi6,ma6;
	minimaxF107AB(mi6,ma6,FAKTORAim*A.x0,FAKTORAim*A.x1);
	NTYP mi7,ma7;
	minimaxF107ABCD(mi7,ma7,A.x0*mi2,A.x0*ma2,A.x1*mi2,A.x1*ma2);
	NTYP mi8,ma8;
	minimaxF107ABCD(mi8,ma8,mi1*A.y0,mi1*A.y1,ma1*A.y0,ma1*A.y1);

	fA.x0=seedC0re+mi3-ma4+x02*A.x0-(3*ma7);
	fA.x1=seedC1re+ma3-mi4+x12*A.x1-(3*mi7);

	fA.y0=mi5+mi6+3*mi8-(y12*A.y1)+seedC0im;
	fA.y1=ma5+ma6+3*ma8-(y02*A.y0)+seedC1im;

	return;
	#endif
	
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x0*A.x0*A.x0-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x1*A.x1*A.x1-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y1*A.y1*A.y1)+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y0*A.y0*A.y0)+seedC1im;
}

// z^4+A*z+c
void getBoundingBoxfA_z4azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	#ifdef _F107
	NTYP x02=A.x0*A.x0;
	NTYP y02=A.y0*A.y0;
	NTYP x12=A.x1*A.x1;
	NTYP y12=A.y1*A.y1;
	// square(A.x0) does not work here

	NTYP x13=x12*A.x1;
	NTYP x03=x02*A.x0;
	NTYP y03=y02*A.y0;
	NTYP y13=y12*A.y1;
	NTYP mi1,ma1;
	minimaxF107ABCD(mi1,ma1,A.x0*(y03),A.x0*(y13),A.x1*(y03),A.x1*(y13));
	NTYP mi2,ma2;
	minimaxF107ABCD(mi2,ma2,(x03)*A.y0,(x03)*A.y1,(x13)*A.y0,(x13)*A.y1);
	NTYP mi3,ma3;
	minimaxF107AB(mi3,ma3,x02,x12);
	NTYP mi4,ma4;
	minimaxF107AB(mi4,ma4,y02,y12);
	NTYP mi5,ma5;
	minimaxF107AB(mi5,ma5,y02*y02,y12*y12);
	NTYP mi6,ma6;
	minimaxF107AB(mi6,ma6,x02*x02,x12*x12);
	NTYP mi7,ma7;
	minimaxF107AB(mi7,ma7,FAKTORAre*A.x0,FAKTORAre*A.x1);
	NTYP mi8,ma8;
	minimaxF107AB(mi8,ma8,FAKTORAim*A.y0,FAKTORAim*A.y1);
	NTYP mi9,ma9;
	minimaxF107AB(mi9,ma9,FAKTORAre*A.y0,FAKTORAre*A.y1);
	NTYP mi10,ma10;
	minimaxF107AB(mi10,ma10,FAKTORAim*A.x0,FAKTORAim*A.x1);
	NTYP mi11,ma11;
	minimaxF107ABCD(mi11,ma11,mi3*mi4,mi3*ma4,ma3*mi4,ma3*ma4);

	fA.x0=seedC0re+mi7-ma8+mi6-(6*ma11)+mi5;
	fA.x1=seedC1re+ma7-mi8+ma6-(6*mi11)+ma5;

	fA.y0=seedC0im+mi9+mi10+4*mi2-(4*ma1);
	fA.y1=seedC1im+ma9+ma10+4*ma2-(4*mi1);
	
	return;
	#endif

	#ifdef _FPA
	FPA x02,x12,y02,y12;
	FPA x04,x14,y04,y14;
	A.x0.squareTo(x02);
	A.x1.squareTo(x12);
	A.y0.squareTo(y02);
	A.y1.squareTo(y12);
	x02.squareTo(x04);
	x12.squareTo(x14);
	y02.squareTo(y04);
	y12.squareTo(y14);
	FPA x03=x02*A.x0,x13=x12*A.x1,y03=y02*A.y0,y13=y12*A.y1;
	FPA arx0=FAKTORAre*A.x0,arx1=FAKTORAre*A.x1;
	FPA aiy0=FAKTORAim*A.y0,aiy1=FAKTORAim*A.y1;
	FPA ary0=FAKTORAre*A.y0,ary1=FAKTORAre*A.y1;
	FPA aix0=FAKTORAim*A.x0,aix1=FAKTORAim*A.x1;
	FPA x03y0=x03*A.y0,x03y1=x03*A.y1,x13y0=x13*A.y0,x13y1=x13*A.y1;
	FPA y03x0=A.x0*y03,y13x0=A.x0*y13,y03x1=A.x1*y03,y13x1=A.x1*y13;

	PFPA mix1,max1,miy1,may1,mix4,max4,miy4,may4;
	PFPA miarx,maarx,miaiy,maaiy,miary,maary,miaix,maaix;
	PFPA mix2,max2;
	PFPA mi3,ma3,mi4,ma4;

	minimaxFPAMIMAAB(mix1,max1,&x02,&x12);
	minimaxFPAMIMAAB(miy1,may1,&y02,&y12);
	minimaxFPAMIMAAB(mix4,max4,&x04,&x14);
	minimaxFPAMIMAAB(miy4,may4,&y04,&y14);
	minimaxFPAMIMAAB(miarx,maarx,&arx0,&arx1);
	minimaxFPAMIMAAB(miaiy,maaiy,&aiy0,&aiy1);
	minimaxFPAMIMAAB(miary,maary,&ary0,&ary1);
	minimaxFPAMIMAAB(miaix,maaix,&aix0,&aix1);
	
	FPA mix1miy1=(*mix1)*(*miy1),mix1may1=(*mix1)*(*may1);
	FPA max1miy1=(*max1)*(*miy1),max1may1=(*max1)*(*may1);
	minimaxFPAMIMAABCD(mix2,max2,&mix1miy1,&mix1may1,&max1miy1,&max1may1);
	
	minimaxFPAMIMAABCD(mi3,ma3,&x03y0,&x03y1,&x13y0,&x13y1);
	minimaxFPAMIMAABCD(mi4,ma4,&y03x0,&y13x0,&y03x1,&y13x1);

	fA.x0=(*miarx)-(*maaiy)+(*mix4)-(6*(*max2))+(*miy4)+seedC0re;
	fA.x1=(*maarx)-(*miaiy)+(*max4)-(6*(*mix2))+(*may4)+seedC1re;
	fA.y0=(*miary)+(*miaix)+4*(*mi3)-(4*(*ma4))+seedC0im;
	fA.y1=(*maary)+(*maaix)+4*(*ma3)-(4*(*mi4))+seedC1im;
	
	return;
	#endif

	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+4*minimumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*maximumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+4*maximumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*minimumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedC1im;
}

// z^5+A*z+c
void getBoundingBoxfA_z5azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	#ifdef _F107
	NTYP x02=A.x0*A.x0;
	NTYP x03=x02*A.x0;
	NTYP x04=x02*x02;
	NTYP y02=A.y0*A.y0;
	NTYP y03=y02*A.y0;
	NTYP y04=y02*y02;
	NTYP x12=A.x1*A.x1;
	NTYP x13=x12*A.x1;
	NTYP x14=x12*x12;
	NTYP y12=A.y1*A.y1;
	NTYP y13=y12*A.y1;
	NTYP y14=y12*y12;
	NTYP mi1,ma1;
	minimaxF107AB(mi1,ma1,x02,x12);
	NTYP mi2,ma2;
	minimaxF107AB(mi2,ma2,y02,y12);
	NTYP mi3,ma3;
	minimaxF107AB(mi3,ma3,x04,x14);
	NTYP mi4,ma4;
	minimaxF107AB(mi4,ma4,y04,y14);
	NTYP mi5,ma5;
	minimaxF107AB(mi5,ma5,FAKTORAre*A.x0,FAKTORAre*A.x1);
	NTYP mi6,ma6;
	minimaxF107AB(mi6,ma6,FAKTORAim*A.y0,FAKTORAim*A.y1);
	NTYP mi7,ma7;
	minimaxF107AB(mi7,ma7,FAKTORAre*A.y0,FAKTORAre*A.y1);
	NTYP mi8,ma8;
	minimaxF107AB(mi8,ma8,FAKTORAim*A.x0,FAKTORAim*A.x1);
	NTYP mi9,ma9;
	minimaxF107ABCD(mi9,ma9,(x03)*mi2,(x03)*ma2,(x13)*mi2,(x13)*ma2);
	NTYP mi10,ma10;
	minimaxF107ABCD(mi10,ma10,A.x0*mi4,A.x0*ma4,A.x1*mi4,A.x1*ma4);
	NTYP mi11,ma11;
	minimaxF107ABCD(mi11,ma11,mi3*A.y0,mi3*A.y1,ma3*A.y0,ma3*A.y1);
	NTYP mi12,ma12;
	minimaxF107ABCD(mi12,ma12,mi1*(y03),mi1*(y13),ma1*(y03),ma1*(y13));

	fA.x0=mi5-ma6+x02*x03-(10*ma9)+5*mi10+seedC0re;
	fA.x1=ma5-mi6+x12*x13-(10*mi9)+5*ma10+seedC1re;

	fA.y0=mi7+mi8+5*mi11-(10*ma12)+y02*y03+seedC0im;
	fA.y1=ma7+ma8+5*ma11-(10*mi12)+y12*y13+seedC1im;
	return;
	#endif

	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x0*A.x0*A.x0*A.x0*A.x0-(2*(5*maximumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*minimumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x1*A.x1*A.x1*A.x1*A.x1-(2*(5*minimumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*maximumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y0*A.y0*A.y0*A.y0*A.y0+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y1*A.y1*A.y1*A.y1*A.y1+seedC1im;
}

// z^6+A*z+c
void getBoundingBoxfA_z6azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	#ifdef _F107
	NTYP x02=A.x0*A.x0;
	NTYP x03=x02*A.x0;
	NTYP x04=x02*x02;
	NTYP x05=x02*x03;
	NTYP y02=A.y0*A.y0;
	NTYP y03=y02*A.y0;
	NTYP y04=y02*y02;
	NTYP y05=y02*y03;
	NTYP x12=A.x1*A.x1;
	NTYP x13=x12*A.x1;
	NTYP x14=x12*x12;
	NTYP x15=x12*x13;
	NTYP y12=A.y1*A.y1;
	NTYP y13=y12*A.y1;
	NTYP y14=y12*y12;
	NTYP y15=y12*y13;
	
	Pf107 mi1,ma1;
	NTYP tmp1=(x05)*A.y0;
	NTYP tmp2=(x05)*A.y1;
	NTYP tmp3=(x15)*A.y0;
	NTYP tmp4=(x15)*A.y1;
	minimaxF107ABCD(mi1,ma1,&tmp1,&tmp2,&tmp3,&tmp4);
	
	Pf107 mi2,ma2;
	NTYP tmp5=(x03)*(y03);
	NTYP tmp6=(x03)*(y13);
	NTYP tmp7=(x13)*(y03);
	NTYP tmp8=(x13)*(y13);
	minimaxF107ABCD(mi2,ma2,&tmp5,&tmp6,&tmp7,&tmp8);
	
	Pf107 mi3,ma3;
	NTYP tmp9=A.x0*(y05);
	NTYP tmp10=A.x0*(y15);
	NTYP tmp11=A.x1*(y05);
	NTYP tmp12=A.x1*(y15);
	minimaxF107ABCD(mi3,ma3,&tmp9,&tmp10,&tmp11,&tmp12);
	
	NTYP mi4,ma4;
	minimaxF107AB(mi4,ma4,x04*x02,x14*x12);
	Pf107 mi5,ma5;
	minimaxF107AB(mi5,ma5,&x02,&x12);
	Pf107 mi6,ma6;
	minimaxF107AB(mi6,ma6,&y02,&y12);
	Pf107 mi7,ma7;
	minimaxF107AB(mi7,ma7,&x04,&x14);
	Pf107 mi8,ma8;
	minimaxF107AB(mi8,ma8,&y04,&y14);
	NTYP mi9,ma9;
	minimaxF107AB(mi9,ma9,y04*y02,y14*y12);
	NTYP mi10,ma10;
	minimaxF107AB(mi10,ma10,FAKTORAre*A.x0,FAKTORAre*A.x1);
	NTYP mi11,ma11;
	minimaxF107AB(mi11,ma11,FAKTORAim*A.y0,FAKTORAim*A.y1);
	NTYP mi12,ma12;
	minimaxF107AB(mi12,ma12,FAKTORAre*A.y0,FAKTORAre*A.y1);
	NTYP mi13,ma13;
	minimaxF107AB(mi13,ma13,FAKTORAim*A.x0,FAKTORAim*A.x1);
	
	Pf107 mi14,ma14;
	NTYP tmp13=(*mi5)*(*mi8);
	NTYP tmp14=(*mi5)*(*ma8);
	NTYP tmp15=(*ma5)*(*mi8);
	NTYP tmp16=(*ma5)*(*ma8);
	minimaxF107ABCD(mi14,ma14,&tmp13,&tmp14,&tmp15,&tmp16);
	
	Pf107 mi15,ma15;
	NTYP tmp17=(*mi7)*(*mi6);
	NTYP tmp18=(*mi7)*(*ma6);
	NTYP tmp19=(*ma7)*(*mi6);
	NTYP tmp20=(*ma7)*(*ma6);
	minimaxF107ABCD(mi15,ma15,&tmp17,&tmp18,&tmp19,&tmp20);

	fA.x0=seedC0re+mi10-ma11+mi4-((15*(*ma15)))+(15*(*mi14))-ma9;
	fA.x1=seedC1re+ma10-mi11+ma4-((15*(*mi15)))+(15*(*ma14))-mi9;

	fA.y0=mi12+mi13+6*(*mi1)-(20*((*ma2)))+6*(*mi3)+seedC0im;
	fA.y1=ma12+ma13+6*(*ma1)-(20*((*mi2)))+6*(*ma3)+seedC1im;
	return;
	#endif

	#ifdef _FPA
	NTYP x02; A.x0.squareTo(x02);
	NTYP x03=x02*A.x0;
	NTYP x04; x02.squareTo(x04);
	NTYP x05=x02*x03;
	NTYP y02; A.y0.squareTo(y02);
	NTYP y03=y02*A.y0;
	NTYP y04; y02.squareTo(y04);
	NTYP y05=y02*y03;
	NTYP x12; A.x1.squareTo(x12);
	NTYP x13=x12*A.x1;
	NTYP x14; x12.squareTo(x14);
	NTYP x15=x12*x13;
	NTYP y12; A.y1.squareTo(y12);
	NTYP y13=y12*A.y1;
	NTYP y14; y12.squareTo(y14);
	NTYP y15=y12*y13;
	PFPA mi1,ma1;
	NTYP tmp1=(x05)*A.y0;
	NTYP tmp2=(x05)*A.y1;
	NTYP tmp3=(x15)*A.y0;
	NTYP tmp4=(x15)*A.y1;
	minimaxFPAMIMAABCD(mi1,ma1,&tmp1,&tmp2,&tmp3,&tmp4);
	
	PFPA mi2,ma2;
	NTYP tmp5=(x03)*(y03);
	NTYP tmp6=(x03)*(y13);
	NTYP tmp7=(x13)*(y03);
	NTYP tmp8=(x13)*(y13);
	minimaxFPAMIMAABCD(mi2,ma2,&tmp5,&tmp6,&tmp7,&tmp8);
	
	PFPA mi3,ma3;
	NTYP tmp9=A.x0*(y05);
	NTYP tmp10=A.x0*(y15);
	NTYP tmp11=A.x1*(y05);
	NTYP tmp12=A.x1*(y15);
	minimaxFPAMIMAABCD(mi3,ma3,&tmp9,&tmp10,&tmp11,&tmp12);
	
	PFPA mi4,ma4;
	NTYP tmp13=x04*x02;
	NTYP tmp14=x14*x12;
	minimaxFPAMIMAAB(mi4,ma4,&tmp13,&tmp14);
	
	PFPA mi5,ma5;
	minimaxFPAMIMAAB(mi5,ma5,&x02,&x12);
	PFPA mi6,ma6;
	minimaxFPAMIMAAB(mi6,ma6,&y02,&y12);
	PFPA mi7,ma7;
	minimaxFPAMIMAAB(mi7,ma7,&x04,&x14);
	PFPA mi8,ma8;
	minimaxFPAMIMAAB(mi8,ma8,&y04,&y14);
	PFPA mi9,ma9;
	NTYP tmpb1=y04*y02;
	NTYP tmpb2=y14*y12;
	minimaxFPAMIMAAB(mi9,ma9,&tmpb1,&tmpb2);
	
	PFPA mi10,ma10;
	NTYP tmp15=FAKTORAre*A.x0;
	NTYP tmp16=FAKTORAre*A.x1;
	minimaxFPAMIMAAB(mi10,ma10,&tmp15,&tmp16);
	
	PFPA mi11,ma11;
	NTYP tmp17=FAKTORAim*A.y0;
	NTYP tmp18=FAKTORAim*A.y1;
	minimaxFPAMIMAAB(mi11,ma11,&tmp17,&tmp18);
	
	PFPA mi12,ma12;
	NTYP tmp19=FAKTORAre*A.y0;
	NTYP tmp20=FAKTORAre*A.y1;
	minimaxFPAMIMAAB(mi12,ma12,&tmp19,&tmp20);

	PFPA mi13,ma13;
	NTYP tmp21=FAKTORAim*A.x0;
	NTYP tmp22=FAKTORAim*A.x1;
	minimaxFPAMIMAAB(mi13,ma13,&tmp21,&tmp22);
	
	PFPA mi14,ma14;
	NTYP tmp23=(*mi5)*(*mi8);
	NTYP tmp24=(*mi5)*(*ma8);
	NTYP tmp25=(*ma5)*(*mi8);
	NTYP tmp26=(*ma5)*(*ma8);
	minimaxFPAMIMAABCD(mi14,ma14,&tmp23,&tmp24,&tmp25,&tmp26);

	PFPA mi15,ma15;
	NTYP tmp27=(*mi7)*(*mi6);
	NTYP tmp28=(*mi7)*(*ma6);
	NTYP tmp29=(*ma7)*(*mi6);
	NTYP tmp30=(*ma7)*(*ma6);
	minimaxFPAMIMAABCD(mi15,ma15,&tmp27,&tmp28,&tmp29,&tmp30);

	fA.x0=seedC0re+(*mi10)-(*ma11)+(*mi4)-((15*(*ma15)))+(15*(*mi14))-(*ma9);
	fA.x1=seedC1re+(*ma10)-(*mi11)+(*ma4)-((15*(*mi15)))+(15*(*ma14))-(*mi9);

	fA.y0=(*mi12)+(*mi13)+6*(*mi1)-(20*(*ma2))+6*(*mi3)+seedC0im;
	fA.y1=(*ma12)+(*ma13)+6*(*ma1)-(20*(*mi2))+6*(*ma3)+seedC1im;
	return;
	#endif

	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(3*(5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+3*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)))-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(3*(5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+3*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)))-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+6*minimumD((A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(4*(5*maximumD((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+6*minimumD(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1))+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+6*maximumD((A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(4*(5*minimumD((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+6*maximumD(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1))+seedC1im;
}


// writing bitmap data to file
void write2(FILE *f,const uint8_t a,const uint8_t b) {
	fwrite(&a,1,sizeof(a),f);
	fwrite(&b,1,sizeof(b),f);
}

void write4(FILE *f,const uint8_t a,const uint8_t b,const uint8_t c,const uint8_t d) {
	fwrite(&a,1,sizeof(a),f);
	fwrite(&b,1,sizeof(b),f);
	fwrite(&c,1,sizeof(c),f);
	fwrite(&d,1,sizeof(d),f);
}

void construct_static_reverse_cellgraph(void) {
	for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
		data5->revcgYX[i].howmany=0;
		data5->revcgYX[i].memused=0;
		data5->revcgYX[i].parent=NULL;
	}
	
	PlaneRect A,bbxfA;
	
	const NTYP DD=REVCGBLOCKWIDTH*scaleRangePerPixel;
	int32_t parentx,parenty;
	
	for(int32_t dl=1;dl<=2;dl++) {
		A.y1=COMPLETE0;
		// first pass: calculuate how many are nedded per square
		// 2nd pass: build cell graph as array
		if (dl==1) printf("counting parents ...\n");
		else printf("setting parents to squares ... ");

		for(int32_t y=0;y<SCREENWIDTH;y+=REVCGBLOCKWIDTH) {
			parenty=(y >> REVCGBITS);
			#ifdef _FPA
			A.y0=A.y1;
			FPA_add_ZAB(A.y1,A.y0,DD);
			#endif
			#ifndef _FPA
			A.y0=y*scaleRangePerPixel + COMPLETE0;
			A.y1=A.y0 + DD;
			#endif
		
			A.x1=COMPLETE0;
			for(int32_t x=0;x<SCREENWIDTH;x+=REVCGBLOCKWIDTH) {
				parentx=(x >> REVCGBITS);
				#ifdef _FPA
				A.x0=A.x1;
				FPA_add_ZAB(A.x1,A.x1,DD);
				#endif
				#ifndef _FPA
				A.x0=x*scaleRangePerPixel + COMPLETE0;
				A.x1=A.x0+DD;
				#endif
			
				getBoundingBoxfA(A,bbxfA);
			
				if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA) > 0) {
					continue;
				}

				ScreenRect scr;
				scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
				scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
				scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
				scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
			
				scr.x0 >>= REVCGBITS;
				scr.x1 >>= REVCGBITS;
				scr.y0 >>= REVCGBITS;
				scr.y1 >>= REVCGBITS;
			
				if (dl==1) {
					// just counting
					for(int32_t by=scr.y0;by<=scr.y1;++by) {
						int64_t yoffset=(int64_t)by*REVCGmaxnumber;
						for(int32_t bx=scr.x0;bx<=scr.x1;++bx) {
							data5->revcgYX[yoffset+bx].memused++;
						}
					}
				} else {
					// setting parents to square
					for(int32_t by=scr.y0;by<=scr.y1;++by) {
						int64_t yoffset=(int64_t)by*REVCGmaxnumber;
						for(int32_t bx=scr.x0;bx<=scr.x1;++bx) {
							data5->revcgYX[yoffset+bx].addParent(parentx,parenty);
						}
					}
				}
			} // x
		} // y
	} // passes
}

void compute(void) {
	// build low-level reverse cell graph with the current seed value
	construct_static_reverse_cellgraph();
	
	propagate_definite();
		
	// propagate potentially white
	// takes most of time - especially if no black
	// will finally emerge
	printf("saving raw data ...\n");
	data5->saveRaw("_temp");
	propagate_potw();
	printf("\ncoloring interior cells ...\n");
	int32_t res=color_changeS32(
		SQUARE_GRAY,
		SQUARE_BLACK,
		SQUARE_GRAY_16_CONSECUTIVE,
		SQUARE_BLACK_16_CONSECUTIVE
	);
	if (res>0) interiorpresent=1;
}

void find_special_exterior_hitting_squares(void) {
	// only definite white is marked here
	// no individual pixels are analyzed
	PlaneRect bbxfA,A16;
	
	encgrayx0=encgrayy0=SCREENWIDTH-16;
	encgrayx1=encgrayy1=0;
	
	// initialising gray enclosement per row
		
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		data5->memgrau[y].g0=0;
		data5->memgrau[y].g1=SCREENWIDTH-1;
		data5->memgrau[y].mem0=0;
		data5->memgrau[y].mem1=(SCREENWIDTH >> 4)-1;
	}
	
	// some step size to output progress
	int32_t still_todo0=(SCREENWIDTH >> 4) >> 2;
	int32_t still_todo=1;
	
	const NTYP DD16SCALE=16.0*scaleRangePerPixel;
	A16.y1=COMPLETE0; 
	for(int32_t y16=0;y16<SCREENWIDTH;y16+=16) {
		A16.y0=A16.y1;
		#ifdef _FPA
		FPA_add_ZAB(A16.y1,A16.y0,DD16SCALE);
		#endif
		#ifndef _FPA
		A16.y1=A16.y1+DD16SCALE;
		#endif
		if ( (--still_todo)==0) {
			// output progress
			printf("%i ",SCREENWIDTH-y16);
			still_todo=still_todo0;
		}
		int32_t gray0=SCREENWIDTH-1;
		int32_t gray1=0;

		A16.x1=COMPLETE0; 
		for(int32_t x16=0;x16<SCREENWIDTH;x16+=16) {
			A16.x0=A16.x1;
			#ifdef _FPA
			FPA_add_ZAB(A16.x1,A16.x0,DD16SCALE);
			#endif
			#ifndef _FPA
			A16.x1=A16.x1+DD16SCALE;
			#endif
					
			getBoundingBoxfA(A16,bbxfA);

			// does the 16x16 square lie completely in the special exterior
			DDBYTE w=SQUARE_GRAY;
			if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA)>0) {
				w=SQUARE_WHITE_16_CONSECUTIVE;
				// yes => all can be colored white
			} else {
				// color large square as gray
				if (x16 < gray0) gray0=x16;
				if ( (x16+15) > gray1) gray1=x16+15;
				
				if (y16 < encgrayy0) encgrayy0=y16;
				if ( (y16+15) > encgrayy1) encgrayy1=y16+15;
			}
			
			int32_t memx0=(x16 >> 4);			
			for(int32_t y2=y16;y2<(y16+16);y2++) {
				SETDATA5BYMEM_MY(memx0,y2,w);
			} // y2
					
		} // x16
		
		if (gray0 < encgrayx0) encgrayx0=gray0;
		if (gray1 > encgrayx1) encgrayx1=gray1;
		
		// adjusting gray enclosement per row
		for(int32_t yy=y16;yy<(y16+16);yy++) {
			data5->memgrau[yy].g0=gray0;
			data5->memgrau[yy].g1=gray1;

			// whole row allocated => mem spans the entire
			// pixel coordinate range
			data5->memgrau[yy].mem0=0;
			data5->memgrau[yy].mem1=(SCREENWIDTH >> 4)-1;
		} // yy
		
	} // y16
	
	// adjusting image gray enclosement
	// one 16er block as buffer
	encgrayx0-=16; 
	if (encgrayx0<0) encgrayx0=0;
	planegrayx0=encgrayx0*scaleRangePerPixel + COMPLETE0;
	encgrayy0-=16; 
	if (encgrayy0<0) encgrayy0=0;
	planegrayy0=encgrayy0*scaleRangePerPixel + COMPLETE0;
	encgrayx1+=16; 
	if (encgrayx1>=SCREENWIDTH) encgrayx1=SCREENWIDTH-16;
	planegrayx1=(encgrayx1+16)*scaleRangePerPixel + COMPLETE0;
	encgrayy1+=16; 
	if (encgrayy1>=SCREENWIDTH) encgrayy1=SCREENWIDTH-16;
	planegrayy1=(encgrayy1+16)*scaleRangePerPixel + COMPLETE0;
	
}

void propagate_definite(void) {
	PlaneRect A,bbxfA;
	ScreenRect scr;
	
	if (data5->readtovisit("_in.def.tovisit") <= 0) {
		// mark every tile as to be visited
		for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
			data5->revcgYX[i].tovisit=1;
		}
	}

	int32_t still_todo0=((SCREENWIDTH >> REVCGBITS) >> 1);
	int32_t still_todo=6;
	
	int32_t changed=1;
	int32_t lastsavetime=0;
	
	int64_t checkclockat=checkclockatbbxcount0;

	while (changed>0) {
		changed=0;
		printf("\npropagating definite color ... ");

		for(int32_t y256=0,YBLOCK=0;y256<SCREENWIDTH;y256+=REVCGBLOCKWIDTH,YBLOCK++) {
			int32_t yrevoffset=YBLOCK*REVCGmaxnumber;
			if ( (--still_todo) <= 0) {
				printf("%i ",SCREENWIDTH-y256);
				still_todo=still_todo0;
			}
			if (ctrbbxfa > checkclockat) {
				checkclockat += checkclockatbbxadd;
				int t2=clock();
				if ((t2-lastsavetime) > CLOCK1HOUR) {
					printf("saving raw data as temporary ... ");
					data5->saveRaw("_temp");
					data5->savetovisit("_temp.def.tovisit");
					lastsavetime=t2;
				}
			}
			// too early => jump further
			if ( (y256+REVCGBLOCKWIDTH) < encgrayy0) continue;
			// outside gray enclosement => jump out
			if (y256 > encgrayy1) break;
		
			for(int32_t x256=0,XBLOCK=0;x256<SCREENWIDTH;x256+=REVCGBLOCKWIDTH,XBLOCK++) {
				if ( (data5->revcgYX[yrevoffset+XBLOCK].tovisit<=0) ) continue;
				
				// block has now been checked
				data5->revcgYX[yrevoffset+XBLOCK].tovisit=0;
		
				const int32_t Y256ENDE=y256+REVCGBLOCKWIDTH;
				#ifdef _FPA
				FPA_mul_ZAuvlong(A.y1,scaleRangePerPixel,y256);
				FPA_add_ZAB(A.y1,A.y1,COMPLETE0);
				#endif
				#ifndef _FPA
				A.y1=y256*scaleRangePerPixel + COMPLETE0;
				#endif
				for(int32_t y=y256;y<Y256ENDE;y++) {
					// am Anfang der Schleife wg. mgl. COntinue
					A.y0=A.y1;
					#ifdef _FPA
					FPA_add_ZAB(A.y1,A.y0,scaleRangePerPixel);
					#endif
					#ifndef _FPA
					A.y1=A.y0+scaleRangePerPixel;
					#endif
					/*A.y0=y*scaleRangePerPixel + COMPLETE0;
					A.y1=A.y0 + scaleRangePerPixel;
					*/

					const int32_t xanf=data5->memgrau[y].g0;
					const int32_t xende=data5->memgrau[y].g1;
					// some rows have no gray
					if (xende < xanf) continue;
			
					// gray in that row lies outside of area to be checked
					if ( 
						(x256 > xende) ||
						( (x256+REVCGBLOCKWIDTH) < xanf) 
					) continue;
			
					int32_t wmem=-1 + (x256 >> 4);
					int32_t GLOBALBY=y >> REVCGBITS;
					int32_t GLOBALBYOFFSET=GLOBALBY*REVCGmaxnumber;
			
					for(int32_t x=x256;x<(x256+REVCGBLOCKWIDTH);x+=16) {
						wmem++;

						//uint32_t w=data5->pixelrow[y][wmem];
						uint32_t w;
						GETDATA5BYMEM_MY(wmem,y,w)
				
						// no gray square in this 32-bit integer
						if (
							(w == SQUARE_WHITE_16_CONSECUTIVE) ||
							(w == SQUARE_BLACK_16_CONSECUTIVE) 
						) continue; 
				
						uint32_t wneu=w;
						int32_t w_changed=0;
						#ifdef _FPA
						FPA_mul_ZAuvlong(A.x1,scaleRangePerPixel,x);
						FPA_add_ZAB(A.x1,A.x1,COMPLETE0);
						#endif
						#ifndef _FPA
						A.x1=x*scaleRangePerPixel + COMPLETE0;
						#endif
						for(int32_t wbith=0;wbith<16;wbith++) {
							// muss HIER stehen wg. CONTINUE unten
							
							A.x0=A.x1;
							//A.x0=(x+wbith)*scaleRangePerPixel + COMPLETE0;
							#ifdef _FPA
							FPA_add_ZAB(A.x1,A.x0,scaleRangePerPixel);
							#endif
							#ifndef _FPA
							A.x1=A.x0+scaleRangePerPixel;
							#endif
							uint32_t globalf=w & 0b11;
							w >>= 2;
					
							// current pixel is already gray
							if (
								(globalf != SQUARE_GRAY) 
								//&& (globalf != SQUARE_GRAY_POTENTIALLY_WHITE) 
							) continue;
									
							getBoundingBoxfA(A,bbxfA);
				
							// bounding box in special exterior
							if ((SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(bbxfA))>0) {
							//if ((SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA))>0) {
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
								w_changed=1;
								// passiert das ?
								continue;
							}
					
							int32_t hits_white=0;
							int32_t hits_black=0;
					
							if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
								// overlaps with white region
								hits_white=1;
							}
							
							scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
							scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
							scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
							scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
					
							for(int32_t ty=scr.y0;ty<=scr.y1;ty++) {
								for(int32_t tx=scr.x0;tx<=scr.x1;tx++) {
									int f;
									GET_SINGLE_CELLCOLOR_XY(tx,ty,f);
									switch (f) {
										case SQUARE_BLACK: hits_black=1; break;
										case SQUARE_WHITE: hits_white=1; break;
										default: hits_black=hits_white=1; break;
									}
									
									if ((hits_white>0)&&(hits_black>0)) break;
								} // tx
								
								if ((hits_white>0)&&(hits_black>0)) break;
							} // ty
					
							if ((hits_white>0) && (hits_black==0) ) {
								// only white pixels in the bounding box
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
								w_changed=1;
							} else if ((hits_black>0) && (hits_white==0) ) {
								// only black cells are intersected
								// start can be colored black as well
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_BLACK[wbith]);
								w_changed=1;
							}
						} // wbith
				
						if (w_changed>0) {
							//data5->pixelrow[y][wmem]=wneu;
							SETDATA5BYMEM_MY(wmem,y,wneu)
							// the parents (preimages) leading
							// to that pixel have to be checked in
							// the next while-loop if the newly generated
							// information can help to judge more pixels
							int32_t BX=x >> REVCGBITS;
							int32_t offset=GLOBALBYOFFSET+BX;
							for(int32_t i=0;i<data5->revcgYX[offset].howmany;i++) {
								data5->revcgYX[
									data5->revcgYX[offset].parent[i].BY*REVCGmaxnumber+
									data5->revcgYX[offset].parent[i].BX
								].tovisit=1;
							}

							changed=1;
						}
					} // x
				} // y
			} // X256
		} // Y256
		
	} // while
}

void propagate_potw(void) {
	PlaneRect A,bbxfA;
	ScreenRect scr;

	if (data5->readtovisit("_in.potw.tovisit") <= 0) {
		// mark every tile as to be visited
		for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
			data5->revcgYX[i].tovisit=1;
		}
	}

	int32_t changed=1;
	int32_t lastsavetime=0;
	
	int64_t checkclockat=checkclockatbbxcount0;
	
	int32_t still_todo0=(SCREENWIDTH >> 4) >> 2;
	int32_t still_todo=1;

	while (changed>0) {
		changed=0;
		printf("\npropagating potentially white ... ");
		for(int32_t y256=0,YBLOCK=0;y256<SCREENWIDTH;y256+=REVCGBLOCKWIDTH,YBLOCK++) {
			int32_t yrevoffset=YBLOCK*REVCGmaxnumber;
			if ( (--still_todo) <= 0) {
				printf("%i ",SCREENWIDTH-y256);
				still_todo=still_todo0;
			}
			if (ctrbbxfa > checkclockat) {
				checkclockat += checkclockatbbxadd;
				int t2=clock();
				if ((t2-lastsavetime) > CLOCK1HOUR) {
					printf("saving raw data as temporary ... ");
					data5->saveRaw("_temp");
					data5->savetovisit("_temp.potw.tovisit");
					lastsavetime=t2;
				}
			}

			// too early => jump further
			if ( (y256+REVCGBLOCKWIDTH) < encgrayy0) continue;
			// outside gray enclosement => jump out
			if (y256 > encgrayy1) break;
		
			for(int32_t x256=0,XBLOCK=0;x256<SCREENWIDTH;x256+=REVCGBLOCKWIDTH,XBLOCK++) {
				if ( (data5->revcgYX[yrevoffset+XBLOCK].tovisit<=0) ) continue;
		
				// block has now been checked
				data5->revcgYX[yrevoffset+XBLOCK].tovisit=0;
		
				const int32_t Y256ENDE=y256+REVCGBLOCKWIDTH;
				#ifdef _FPA
				FPA_mul_ZAuvlong(A.y1,scaleRangePerPixel,y256);
				FPA_add_ZAB(A.y1,A.y1,COMPLETE0);
				#endif
				#ifndef _FPA
				A.y1=y256*scaleRangePerPixel + COMPLETE0;
				#endif
				for(int32_t y=y256;y<Y256ENDE;y++) {
					// am Anfang der Schleife wg. mgl. COntinue
					A.y0=A.y1;
					#ifdef _FPA
					FPA_add_ZAB(A.y1,A.y0,scaleRangePerPixel);
					#endif
					#ifndef _FPA
					A.y1=A.y0+scaleRangePerPixel;
					#endif
					/*A.y0=y*scaleRangePerPixel + COMPLETE0;
					A.y1=A.y0 + scaleRangePerPixel;
					*/

					const int32_t xanf=data5->memgrau[y].g0;
					const int32_t xende=data5->memgrau[y].g1;
					// some rows have no gray
					if (xende < xanf) continue;
			
					// gray in that row lies outside of area to be checked
					if ( 
						(x256 > xende) ||
						( (x256+REVCGBLOCKWIDTH) < xanf) 
					) continue;
			
					int32_t wmem=-1 + (x256 >> 4);
					int32_t GLOBALBY=y >> REVCGBITS;
					int32_t GLOBALBYOFFSET=GLOBALBY*REVCGmaxnumber;
			
					for(int32_t x=x256;x<(x256+REVCGBLOCKWIDTH);x+=16) {
						wmem++;

						//uint32_t w=data5->pixelrow[y][wmem];
						uint32_t w;
						GETDATA5BYMEM_MY(wmem,y,w)
				
						// no gray square in this 32-bit integer
						if (
							(w == SQUARE_WHITE_16_CONSECUTIVE) ||
							(w == SQUARE_BLACK_16_CONSECUTIVE) 
						) continue; 
				
						uint32_t wneu=w;
						int32_t w_changed=0;
						#ifdef _FPA
						FPA_mul_ZAuvlong(A.x1,scaleRangePerPixel,x);
						FPA_add_ZAB(A.x1,A.x1,COMPLETE0);
						#endif
						#ifndef _FPA
						A.x1=x*scaleRangePerPixel + COMPLETE0;
						#endif
						for(int32_t wbith=0;wbith<16;wbith++) {
							// muss HIER stehen wg. CONTINUE unten
							
							A.x0=A.x1;
							//A.x0=(x+wbith)*scaleRangePerPixel + COMPLETE0;
							#ifdef _FPA
							FPA_add_ZAB(A.x1,A.x0,scaleRangePerPixel);
							#endif
							#ifndef _FPA
							A.x1=A.x0+scaleRangePerPixel;
							#endif
							uint32_t globalf=w & 0b11;
							w >>= 2;
					
							// current pixel is already gray
							if (
								(globalf != SQUARE_GRAY) 
								//&& (globalf != SQUARE_GRAY_POTENTIALLY_WHITE) 
							) continue;
									
							getBoundingBoxfA(A,bbxfA);
				
							// bounding box in special exterior
							if ((SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(bbxfA))>0) {
							//if ((SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA))>0) {
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
								w_changed=1;
								// passiert das ?
								continue;
							}
					
							int32_t hits_graypotw=0;
							int32_t hits_white=0;
							int32_t hits_gray=0;
							int32_t hits_black=0;
					
							if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
								// overlaps with white region
								hits_white=1;
							}
							
							scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
							scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
							scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
							scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
					
							for(int32_t ty=scr.y0;ty<=scr.y1;ty++) {
								for(int32_t tx=scr.x0;tx<=scr.x1;tx++) {
									int f;
									GET_SINGLE_CELLCOLOR_XY(tx,ty,f);
									switch (f) {
										case SQUARE_BLACK: hits_black=1; break;
										case SQUARE_WHITE: hits_white=1; break;
										case SQUARE_GRAY: hits_gray=1; break;
										case SQUARE_GRAY_POTENTIALLY_WHITE: hits_graypotw=1; break;
									}
									
									if (
										(hits_white>0) &&
										( (hits_black>0) || (hits_gray>0) ) 
									) hits_graypotw=1;
									
									if (hits_graypotw>0) break;
								} // tx
								
								if (hits_graypotw>0) break;
							} // ty
					
							if (hits_graypotw>0) {
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_GRAYPOTW[wbith]);
								w_changed=1;
							}
						} // wbith
				
						if (w_changed>0) {
							//data5->pixelrow[y][wmem]=wneu;
							SETDATA5BYMEM_MY(wmem,y,wneu)
							// the parents (preimages) leading
							// to that pixel have to be checked in
							// the next while-loop if the newly generated
							// information can help to judge more pixels
							int32_t BX=x >> REVCGBITS;
							int32_t offset=GLOBALBYOFFSET+BX;
							for(int32_t i=0;i<data5->revcgYX[offset].howmany;i++) {
								data5->revcgYX[
									data5->revcgYX[offset].parent[i].BY*REVCGmaxnumber+
									data5->revcgYX[offset].parent[i].BX
								].tovisit=1;
							}

							changed=1;
						}
					} // x
				} // y
			} // X256
		} // Y256
		
	} // while
}

int color_changeS32(
	const DDBYTE source1,const DDBYTE target1,
	const DDBYTE source16,const DDBYTE target16
) {
	// swap colors
	int32_t res=0;
	int noch0=SCREENWIDTH >> 3;
	int noch=1;
	
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",y);
			noch=noch0;
		}
		const int32_t xanf=data5->memgrau[y].g0;
		const int32_t xende=data5->memgrau[y].g1;
		int32_t mem=-1 + (xanf >> 4);

		for(int32_t x=xanf;x<=xende;x+=16) {
			mem++; 
			uint32_t w;
			GETDATA5BYMEM_MY(mem,y,w)
			uint32_t wneu=w;
			uint8_t w_changed=0;
			
			if (w == source16) {
				wneu=target16;
				w_changed=1;
				res=1;
			} else {
				for(int32_t b=0;b<16;b++) {
					int32_t f=w & 0b11;
					w >>= 2;
					if (f == source1) {
						w_changed=1;
						switch (target1) {
							case SQUARE_BLACK:
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[b],ARRAY_SQUARE_BLACK[b]);
								break;
							case SQUARE_GRAY_POTENTIALLY_WHITE:
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[b],ARRAY_SQUARE_GRAYPOTW[b]);
								break;
							case SQUARE_WHITE:
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[b],ARRAY_SQUARE_WHITE[b]);
								break;
							default:
								// setting a cell to GRAY
								// is never an error
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[b],ARRAY_SQUARE_GRAY[b]);
								break;
						}
					
						res=1;
					} else if (f == SQUARE_BLACK) {
						// counting interior
						res=1;
					}
				} // b
			} // single
			
			if (w_changed>0) {
				SETDATA5BYMEM_MY(mem,y,wneu)
			}
		} // x
		
	} // y
	
	return res;
}

ParentManager::ParentManager() {
	lastallocated=NULL;
	memused=0;
	freefrom=0;
	anzptr=0;
}

ParentManager::~ParentManager() {
	for(int32_t i=0;i<anzptr;i++) {
		if (ptr[i]) delete[] ptr[i];
	}
	anzptr=0;
}

Parent* ParentManager::getParentSpace(const int32_t aneeded) {
	if (
		(!lastallocated) ||
		( (memused-freefrom) < aneeded)
	) {
		// allocate new memory in 1 GB chunks
		// deleting is done in a dirty manner when
		// program ends
		if (anzptr >= MAXPTR) {
			LOGMSG("Error. Memory parentManager.\n");
			exit(99);
		}
		printf("x");
		int64_t all=CHUNKSIZE / sizeof(Parent);
		ptr[anzptr]=lastallocated=new Parent[all];
		anzptr++;
		if (!lastallocated) {
			LOGMSG("Error/2. Memory parentManager.\n");
			exit(99);
		}
		memused=all;
		freefrom=0;
	}
	
	int32_t ret=freefrom;
	freefrom += aneeded;
	return (&lastallocated[ret]);
}

char* seedCstr(char* erg) {
	#ifdef _FPA
	sprintf(erg,"c_ia_%.20lg_%.20lg_x_%.20lg_%.20lg",
		seedC0re.convert_to_double(),
		seedC1re.convert_to_double(),
		seedC0im.convert_to_double(),
		seedC1im.convert_to_double()
	);
	#endif
	#ifndef _FPA
	sprintf(erg,"c_ia_%.20lg_%.20lg_x_%.20lg_%.20lg",
		(double)seedC0re,(double)seedC1re,(double)seedC0im,(double)seedC1im);
	#endif
	
	return erg;
}

char* FAKTORAstr(char* erg) {
	#ifdef _FPA
	sprintf(erg,"A_%.20lg_%.20lg",
		FAKTORAre.convert_to_double(),
		FAKTORAim.convert_to_double()
	);
	#endif
	#ifndef _FPA
	sprintf(erg,"A_%.20lg_%.20lg",
		(double)FAKTORAre,(double)FAKTORAim);
	#endif
	return erg;
}

char* seedCstr225(char* erg) {
	#ifdef _FPA
	char t1[1024],t2[1024],t3[1024],t4[1024];
	sprintf(erg,"c_ia_%I64d_%I64d_x_%I64d_%I64d",
		seedC0re.str(t1),
		seedC1re.str(t2),
		seedC0im.str(t3),
		seedC1im.str(t4)
	);
	#endif
	#ifndef _FPA
	sprintf(erg,"c_ia_%I64d_%I64d_x_%I64d_%I64d",
		(int64_t)floor(DENOM225*(double)seedC0re),
		(int64_t)floor(DENOM225*(double)seedC1re),
		(int64_t)floor(DENOM225*(double)seedC0im),
		(int64_t)floor(DENOM225*(double)seedC1im)
	);
	#endif
	
	return erg;
}

char* FAKTORAstr225(char* erg) {
	#ifdef _FPA
	char t1[1024],t2[1024];
	sprintf(erg,"A_%I64d_%I64d",
		FAKTORAre.str(t1),
		FAKTORAim.str(t2)
	);
	#endif
	#ifndef _FPA
	sprintf(erg,"A_%I64d_%I64d",
		(int64_t)floor(DENOM225*(double)FAKTORAre),
		(int64_t)floor(DENOM225*(double)FAKTORAim)
	);
	#endif
		
	return erg;
}

void periodicity(char* afn) {
	PDBYTE *dbY=new PDBYTE[SCREENWIDTH];
	DByteMgr *mgr=new DByteMgr; // initialisiert
	
	#define SETDBY(XX,YY,FF) \
	{\
		if (dbY[YY]) {\
			if (\
				(XX >= data5->memgrau[YY].g0) &&\
				(XX <= data5->memgrau[YY].g1) \
			) {\
				dbY[YY][XX-data5->memgrau[YY].g0]=FF;\
			}\
		}\
	}

	int64_t memused=0;
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if (data5->memgrau[y].g1 < data5->memgrau[y].g0) {
			dbY[y]=NULL;
		} else {
			// mem und g muessen NICHT identisch sein
			// hier wird g0/g1 auf durch 16 gebracht
			int m0=(data5->memgrau[y].g0 >> 4);
			int m1=( (data5->memgrau[y].g1 >> 4) + 1);
			data5->memgrau[y].g0 = (m0 << 4);
			data5->memgrau[y].g1 = ( ((m1+1) << 4) - 1);
			
			int32_t xlen=data5->memgrau[y].g1-data5->memgrau[y].g0+1;
			memused += (xlen*sizeof(DBYTE));

			if (xlen>0) {
				dbY[y]=mgr->getMemory(xlen);
				if (!dbY[y]) {
					printf("Fehler. MgrMemoryPeriod ohne Speicher.\n");
					exit(99);
				}
				int32_t dbx=data5->memgrau[y].g0;
				for(int32_t m=m0;m<=m1;m++) {
					DDBYTE w;
					GETDATA5BYMEM_MY(m,y,w)
					for(int32_t b=0;b<16;b++) {
						int32_t f=w & 0b11;
						w >>= 2;
						if (
							(f==SQUARE_GRAY) ||
							(f==SQUARE_GRAY_POTENTIALLY_WHITE)
						) {
							SETDBY(dbx,y,SQUARE_GRAY);
						} else if (f==SQUARE_BLACK) {
							SETDBY(dbx,y,SQUARE_BLACK);
						} else {
							SETDBY(dbx,y,f);
						}
						
						dbx++;
					}
				} 
			} else {
				dbY[y]=NULL;
			}
		}
	}
	printf("periodicity memory used %I64d GB\n",1+(memused >> 30));
	
	FatouComponent *oneorbit=new FatouComponent[MAXFATOUCOMPONENTS];
	int32_t anzfatouinorbit=0;
	// routine-global Components of immediate basins
	ibfcomponents=new FatouComponent[MAXFATOUCOMPONENTS];
	anzibf=0;
	cycles=new Cycle[MAXCYCLES];
	anzcycles=0;
	
	const int32_t MINTEMPCOLOR=256;
	int32_t orbitfcnbr=MINTEMPCOLOR; // >= 256
	int32_t cyclesetnbrimmediate=FATOUCOMPONENTCOLOROFFSET;
	int32_t cyclesetnbrattraction=FATOUCOMPONENTCOLOROFFSET+1;
	DBYTE blobaktiv=FATOUCOMPONENTCOLOROFFSET-1;
	
	int32_t noch0=SCREENWIDTH >> 2;
	int32_t noch=1;
	const VLONG MAXINLISTE=CHUNKSIZE / sizeof(Int2); // a 1 GB Speicher
	Int2 *liste=new Int2[MAXINLISTE];
	VLONG anzliste=0;
	printf("searching for cycles ");
	if (!liste) {
		printf("\n  periodicity search using cube only (slow routine)\n");
	} 
	
	Palette4 periodpal;
	periodpal.setPaletteRGB(SQUARE_BLACK,0,0,0);
	periodpal.setPaletteRGB(SQUARE_WHITE,255,255,255);
	periodpal.setPaletteRGB(SQUARE_GRAY,127,127,127);
	periodpal.setPaletteRGB(SQUARE_GRAY_POTENTIALLY_WHITE,255,0,0);

	// some shuffled version of a heat-map
	double d=0.0;
	double dst=0.19;
	for(int32_t i=FATOUCOMPONENTCOLOROFFSET;i<=255;i+=2) {
		int32_t r,g,b;
		basinpal.getColor(d,r,g,b);
		periodpal.setPaletteRGB(i,r,g,b);
		periodpal.setPaletteRGB(i+1,0.67*r,0.67*g,0.67*b);
		d += dst;
		while (d >= 1.0) d -= 1.0;
	}

	#define SETCYCLE(NR,RR,GG,BB) \
	{\
		periodpal.setPaletteRGB(2*(NR)+FATOUCOMPONENTCOLOROFFSET,RR,GG,BB);\
		periodpal.setPaletteRGB(1+2*(NR)+FATOUCOMPONENTCOLOROFFSET,0.67*RR,0.67*GG,0.67*BB);\
	}
	
	// first cycle is "a bit" random in color
	int32_t rd[7]={0,1,2,3,4,5,6};
	int32_t rand1=clock()%7;
	if (rand1>0) {
		int32_t c=rd[0];
		rd[0]=rd[rand1];
		rd[rand1]=c;
	}
	SETCYCLE(rd[0],0,255,255)
	SETCYCLE(rd[1],255,0,255)
	SETCYCLE(rd[2],255,0,0)
	SETCYCLE(rd[3],0,255,0)
	SETCYCLE(rd[4],255,255,0)
	SETCYCLE(rd[5],193,193,255)
	SETCYCLE(rd[6],193,63,255)
	
	#define ORDBY(XX,YY,FF) \
	{\
		if (dbY[YY]) {\
			if (\
				(XX >= data5->memgrau[YY].g0) &&\
				(XX <= data5->memgrau[YY].g1) \
			) {\
				dbY[YY][XX-data5->memgrau[YY].g0] |= FF;\
			}\
		}\
	}
	
	#define ANDDBY(XX,YY,FF) \
	{\
		if (dbY[YY]) {\
			if (\
				(XX >= data5->memgrau[YY].g0) &&\
				(XX <= data5->memgrau[YY].g1) \
			) {\
				dbY[YY][XX-data5->memgrau[YY].g0] &= FF;\
			}\
		}\
	}

	#define GETDBY(XX,YY,ERG) \
	{\
		if (dbY[YY]) {\
			if (\
				(XX < data5->memgrau[YY].g0) ||\
				(XX > data5->memgrau[YY].g1) \
			) {\
				ERG=SQUARE_WHITE;\
			} else {\
				ERG=dbY[YY][XX-data5->memgrau[YY].g0];\
			}\
		} else ERG=SQUARE_WHITE;\
	}

	#define ADDLISTE(XX,YY) \
	{\
		if ( (liste) && (anzliste < (MAXINLISTE-8)) ) {\
			liste[anzliste].x=XX;\
			liste[anzliste].y=YY;\
			anzliste++;\
		}\
	}
	
	#define STRAHLEN(BX,BY) \
	{\
		for(int xx=(BX-1);xx>=0;xx--) {\
			int ff;\
			GETDBY(xx,BY,ff);\
			if (ff == SQUARE_BLACK) {\
				ADDLISTE(xx,BY)\
				if (xx < bx0) bx0=xx;\
				if (xx > bx1) bx1=xx;\
				SETDBY(xx,BY,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t xx=(BX+1);xx<SCREENWIDTH;xx++) {\
			int32_t ff;\
			GETDBY(xx,BY,ff);\
			if (ff == SQUARE_BLACK) {\
				ADDLISTE(xx,BY)\
				if (xx < bx0) bx0=xx;\
				if (xx > bx1) bx1=xx;\
				SETDBY(xx,BY,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t yy=(BY+1);yy<SCREENWIDTH;yy++) {\
			int32_t ff;\
			GETDBY(BX,yy,ff);\
			if (ff == SQUARE_BLACK) {\
				ADDLISTE(BX,yy)\
				if (yy < by0) by0=yy;\
				if (yy > by1) by1=yy;\
				SETDBY(BX,yy,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t yy=(BY-1);yy>=0;yy--) {\
			int32_t ff;\
			GETDBY(BX,yy,ff);\
			if (ff == SQUARE_BLACK) {\
				ADDLISTE(BX,yy)\
				if (yy < by0) by0=yy;\
				if (yy > by1) by1=yy;\
				SETDBY(BX,yy,blobaktiv);\
			} else break;\
		}\
		\
	}
	
	// exchange colors from temporary to final
	#define COLORFC(FC,QF,ZF) \
	{\
		for(int32_t yy=FC.scrc.y0;yy<=FC.scrc.y1;yy++) {\
			for(int32_t xx=FC.scrc.x0;xx<=FC.scrc.x1;xx++) {\
				int32_t ff;\
				GETDBY(xx,yy,ff);\
				if (ff == QF) {\
					SETDBY(xx,yy,ZF);\
				}\
			}\
		}\
	}
	
	int32_t maxorbitlen=0;
	
	for(int32_t yb=0;yb<SCREENWIDTH;yb++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-yb);
			noch=noch0;
		}
		for(int32_t xb=0;xb<SCREENWIDTH;xb++) {
			int32_t ff;
			GETDBY(xb,yb,ff);
			if (ff != SQUARE_BLACK) continue;
				
			int32_t x=xb,y=yb;
			orbitfcnbr=MINTEMPCOLOR;
			anzfatouinorbit=0; // starten
				
			// no Voxel: color >= MINTEMPCOLOR
			// no voxel: color blobaktiv
				
			while (1) {
				// current orbit ongoing
				// starts at x,y
				
				int32_t bx0=x,bx1=x,by0=y,by1=y;
				int32_t currorbitidx=anzfatouinorbit;
				oneorbit[currorbitidx].currentOrbitColorIdxTemp=orbitfcnbr;
				oneorbit[currorbitidx].inCycleNbr=-1; 
				oneorbit[currorbitidx].isimmediate=0;
				anzfatouinorbit++;
			
				if (anzliste != 0) {
					LOGMSG3("Implementation error. new blob#%i, but list with %I64d elements\n",anzfatouinorbit,anzliste);
					exit(99);
				}
				ADDLISTE(x,y)
				bx0=bx1=x;
				by0=by1=y;
				SETDBY(x,y,blobaktiv);
				
				// Floodfill 
				int32_t changed=1;
				while (changed>0) {
					changed=0;
					
					int32_t ey0=by0,ey1=by1,ex0=bx0,ex1=bx1;
					
					if (anzliste>0) {
						int32_t lx=liste[anzliste-1].x;
						int32_t ly=liste[anzliste-1].y;
						if (lx < bx0) bx0=lx;
						if (lx > bx1) bx1=lx;
						if (ly < by0) by0=ly;
						if (ly > by1) by1=ly;
						anzliste--;
						int32_t lff;
						GETDBY(lx,ly,lff);
						if (lff == blobaktiv) {
							SETDBY(lx,ly,orbitfcnbr);
							STRAHLEN(lx,ly)
						}
						
						// go on till list is empty
						changed=1;
						continue;
					};
					
					// if list is empty or has changed the cube
					// check for voxels to follow
					if ((anzliste>0) || (changed>0)) {
						LOGMSG2("Implementation error. list mistake %I64d elements, 0 expected\n",anzliste);
						exit(99);
					} else {
						// per cube
						for(int32_t by=ey0;by<=ey1;by++) {
							for(int32_t bx=ex0;bx<=ex1;bx++) {
								int32_t bff;
								GETDBY(bx,by,bff);
								if (bff != blobaktiv) continue;
									SETDBY(bx,by,orbitfcnbr);
									changed=1;
									if (bx < bx0) bx0=bx;
									if (bx > bx1) bx1=bx;
									if (by < by0) by0=by;
									if (by > by1) by1=by;
								
									STRAHLEN(bx,by)
								} // bx
							} // by
						} // per
					} // while
			
					// enclosement cube of current Fatou component
					// SETCOL
					oneorbit[currorbitidx].scrc.x0=bx0;
					oneorbit[currorbitidx].scrc.x1=bx1;
					oneorbit[currorbitidx].scrc.y0=by0;
					oneorbit[currorbitidx].scrc.y1=by1;
					
					// follow ONE arbitrary voxel to
					// find the target Fatou component
					PlaneRect A,bbxfA;
					A.x0=x*scaleRangePerPixel+COMPLETE0;
					A.x1=A.x0+scaleRangePerPixel;
					A.y0=y*scaleRangePerPixel+COMPLETE0;
					A.y1=A.y0+scaleRangePerPixel;
					
					getBoundingBoxfA(A,bbxfA);
					
					if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA) > 0) {
						LOGMSG("Implementation error. No target Fatou component.\n");
						exit(99);
					}
				
					ScreenRect scr;
					scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
					scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
					scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
					scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
					
					if (
						(scr.x0<0) || (scr.x0 >= SCREENWIDTH) ||
						(scr.y0<0) || (scr.y0 >= SCREENWIDTH) 
					) {
						LOGMSG("Implementation error. BbxfA is inconsistent.\n");
						exit(99);
					}
					
					int32_t vf;
					GETDBY(scr.x0,scr.y0,vf);
					
					if (vf == SQUARE_BLACK) {
						// new Fatou component found
						x=scr.x0;
						y=scr.y0;
						orbitfcnbr++;
						
						continue;
					} else if (vf >= MINTEMPCOLOR) {
						// new cycle found
						LOGMSG("\n  Cycle ");
						
						int32_t o0=-1;
						for(int32_t oi=0;oi<=currorbitidx;oi++) {
							if (oneorbit[oi].currentOrbitColorIdxTemp == vf) {
								o0=oi;
								break;
							}
						} // oi
						if (o0 < 0) {
							LOGMSG("Implementation error. No cycle in orbit\n");
							LOGMSG4("anzcycles %i currorbitidx %i currentOrbitColorIdxTemp %i\n",
								anzcycles,currorbitidx,oneorbit[currorbitidx].currentOrbitColorIdxTemp);
							exit(99);
						}
						
						cycles[anzcycles].len=currorbitidx-o0+1;

						LOGMSG2("len %i found  ",cycles[anzcycles].len);
						// up to first cfyclic point => attraction basin
						for(int32_t oi=0;oi<o0;oi++) {
							COLORFC(oneorbit[oi],oneorbit[oi].currentOrbitColorIdxTemp,cyclesetnbrattraction);
						}

						cycles[anzcycles].fatouidx0=anzibf;

						for(int oi=o0;oi<=currorbitidx;oi++) {
							// immediate basins
							ibfcomponents[anzibf].scrc.x0=oneorbit[oi].scrc.x0;
							ibfcomponents[anzibf].scrc.x1=oneorbit[oi].scrc.x1;
							ibfcomponents[anzibf].scrc.y0=oneorbit[oi].scrc.y0;
							ibfcomponents[anzibf].scrc.y1=oneorbit[oi].scrc.y1;
							ibfcomponents[anzibf].currentOrbitColorIdxTemp=cyclesetnbrimmediate;
							ibfcomponents[anzibf].inCycleNbr=anzcycles;
							ibfcomponents[anzibf].isimmediate=1;
							anzibf++;
							COLORFC(oneorbit[oi],oneorbit[oi].currentOrbitColorIdxTemp,cyclesetnbrimmediate);
						} // copy
						
						cycles[anzcycles].fatouidx1=-1 + anzibf;
						cycles[anzcycles].immediateBasinColorIdx=cyclesetnbrimmediate;
						cycles[anzcycles].attractionBasinColorIdx=cyclesetnbrattraction;
						anzcycles++;
						if (anzcycles >= (MAXCYCLES-2)) {
							LOGMSG("Not possible: Too many cycles (yet to be implemented).\n");
							exit(99);
						}
						
						cyclesetnbrimmediate+=2;;
						cyclesetnbrattraction+=2;
						if (cyclesetnbrattraction >= MINTEMPCOLOR) {
							LOGMSG("Implementation error. Too many cycles.\n");
							exit(99);
						}
						
						if (currorbitidx > maxorbitlen) maxorbitlen=currorbitidx;
						break; // beendet
					} else {
						// orbit lands in already found cycle
						int32_t zyklus=-1;
						for(int32_t cyc=0;cyc<anzcycles;cyc++) {
							if (
								(cycles[cyc].immediateBasinColorIdx == vf) ||
								(cycles[cyc].attractionBasinColorIdx == vf)
							) {
								zyklus=cyc;
								break;
							}
						} // cyc
						
						if (zyklus<0) {
							LOGMSG("Implementation error. Found cycle not detected in orbit\n");
							exit(99);
						}
						
						// whole orbit is attraction basin
						for(int32_t oi=0;oi<=currorbitidx;oi++) {
							COLORFC(oneorbit[oi],oneorbit[oi].currentOrbitColorIdxTemp,cycles[zyklus].attractionBasinColorIdx);
						}
						
					break;
				}
			} // while in Orbit
		} // xb
	} // yb
	
	LOGMSG3("\n%i cycles (max. orbit length %i)\n",anzcycles,maxorbitlen);
	for(int32_t i=0;i<anzcycles;i++) {
		LOGMSG2("  Cycle #%i: ",i);
		LOGMSG5("len=%i immediate RGB(%i,%i,%i) ",cycles[i].len,
			periodpal.rgbs[cycles[i].immediateBasinColorIdx].R,
			periodpal.rgbs[cycles[i].immediateBasinColorIdx].G,
			periodpal.rgbs[cycles[i].immediateBasinColorIdx].B);
		LOGMSG4("attraction RGB(%i,%i,%i)\n",
			periodpal.rgbs[cycles[i].attractionBasinColorIdx].R,
			periodpal.rgbs[cycles[i].attractionBasinColorIdx].G,
			periodpal.rgbs[cycles[i].attractionBasinColorIdx].B);
	}
	
	// save as (downscaled) image, at most 2^16 x 2^16 pixels
	char tt[1024];
	int32_t DSEXPONENT=0;
	VLONG totalsize=(VLONG)SCREENWIDTH * SCREENWIDTH;
	VLONG fourgb=( (VLONG)1 << 32 );
	//fourgb=(1 << 20);
	while (totalsize > fourgb) {
		DSEXPONENT++;
		totalsize >>= 1;
	}
	int32_t TWDSTEP=(1 << DSEXPONENT); // n x n pixels are trustworthily downscaled
	int32_t dslen=(SCREENWIDTH >> DSEXPONENT);

	// so that the final image has at max 4 GB
	int32_t ybytes=dslen; // 1 Pixel pro Byte
	ybytes=(int)(4*ceil(ybytes*0.25));
	char* rgbz=new char[dslen+16];

	if (TWDSTEP>1) {
		sprintf(tt,"%s_period_twd_%i_fold.bmp",afn,DSEXPONENT);
	} else {
		sprintf(tt,"%s_period.bmp",afn);
	}
	
	if (DSEXPONENT>0) printf("  trustworthily downscaled 2^%i-fold\n",DSEXPONENT);

	unsigned int off
		=		14 // FILEHeader
			+	40 // Bitmapheader
			+	256*4 // ColorPalette
		;
	unsigned int filelen
			=	off
			+	(ybytes*dslen);
		;

	FILE *fbmp=fopen(tt,"wb");
	
	#define BMPHEADER1(FBMP) \
	{\
		write2(FBMP,66,77); \
		fwrite(&filelen,1,sizeof(filelen),FBMP);\
		write4(FBMP,0,0,0,0);\
		fwrite(&off,1,sizeof(off),FBMP); \
		write4(FBMP,40,0,0,0);\
		unsigned int w = dslen;\
		fwrite(&w,sizeof(w),1,FBMP);\
		w = dslen;\
		fwrite(&w,sizeof(w),1,FBMP);\
		write2(FBMP,1,0);\
		write2(FBMP,8,0);\
		write4(FBMP,0,0,0,0);\
		write4(FBMP,0,0,0,0);\
		write4(FBMP,19,10,0,0);\
		write4(FBMP,19,10,0,0);\
		write4(FBMP,0,1,0,0);\
		write4(FBMP,0,0,0,0);\
		fwrite(periodpal.rgbs,sizeof(RGB4),256,FBMP);\
	}
	
	BMPHEADER1(fbmp)
	
	for(int32_t y=0;y<SCREENWIDTH;y+=TWDSTEP) {
		int32_t setx=-1;
		for(int32_t x=0;x<SCREENWIDTH;x+=TWDSTEP) {
			setx++;
			int32_t f=-1;
			for(int32_t dy=0;dy<TWDSTEP;dy++) {
				for(int32_t dx=0;dx<TWDSTEP;dx++) {
					int32_t tmpf;
					GETDBY(x+dx,y+dy,tmpf);
					if (tmpf == SQUARE_GRAY_POTENTIALLY_WHITE) tmpf=SQUARE_GRAY;
					if (f<0) f=tmpf;
					else if (f != tmpf) {
						f=SQUARE_GRAY;
						break;
					}
					if (f == SQUARE_GRAY) break;
				}
				if (f == SQUARE_GRAY) break;				
			}
			if ( (f<0) || (f>=256) ) {
				LOGMSG2("Periodicity. Farbfehler %i\n",f);
				exit(99);
			}
			rgbz[setx]=f;
		}
		fwrite(rgbz,sizeof(char),dslen,fbmp);
	}
		
	fclose(fbmp);
	
	delete[] oneorbit;
	
	#define ENDROUTINE \
	{\
		delete[] rgbz;\
		delete mgr; \
		mgr=NULL;\
		delete[] dbY;\
	}

	if (_PERIODICPOINTS <= 0) {
		ENDROUTINE
		return;
	}
	
	// ibfcomponents und cycles BLEIBEN
	// für äußere Anwendung
	
	// color value in dbY reflects CYCLE, not the individual
	// fatou component id
	
	const DBYTE PP_FATOUOHNEFLAGS=0b0001111111111111;
	const DBYTE PP_FLAGS         =0b1110000000000000;
	const DBYTE PP_TOVISIT       =0b1000000000000000;
	const DBYTE PP_UN_TOVISIT    =0b0111111111111111;
	const DBYTE PP_VISITED       =0b0100000000000000;
	const DBYTE PP_UN_VISITED    =0b1011111111111111;
	const DBYTE PP_POSSIBLEPER   =0b0010000000000000;
	const DBYTE PP_UN_POSSIBLEPER=0b1101111111111111;
	
	printf("following interior points ... ");
	
	#define DBYIMG(FN) \
	{\
		FILE *f=fopen(FN,"wb");\
		BMPHEADER1(f)\
		for(int32_t y=0;y<SCREENWIDTH;y++) {\
			for(int32_t x=0;x<SCREENWIDTH;x++) {\
				int32_t ff;\
				GETDBY(x,y,ff);\
				if (\
					(ff == SQUARE_WHITE) || \
					(ff == SQUARE_GRAY) ||\
					(ff == SQUARE_GRAY_POTENTIALLY_WHITE)\
				) {\
					rgbz[x]=ff;\
				} else {\
					if ( (ff & PP_POSSIBLEPER) != 0) {\
						rgbz[x]=SQUARE_BLACK;\
					} else {\
						int32_t f=ff & PP_FATOUOHNEFLAGS;\
						if (f <= 255) {\
							rgbz[x]=f;\
						} else {\
							LOGMSG("Farbfehler dbyimg\n");\
							exit(99);\
						}\
					}\
				}\
			}\
			fwrite(rgbz,sizeof(char),ybytes,f);\
		}\
		fclose(f);\
	}
	
	delete[] liste;
	liste=NULL;
	anzliste=0;

	anzliste=0;
	ListeDFS orbit;
	ListeFIFO possibleper; 
	int32_t ppx,ppy;
	ScreenRect* ppscr=new ScreenRect[MAXPERIODICPOINTS];
	int32_t rausausfc=0;
					
	for(int cyc=0;cyc<anzcycles;cyc++) {
		printf("\ncycle #%i periodic start ",cyc);
		int PRELEN=-1 + cycles[cyc].len;
		ppx=ppy=-1;

		int32_t fc=cycles[cyc].fatouidx0;
		int64_t minarea=(int64_t)(ibfcomponents[fc].scrc.x1-ibfcomponents[fc].scrc.x0) * (ibfcomponents[fc].scrc.y1-ibfcomponents[fc].scrc.y0);
		for(int32_t fctest=(cycles[cyc].fatouidx0+1);fctest<=cycles[cyc].fatouidx1;fctest++) {
			int64_t area=(int64_t)(ibfcomponents[fctest].scrc.x1-ibfcomponents[fctest].scrc.x0) * (ibfcomponents[fctest].scrc.y1-ibfcomponents[fctest].scrc.y0);
			if (area < minarea) {
				minarea=area;
				fc=fctest;
			}
		}
		
		rausausfc=0;
		int32_t noch=1;
		int32_t noch0=4;
		// as periodic points often lie deep within a fatou component
		// I start at the geometric center and alternate between +-d
		int32_t my=((int64_t)ibfcomponents[fc].scrc.y1+(int64_t)ibfcomponents[fc].scrc.y0) >> 1;
		int32_t deltay=0;
			
		#define NEXTDELTAY \
		{\
			if (deltay==0) deltay=1;\
			else if (deltay>0) deltay=-deltay;\
			else deltay=-deltay + 1;\
		}
			
		#define NEXTDELTAX \
		{\
			if (deltax==0) deltax=1;\
			else if (deltax>0) deltax=-deltax;\
			else deltax=-deltax + 1;\
		}

		for(int32_t vy=ibfcomponents[fc].scrc.y0;vy<=ibfcomponents[fc].scrc.y1;vy++) {
			int32_t iby=my+deltay;
			if (
				(iby<ibfcomponents[fc].scrc.y0) ||
				(iby>ibfcomponents[fc].scrc.y1)
			) {
				NEXTDELTAY
				continue;
			}
			
			// naechstes Delta
			NEXTDELTAY
			
			if ((--noch)<=0) {
				noch=noch0;
				printf(".");
			}
			
			int32_t mx=((int64_t)ibfcomponents[fc].scrc.x1+(int64_t)ibfcomponents[fc].scrc.x0) >> 1;
			int32_t deltax=0;

			for(int32_t vx=ibfcomponents[fc].scrc.x0;vx<=ibfcomponents[fc].scrc.x1;vx++) {
				int32_t ibx=mx+deltax;
				if (
					(ibx<ibfcomponents[fc].scrc.x0) ||
					(ibx>ibfcomponents[fc].scrc.x1)
				) {
					NEXTDELTAX
				}
			
				NEXTDELTAX

				int32_t iff;
				GETDBY(ibx,iby,iff);
				if (
					(iff & PP_FATOUOHNEFLAGS) != cycles[cyc].immediateBasinColorIdx
				) continue;
				
				orbit.start();
				orbit.write(ibx,iby,0);
					
				int32_t ox,oy;
				DBYTE ot;
				while (orbit.read(ox,oy,ot) > 0) {
					if (ot >= cycles[cyc].len) continue;
					
					PlaneRect A,bbxfA;
					A.x0=ox*scaleRangePerPixel + COMPLETE0;
					A.x1=A.x0+scaleRangePerPixel;
					A.y0=oy*scaleRangePerPixel + COMPLETE0;
					A.y1=A.y0+scaleRangePerPixel;
					
					getBoundingBoxfA(A,bbxfA);
					// bxfA in complete
					
					ScreenRect scr;
					scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
					scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
					scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
					scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
					
		#define PPUEBERLAPPTMITIB(IBXX,IBYY,BX,BY,ERG) \
		{\
			ERG=0;\
			int32_t dx=(IBXX)-(BX);\
			int32_t dy=(IBYY)-(BY);\
			if (\
				(dx>=-1) && (dx<=1) &&\
				(dy>=-1) && (dy<=1)\
			) {\
				ERG=1;\
			}\
		}
					rausausfc=0;
					for(int32_t by=scr.y0;by<=scr.y1;by++) {
						for(int32_t bx=scr.x0;bx<=scr.x1;bx++) {
							int32_t erg;
							if (ot == PRELEN) {
								PPUEBERLAPPTMITIB(ibx,iby,bx,by,erg);
								if (erg>0) {
									ppx=ibx;
									ppy=iby;
									rausausfc=1;
									break;
								}
							} else {
								orbit.write(bx,by,ot+1);
							}
								
						} // bx
						
						if (rausausfc>0) break;
					} // by
						
					if (rausausfc>0) break;
					
				} // orbit while
					
				if (rausausfc>0) break;
			} // ibx
			if (rausausfc>0) break;
		} // iby
		
		possibleper.start();
		possibleper.write(ppx,ppy);

		int32_t ff;
		GETDBY(ppx,ppy,ff);
		ff |= PP_TOVISIT;
		ff |= PP_POSSIBLEPER;
		SETDBY(ppx,ppy,ff);
		
		int32_t wx,wy;
		while (possibleper.read(wx,wy) > 0) {
			// bereits VISITED
			int32_t ff;
			GETDBY(wx,wy,ff);
			if (
				( (ff & PP_TOVISIT) == 0)
			) continue;
			
			ff &= PP_UN_TOVISIT;
			ff |= PP_VISITED;
			SETDBY(wx,wy,ff)
			
			PlaneRect A,bbxfA;
			A.x0=wx*scaleRangePerPixel + COMPLETE0;
			A.x1=A.x0+scaleRangePerPixel;
			A.y0=wy*scaleRangePerPixel + COMPLETE0;
			A.y1=A.y0+scaleRangePerPixel;
						
			getBoundingBoxfA(A,bbxfA);
			ScreenRect scr;
			scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
			scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
			scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
			scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
						
			for(int32_t by=scr.y0;by<=scr.y1;by++) {
				for(int32_t bx=scr.x0;bx<=scr.x1;bx++) {
					int32_t bff;
					GETDBY(bx,by,bff);
					if (
						((bff & PP_TOVISIT) != 0) ||
						((bff & PP_VISITED) != 0)
					) continue;
					
					bff |= (PP_TOVISIT | PP_POSSIBLEPER);
					SETDBY(bx,by,bff)
					possibleper.write(bx,by);
				} // bx
			} // by
			
		} // while

		// look for adjacent components of identified
		// ppscr
		int32_t anzpp=0;
		for(int32_t i=0;i<possibleper.next_writepos;i++) {
			int32_t idx=-1;
			int32_t sx=possibleper.werte[i].x;
			int32_t sy=possibleper.werte[i].y;

			for(int k=0;k<anzpp;k++) {
				if (
					(ppscr[k].x0 <= sx) &&
					(sx <= ppscr[k].x1) &&
					(ppscr[k].y0 <= sy) &&
					(sy <= ppscr[k].y1)
				) {
					idx=k;
					break;
				}
					
				// adjacent
				
				if (
					(ppscr[k].x0 <= sx) &&
					(sx <= ppscr[k].x1)
				) {
					if (sy == (ppscr[k].y0-1)) {
						ppscr[k].y0--;
						idx=k;
						break;
					} else
					if (sy == (ppscr[k].y1+1)) {
						ppscr[k].y1++;
						idx=k;
						break;
					}
				}
				
				if (
					(ppscr[k].y0 <= sy) &&
					(sy <= ppscr[k].y1)
				) {
					if (sx == (ppscr[k].x0-1)) {
						ppscr[k].x0--;
						idx=k;
						break;
					} else
					if (sx == (ppscr[k].x1+1)) {
						ppscr[k].x1++;
						idx=k;
						break;
					}
				}

			}
			
			if (idx>=0) continue;
			
			// go left, right etc as long as in the region
			int32_t ccx0=sx;
			while (ccx0>=0) {
				int32_t bff;
				GETDBY(ccx0,sy,bff);
				if (((bff & PP_POSSIBLEPER) == 0)) {
					ccx0++;
					break;
				}
				ccx0--;
			} // while
			
			int32_t ccx1=sx;
			while (ccx1<SCREENWIDTH) {
				int32_t bff;
				GETDBY(ccx1,sy,bff);
				if (((bff & PP_POSSIBLEPER) == 0)) {
					ccx1--;
					break;
				}
				ccx1++;
			} // while

			int32_t ccy0=sy;
			while (ccy0>=0) {
				int32_t bff;
				GETDBY(sx,ccy0,bff);
				if (((bff & PP_POSSIBLEPER) == 0)) {
					ccy0++;
					break;
				}
				ccy0--;
			} // while
			
			int32_t ccy1=sy;
			while (ccy1<SCREENWIDTH) {
				int32_t bff;
				GETDBY(sx,ccy1,bff);
				if (((bff & PP_POSSIBLEPER) == 0)) {
					ccy1--;
					break;
				}
				ccy1++;
			} // while

			if (anzpp >= (MAXPERIODICPOINTS-8)) {
				LOGMSG("Error. Too many periodic point regions.\n");
				return;
			}
			
			ppscr[anzpp].x0=ccx0;
			ppscr[anzpp].x1=ccx1;
			ppscr[anzpp].y0=ccy0;
			ppscr[anzpp].y1=ccy1;
			anzpp++;
			
			int32_t changed=1;
			if (anzpp==1) changed=0; // nix zu tun
			while (changed>0) {
				changed=0;
				for(int32_t p0=0;p0<=(anzpp-2);p0++) {
					// uberlappt bzw. grenzt an p0 mit dem LETZTEN
					// ja => fusionieren und in p0 setzen
					if (
						(ppscr[p0].x1 < (ppscr[anzpp-1].x0-1)) ||
						(ppscr[p0].x0 > (ppscr[anzpp-1].x1+1)) ||
						(ppscr[p0].y1 < (ppscr[anzpp-1].y0-1)) ||
						(ppscr[p0].y0 > (ppscr[anzpp-1].y1+1)) ||
						(ppscr[anzpp-1].x1 < (ppscr[p0].x0-1)) ||
						(ppscr[anzpp-1].x0 > (ppscr[p0].x1+1)) ||
						(ppscr[anzpp-1].y1 < (ppscr[p0].y0-1)) ||
						(ppscr[anzpp-1].y0 > (ppscr[p0].y1+1))
					) {
						continue;
					}
					
					ppscr[p0].x0=minimumI(ppscr[p0].x0,ppscr[anzpp-1].x0);
					ppscr[p0].x1=maximumI(ppscr[p0].x1,ppscr[anzpp-1].x1);
					ppscr[p0].y0=minimumI(ppscr[p0].y0,ppscr[anzpp-1].y0);
					ppscr[p0].y1=maximumI(ppscr[p0].y1,ppscr[anzpp-1].y1);
					changed=1;
					anzpp--;
					break;
				} // p0
			} // changed
			
		} // i
		
		LOGMSG2("\n%i possible periodic regions\n",anzpp);
		for(int32_t p=0;p<anzpp;p++) {
			#ifdef _FPA
			NTYP w1=(ppscr[p].x0*scaleRangePerPixel+COMPLETE0);
			NTYP w2=((ppscr[p].x1+1)*scaleRangePerPixel+COMPLETE0);
			NTYP w3=(ppscr[p].y0*scaleRangePerPixel+COMPLETE0);
			NTYP w4=((ppscr[p].y1+1)*scaleRangePerPixel+COMPLETE0);

			LOGMSG6("#%i: [%.20lg..%.20lg] x [%.20lg..%.20lg]\n",
				p,
				w1.convert_to_double(),
				w2.convert_to_double(),
				w3.convert_to_double(),
				w4.convert_to_double()
			);
			#endif
			#ifndef _FPA
			LOGMSG6("#%i: [%.20lg..%.20lg] x [%.20lg..%.20lg]\n",
				p,
				(double)(ppscr[p].x0*scaleRangePerPixel+COMPLETE0),
				(double)((ppscr[p].x1+1)*scaleRangePerPixel+COMPLETE0),
				(double)(ppscr[p].y0*scaleRangePerPixel+COMPLETE0),
				(double)((ppscr[p].y1+1)*scaleRangePerPixel+COMPLETE0)
			);
			#endif
		} // p
				
	} // cyc
	
	sprintf(tt,"%s_periodic_points.bmp",afn);
	DBYIMG(tt)
	
	delete[] ppscr;
	ENDROUTINE
}

// struct ColorPalette

ColorPalette::ColorPalette() {
	anz=0;
	rgbs=NULL;
}

ColorPalette::~ColorPalette() {
	if ((rgbs)&&(anz>0)) delete[] rgbs;
	rgbs=NULL;
	anz=0;
}

void ColorPalette::setlen(const int al) {
	if (rgbs) delete[] rgbs;
	rgbs=new RGB4[al];
	anz=al;
}

void ColorPalette::setInterval(
	const double p0,const double p1,
	const int32_t ar,const int32_t ag,const int32_t ab,
	const int32_t br,const int32_t bg,const int32_t bb) {

	int32_t i0=(int32_t)floor(p0*anz); 
	if (i0<0) i0=0;
	int32_t i1=(int32_t)floor(p1*anz); 
	if (i1>= anz) i1=anz-1;

	double dr=br-ar; dr /= (i1-i0);
	double dg=bg-ag; dg /= (i1-i0);
	double db=bb-ab; db /= (i1-i0);

	for(int32_t i=i0;i<=i1;i++) {
		rgbs[i].R=ar + (int32_t)floor(dr*(i-i0));
		rgbs[i].G=ag + (int32_t)floor(dg*(i-i0));
		rgbs[i].B=ab + (int32_t)floor(db*(i-i0));
	}
}

void ColorPalette::getColor(const double p,int32_t& fr,int32_t& fg,int32_t& fb) {
	int32_t idx=(int)floor(p*anz);
	if (idx<0) idx=0;
	if (idx>=anz) idx=anz-1;
	if (!rgbs) {
		fr=fg=fb=0;
	} else {
		fr=rgbs[idx].R;
		fg=rgbs[idx].G;
		fb=rgbs[idx].B;
	}
};

// struct memory Manager

DByteMgr::DByteMgr() {
	current=NULL;
	allocatedIdx=0;
	freeFromIdx=-1;
	// in 1 GB chuncks
	double d=CHUNKSIZE; 
	d /= sizeof(DBYTE);
	allocatePerBlockIdx=(int32_t)floor(d);
	anzptr=0;
}

DBYTE* DByteMgr::getMemory(const int32_t aanz) {
	if (
		(!current) ||
		((freeFromIdx + aanz) >= allocatedIdx)
	) {
		if (anzptr >= MAXPTR) {
			LOGMSG("Error. Memory-DByteMgr.\n");
			exit(99);
		}
		printf("x");
		ptr[anzptr]=current=new DBYTE[allocatePerBlockIdx];
		anzptr++;
		if (!current) {
			LOGMSG("Error/2. Memory-DByteMgr.\n");
			exit(99);
		}
		freeFromIdx=0;
		allocatedIdx=allocatePerBlockIdx;
	}
	
	DBYTE *p=&current[freeFromIdx];
	freeFromIdx += aanz;
	
	return p;
}

DByteMgr::~DByteMgr() {
	for(int32_t i=0;i<anzptr;i++) {
		if (ptr[i]) delete[] ptr[i];
	}
	anzptr=0;
}

// struct Palette4

void Palette4::setPaletteRGB(const int aidx,const int ar,const int ag,const int ab) {
	if ( (aidx>=0) && (aidx < 256) ) {
		rgbs[aidx].R=ar;
		rgbs[aidx].G=ag;
		rgbs[aidx].B=ab;
		rgbs[aidx].alpha=0;
	}
}

int getfuncidx(const char* s) {
	for(int32_t i=0;i<FUNCANZ;i++) {
		if (!strcmp(s,funcname[i])) return i;
	}
	
	return -1;
}

int32_t bitsSufficient(
	const char* as,
	const int32_t aRANGE,
	const int32_t aREFINEMENT,
	const char* aNTS
) {
	char rl[128];
	char *test=new char[strlen(as)+16];
	strcpy(test,as);
	upper(test);
	sprintf(rl,";R%iL%i,",aRANGE,aREFINEMENT);
	upper(rl);
	char *p=strstr(test,rl);
	if (!p) {
		delete[] test;
		return 0;
	}
	// ;R2L8,A,D,LD,F1,QD,FP,;
	char *p2=strstr(p+1,",;");
	if (!p2) {
		delete[] test;
		return 0;
	}
	p2[1]=0;
	sprintf(rl,",%s,",aNTS);
	delete[] test;
	if (strstr(p,rl)) return 1;

	return 0;
}

int32_t testA(void) {
	#ifdef _FPA
	double c1=FAKTORAre.convert_to_double();
	double c2=FAKTORAim.convert_to_double();
	if (
		( fabs(c1) > 2.0 ) ||
		( fabs(c2) > 2.0 )
		) return 0;
	#else
	if (\
		( fabs( (double)FAKTORAre ) > 2.0 ) ||\
		( fabs( (double)FAKTORAim ) > 2.0 )\
	) return 0;
	#endif
	
	return 1;
}


void setfunc_and_bitprecision(const int afunc,char* afn) {
	char tmp2[1024],tmp3[1024];
	
	int8_t bitprecision=-1; // -1 = not tested
	
	// |Creal|<=2, |Cimag|<=2
	// und wenn nötig: |Areal|<=2, |Aimag|<=2
	#ifdef _FPA
	double c1=seedC0re.convert_to_double();
	double c2=seedC0im.convert_to_double();
	double c3=seedC1re.convert_to_double();
	double c4=seedC1im.convert_to_double();
	if (
		(fabs(c1) > 2.0) ||
		(fabs(c2) > 2.0) ||
		(fabs(c3) > 2.0) ||
		(fabs(c4) > 2.0)
	) bitprecision=0;
	#else
	if (
		( fabs( (double)seedC0re ) > 2.0 ) || 
		( fabs( (double)seedC0im ) > 2.0 ) || 
		( fabs( (double)seedC1re ) > 2.0 ) || 
		( fabs( (double)seedC1im ) > 2.0 ) 
	) bitprecision=0;
	#endif
	
	switch (afunc) {
		case FUNC_Z3AZC: {
			checkclockatbbxadd >>= 1;
			sprintf(afn,"_L%02i_%s_z3azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z3azc;
			
			bitprecision=1;
			if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,D,LD,F1,QD,FP,;R2L15,A,D,LD,F1,QD,FP,;R2L16,A,D,LD,F1,QD,FP,;R2L17,A,D,LD,F1,QD,FP,;R2L18,A,LD,F1,QD,FP,;R2L19,A,LD,F1,QD,FP,;R2L20,A,LD,F1,QD,FP,;R2L21,A,F1,QD,FP,;R2L22,A,F1,QD,FP,;R2L23,A,F1,QD,FP,;R2L24,A,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,D,LD,F1,QD,FP,;R4L15,A,D,LD,F1,QD,FP,;R4L16,A,D,LD,F1,QD,FP,;R4L17,A,D,LD,F1,QD,FP,;R4L18,A,LD,F1,QD,FP,;R4L19,A,LD,F1,QD,FP,;R4L20,A,LD,F1,QD,FP,;R4L21,A,F1,QD,FP,;R4L22,A,F1,QD,FP,;R4L23,A,F1,QD,FP,;R4L24,A,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,D,LD,F1,QD,FP,;R8L15,A,D,LD,F1,QD,FP,;R8L16,A,D,LD,F1,QD,FP,;R8L17,A,D,LD,F1,QD,FP,;R8L18,A,LD,F1,QD,FP,;R8L19,A,LD,F1,QD,FP,;R8L20,A,LD,F1,QD,FP,;R8L21,A,F1,QD,FP,;R8L22,A,F1,QD,FP,;R8L23,A,F1,QD,FP,;R8L24,A,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,D,LD,F1,QD,FP,;R2L15,A,D,LD,F1,QD,FP,;R2L16,A,D,LD,F1,QD,FP,;R2L17,A,D,LD,F1,QD,FP,;R2L18,A,LD,F1,QD,FP,;R2L19,A,LD,F1,QD,FP,;R2L20,A,LD,F1,QD,FP,;R2L21,A,F1,QD,FP,;R2L22,A,F1,QD,FP,;R2L23,A,F1,QD,FP,;R2L24,A,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,D,LD,F1,QD,FP,;R4L15,A,D,LD,F1,QD,FP,;R4L16,A,D,LD,F1,QD,FP,;R4L17,A,D,LD,F1,QD,FP,;R4L18,A,LD,F1,QD,FP,;R4L19,A,LD,F1,QD,FP,;R4L20,A,LD,F1,QD,FP,;R4L21,A,F1,QD,FP,;R4L22,A,F1,QD,FP,;R4L23,A,F1,QD,FP,;R4L24,A,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,D,LD,F1,QD,FP,;R8L15,A,D,LD,F1,QD,FP,;R8L16,A,D,LD,F1,QD,FP,;R8L17,A,D,LD,F1,QD,FP,;R8L18,A,LD,F1,QD,FP,;R8L19,A,LD,F1,QD,FP,;R8L20,A,LD,F1,QD,FP,;R8L21,A,F1,QD,FP,;R8L22,A,F1,QD,FP,;R8L23,A,F1,QD,FP,;R8L24,A,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			break;
		} 
		case FUNC_Z4AZC: {
			checkclockatbbxadd >>= 2;
			sprintf(afn,"_L%02i_%s_z4azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z4azc;
			bitprecision=1;
			
			if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,LD,F1,QD,FP,;R2L15,A,LD,F1,QD,FP,;R2L16,A,F1,QD,FP,;R2L17,A,F1,QD,FP,;R2L18,A,F1,QD,FP,;R2L19,A,F1,QD,FP,;R2L20,A,F1,QD,FP,;R2L21,A,F1,QD,FP,;R2L22,A,F1,QD,FP,;R2L23,A,F1,QD,FP,;R2L24,A,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,LD,F1,QD,FP,;R4L15,A,LD,F1,QD,FP,;R4L16,A,F1,QD,FP,;R4L17,A,F1,QD,FP,;R4L18,A,F1,QD,FP,;R4L19,A,F1,QD,FP,;R4L20,A,F1,QD,FP,;R4L21,A,F1,QD,FP,;R4L22,A,F1,QD,FP,;R4L23,A,F1,QD,FP,;R4L24,A,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,LD,F1,QD,FP,;R8L15,A,LD,F1,QD,FP,;R8L16,A,F1,QD,FP,;R8L17,A,F1,QD,FP,;R8L18,A,F1,QD,FP,;R8L19,A,F1,QD,FP,;R8L20,A,F1,QD,FP,;R8L21,A,F1,QD,FP,;R8L22,A,F1,QD,FP,;R8L23,A,F1,QD,FP,;R8L24,A,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,LD,F1,QD,FP,;R2L15,A,LD,F1,QD,FP,;R2L16,A,F1,QD,FP,;R2L17,A,F1,QD,FP,;R2L18,A,F1,QD,FP,;R2L19,A,F1,QD,FP,;R2L20,A,F1,QD,FP,;R2L21,A,F1,QD,FP,;R2L22,A,F1,QD,FP,;R2L23,A,F1,QD,FP,;R2L24,A,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,LD,F1,QD,FP,;R4L15,A,LD,F1,QD,FP,;R4L16,A,F1,QD,FP,;R4L17,A,F1,QD,FP,;R4L18,A,F1,QD,FP,;R4L19,A,F1,QD,FP,;R4L20,A,F1,QD,FP,;R4L21,A,F1,QD,FP,;R4L22,A,F1,QD,FP,;R4L23,A,F1,QD,FP,;R4L24,A,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,LD,F1,QD,FP,;R8L15,A,LD,F1,QD,FP,;R8L16,A,F1,QD,FP,;R8L17,A,F1,QD,FP,;R8L18,A,F1,QD,FP,;R8L19,A,F1,QD,FP,;R8L20,A,F1,QD,FP,;R8L21,A,F1,QD,FP,;R8L22,A,F1,QD,FP,;R8L23,A,F1,QD,FP,;R8L24,A,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
		case FUNC_Z5AZC: {
			checkclockatbbxadd >>= 4;
			sprintf(afn,"_L%02i_%s_z5azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z5azc;
			bitprecision=1;
			
			if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,LD,F1,QD,FP,;R2L12,A,LD,F1,QD,FP,;R2L13,A,F1,QD,FP,;R2L14,A,F1,QD,FP,;R2L15,A,F1,QD,FP,;R2L16,A,F1,QD,FP,;R2L17,A,F1,QD,FP,;R2L18,A,F1,QD,FP,;R2L19,A,F1,QD,FP,;R2L20,A,F1,QD,FP,;R2L21,A,F1,QD,FP,;R2L22,A,QD,;R2L23,A,QD,;R2L24,A,QD,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,LD,F1,QD,FP,;R4L12,A,LD,F1,QD,FP,;R4L13,A,F1,QD,FP,;R4L14,A,F1,QD,FP,;R4L15,A,F1,QD,FP,;R4L16,A,F1,QD,FP,;R4L17,A,F1,QD,FP,;R4L18,A,F1,QD,FP,;R4L19,A,F1,QD,FP,;R4L20,A,F1,QD,FP,;R4L21,A,F1,QD,FP,;R4L22,A,QD,FP,;R4L23,A,QD,;R4L24,A,QD,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,LD,F1,QD,FP,;R8L12,A,LD,F1,QD,FP,;R8L13,A,F1,QD,FP,;R8L14,A,F1,QD,FP,;R8L15,A,F1,QD,FP,;R8L16,A,F1,QD,FP,;R8L17,A,F1,QD,FP,;R8L18,A,F1,QD,FP,;R8L19,A,F1,QD,FP,;R8L20,A,F1,QD,FP,;R8L21,A,F1,QD,FP,;R8L22,A,QD,FP,;R8L23,A,QD,FP,;R8L24,A,QD,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,LD,F1,QD,FP,;R2L12,A,LD,F1,QD,FP,;R2L13,A,F1,QD,FP,;R2L14,A,F1,QD,FP,;R2L15,A,F1,QD,FP,;R2L16,A,F1,QD,FP,;R2L17,A,F1,QD,FP,;R2L18,A,F1,QD,FP,;R2L19,A,F1,QD,FP,;R2L20,A,F1,QD,FP,;R2L21,A,F1,QD,FP,;R2L22,A,QD,;R2L23,A,QD,;R2L24,A,QD,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,LD,F1,QD,FP,;R4L12,A,LD,F1,QD,FP,;R4L13,A,F1,QD,FP,;R4L14,A,F1,QD,FP,;R4L15,A,F1,QD,FP,;R4L16,A,F1,QD,FP,;R4L17,A,F1,QD,FP,;R4L18,A,F1,QD,FP,;R4L19,A,F1,QD,FP,;R4L20,A,F1,QD,FP,;R4L21,A,F1,QD,FP,;R4L22,A,QD,FP,;R4L23,A,QD,;R4L24,A,QD,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,LD,F1,QD,FP,;R8L12,A,LD,F1,QD,FP,;R8L13,A,F1,QD,FP,;R8L14,A,F1,QD,FP,;R8L15,A,F1,QD,FP,;R8L16,A,F1,QD,FP,;R8L17,A,F1,QD,FP,;R8L18,A,F1,QD,FP,;R8L19,A,F1,QD,FP,;R8L20,A,F1,QD,FP,;R8L21,A,F1,QD,FP,;R8L22,A,QD,FP,;R8L23,A,QD,FP,;R8L24,A,QD,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
		case FUNC_Z6AZC: {
			checkclockatbbxadd >>= 5;
			sprintf(afn,"_L%02i_%s_z6azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z6azc;
			
			bitprecision=1;
			if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,LD,F1,QD,FP,;R2L10,A,LD,F1,QD,FP,;R2L11,A,F1,QD,FP,;R2L12,A,F1,QD,FP,;R2L13,A,F1,QD,FP,;R2L14,A,F1,QD,FP,;R2L15,A,F1,QD,FP,;R2L16,A,F1,QD,FP,;R2L17,A,F1,QD,FP,;R2L18,A,QD,;R2L19,A,QD,;R2L20,A,QD,;R2L21,A,QD,;R2L22,A,;R2L23,A,;R2L24,A,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,LD,F1,QD,FP,;R4L10,A,LD,F1,QD,FP,;R4L11,A,F1,QD,FP,;R4L12,A,F1,QD,FP,;R4L13,A,F1,QD,FP,;R4L14,A,F1,QD,FP,;R4L15,A,F1,QD,FP,;R4L16,A,F1,QD,FP,;R4L17,A,F1,QD,FP,;R4L18,A,QD,FP,;R4L19,A,QD,;R4L20,A,QD,;R4L21,A,QD,;R4L22,A,;R4L23,A,;R4L24,A,;;R8L8,A,LD,F1,QD,FP,;R8L9,A,LD,F1,QD,FP,;R8L10,A,LD,F1,QD,FP,;R8L11,A,F1,QD,FP,;R8L12,A,F1,QD,FP,;R8L13,A,F1,QD,FP,;R8L14,A,F1,QD,FP,;R8L15,A,F1,QD,FP,;R8L16,A,F1,QD,FP,;R8L17,A,F1,QD,FP,;R8L18,A,QD,FP,;R8L19,A,QD,FP,;R8L20,A,QD,;R8L21,A,QD,;R8L22,A,;R8L23,A,;R8L24,A,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,LD,F1,QD,FP,;R2L10,A,LD,F1,QD,FP,;R2L11,A,F1,QD,FP,;R2L12,A,F1,QD,FP,;R2L13,A,F1,QD,FP,;R2L14,A,F1,QD,FP,;R2L15,A,F1,QD,FP,;R2L16,A,F1,QD,FP,;R2L17,A,F1,QD,FP,;R2L18,A,QD,;R2L19,A,QD,;R2L20,A,QD,;R2L21,A,QD,;R2L22,A,;R2L23,A,;R2L24,A,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,LD,F1,QD,FP,;R4L10,A,LD,F1,QD,FP,;R4L11,A,F1,QD,FP,;R4L12,A,F1,QD,FP,;R4L13,A,F1,QD,FP,;R4L14,A,F1,QD,FP,;R4L15,A,F1,QD,FP,;R4L16,A,F1,QD,FP,;R4L17,A,F1,QD,FP,;R4L18,A,QD,FP,;R4L19,A,QD,;R4L20,A,QD,;R4L21,A,QD,;R4L22,A,;R4L23,A,;R4L24,A,;;R8L8,A,LD,F1,QD,FP,;R8L9,A,LD,F1,QD,FP,;R8L10,A,LD,F1,QD,FP,;R8L11,A,F1,QD,FP,;R8L12,A,F1,QD,FP,;R8L13,A,F1,QD,FP,;R8L14,A,F1,QD,FP,;R8L15,A,F1,QD,FP,;R8L16,A,F1,QD,FP,;R8L17,A,F1,QD,FP,;R8L18,A,QD,FP,;R8L19,A,QD,FP,;R8L20,A,QD,;R8L21,A,QD,;R8L22,A,;R8L23,A,;R8L24,A,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
		case FUNC_Z2AZC: {
			sprintf(afn,"_L%02i_%s_z2azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z2azc;
			
			bitprecision=1;
			if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,D,LD,F1,QD,FP,;R2L15,A,D,LD,F1,QD,FP,;R2L16,A,D,LD,F1,QD,FP,;R2L17,A,D,LD,F1,QD,FP,;R2L18,A,D,LD,F1,QD,FP,;R2L19,A,D,LD,F1,QD,FP,;R2L20,A,D,LD,F1,QD,FP,;R2L21,A,D,LD,F1,QD,FP,;R2L22,A,D,LD,F1,QD,FP,;R2L23,A,D,LD,F1,QD,FP,;R2L24,A,D,LD,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,D,LD,F1,QD,FP,;R4L15,A,D,LD,F1,QD,FP,;R4L16,A,D,LD,F1,QD,FP,;R4L17,A,D,LD,F1,QD,FP,;R4L18,A,D,LD,F1,QD,FP,;R4L19,A,D,LD,F1,QD,FP,;R4L20,A,D,LD,F1,QD,FP,;R4L21,A,D,LD,F1,QD,FP,;R4L22,A,D,LD,F1,QD,FP,;R4L23,A,D,LD,F1,QD,FP,;R4L24,A,D,LD,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,D,LD,F1,QD,FP,;R8L15,A,D,LD,F1,QD,FP,;R8L16,A,D,LD,F1,QD,FP,;R8L17,A,D,LD,F1,QD,FP,;R8L18,A,D,LD,F1,QD,FP,;R8L19,A,D,LD,F1,QD,FP,;R8L20,A,D,LD,F1,QD,FP,;R8L21,A,D,LD,F1,QD,FP,;R8L22,A,D,LD,F1,QD,FP,;R8L23,A,D,LD,F1,QD,FP,;R8L24,A,LD,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,D,LD,F1,QD,FP,;R2L15,A,D,LD,F1,QD,FP,;R2L16,A,D,LD,F1,QD,FP,;R2L17,A,D,LD,F1,QD,FP,;R2L18,A,D,LD,F1,QD,FP,;R2L19,A,D,LD,F1,QD,FP,;R2L20,A,D,LD,F1,QD,FP,;R2L21,A,D,LD,F1,QD,FP,;R2L22,A,D,LD,F1,QD,FP,;R2L23,A,D,LD,F1,QD,FP,;R2L24,A,D,LD,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,D,LD,F1,QD,FP,;R4L15,A,D,LD,F1,QD,FP,;R4L16,A,D,LD,F1,QD,FP,;R4L17,A,D,LD,F1,QD,FP,;R4L18,A,D,LD,F1,QD,FP,;R4L19,A,D,LD,F1,QD,FP,;R4L20,A,D,LD,F1,QD,FP,;R4L21,A,D,LD,F1,QD,FP,;R4L22,A,D,LD,F1,QD,FP,;R4L23,A,D,LD,F1,QD,FP,;R4L24,A,D,LD,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,D,LD,F1,QD,FP,;R8L15,A,D,LD,F1,QD,FP,;R8L16,A,D,LD,F1,QD,FP,;R8L17,A,D,LD,F1,QD,FP,;R8L18,A,D,LD,F1,QD,FP,;R8L19,A,D,LD,F1,QD,FP,;R8L20,A,D,LD,F1,QD,FP,;R8L21,A,D,LD,F1,QD,FP,;R8L22,A,D,LD,F1,QD,FP,;R8L23,A,D,LD,F1,QD,FP,;R8L24,A,LD,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			break;
		}
		default: {
			getBoundingBoxfA=getBoundingBoxfA_z2c;
			sprintf(afn,"_L%02i%s_z2c_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2));

			bitprecision=1;
			
			// real part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,D,LD,F1,QD,FP,;R2L15,A,D,LD,F1,QD,FP,;R2L16,A,D,LD,F1,QD,FP,;R2L17,A,D,LD,F1,QD,FP,;R2L18,A,D,LD,F1,QD,FP,;R2L19,A,D,LD,F1,QD,FP,;R2L20,A,D,LD,F1,QD,FP,;R2L21,A,D,LD,F1,QD,FP,;R2L22,A,D,LD,F1,QD,FP,;R2L23,A,D,LD,F1,QD,FP,;R2L24,A,D,LD,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,D,LD,F1,QD,FP,;R4L15,A,D,LD,F1,QD,FP,;R4L16,A,D,LD,F1,QD,FP,;R4L17,A,D,LD,F1,QD,FP,;R4L18,A,D,LD,F1,QD,FP,;R4L19,A,D,LD,F1,QD,FP,;R4L20,A,D,LD,F1,QD,FP,;R4L21,A,D,LD,F1,QD,FP,;R4L22,A,D,LD,F1,QD,FP,;R4L23,A,D,LD,F1,QD,FP,;R4L24,A,D,LD,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,D,LD,F1,QD,FP,;R8L15,A,D,LD,F1,QD,FP,;R8L16,A,D,LD,F1,QD,FP,;R8L17,A,D,LD,F1,QD,FP,;R8L18,A,D,LD,F1,QD,FP,;R8L19,A,D,LD,F1,QD,FP,;R8L20,A,D,LD,F1,QD,FP,;R8L21,A,D,LD,F1,QD,FP,;R8L22,A,D,LD,F1,QD,FP,;R8L23,A,D,LD,F1,QD,FP,;R8L24,A,D,LD,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,FP,;R2L9,A,D,LD,F1,QD,FP,;R2L10,A,D,LD,F1,QD,FP,;R2L11,A,D,LD,F1,QD,FP,;R2L12,A,D,LD,F1,QD,FP,;R2L13,A,D,LD,F1,QD,FP,;R2L14,A,D,LD,F1,QD,FP,;R2L15,A,D,LD,F1,QD,FP,;R2L16,A,D,LD,F1,QD,FP,;R2L17,A,D,LD,F1,QD,FP,;R2L18,A,D,LD,F1,QD,FP,;R2L19,A,D,LD,F1,QD,FP,;R2L20,A,D,LD,F1,QD,FP,;R2L21,A,D,LD,F1,QD,FP,;R2L22,A,D,LD,F1,QD,FP,;R2L23,A,D,LD,F1,QD,FP,;R2L24,A,D,LD,F1,QD,FP,;;R4L8,A,D,LD,F1,QD,FP,;R4L9,A,D,LD,F1,QD,FP,;R4L10,A,D,LD,F1,QD,FP,;R4L11,A,D,LD,F1,QD,FP,;R4L12,A,D,LD,F1,QD,FP,;R4L13,A,D,LD,F1,QD,FP,;R4L14,A,D,LD,F1,QD,FP,;R4L15,A,D,LD,F1,QD,FP,;R4L16,A,D,LD,F1,QD,FP,;R4L17,A,D,LD,F1,QD,FP,;R4L18,A,D,LD,F1,QD,FP,;R4L19,A,D,LD,F1,QD,FP,;R4L20,A,D,LD,F1,QD,FP,;R4L21,A,D,LD,F1,QD,FP,;R4L22,A,D,LD,F1,QD,FP,;R4L23,A,D,LD,F1,QD,FP,;R4L24,A,D,LD,F1,QD,FP,;;R8L8,A,D,LD,F1,QD,FP,;R8L9,A,D,LD,F1,QD,FP,;R8L10,A,D,LD,F1,QD,FP,;R8L11,A,D,LD,F1,QD,FP,;R8L12,A,D,LD,F1,QD,FP,;R8L13,A,D,LD,F1,QD,FP,;R8L14,A,D,LD,F1,QD,FP,;R8L15,A,D,LD,F1,QD,FP,;R8L16,A,D,LD,F1,QD,FP,;R8L17,A,D,LD,F1,QD,FP,;R8L18,A,D,LD,F1,QD,FP,;R8L19,A,D,LD,F1,QD,FP,;R8L20,A,D,LD,F1,QD,FP,;R8L21,A,D,LD,F1,QD,FP,;R8L22,A,D,LD,F1,QD,FP,;R8L23,A,D,LD,F1,QD,FP,;R8L24,A,D,LD,F1,QD,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			break;
		}
	} // switch
	
	if (bitprecision<0) {
		LOGMSG("Not able to check bit precision requirements.\n");
	} else if (bitprecision==0) {
		LOGMSG("Datatype does NOT provide sufficient precision.\n");
	} else {
		LOGMSG("Datatype provides sufficient precision.\n");
	}
}

// struct ListeFIFO
ListeFIFO::ListeFIFO() {
	next_readpos=next_writepos=0;
	werte=NULL;
	allokiereAmStueck=((CHUNKSIZE>>1) / sizeof(Int2));
}

ListeFIFO::~ListeFIFO() {
	if (werte) delete[] werte;
}

void ListeFIFO::start(void) {
	next_readpos=next_writepos=0;
}

char ListeFIFO::read(int32_t& x,int32_t& y) {
	if (!werte) return 0;
	if (next_readpos >= next_writepos) return 0;
	x=werte[next_readpos].x;
	y=werte[next_readpos].y;
	next_readpos++;
	
	return 1;
}

char ListeFIFO::write(const int32_t x,const int32_t y) {
	if (!werte) {
		// allokieren
		werte=new Int2[allokiereAmStueck];
		if (!werte) {
			LOGMSG("ListeFIFO: Kein Speicher\n");
			exit(99);
		}
		next_writepos=0;
	}
	
	if (next_writepos >= allokiereAmStueck) {
		// kein Platz mehr
		return 0;
	}
	
	werte[next_writepos].x=x;
	werte[next_writepos].y=y;
	next_writepos++;
	
	return 1;
}

// struct ListeDFS
ListeDFS::ListeDFS() {
	anz=0;
	werte=NULL;
	allokiereAmStueck=(CHUNKSIZE / sizeof(DFSPunkt));
}

ListeDFS::~ListeDFS() {
	if (werte) delete[] werte;
}

void ListeDFS::start(void) {
	anz=0;
}

char ListeDFS::read(int32_t& x,int32_t& y,DBYTE& t) {
	if (!werte) return 0;
	if (anz <= 0) return 0;
	x=werte[anz-1].x;
	y=werte[anz-1].y;
	t=werte[anz-1].tiefe;
	anz--;
	
	return 1;
}

char ListeDFS::write(const int32_t x,const int32_t y,const DBYTE t) {
	if (!werte) {
		// allokieren
		werte=new DFSPunkt[allokiereAmStueck];
		if (!werte) {
			LOGMSG("ListeFIFO: Kein Speicher\n");
			exit(99);
		}
		anz=0;
	}
	
	if (anz >= allokiereAmStueck) {
		return 0;
	}
	
	werte[anz].x=x;
	werte[anz].y=y;
	werte[anz].tiefe=t;
	anz++;
	
	return 1;
}

// struct ArrayDDByteManager

ArrayDDByteManager::ArrayDDByteManager() {
	current=NULL;
	allokierteIdx=0;
	freiAbIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(DDBYTE);
	allokierePerBlockIdx=(int32_t)floor(d);
}

ArrayDDByteManager::~ArrayDDByteManager() {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
}

PDDBYTE ArrayDDByteManager::getMemory(const int aanz) {
	// lieber etwas Puffer, daher +2
	if (anzptr >= (MAXPTR-8)) {
		LOGMSG("ArrayDDByteManager:: Zu wenig Speicher\n");
		exit(99);
	}
	if (
		(!current) ||
		((freiAbIdx + aanz + 2) >= allokierteIdx)
	) {
		printf("x");
		ptr[anzptr]=current=new DDBYTE[allokierePerBlockIdx];
		anzptr++;
		if (!current) {
			printf("Memory-Fehler. ArrayDDByteManager.\n");
			exit(99);
		}
		freiAbIdx=0;
		allokierteIdx=allokierePerBlockIdx;
	}
	
	// nun ist current allokiert udn hat genügend Platz
	PDDBYTE p=&current[freiAbIdx];
	freiAbIdx += aanz;
	return p;
}

void convert_raw_into_newstructure(void) {
	// convert old raw-file structure into new
	// dynamical one
	// input files: _in.raw_header, _in.raw_0001 etc.
	// output file: _2d.raw
	
	LOGMSG("converting old file structure ... ");
	
	char fn[1024];
	sprintf(fn,"_in.raw_header");
	FILE *f=fopen(fn,"rb");
	int32_t scrb;
	fread(&scrb,sizeof(scrb),1,f);
	fclose(f);
	
	int32_t memwidth=scrb >> 4;
	DDBYTE *one=new DDBYTE[memwidth];
	int32_t fctr=0;
	// _in.raw_0001
	FILE *fout=fopen("_2d.raw","wb");
	fwrite(&scrb,1,sizeof(DDBYTE),fout);
	
	while (1) {
		fctr++;
		sprintf(fn,"_in.raw_%04i",fctr);
		f=fopen(fn,"rb");
		if (!f) break;
		printf("%s\n",fn);
		int32_t spla;
		fread(&spla,sizeof(spla),1,f);
		for(int32_t yread=0;yread<spla;yread++) {
			fread(one,sizeof(DDBYTE),memwidth,f);
			int32_t m0=memwidth,m1=0;
			for(int32_t x=0;x<memwidth;x++) {
				if (one[x] != SQUARE_WHITE_16_CONSECUTIVE) {
					if (x < m0) m0=x;
					if (x > m1) m1=x;
				}
			}

			int32_t laenge=m1-m0+1;
			if (laenge > 0) {
				fwrite(&m0,1,sizeof(m0),fout);
				fwrite(&laenge,1,sizeof(laenge),fout);
				fwrite(
					&one[m0],
					laenge,
					sizeof(DDBYTE),
					fout);
			} else {
				int32_t start=0;
				laenge=0;
				fwrite(&start,1,sizeof(start),fout);
				fwrite(&laenge,1,sizeof(laenge),fout);
			}
		}
		fclose(f);
	} // file
	
	fclose(fout);
	delete[] one;
	
	LOGMSG("done\n");
}

#ifdef _FPA
// struct FPA

inline FPA& FPA::operator=(const FPA avalue) {
	if (this != &avalue) {
		a=avalue.a;
		b=avalue.b;
		c=avalue.c;
		d=avalue.d;
		vorz=avalue.vorz;
	}
	
	return *this;
}

inline bool operator!=(const FPA a,const FPA b) {
	int vgl=FPA_vgl(a,b);
	if (vgl != 0) return 0;
	return 1;
}

inline bool operator>(const FPA a,const FPA b) {
	int vgl=FPA_vgl(a,b);
	if (vgl > 0) return 1;
	return 0;
}

inline bool operator<(const FPA a,const FPA b) {
	int vgl=FPA_vgl(a,b);
	if (vgl < 0) return 1;
	return 0;
}

void minimaxFPAMIMAAB(FPA& mi,FPA& ma,const FPA term1,const FPA term2) {
	int vgl=FPA_vgl(term1,term2);
	if (vgl>=0) { mi=term2; ma=term1; }
	else { mi=term1; ma=term2; }
}

void minimaxFPAMIMAAB(PFPA& mi,PFPA& ma,PFPA term1,PFPA term2) {
	int vgl=FPA_vgl(term1,term2);
	if (vgl>=0) { mi=term2; ma=term1; }
	else { mi=term1; ma=term2; }
}

void minimaxFPAMIMAABCD(FPA& mi,FPA& ma,const FPA term1,const FPA term2,const FPA term3,const FPA term4) {
	FPA miab,maab,micd,macd;

	minimaxFPAMIMAAB(miab,maab,term1,term2);	
	minimaxFPAMIMAAB(micd,macd,term3,term4);	
	if (FPA_vgl(miab,micd) <= 0) mi=miab; else mi=micd;
	if (FPA_vgl(maab,macd) >= 0) ma=maab; else ma=macd;
}

void minimaxFPAMIMAABCD(PFPA& mi,PFPA& ma,PFPA term1,PFPA term2,PFPA term3,PFPA term4) {
	PFPA miab,maab,micd,macd;

	minimaxFPAMIMAAB(miab,maab,term1,term2);	
	minimaxFPAMIMAAB(micd,macd,term3,term4);	
	if (FPA_vgl(miab,micd) <= 0) mi=miab; else mi=micd;
	if (FPA_vgl(maab,macd) >= 0) ma=maab; else ma=macd;
}

inline FPA operator+(const FPA a,const FPA b) {
    FPA ret;
    FPA_add_ZAB(ret,a,b);
    
    return ret;
}

inline FPA operator*(const FPA a,const FPA b) {
    FPA ret;
    FPA_mul_ZAB(ret,a,b);
    
    return ret;
}

inline FPA operator*(const int32_t a,const FPA b) {
    FPA ret;
    if (a<0) {
		FPA_mul_ZAuvlong(ret,b,-a);
		ret.vorz *= -1;
    } else if (a==0) {
		ret.setNull();
    } else {
		FPA_mul_ZAuvlong(ret,b,a);
	}
    
    return ret;
}

inline FPA operator-(const FPA a,const FPA b) {
    FPA ret;
    FPA_sub_ZAB(ret,a,b);
    
    return ret;
}

inline int64_t floorFPA(const FPA term) {
	if (term.vorz==0) return 0;
	
	if (term.vorz>0) {
		return (VLONG)term.a;
	}
	
	// negative
	// add one
	if (
		(term.b != 0) ||
		(term.c != 0) ||
		(term.d != 0)
	) return (-(VLONG)term.a-(VLONG)1);
	
	return -(VLONG)term.a;
}

inline void FPA_mul_ZAuvlong(FPA& erg,const FPA term,const UVLONG intmul) {
	if ( (intmul==0) || (term.vorz==0) ) {
		erg.setNull();
		return;
	}
	
	UVLONG w=(UVLONG)term.d * (UVLONG)intmul;
	erg.d=(w & MAXDDBYTE);
	w=(UVLONG)term.c * (UVLONG)intmul + (w >> 32);
	erg.c=(w & MAXDDBYTE);
	w=(UVLONG)term.b * (UVLONG)intmul + (w >> 32);
	erg.b=(w & MAXDDBYTE);
	w=(UVLONG)term.a * (UVLONG)intmul + (w >> 32);
	
	if (w > MAXDDBYTE) {
		LOGMSG("mul_zap out of range\n");
		exit(99);
	}
	erg.a=(w & MAXDDBYTE);
	erg.vorz=term.vorz;
}

inline void FPA_mul_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if ( (term1.vorz==0) || (term2.vorz==0) ) {
		erg.setNull();
		return;
	}
	
	// beide ungleich null
	
	// R=2^-32
	// Zahl 1: a + b*R + c*R^2 + d*R^3
	// Zahl 2: e + f*R + g*R^2 + h*R^3
	/*
		a*e
		a*f*R + b*e*R
		a*g*R^2 + b*f*R^2 + c*e*R^2
		a*h*R^3 + b*g*R^3 + c*f*R^3 + d*e*R^3
		
		higher order terms: must result in 0, but can
		create a carry-over into an R^3 term
				
		b*h*R^4 + c*g*R^4 + d*f*R^4
		c*h*R^5 + d*g*R^5
		d*h*R^6
	*/
	
	FPA w;
	UVLONG tmp;
	erg.setNull();
	w.a=w.b=w.c=w.d=0;
	w.vorz=1; // not 0, because later digits are set into

	// R^6
	FPA testerg;
	testerg.a=testerg.b=testerg.c=testerg.d=0;
	testerg.vorz=1;
	
	if ((term1.d!=0)&&(term2.d!=0)) {
		//d*h*R^6 => R^3
		tmp=(UVLONG)term1.d*(UVLONG)term2.d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	w.d=0;
	if ((term1.c!=0)&&(term2.d!=0)) {
		//c*h*R^5 + 
		tmp=(UVLONG)term1.c*(UVLONG)term2.d;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	if ((term1.d!=0)&&(term2.c!=0)) {
		//d*g*R^5
		tmp=(UVLONG)term1.d*(UVLONG)term2.c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	w.c=0;
	if ((term1.b!=0)&&(term2.d!=0)) {
		// b*h*R^4 
		tmp=(UVLONG)term1.b*(UVLONG)term2.d;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if ((term1.c!=0)&&(term2.c!=0)) {
		// + c*g*R^4 
		tmp=(UVLONG)term1.c*(UVLONG)term2.c;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if ((term1.d!=0)&&(term2.b!=0)) {
		// + d*f*R^4
		tmp=(UVLONG)term1.d*(UVLONG)term2.b;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	#define FEHLER(TT) \
	{\
		LOGMSG("Implementation error. Precision FPA not sufficient in mul\n");\
		LOGMSG2("%s",TT);\
		exit(99);\
	}

	// w.a can != 0 
	// w.b==w.c==w.d==0
	
	if (
		(testerg.b != 0) ||
		(testerg.c != 0) ||
		(testerg.d != 0)
	) {
		FEHLER("Out-of-range R4-R6 mul");
	}
	
	erg.d=testerg.a; // possible carry over
	
	w.a=w.b=w.c=w.d=0;
	// R^3
	// d*e*R^3
	if ( (term1.d!=0) && (term2.a!=0) ) {
		tmp=(UVLONG)term1.d*(UVLONG)term2.a;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	
	// a*h*R^3 + 
	if ( (term1.a!=0) && (term2.d!=0) ) {
		tmp=(UVLONG)term1.a*(UVLONG)term2.d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// b*g*R^3 + 
	if ( (term1.b!=0) && (term2.c!=0) ) {
		tmp=(UVLONG)term1.b*(UVLONG)term2.c;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// c*f*R^3 + 
	if ( (term1.c!=0) && (term2.b!=0) ) {
		tmp=(UVLONG)term1.c*(UVLONG)term2.b;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R^2
	w.d=0;
	// a*g*R^2 + 
	if ( (term1.a!=0) && (term2.c!=0) ) {
		tmp=(UVLONG)term1.a*(UVLONG)term2.c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// b*f*R^2 + 
	if ( (term1.b!=0) && (term2.b!=0) ) {
		tmp=(UVLONG)term1.b*(UVLONG)term2.b;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	// c*e*R^2
	if ( (term1.c!=0) && (term2.a!=0) ) {
		tmp=(UVLONG)term1.c*(UVLONG)term2.a;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R
	w.c=0;
	// a*f*R 
	if ( (term1.a!=0) && (term2.b!=0) ) {
		tmp=(UVLONG)term1.a*(UVLONG)term2.b;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	// + b*e*R
	if ( (term1.b!=0) && (term2.a!=0) ) {
		tmp=(UVLONG)term1.b*(UVLONG)term2.a;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
		
	// integer-Part
	w.b=0;
	// a*e
	if ( (term1.a!=0) && (term2.a!=0) ) {
		tmp=(UVLONG)term1.a*(UVLONG)term2.a;
		w.a=tmp & MAXDDBYTE;
		if ( (tmp >> 32) != 0) {
			LOGMSG("Overflow FPA mul.\n");
			exit(99);
		}
		FPA_add_ZAB(erg,erg,w);
	}
				
	if (
		(erg.b==0) &&
		(erg.a==0) &&
		(erg.c==0) &&
		(erg.d==0)
	) erg.vorz=0;
	else if (term1.vorz==term2.vorz) erg.vorz=1;
	else erg.vorz=-1;
}

inline void FPA::squareTo(FPA& erg) {
	// works faster than simple multipliucation
	if (vorz==0) {
		erg.setNull();
		return;
	}
	
	// R=2^-32
	// Zahl: a + b*R + c*R^2 + d*R^3
	/*
		a^2
		2abR
		2acR^2 + b^2*R^2
		2adR^3 + 2bc*R^3

		2bdR^4 + c^2*R^4 
		2cdR^5
		d^2R^6
		
	*/
	
	#define FEHLER(TT) \
	{\
		LOGMSG("Implementation error. Precision FPA not sufficient in mul\n");\
		LOGMSG2("%s",TT);\
		exit(99);\
	}
	
	FPA w;
	UVLONG tmp;
	erg.setNull();
	w.a=w.b=w.c=w.d=0;
	w.vorz=1; // not 0, as values are set into w later

	// R^6
	FPA testerg;
	testerg.a=testerg.b=testerg.c=testerg.d=0;
	testerg.vorz=1;
	
	if (d!=0) {
		// d^2R^6 => R^3 locally here in testerg
		tmp=(UVLONG)d*(UVLONG)d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	w.d=0;
	if ((c!=0)&&(d!=0)) {
		//2cdR^5: c < 2^32, d <2^32
		tmp=(UVLONG)c*(UVLONG)d;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		w.shiftLeft(1);
		FPA_add_ZAB(testerg,testerg,w);
		// maybe a carry-over was shifted into w.a
		// so has to be reset to 0
		w.a=0; 
	}
	
	w.c=0;
	if ((b!=0)&&(d!=0)) {
		// 2bdR^4
		tmp=(UVLONG)b*(UVLONG)d;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		w.shiftLeft(1);
		FPA_add_ZAB(testerg,testerg,w);
		// only a and b can be set
	}

	if (c!=0) {
		// + c^2*R^4 
		tmp=(UVLONG)c*(UVLONG)c;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if (
		(testerg.b != 0) ||
		(testerg.c != 0) ||
		(testerg.d != 0)
	) {
		FEHLER("out-of-range fpa sqzare -96");
	}
	
	erg.d=testerg.a; // Übertrag
	w.a=w.b=w.c=w.d=0;

	// R^3
	if ((d!=0)&&(a!=0)) {
		// 2adR^3
		tmp=2*(UVLONG)d*(UVLONG)a;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	
	if ( (b!=0) && (c!=0) ) {
		// + 2bc*R^3
		tmp=2*(UVLONG)b*(UVLONG)c;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R^2
	w.d=0;
	if ( (a!=0) && (c!=0) ) {
		//2acR^2 + 
		tmp=2*(UVLONG)a*(UVLONG)c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	if (b!=0) {
		// b^2*R^2
		tmp=(UVLONG)b*(UVLONG)b;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R
	w.c=0;
	if ( (a!=0) && (b!=0) ) {
		// 2abR
		tmp=2*(UVLONG)a*(UVLONG)b;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
		
	// integer-Part
	w.b=0;
	if (a!=0) {
		//a^2
		tmp=(UVLONG)a*(UVLONG)a;
		w.a=tmp & MAXDDBYTE;
		if ( (tmp >> 32) != 0) {
			FEHLER("a overflow");
		}
		FPA_add_ZAB(erg,erg,w);
	}
				
	if (
		(erg.b==0) &&
		(erg.a==0) &&
		(erg.c==0) &&
		(erg.d==0)
	) erg.vorz=0;
	else erg.vorz=1; // immer positiv
}

FPA::FPA(const VLONG a) {
	set_vlong(a);
}

FPA::FPA(const double a) {
	set_double32(a);
}

void FPA::shiftLeft(const int32_t bits) {
	// especially for multiplication in scrcoord_as_lowerleft
	// as scalePixelPerRange is a power of 2
	// sign is unaffected
	UVLONG tmpd=( (UVLONG)d << bits);
	d=tmpd & MAXDDBYTE;
	UVLONG tmpc=(((UVLONG)c << bits) | (tmpd >> 32));
	c=tmpc & MAXDDBYTE;
	UVLONG tmpb=(((UVLONG)b << bits) | (tmpc >> 32));
	b=tmpb & MAXDDBYTE;
	UVLONG tmpa=(((UVLONG)a << bits) | (tmpb >> 32));
	a=tmpa & MAXDDBYTE;

	if ((tmpa >> 32) != 0) {
		// overflow
		LOGMSG("Implementation error. Shift-FPA overflow.\n");
		exit(99);
	}

}

double FPA::convert_to_double(void) {
	if (vorz==0) return 0.0;
	
	double erg=a;
	if (b!=0) erg+=b*exp(log(2.0)*-32); 
	if (c!=0) erg+=c*exp(log(2.0)*-64);
	if (d!=0) erg+=d*exp(log(2.0)*-96);
	
	if (vorz<0) return (-erg);
	
	return erg;
}

inline void FPA::copyFrom(const FPA w) {
	vorz=w.vorz;
	a=w.a;
	b=w.b;
	c=w.c;
	d=w.d;
}

inline void FPA::checkNull(void) {
	if (
		(b == 0) &&
		(c == 0) &&
		(a == 0) &&
		(d == 0)
	) vorz=0;
}

inline void FPA::setNull(void) {
	vorz=0;
	a=b=c=d=0;
}

FPA::FPA() {
	vorz=0;
}

char* FPA::str(char* erg) {
	if (vorz>0) strcpy(erg,"+");
	else if (vorz<0) strcpy(erg,"-(");
	else {
		strcpy(erg,"0");
		return erg;
	}
	
	sprintf(&erg[strlen(erg)],"%u + %u*2^-32 + %u*2^-64 + %u*2^-96",
		a,b,c,d);
	
	return erg;
}

void FPA::set_double32(const double avalue) {
	// integerpart
	double w;
	if (avalue < 0) {
		vorz=-1;
		w=-avalue;
	} else {
		vorz=1;
		w=avalue;
	}
	
	VLONG vl=(VLONG)floor(w);
	if (
		(vl > (VLONG)MAXDDBYTE) ||
		(vl < (-(VLONG)MAXDDBYTE) )
	) {
		printf("setdouble32 out of range. FAP %I64d\n",vl);
		exit(99);
	}
	
	a=(vl & MAXDDBYTE);
	w -= vl;
	w *= ZWEIHOCH32;
	
	vl=(VLONG)floor(w);
	b=(vl & MAXDDBYTE);
	w -= vl;
	w *= ZWEIHOCH32;

	vl=(VLONG)floor(w);
	c=(vl & MAXDDBYTE);

	d=0; // sonst sind's zuviele Nachkommastellen
	
	// d ist ja 0
	if (
		(b==0) &&
		(a==0) &&
		(c==0)
	) vorz=0;
}

void FPA::set_vlong(const VLONG avalue) {
	UVLONG wert;
	if (avalue == 0) {
		vorz=0;
	} else {
		if (avalue < 0) {
			vorz=-1;
			wert=-avalue;
		} else {
			vorz=+1;
			wert=avalue;
		}
		
		if (wert > MAXDDBYTE) {
			printf("out of range. setvlong\n");
			exit(99);
		}
		
		a=wert & MAXDDBYTE;
		b=c=d=0;
	}
}

inline void FPA_add_abs_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if (term1.vorz==0) {
		if (term2.vorz==0) erg.setNull();
		else erg.copyFrom(term2);
	}
	else if (term2.vorz==0) erg.copyFrom(term1);
	else {
		erg.vorz=1; // null ist nicht möglich, da kein
		UVLONG sumdd=(UVLONG)term1.d+(UVLONG)term2.d;
		UVLONG sumcc=(UVLONG)term1.c+(UVLONG)term2.c;
		UVLONG sumbb=(UVLONG)term1.b+(UVLONG)term2.b;
		UVLONG sumaa=(UVLONG)term1.a+(UVLONG)term2.a;
	
		erg.d=(sumdd & MAXDDBYTE);
		sumcc += (sumdd >> 32);
		erg.c=(sumcc & MAXDDBYTE);
		sumbb += (sumcc >> 32);
		erg.b=(sumbb & MAXDDBYTE);
		sumaa += (sumbb >> 32);
		erg.a=(sumaa & MAXDDBYTE);
		if ((sumaa >> 32) != 0) {
			LOGMSG("Overflow add_abs_ZAB\n");
			exit(99);
		}
	}
}

inline void FPA_sub_abs_ovgl_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	// only absolute part, |a| >= |b|
	if (term1.vorz==0) {
		if (term2.vorz==0) {
			erg.setNull();
			return;
		}
		else {
			// Fehler
			LOGMSG("sub_abs_vgl out of range\n");
			exit(99);
		}
	} else if (term2.vorz==0) {
		erg.copyFrom(term1);
		return;
	}
	
	VLONG w=(VLONG)term1.d-(VLONG)term2.d;
	
	if (w < 0) {
		w += (VLONG)ZWEIHOCH32;
		erg.d=(w & MAXDDBYTE);	
		// carry over
		w=(VLONG)term1.c-((VLONG)term2.c+(VLONG)1);
	} else {
		erg.d=(w & MAXDDBYTE);	
		w=(VLONG)term1.c-(VLONG)term2.c;
	}

	if (w < 0) {
		w += (VLONG)ZWEIHOCH32;
		erg.c=(w & MAXDDBYTE);	
		// 1 Uebertrag
		w=(VLONG)term1.b-((VLONG)term2.b+(VLONG)1);
	} else {
		erg.c=(w & MAXDDBYTE);	
		w=(VLONG)term1.b-(VLONG)term2.b;
	}

	if (w < 0) {
		w += (VLONG)ZWEIHOCH32;
		erg.b=(w & MAXDDBYTE);	
		// carry over
		w=(VLONG)term1.a-((VLONG)term2.a+(VLONG)1);
	} else {
		erg.b=(w & MAXDDBYTE);	
		w=(VLONG)term1.a-(VLONG)term2.a;
	}

	if (w < 0) {
		printf("sub_ovgl: out of rnage\n");
		exit(99);
		erg.a=(w & MAXDDBYTE);	
	} else {
		erg.a=(w & MAXDDBYTE);	
	}
	
	if (
		(erg.b == 0) &&
		(erg.c == 0) &&
		(erg.a == 0) &&
		(erg.d == 0)
	) erg.vorz=0;
	else erg.vorz=1;
}

inline void FPA_sub_abs_ZAB(FPA& erg,const FPA a,const FPA b) {
	// only absolute value
	int32_t vgl=FPA_vgl_abs(a,b);
	if (vgl==0) {
		erg.setNull();
		return;
	}
	
	if (vgl>0) {
		// |a| > |b|
		FPA_sub_abs_ovgl_ZAB(erg,a,b);
	} else {
		// |a| < |b|
		FPA_sub_abs_ovgl_ZAB(erg,b,a);
		erg.vorz=-1;
	}
	
	erg.checkNull();
}

inline void FPA_sub_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if (term1.vorz == 0) {
		if (term2.vorz == 0) {
			erg.setNull();
		} else {
			erg.copyFrom(term2);
			erg.vorz*=-1;
			erg.checkNull();
		}
		return;
	} else if (term2.vorz == 0) {
		erg.copyFrom(term1);
		erg.checkNull();
		return;
	}

	erg.setNull();
	if (term1.vorz>0) {
		// a>=0
		if (term2.vorz>0) {
			// a>=0, b>0
			FPA_sub_abs_ZAB(erg,term1,term2);
		} else {
			// a>=0, b<0: a+|b|
			FPA_add_abs_ZAB(erg,term1,term2);
		}
	} else {
		// a<0
		if (term2.vorz>0) {
			// a<0,b>=0:  -|a|-|b| = -(|a|+|b|)
			FPA_add_abs_ZAB(erg,term1,term2);
			erg.vorz=-1;
		} else {
			// a<0,b<0: -|a| - -|b| = |b|-|a|
			FPA_sub_abs_ZAB(erg,term2,term1);
		}
	}
	
	erg.checkNull();
}

inline void FPA_sub2tpos_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	// b is assumed positive (sign ignored=

	if (term1.vorz == 0) {
		if (term2.vorz == 0) {
			erg.setNull();
		} else {
			erg.copyFrom(term2);
			erg.vorz=-1;
			erg.checkNull();
		}
		return;
	} else if (term2.vorz == 0) {
		erg.copyFrom(term1);
		erg.checkNull();
		return;
	}

	if (term1.vorz>0) {
		// a>=0
		FPA_sub_abs_ZAB(erg,term1,term2);
	} else {
		// a<0
		FPA_add_abs_ZAB(erg,term1,term2);
		erg.vorz=-1;
	}
	
	erg.checkNull();
}

inline void FPA_add_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if (term1.vorz==0) {
		erg.copyFrom(term2);
		erg.checkNull();
		return;
	} else 
	if (term2.vorz==0) {
		erg.copyFrom(term1);
		erg.checkNull();
		return;
	}
	
	erg.setNull();
	// check sign
	if (term1.vorz == term2.vorz) {
		FPA_add_abs_ZAB(erg,term1,term2);
		erg.vorz=term1.vorz;
		return; 
	} else {
		if ( (term1.vorz > 0) && (term2.vorz < 0) ) {
			// erg=a-|b|
			FPA_sub2tpos_ZAB(erg,term1,term2);
		} else {
			// a < 0, b > 0
			FPA_sub2tpos_ZAB(erg,term2,term1);
		}
	}

	erg.checkNull();
}

inline int FPA_vgl(const FPA aw,const FPA bw) {
	if (aw.vorz > 0) {
		if (bw.vorz <= 0) return +1;
		// a > 0
		return FPA_vgl_abs(aw,bw);
	} else if (aw.vorz == 0) {
		if (bw.vorz==0) return 0;
		else if (bw.vorz > 0) return -1;
		else if (bw.vorz < 0) return +1;
	} else {
		if (bw.vorz >= 0) return -1;
		return (-FPA_vgl_abs(aw,bw));
	}
	
	return 0;
}

inline int FPA_vgl_abs(const FPA term1,const FPA term2) {
	if (term1.a > term2.a) return +1;
	if (term1.a < term2.a) return -1;

	if (term1.b > term2.b) return +1;
	if (term1.b < term2.b) return -1;

	if (term1.c > term2.c) return +1;
	if (term1.c < term2.c) return -1;

	if (term1.d > term2.d) return +1;
	if (term1.d < term2.d) return -1;

	return 0;
}

inline int FPA_vgl(PFPA aw,PFPA bw) {
	if (aw->vorz > 0) {
		if (bw->vorz <= 0) return +1;
		// a > 0
		return FPA_vgl_abs(aw,bw);
	} else if (aw->vorz == 0) {
		if (bw->vorz==0) return 0;
		else if (bw->vorz > 0) return -1;
		else if (bw->vorz < 0) return +1;
	} else {
		if (bw->vorz >= 0) return -1;
		return (-FPA_vgl_abs(aw,bw));
	}
	
	return 0;
}

inline int FPA_vgl_abs(PFPA term1,PFPA term2) {
	if (term1->a > term2->a) return +1;
	if (term1->a < term2->a) return -1;

	if (term1->b > term2->b) return +1;
	if (term1->b < term2->b) return -1;

	if (term1->c > term2->c) return +1;
	if (term1->c < term2->c) return -1;

	if (term1->d > term2->d) return +1;
	if (term1->d < term2->d) return -1;

	return 0;
}

inline FPA maximumD(const FPA a,const FPA b) {
	FPA erg;
	if (FPA_vgl(a,b) < 0) {
		erg.copyFrom(b);
	} else {
		erg.copyFrom(a);
	}
	
	return erg;
}

inline FPA maximumD(const FPA a,const FPA b,const FPA c,const FPA d) {
	FPA ab=maximumD(a,b);
	FPA cd=maximumD(c,d);
	
	return maximumD(ab,cd);
}

inline FPA minimumD(const FPA a,const FPA b) {
	FPA erg;
	if (FPA_vgl(a,b) > 0) {
		erg.copyFrom(b);
	} else {
		erg.copyFrom(a);
	}
	
	return erg;
}

inline FPA minimumD(const FPA a,const FPA b,const FPA c,const FPA d) {
	FPA ab=minimumD(a,b);
	FPA cd=minimumD(c,d);

	return minimumD(ab,cd);
}
#endif

#ifdef _F107
f107_o operator*(int a,const f107_o b) {
	f107_o w=a;
	f107_o erg=w*b;
	
	return erg;
}

inline void minimaxF107AB(f107_o& mi,f107_o& ma,const f107_o a,const f107_o b) {
	if (a < b) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

inline void minimaxF107ABCD(f107_o& mi,f107_o& ma,
	const f107_o a,const f107_o b,
	const f107_o c,const f107_o d
) {
	f107_o miab,micd,maab,macd;
	minimaxF107AB(miab,maab,a,b);
	minimaxF107AB(micd,macd,c,d);
	
	if (miab < micd) {
		mi=miab;
	} else {
		mi=micd;
	}

	if (maab > macd) {
		ma=maab;
	} else {
		ma=macd;
	}
}

inline void minimaxF107AB(Pf107& mi,Pf107& ma,f107_o* a,f107_o* b) {
	if ( (*a) < (*b) ) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

inline void minimaxF107ABCD(Pf107& mi,Pf107& ma,
	f107_o* a,f107_o* b,
	f107_o* c,f107_o* d
) {
	Pf107 miab,micd,maab,macd;
	minimaxF107AB(miab,maab,a,b);
	minimaxF107AB(micd,macd,c,d);
	
	if ((*miab) < (*micd)) {
		mi=miab;
	} else {
		mi=micd;
	}

	if ((*maab) > (*macd)) {
		ma=maab;
	} else {
		ma=macd;
	}
}
#endif

int32_t getPower2Exponent(const uint64_t aw) {
	int32_t exponent=-1;
	uint64_t w=aw;
	for(int32_t i=0;i<64;i++) {
		if ( (w & 0b1) != 0) {
			if (exponent<0) exponent=i;
			else {
				// more than one => not a power of 2
				exponent=-1;
				break;
			}
		}
		
		w >>= 1;
	}
	
	if (exponent<0) {
		LOGMSG("Error. Range must be a power of 2.\n");
		exit(99);
	}
	
	return exponent;
}

int32_t makePowerOf2(int32_t& avalue) {
	return (1 << ( (int32_t)ceil(log(avalue)/log(2.0)) ) );
}

void freeRevCGMem(void) {
	if (parentmgr) {
		delete parentmgr;
		parentmgr=NULL;
	}
	
	if (data5->revcgYX) {
		delete[] data5->revcgYX;
		data5->revcgYX=NULL;
	}
}


// main

int32_t main(int32_t argc,char** argv) {
	int32_t c0=clock();
	
	flog=fopen("juliatsacoredyn.log.txt","at");
	fprintf(flog,"\n-----------------\n");
	printf("\n  FUNC=string\n  CMD=calc\n  or CMD=period[,PP,M2]\n  LEN=n / c=re,im or re,re,im,im / A=re,im / REVCG=n\n  NOPOTW\n\n");

	parentmgr=new ParentManager;
	
	int32_t cmd=CMD_CALC;
	
	// standard
	getBoundingBoxfA=getBoundingBoxfA_z2c;
	_FUNC=FUNC_Z2C;
	RANGE0=-2;
	RANGE1=2;
	seedC0re=seedC1re=-1.0;
	seedC0im=seedC1im=0.0;
	FAKTORAre=FAKTORAim=0.0;
	REVCGBITS=4;
	SCREENWIDTH=(1 << 10);
	
	for(int32_t i=1;i<argc;i++) {
		upper(argv[i]);
		if (strstr(argv[i],"FUNC=")==argv[i]) {
			_FUNC=getfuncidx(&argv[i][5]);
		} else
		if (strstr(argv[i],"CMD=")==argv[i]) {
			if (strstr(&argv[i][4],"PERIOD")==&argv[i][4]) {
				cmd=CMD_PERIOD;
				if (strstr(argv[i],",PP")) _PERIODICPOINTS=1;
			} else
			if (strstr(&argv[i][4],"CONVERT")==&argv[i][4]) {
				convert_raw_into_newstructure();
				fclose(flog);
				return 0;
			} 
		} else if (strstr(argv[i],"C=")==argv[i]) {
			double r0,r1,i0,i1; // not NTYP
			// command line parameters are always considered double no matter the datatype used
			if (sscanf(&argv[i][2],"%lf,%lf,%lf,%lf",&r0,&r1,&i0,&i1) == 4) {
				double w1=floor(r0*DENOM225)/DENOM225;
				double w2=floor(r1*DENOM225)/DENOM225;
				if (w1 > w2) {
					seedC0re=w2;
					seedC1re=w1;
				} else {
					seedC0re=w1;
					seedC1re=w2;
				}
				w1=floor(i0*DENOM225)/DENOM225;
				w2=floor(i1*DENOM225)/DENOM225;
				if (w1 > w2) {
					seedC0im=w2;
					seedC1im=w1;
				} else {
					seedC0im=w1;
					seedC1im=w2;
				}
			} else
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				seedC0re=seedC1re=floor(r0*DENOM225)/DENOM225;
				seedC0im=seedC1im=floor(i0*DENOM225)/DENOM225;
			}
		} else if (strstr(argv[i],"A=")==argv[i]) {
			double r0,i0;
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				FAKTORAre=floor(r0*DENOM225)/DENOM225;
				FAKTORAim=floor(i0*DENOM225)/DENOM225;
			}
		} 
		if (strstr(argv[i],"LEN=")==argv[i]) {
			int a;
			if (sscanf(&argv[i][4],"%i",&a) == 1) {
				if (a < 8) a=8;
				if (a > 31) a=31;
				SCREENWIDTH=(1 << a);
			}
		} else
		if (strstr(argv[i],"REVCG=")==argv[i]) {
			int a;
			if (sscanf(&argv[i][6],"%i",&a) == 1) REVCGBITS=a;
		} else
		if (strstr(argv[i],"RANGE=")==argv[i]) {
			int a;
			if (sscanf(&argv[i][6],"%i",&a) == 1) {
				if (a<0) a=-a;
				RANGE1=makePowerOf2(a);
				if (RANGE1 != a) LOGMSG2("range adjusted to next-bigger power of 2: %i\n",RANGE1);
				RANGE0=-RANGE1;
			}
		}
	} // i
	
	#ifdef _FPA
	COMPLETE0.set_vlong(RANGE0);
	COMPLETE1.set_vlong(RANGE1);
	#endif
	#ifndef _FPA
	COMPLETE0=RANGE0;
	COMPLETE1=RANGE1;
	#endif
	
	if (SCREENWIDTH < 256) SCREENWIDTH=256;
	
	REFINEMENTLEVEL=(int)ceil(log(SCREENWIDTH)/log(2.0));

	if (REVCGBITS < 4) REVCGBITS=4;
	
	// SCREENWIDTH / (2^REVCGBITS) <= 2^15
	// REVCG allocated en bloc => adjust REVCGBITS
	// revcg-Parents numbered with 16 bit
	while ( (SCREENWIDTH >> REVCGBITS) > (1 << 15)) REVCGBITS++;

	REVCGBLOCKWIDTH=(1 << REVCGBITS);
	if (SCREENWIDTH >= REVCGBLOCKWIDTH) {
		REVCGmaxnumber=SCREENWIDTH >> REVCGBITS;
	} else REVCGmaxnumber=1;
	REVCGmaxnumberQ=REVCGmaxnumber*REVCGmaxnumber;
	
	char fn[1024];
	fn[0]=0;
	setfunc_and_bitprecision(_FUNC,fn);

	if (fn[0]<=0) {
		LOGMSG("Error. Name of function not defined.\n");
		exit(99);
	}
	
	LOGMSG2("file principal part %s\n",fn);
	
	data5=new Data5;

	double w=(RANGE1-RANGE0) / (double)SCREENWIDTH;
	scaleRangePerPixel=w;

	w=(double)SCREENWIDTH / (RANGE1-RANGE0);
	scalePixelPerRange=w;
	// w is a power of 2
	scalePixelPerRangeExponent=getPower2Exponent( (uint64_t)w );
	
	// heat-map color table
	basinpal.setlen(800);
	basinpal.setInterval(0.0/8.0,1.0/8.0,255,255,0,255,0,0);
	basinpal.setInterval(1.0/8.0,2.0/8.0,255,255,0,0,255,0);
	basinpal.setInterval(2.0/8.0,3.0/8.0,0,255,0,0,255,255);
	basinpal.setInterval(3.0/8.0,4.0/8.0,0,255,255,0,0,255);
	basinpal.setInterval(4.0/8.0,5.0/8.0,0,0,255,255,0,255);
	basinpal.setInterval(5.0/8.0,6.0/8.0,255,0,255,255,127,0);
	basinpal.setInterval(6.0/8.0,7.0/8.0,255,127,0,127,127,255);
	basinpal.setInterval(7.0/8.0,8.0/8.0,127,127,255,255,255,127);
	
	PlaneRect plane;
	
	// if raw data file exists: read it and if necessary blow up the pixels 2fold
	if (data5->readRawBlowUp() <= 0) {
		// data5 object is - no matter what data it holds - considered uninitialised
		printf("searching for special exterior ... ");

		// at the start: everything is considered gray
		encgrayx0=encgrayy0=0;
		encgrayx1=encgrayy1=SCREENWIDTH-16;
		planegrayx0=planegrayy0=COMPLETE0;
		planegrayx1=planegrayy1=COMPLETE1;
		plane.x0=COMPLETE0;
		plane.x1=COMPLETE1;
		plane.y0=COMPLETE0;
		plane.y1=COMPLETE1;

		int32_t MEMWIDTH=(SCREENWIDTH >> 4);
		for(int32_t y=0;y<SCREENWIDTH;y++) {
			data5->zeilen[y]=data5->datamgr->getMemory(MEMWIDTH);
		}

		// squares(pixels) whose bounding box lies completely in the special exterior
		find_special_exterior_hitting_squares();
	} else {
		// add/substract 16 for buffer reason's
		plane.x0=(encgrayx0-16) * scaleRangePerPixel + COMPLETE0;
		plane.x1=(encgrayx1+16) * scaleRangePerPixel + COMPLETE0;
		plane.y0=(encgrayy0-16) * scaleRangePerPixel + COMPLETE0;
		plane.y1=(encgrayy1+16) * scaleRangePerPixel + COMPLETE0;
	}
	
	LOGMSG5("\ngray in pixel region [%i..%i] x [%i..%i]\n",
		encgrayx0,encgrayx1,encgrayy0,encgrayy1);
		
	#ifdef _FPA
	LOGMSG5("  roughly %.5lg..%.5lg x %.5lg..%.5lg used\n",
		plane.x0.convert_to_double(),
		plane.x1.convert_to_double(),
		plane.y0.convert_to_double(),
		plane.y1.convert_to_double());
	#endif
	#ifndef _FPA
	LOGMSG5("  roughly %.5lg..%.5lg x %.5lg..%.5lg used\n",
		(double)plane.x0,(double)plane.x1,(double)plane.y0,(double)plane.y1);
	#endif
	
	if (
		(encgrayx0 > (SCREENWIDTH >> 2) ) &&
		(encgrayx1 < (3*(SCREENWIDTH >> 2) ) ) &&
		(encgrayy0 > (SCREENWIDTH >> 2) ) &&
		(encgrayy1 < (3*(SCREENWIDTH >> 2) ) ) 
	) {
		LOGMSG("  range could be adjusted (half is enough)\n");
	}
	
	char tmp[4096];
	fprintf(flog,"%s * 2^-%i\n",seedCstr225(tmp),BASEDENOMINATOR);
	fprintf(flog,"(if needed): %s * 2^-%i\n",FAKTORAstr225(tmp),BASEDENOMINATOR);
	
	// main routine
	// has to be done for CALC, PERIOD
	compute(); 

	// storing image(s)
	if (
		(SAVE128KIMAGE>0) ||
		(SCREENWIDTH <= 65536)
	) {
		printf("saving image ...\n");
		data5->saveBitmap4(fn);
	}
	// storing a trustworthily downscaled image
	if (SCREENWIDTH > 65536) {
		printf("downscaling in a trustworthy manner ...\n");
		data5->saveBitmap4_twd(fn,-1);
	}
	// storing raw data
	printf("saving raw data ...\n");
	data5->saveRaw(fn);
	
	// free memory to get enough to allocate for the periodicty check
	printf("freeing reverse cell graph memory ...\n");
	freeRevCGMem();

	// data is now computed or loaded
	if (cmd==CMD_PERIOD) {
		if (interiorpresent>0) {
			periodicity(fn);
		} else {
			LOGMSG("No interior present. Periodicity check skipped.\n");
		}
	} 
	
	delete data5;
	
	int32_t c1=clock();
	LOGMSG2("\nduration %.0lf sec\n",(double)(c1-c0)/CLOCKS_PER_SEC);
	
	LOGMSG2("%I64d bounding boxes calculated\n",ctrbbxfa);

	fclose(flog);
	
    return 0;
}

