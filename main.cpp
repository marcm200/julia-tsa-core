/*

	Marc Meidlinger, July 2019-March 2020
	Implementation of the ingenious algorithm from:
	
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

#define _DOUBLE
//#define _LONGDOUBLE
//#define _QUADMATH
//#define _FPA
//#define _F107
//#define _F161


// typedefs

typedef uint8_t BYTE;
typedef BYTE* PBYTE;
typedef uint16_t DBYTE;
typedef DBYTE *PDBYTE;
typedef uint32_t DDBYTE;
typedef DDBYTE *PDDBYTE;

// chunk size depe4ndend on operating system
// win32 => 512 MB
//#define _CHUNK512

#ifdef _CHUNK512
const uint64_t CHUNKSIZE=( (uint64_t)1 << 27 );
#else
// win64 => 1 GB
const uint64_t CHUNKSIZE=( (uint64_t)1 << 30 );
#endif

#ifdef _QUADMATH
#include "quadmath.h"
typedef __float128 NTYP;
typedef __float128 *Pfloat128;
const char NNTYPSTR[]="f128_";
inline void	minimaxQDAB(__float128&,__float128&,const __float128,const __float128);
inline void	minimaxQDABCD(__float128&,__float128&,const __float128,const __float128,const __float128,const __float128);
const char NTS[]="QD";
#endif

#ifdef  _LONGDOUBLE
typedef long double NTYP;
const char NNTYPSTR[]="ld_";
const char NTS[]="LD";
inline void minimaxldAB(long double&,long double&,const long double,const long double);
inline void minimaxldABCD(long double&,long double&,const long double,const long double,const long double,const long double);
#endif

#ifdef _DOUBLE
typedef double NTYP;
const char NNTYPSTR[]="";
const char NTS[]="D";
#endif

// those are generally eecessary for dtcheck
inline void	minimaxdAB(double&,double&,const double,const double);
inline void	minimaxdABCD(double&,double&,const double,const double,const double,const double);

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

#ifdef _F161
// https://mrob.com/pub/math/f161.html#source
// f107 ist auch nötig
#include "f107_o.cpp"
#include "f161_o.cpp"
typedef f161_o NTYP;
typedef f161_o *Pf161;
inline void	minimaxF161AB(f161_o&,f161_o&,const f161_o,const f161_o);
inline void minimaxF161ABCD(f161_o&,f161_o&,const f161_o,const f161_o,const f161_o,const f161_o);
inline void minimaxF161AB(Pf161&,Pf161&,f161_o*,f161_o*);
inline void minimaxF161ABCD(Pf161&,Pf161&,f161_o*,f161_o*,f161_o*,f161_o*);
const char NNTYPSTR[]="f161_";
const char NTS[]="F6";
f161_o operator*(int,const f161_o);
double specialFloor161(const f161_o);
#endif

#ifdef _FPA
struct FPA {
	// UVLONG or DDBYTE not performance-relevant
	int8_t vorz; // -1,0,1
	DDBYTE a,b,c,d; // a + b*2^-32 + c*2^-64 + d*2^-96
	
	FPA();
	FPA(const double);
	FPA(const int64_t);
	
	void set_vlong(const int64_t);
	void setNull(void);
	void squareTo(FPA&);
	void checkNull(void);
	void copyFrom(const FPA);
	void copyFrom(FPA*);
	void set_double32(const double);
	char* str(char*);
	double convert_to_double(void);

	void shiftLeft(const int32_t); // Multiplication
	void shiftRight(const int32_t); // Division
	FPA& operator=(const FPA);
};

typedef FPA NTYP;
typedef FPA *PFPA;
const char NNTYPSTR[]="fpa_";
const char NTS[]="FP";

static inline int32_t scrcoord_as_lowerleft(FPA*);
bool operator<(const FPA,const FPA);
bool operator>(const FPA,const FPA);
bool operator<=(const FPA,const FPA);
bool operator>=(const FPA,const FPA);
bool operator!=(const FPA,const FPA);
bool operator==(const FPA,const FPA);
FPA operator-(const FPA,const FPA);
FPA operator+(const FPA,const FPA);
FPA operator*(const FPA,const FPA);
FPA operator*(const int,const FPA);
inline void FPA_sub_abs_ZAB(FPA&,const FPA,const FPA);
inline void FPA_add_abs_ZAB(FPA&,const FPA,const FPA);
inline void FPA_add_ZAB(FPA&,const FPA,const FPA);
inline void FPA_sub_ZAB(FPA&,const FPA,const FPA);
inline void FPA_sub2tpos_ZAB(FPA&,const FPA,const FPA);
inline void FPA_sub_abs_ovgl_ZAB(FPA&,const FPA,const FPA);

// addition with passing pointers instead of const objects
inline void FPA_sub_abs_ZAB(FPA&,FPA*,FPA*);
inline void FPA_add_abs_ZAB(FPA&,FPA*,FPA*);
inline void FPA_add_ZAB(FPA&,FPA*,FPA*);
inline void FPA_sub_ZAB(FPA&,FPA*,FPA*);
inline void FPA_sub2tpos_ZAB(FPA&,FPA*,FPA*);
inline void FPA_sub_abs_ovgl_ZAB(FPA&,FPA*,FPA*);

inline int FPA_vgl_abs(const FPA,const FPA);
inline int FPA_vgl(const FPA,const FPA);
inline int FPA_vgl_abs(PFPA,PFPA);
inline int FPA_vgl(PFPA,PFPA);
inline void FPA_mul_ZAuvlong(FPA&,const FPA,const uint64_t);
inline void FPA_mul_ZAB(FPA&,const FPA,const FPA);
inline void FPA_mul_ZAB(FPA&,FPA*,FPA*);
inline void minimaxFPAAB(FPA&,FPA&,const FPA,const FPA);
inline void minimaxFPAAB(PFPA&,PFPA&,PFPA,PFPA);
inline void minimaxFPAABCD(FPA&,FPA&,const FPA,const FPA,const FPA,const FPA);
// array of 4 elements
inline void minimaxFPAABCD(PFPA&,PFPA&,FPA*);
inline void minimaxFPAABCD(PFPA&,PFPA&,PFPA,PFPA,PFPA,PFPA);
inline int64_t floorFPA(const FPA);
#endif


// const definitions

const int64_t UINT32MAX=(((int64_t)1 << 32) - 1);
const int64_t MAXDDBYTE=UINT32MAX;
const int64_t ZWEIHOCH32=((int64_t)1 << 32);
const uint32_t SQUARE_GRAY=0b00;
const uint32_t SQUARE_WHITE=0b01;
const uint32_t SQUARE_BLACK=0b10;
const uint32_t SQUARE_GRAY_POTENTIALLY_WHITE=0b11;
const uint32_t COLORRED=4;

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

const int32_t MAXPTR=2048;

const uint32_t power2[]={
	(uint32_t)(((uint32_t)1 << 1)-1),
	(uint32_t)(((uint32_t)1 << 2)-1),
	(uint32_t)(((uint32_t)1 << 3)-1),
	(uint32_t)(((uint32_t)1 << 4)-1),
	(uint32_t)(((uint32_t)1 << 5)-1),
	(uint32_t)(((uint32_t)1 << 6)-1),
	(uint32_t)(((uint32_t)1 << 8)-1),
	(uint32_t)(((uint32_t)1 << 9)-1),
	(uint32_t)(((uint32_t)1 << 0)-1),
	(uint32_t)(((uint32_t)1 << 10)-1),
	(uint32_t)(((uint32_t)1 << 11)-1),
	(uint32_t)(((uint32_t)1 << 12)-1),
	(uint32_t)(((uint32_t)1 << 13)-1),
	(uint32_t)(((uint32_t)1 << 14)-1),
	(uint32_t)(((uint32_t)1 << 15)-1),
	(uint32_t)(((uint32_t)1 << 16)-1),
	(uint32_t)(((uint32_t)1 << 17)-1),
	(uint32_t)(((uint32_t)1 << 18)-1),
	(uint32_t)(((uint32_t)1 << 19)-1),
	(uint32_t)(((uint32_t)1 << 20)-1),
	(uint32_t)(((uint32_t)1 << 21)-1),
	(uint32_t)(((uint32_t)1 << 22)-1),
	(uint32_t)(((uint32_t)1 << 23)-1),
	(uint32_t)(((uint32_t)1 << 24)-1),
	(uint32_t)(((uint32_t)1 << 25)-1),
	(uint32_t)(((uint32_t)1 << 26)-1),
	(uint32_t)(((uint32_t)1 << 27)-1),
	(uint32_t)(((uint32_t)1 << 28)-1),
	(uint32_t)(((uint32_t)1 << 29)-1),
	(uint32_t)(((uint32_t)1 << 30)-1)
};

enum { 
	CMD_CALC=1,
	CMD_PERIOD,
	CMD_FASTDTCHECK
};

// not all are implemented, but values need be the same as
// in the TSApredictor
enum { 
	FUNC_Z2C=0,FUNC_Z2AZC=1,FUNC_Z3AZC=2,
	FUNC_Z4AZC=3,FUNC_Z5AZC=4,FUNC_Z6AZC=5,
	FUNC_AZ2ZC,FUNC_Z7AZC,FUNC_Z8AZC,
	FUNC_2ITZ2C,FUNC_BZ2AZC,FUNC_BZ5AZC,
	FUNC_BZ3AZC,
	
	FUNCANZ
};

const char funcname[][32] = {
	"Z2C","Z2AZC","Z3AZC","Z4AZC","Z5AZC","Z6AZC",
	"AZ2ZC","Z7AZC","Z8AZC","2ITZ2C","BZ2AZC",
	"BZ5AZC","BZ3AZC"
};

// clearing a specific pixel's color via bitwise AND-mask
const uint32_t COLOR_CLEARMASK[] = {
	(uint32_t)(UINT32MAX - (0b11)),
	(uint32_t)(UINT32MAX - (0b11 << 2)),
	(uint32_t)(UINT32MAX - (0b11 << 4)),
	(uint32_t)(UINT32MAX - (0b11 << 6)),
	(uint32_t)(UINT32MAX - (0b11 << 8)),
	(uint32_t)(UINT32MAX - (0b11 << 10)),
	(uint32_t)(UINT32MAX - (0b11 << 12)),
	(uint32_t)(UINT32MAX - (0b11 << 14)),
	(uint32_t)(UINT32MAX - (0b11 << 16)),
	(uint32_t)(UINT32MAX - (0b11 << 18)),
	(uint32_t)(UINT32MAX - (0b11 << 20)),
	(uint32_t)(UINT32MAX - (0b11 << 22)),
	(uint32_t)(UINT32MAX - (0b11 << 24)),
	(uint32_t)(UINT32MAX - (0b11 << 26)),
	(uint32_t)(UINT32MAX - (0b11 << 28)),
	(uint32_t)(UINT32MAX - (0b11 << 30))
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
const int32_t MAXFATOUCOMPONENTS=8192;
#else
const int32_t MAXFATOUCOMPONENTS=65500;
#endif
const int32_t MAXCYCLES=110;
const int32_t MAXPERIODICPOINTS=1024;
const int32_t FATOUCOMPONENTCOLOROFFSET=24;

const int32_t M3MAXCYCLES=1024;
const int32_t M3MAXORBITLEN=(1 << 20);


// structs

struct RGB4 {
	uint8_t B,G,R,alpha; // in the order of the bitmap file structure
};

struct Palette4 {
	RGB4 rgbs[256]; 
	
	void setPaletteRGB(const int,const int,const int,const int);
};

struct PlaneRect {
	NTYP x0,x1,y0,y1;
};

// for use in fastdtcheck
struct PlaneRect_double {
	double x0,x1,y0,y1;
};

// for the reverse cell graph: position of one tile
struct Parent {
	uint16_t BX,BY; 
};

typedef Parent *PParent;

struct Int2 {
	int32_t x,y;
};

typedef Int2 *PInt2;

struct DFSPunkt {
	int32_t x,y;
	DBYTE tiefe;
};

struct ListeFIFO {
	int32_t next_readpos,next_writepos;
	int32_t allokiereAmStueck;
	Int2* werte;
	
	ListeFIFO();
	virtual ~ListeFIFO();
	void start(void);
	char read(int32_t&,int32_t&);
	char write(const int32_t,const int32_t);
};

struct ListeDFS {
	int32_t anz;
	int32_t allokiereAmStueck;
	DFSPunkt* werte;
	
	ListeDFS();
	virtual ~ListeDFS();
	void start(void);
	char read(int32_t&,int32_t&,DBYTE&);
	char write(const int32_t,const int32_t,const DBYTE);
};

struct Cycle {
	int32_t len;
	DBYTE immediateBasinColorIdx;
	DBYTE attractionBasinColorIdx;
	DBYTE fatouidx0,fatouidx1; // Index auf ibfcomponents = immedaite basin fatou components
};

struct CycleM3 {
	uint8_t color;
	int32_t len;
	int32_t *perblobs;
};

// pixel-coordinate rectangle
struct ScreenRect {
	int32_t x0,x1,y0,y1;
};

typedef ScreenRect *PScreenRect;

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
	// how many parents are needed here? First pass result
	int32_t howmany; 
	// flag, whether this tile has to be checked (again) for gray pixels and the bounding
	// boxes' hits after one iteration
	int8_t tovisit;
	int8_t containsgray;
	int32_t memused;
	// parents as variable size array
	Parent* parent;
		
	RevCGBlock();
	void addParent(const int32_t,const int32_t);
};

struct Int2Manager {
	Int2* current;
	int32_t allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PInt2 ptr[MAXPTR];
	int32_t anzptr;
	
	Int2Manager();
	virtual ~Int2Manager();
	Int2* getMemory(const int32_t);
};

struct ScreenRectManager {
	ScreenRect* current;
	int32_t allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PScreenRect ptr[MAXPTR];
	int32_t anzptr;
	
	ScreenRectManager();
	virtual ~ScreenRectManager();
	ScreenRect* getMemory(const int32_t);
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

typedef int32_t *PInt;

struct IntManager {
	int32_t* current;
	int32_t allocatedIdx,freeFromIdx,allocatePerBlockIdx;
	PInt ptr[MAXPTR];
	int32_t anzptr;
	
	IntManager();
	virtual ~IntManager();
	int32_t* getMemory(const int32_t);
};

struct ByteManager {
	BYTE* current;
	int32_t allocatedIdx,freeFromIdx,allocatePerBlockIdx;
	PBYTE ptr[MAXPTR];
	int32_t anzptr;
	
	ByteManager();
	virtual ~ByteManager();
	BYTE* getMemory(const int32_t);
};

struct ArrayDDByteManager {
	DDBYTE* current;
	int allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PDDBYTE ptr[MAXPTR];
	int anzptr;
	
	ArrayDDByteManager ();
	virtual ~ArrayDDByteManager ();
	PDDBYTE getMemory(const int32_t);
};

struct VGridRow {
	int32_t x0,x1; // offset
	BYTE *shrinkptr;
};

// main object
struct Data5 {
	uint32_t** zeilen;
	Gray_in_row* memgrau;
	uint8_t* graudensity;
	RevCGBlock* revcgYX;
	ArrayDDByteManager* datamgr;
	VGridRow *vgridYX;
	PScreenRect *pcscr;
	ScreenRectManager *pcscrmgr;
		
	Data5();
	virtual ~Data5();
		
	void saveBitmap4_twd(const char*,const int32_t);
	void saveRaw(const char*);
	int32_t readRawBlowUp(void);
	
	void precomputeScreenRect(void);
};

struct RefPoint {
	int32_t x;
	int32_t blobid;
};

typedef RefPoint *PRefPoint;

struct Streak {
	int32_t x0,x1,y;
};

typedef Streak *PStreak;

const int32_t STREAK_MAXBLOCKS=1024;
const int32_t STREAK_PERBLOCK=((int32_t)1 << 28);

// includes StreakManager
struct StreakArray {
	PStreak ptr[STREAK_MAXBLOCKS];
	int32_t anz[STREAK_MAXBLOCKS];
	int32_t writeidx,writenextpos;
	
	StreakArray();
	virtual ~StreakArray();
	
	void fastEmpty(void);
	void pushStreak(const int32_t,const int32_t,const int32_t);
	void popStreak(Streak&);
};

struct RefPointManager {
	RefPoint* current;
	int32_t allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PRefPoint ptr[MAXPTR];
	int32_t anzptr;
	
	RefPointManager();
	virtual ~RefPointManager();
	RefPoint* getMemory(const int32_t);
};

const int32_t ANZREFPOINTSINBLOCK=16384;

struct RefList {
	RefPoint *points;
	int32_t anz;
	int32_t memused;
	
	RefList();
	
	// garbage collection is done outside in manager
	void addXB(const int32_t,const int32_t);
	RefPoint* getRefPtr(const int32_t);
};

typedef RefList *PRefList;

struct RefPointArray {
	RefList *listY;
	
	RefPointArray();
	virtual ~RefPointArray();
	
	void addRefPoint(const int32_t,const int32_t,const int32_t);
	RefPoint* getRefPtr(const int32_t,const int32_t);
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

const int32_t MAXHELPERVALUES=32;
// adjust HELPERPERBLOCKMODULO if BITS are changed
const int32_t HELPERPERBLOCKBITS=16;
const int32_t HELPERPERBLOCKMODULO=0b1111111111111111;
const int32_t MAXHELPERPERBLOCK=(1 << HELPERPERBLOCKBITS);
const int32_t MAXHELPERBLOCKS=(1 << 16); // i.e. max refinement 32
const int32_t DIRECTIONX=1;
const int32_t DIRECTIONY=2;

struct Helper {
	// values indexed by constants of the form
	// _Z6AZC_X_FAKTORAREmulAx0 etc.
	NTYP val[MAXHELPERVALUES];
};

struct Helper_double {
	double val[MAXHELPERVALUES];
};

typedef Helper *PHelper;
typedef Helper_double *PHelper_double;

struct HelperManager {
	Helper* current;
	int32_t allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PHelper ptr[MAXPTR];
	int32_t anzptr;
	
	HelperManager();
	virtual ~HelperManager();
	PHelper getMemory(const int32_t);
};

struct Helper_doubleManager {
	Helper_double* current;
	int32_t allokierteIdx,freiAbIdx,allokierePerBlockIdx;
	PHelper_double ptr[MAXPTR];
	int32_t anzptr;
	
	Helper_doubleManager();
	virtual ~Helper_doubleManager();
	PHelper_double getMemory(const int32_t);
};

// one object for x-independent IA-expressions
// one object for y-independent
// accessed per screen coordinate of lower left
// assuming size of the interval being exact
// scaleRangePerPixel
// CAVE: Helper cannot be used to find_special_exterior
// as there the "tiles" have wider width.
struct HelperAccess {
	// access by e.g. y-coordinate works as
	// helperblock[y >> HELPERPERBLOCKBITS][y & HELPERPERBLOCKMODULO];
	// to get the helper-object
	// Helpers are grouped in chunks of MAXHELPERPERBLOCK
	// to reduce fragmentation of memory while still allowing
	// in principle arbitrary many rows and refinement levels
	PHelper *helperblocks;
	int32_t blockanz;
	
	HelperAccess();
	virtual ~HelperAccess();
	
	void initMemory(void);
	inline PHelper getHelper(const int32_t);
	void precompute(const int32_t);
};

struct HelperAccess_double {
	PHelper_double *helperblocks;
	int32_t blockanz;
	
	HelperAccess_double();
	virtual ~HelperAccess_double();
	
	void initMemory(void);
	inline PHelper_double getHelper(const int32_t);
	void precompute(const int32_t);
};

// HelperObject-constants

enum {
	_HELPER_Z8AZC_Ydep_mi4=0,
	_HELPER_Z8AZC_Ydep_ma4,
	_HELPER_Z8AZC_Ydep_mi228,
	_HELPER_Z8AZC_Ydep_ma228,
	_HELPER_Z8AZC_Ydep_mi13,
	_HELPER_Z8AZC_Ydep_ma13,
	_HELPER_Z8AZC_Ydep_y03,
	_HELPER_Z8AZC_Ydep_y13,
	_HELPER_Z8AZC_Ydep_y05,
	_HELPER_Z8AZC_Ydep_y15,
	_HELPER_Z8AZC_Ydep_t4a,
	_HELPER_Z8AZC_Ydep_t4b,
	_HELPER_Z8AZC_Ydep_mi7,
	_HELPER_Z8AZC_Ydep_ma7,
	_HELPER_Z8AZC_Ydep_ma15_minus_mi6,
	_HELPER_Z8AZC_Ydep_mi15_minus_ma6,

	_HELPER_Z8AZC_Ydep_ANZ
};

enum {
	_HELPER_Z8AZC_Xdep_mi370=0,
	_HELPER_Z8AZC_Xdep_ma370,
	_HELPER_Z8AZC_Xdep_mi11,
	_HELPER_Z8AZC_Xdep_ma11,
	_HELPER_Z8AZC_Xdep_mi128,
	_HELPER_Z8AZC_Xdep_ma128,
	_HELPER_Z8AZC_Xdep_t1a,
	_HELPER_Z8AZC_Xdep_t1b,
	_HELPER_Z8AZC_Xdep_t2a,
	_HELPER_Z8AZC_Xdep_t2b,
	_HELPER_Z8AZC_Xdep_t3a,
	_HELPER_Z8AZC_Xdep_t3b,
	_HELPER_Z8AZC_Xdep_C0re_plus_mi5_plus_mi10,
	_HELPER_Z8AZC_Xdep_C1re_plus_ma5_plus_ma10,
	_HELPER_Z8AZC_Xdep_C0im_plus_mi8,
	_HELPER_Z8AZC_Xdep_C1im_plus_ma8,

	_HELPER_Z8AZC_Xdep_ANZ
};

enum {
	_HELPER_Z7AZC_Xdep_C0im_plus_mi4=0,
	_HELPER_Z7AZC_Xdep_C1im_plus_ma4,
	_HELPER_Z7AZC_Xdep_mi1,
	_HELPER_Z7AZC_Xdep_ma1,
	_HELPER_Z7AZC_Xdep_21mi1,
	_HELPER_Z7AZC_Xdep_21ma1,
	_HELPER_Z7AZC_Xdep_35mi1,
	_HELPER_Z7AZC_Xdep_35ma1,
	_HELPER_Z7AZC_Xdep_mi6,
	_HELPER_Z7AZC_Xdep_ma6,
	_HELPER_Z7AZC_Xdep_7mi6,
	_HELPER_Z7AZC_Xdep_7ma6,
	_HELPER_Z7AZC_Xdep_Are_plus_mi6,
	_HELPER_Z7AZC_Xdep_Are_plus_ma6,
	
	_HELPER_Z7AZC_Xdep_ANZ
};

enum {
	_HELPER_Z7AZC_Ydep_C0re_minus_ma3=0,
	_HELPER_Z7AZC_Ydep_C1re_minus_mi3,
	_HELPER_Z7AZC_Ydep_mi2,
	_HELPER_Z7AZC_Ydep_ma2,
	_HELPER_Z7AZC_Ydep_21mi2,
	_HELPER_Z7AZC_Ydep_21ma2,
	_HELPER_Z7AZC_Ydep_35mi2,
	_HELPER_Z7AZC_Ydep_35ma2,
	_HELPER_Z7AZC_Ydep_mi7,
	_HELPER_Z7AZC_Ydep_ma7,
	_HELPER_Z7AZC_Ydep_7mi7,
	_HELPER_Z7AZC_Ydep_7ma7,
	_HELPER_Z7AZC_Ydep_Are_minus_mi7,
	_HELPER_Z7AZC_Ydep_Are_minus_ma7,
	
	_HELPER_Z7AZC_Ydep_ANZ
};

enum {
	_HELPER_Z6AZC_Xdep_mi6=0,
	_HELPER_Z6AZC_Xdep_ma6,
	_HELPER_Z6AZC_Xdep_tmp1,
	_HELPER_Z6AZC_Xdep_tmp2,
	_HELPER_Z6AZC_Xdep_mi1,
	_HELPER_Z6AZC_Xdep_ma1,
	_HELPER_Z6AZC_Xdep_mi3,
	_HELPER_Z6AZC_Xdep_ma3,
	_HELPER_Z6AZC_Xdep_mi4,
	_HELPER_Z6AZC_Xdep_ma4,
	_HELPER_Z6AZC_Xdep_20_mul_ma1,
	_HELPER_Z6AZC_Xdep_20_mul_mi1,
	_HELPER_Z6AZC_Xdep_C0im_plus_mi6,
	_HELPER_Z6AZC_Xdep_C1im_plus_ma6,

	_HELPER_Z6AZC_Xdep_ANZ
};

enum {
	_HELPER_Z6AZC_Ydep_mi5=0,
	_HELPER_Z6AZC_Ydep_ma5,
	_HELPER_Z6AZC_Ydep_y03,
	_HELPER_Z6AZC_Ydep_y13,
	_HELPER_Z6AZC_Ydep_mi2,
	_HELPER_Z6AZC_Ydep_ma2,
	_HELPER_Z6AZC_Ydep_mi7,
	_HELPER_Z6AZC_Ydep_ma7,
	_HELPER_Z6AZC_Ydep_15_mul_mi2,
	_HELPER_Z6AZC_Ydep_15_mul_ma2,
	_HELPER_Z6AZC_Ydep_6_mul_mi2,
	_HELPER_Z6AZC_Ydep_6_mul_ma2,
	_HELPER_Z6AZC_Ydep_C0re_minus_ma5_minus_ma7,
	_HELPER_Z6AZC_Ydep_C1re_minus_mi5_minus_mi7,

	_HELPER_Z6AZC_Ydep_ANZ
};

enum {
	_HELPER_Z5AZC_Xdep_mi1=0,
	_HELPER_Z5AZC_Xdep_ma1,
	_HELPER_Z5AZC_Xdep_mi3,
	_HELPER_Z5AZC_Xdep_ma3,
	_HELPER_Z5AZC_Xdep_mi6,
	_HELPER_Z5AZC_Xdep_ma6,
	_HELPER_Z5AZC_Xdep_10_mul_ma1,
	_HELPER_Z5AZC_Xdep_10_mul_mi1,
	_HELPER_Z5AZC_Xdep_tmp5,
	_HELPER_Z5AZC_Xdep_tmp6,
	_HELPER_Z5AZC_Xdep_C1im_plus_ma6,
	_HELPER_Z5AZC_Xdep_C0im_plus_mi6,
	_HELPER_Z5AZC_Xdep_x03,
	_HELPER_Z5AZC_Xdep_x13,

	_HELPER_Z5AZC_Xdep_ANZ
};

enum {
	_HELPER_Z5AZC_Ydep_mi2=0,
	_HELPER_Z5AZC_Ydep_ma2,
	_HELPER_Z5AZC_Ydep_mi4,
	_HELPER_Z5AZC_Ydep_ma4,
	_HELPER_Z5AZC_Ydep_mi5,
	_HELPER_Z5AZC_Ydep_ma5,
	_HELPER_Z5AZC_Ydep_tmp1,
	_HELPER_Z5AZC_Ydep_tmp2,
	_HELPER_Z5AZC_Ydep_10_mul_ma2,
	_HELPER_Z5AZC_Ydep_10_mul_mi2,
	_HELPER_Z5AZC_Ydep_C0re_minus_ma5,
	_HELPER_Z5AZC_Ydep_C1re_minus_mi5,
	_HELPER_Z5AZC_Ydep_y03,
	_HELPER_Z5AZC_Ydep_y13,

	_HELPER_Z5AZC_Ydep_ANZ
};

enum {
	_HELPER_Z4AZC_Xdep_mi1=0,
	_HELPER_Z4AZC_Xdep_ma1,
	_HELPER_Z4AZC_Xdep_mi3,
	_HELPER_Z4AZC_Xdep_ma3,
	_HELPER_Z4AZC_Xdep_C0re_plus_mi3,
	_HELPER_Z4AZC_Xdep_C1re_plus_ma3,

	_HELPER_Z4AZC_Xdep_ANZ
};

enum {
	_HELPER_Z4AZC_Ydep_mi2=0,
	_HELPER_Z4AZC_Ydep_ma2,
	_HELPER_Z4AZC_Ydep_6mi2,
	_HELPER_Z4AZC_Ydep_6ma2,
	_HELPER_Z4AZC_Ydep_mi4,
	_HELPER_Z4AZC_Ydep_ma4,
	_HELPER_Z4AZC_Ydep_mi5,
	_HELPER_Z4AZC_Ydep_ma5,
	_HELPER_Z4AZC_Ydep_4y0,
	_HELPER_Z4AZC_Ydep_4y1,
	_HELPER_Z4AZC_Ydep_C0im_plus_mi4,
	_HELPER_Z4AZC_Ydep_C1im_plus_ma4,
	
	_HELPER_Z4AZC_Ydep_ANZ
};

enum {
	_HELPER_2ITZ2C_Xdep_6mi3=0,
	_HELPER_2ITZ2C_Xdep_6ma3,
	_HELPER_2ITZ2C_Xdep_mi5,
	_HELPER_2ITZ2C_Xdep_ma5,
	_HELPER_2ITZ2C_Xdep_mi8,
	_HELPER_2ITZ2C_Xdep_ma8,
	_HELPER_2ITZ2C_Xdep_2mi9,
	_HELPER_2ITZ2C_Xdep_2ma9,
	_HELPER_2ITZ2C_Xdep_mi13,
	_HELPER_2ITZ2C_Xdep_ma13,
	_HELPER_2ITZ2C_Xdep_2mi15,
	_HELPER_2ITZ2C_Xdep_2ma15,
	_HELPER_2ITZ2C_Xdep_2mi16,
	_HELPER_2ITZ2C_Xdep_2ma16,
	_HELPER_2ITZ2C_Xdep_x03,
	_HELPER_2ITZ2C_Xdep_x13,
	_HELPER_2ITZ2C_Xdep_mi8_plus_2mi9_plus_C0RE_plus_mi5_minus_ma13,
	_HELPER_2ITZ2C_Xdep_ma8_plus_2ma9_plus_C1RE_plus_ma5_minus_mi13,
	_HELPER_2ITZ2C_Xdep_ANZ
};

enum {
	_HELPER_2ITZ2C_Ydep_y03=0,
	_HELPER_2ITZ2C_Ydep_y13,
	_HELPER_2ITZ2C_Ydep_mi4,
	_HELPER_2ITZ2C_Ydep_ma4,
	_HELPER_2ITZ2C_Ydep_mi6,
	_HELPER_2ITZ2C_Ydep_ma6,
	_HELPER_2ITZ2C_Ydep_2mi10,
	_HELPER_2ITZ2C_Ydep_2ma10,
	_HELPER_2ITZ2C_Ydep_2mi17,
	_HELPER_2ITZ2C_Ydep_2ma17,

	_HELPER_2ITZ2C_Ydep_ANZ
};

enum {
	_HELPER_Z3AZC_Xdep_x02=0,
	_HELPER_Z3AZC_Xdep_x12,
	_HELPER_Z3AZC_Xdep_mi1,
	_HELPER_Z3AZC_Xdep_ma1,
	_HELPER_Z3AZC_Xdep_mi4,
	_HELPER_Z3AZC_Xdep_ma4,
	_HELPER_Z3AZC_Xdep_3_mul_mi1,
	_HELPER_Z3AZC_Xdep_3_mul_ma1,
	_HELPER_Z3AZC_Xdep_seedC0im_plus_mi4,
	_HELPER_Z3AZC_Xdep_seedC1im_plus_ma4,
	
	_HELPER_Z3AZC_Xdep_ANZ
};

enum {
	_HELPER_Z3AZC_Ydep_y02=0,
	_HELPER_Z3AZC_Ydep_y12,
	_HELPER_Z3AZC_Ydep_mi2,
	_HELPER_Z3AZC_Ydep_ma2,
	_HELPER_Z3AZC_Ydep_mi3,
	_HELPER_Z3AZC_Ydep_ma3,
	_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_ma2,
	_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_mi2,
	_HELPER_Z3AZC_Ydep_FAKTORAre_minus_ma2,
	_HELPER_Z3AZC_Ydep_FAKTORAre_minus_mi2,
	_HELPER_Z3AZC_Ydep_seedC0re_minus_ma3,
	_HELPER_Z3AZC_Ydep_seedC1re_minus_mi3,
			
	_HELPER_Z3AZC_Ydep_ANZ
};

enum {
	_HELPER_Z2C_Xdep_mi1=0,
	_HELPER_Z2C_Xdep_ma1,
	
	_HELPER_Z2C_Xdep_ANZ
};

enum {
	_HELPER_Z2C_Ydep_C0re_minus_ma2=0,
	_HELPER_Z2C_Ydep_C1re_minus_mi2,
	
	_HELPER_Z2C_Ydep_ANZ
};


// globals

int8_t SAVEIMAGE=1;
int64_t ctrbbxfa=0;
int8_t _RESETPOTW=0;
int8_t _PRECOMPUTEBBXMEMORYGB=0;
ByteManager vgridmgr;
int32_t _SHRINKAGEUSABLE=0;
int32_t _VIRTUALGRIDBITS=0;
int64_t _VIRTUALGRIDMEMORY=0;
// pointer variables need be always declared
HelperAccess *helperYdep=NULL;
HelperAccess *helperXdep=NULL;
HelperAccess_double *helperYdep_double=NULL;
HelperAccess_double *helperXdep_double=NULL;
HelperManager *helpermgr;
Helper_doubleManager *helper_doublemgr;
int8_t interiorpresent=0;
int64_t checkclockatbbxcount0=10000000;
int64_t checkclockatbbxadd=(1 << 26);
int8_t _PERIODICITYMETHOD=1;
int32_t CLOCKHOURSTOSAVE=CLOCKS_PER_SEC*3600*2;
int8_t _PERIODICPOINTS=0;
int8_t _PROPAGATEDEF=1;
int8_t _PROPAGATEPOTW=1;
FILE *flog=NULL;
Cycle* cycles=NULL;
FatouComponent* ibfcomponents=NULL;
int32_t anzibf=0;
int32_t anzcycles=0;
ColorPalette basinpal;
void (*getBoundingBoxfA)(PlaneRect&,PlaneRect&) = NULL;
void (*getBoundingBoxfA_double)(PlaneRect_double&,PlaneRect_double&,Helper_double*,Helper_double*) = NULL;
void (*getBoundingBoxfA_double_oh)(PlaneRect_double&,PlaneRect_double&) = NULL;
void (*getBoundingBoxfA_helper)(PlaneRect&,PlaneRect&,Helper*,Helper*) = NULL;
void (*precompute_helperYdep)(PlaneRect&,Helper*) = NULL;
void (*precompute_helperXdep)(PlaneRect&,Helper*) = NULL;
void (*precompute_helperYdep_double)(PlaneRect_double&,Helper_double*) = NULL;
void (*precompute_helperXdep_double)(PlaneRect_double&,Helper_double*) = NULL;
int32_t _FUNC;
ParentManager* parentmgr=NULL;
int32_t REFINEMENTLEVEL=0;
NTYP seedC0re,seedC1re,seedC0im,seedC1im; 
NTYP FAKTORAre,FAKTORAim;
Data5 *data5;
NTYP scaleRangePerPixel,scalePixelPerRange;

double seedC0re_double,seedC1re_double,seedC0im_double,seedC1im_double; 
double FAKTORAre_double,FAKTORAim_double;
double FAKTORBre_double,FAKTORBim_double;
double scaleRangePerPixel_double,scalePixelPerRange_double;

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
double RANGE0=-2,RANGE1=2; // no longer int, as ranges 0.5 etc shall be possible
NTYP COMPLETE0,COMPLETE1;
double COMPLETE0_double,COMPLETE1_double;
// for the low-resolution reverse cell graph working on (usually) 64x64 pixel squares or bigger
int32_t REVCGBITS,REVCGBLOCKWIDTH;
int32_t REVCGmaxnumber,REVCGmaxnumberQ;
// variables for bit precision checks


// forward declarations

void compute(void);
void freeRevCGMem(void);
void construct_static_reverse_cellgraph(void);
void find_special_exterior_hitting_squares(void);
void propagate_definite(void);
void propagate_potw(void);
int32_t color_changeS32(const DDBYTE,const DDBYTE,const DDBYTE,const DDBYTE);
void copy_pixel_to_2x2grid(const uint32_t,uint32_t*);

static inline int32_t scrcoord_as_lowerleft(const NTYP&);
// used for fastdtcheck
static inline int32_t scrcoord_as_lowerleft_double(const double&);

void write2(FILE*,const uint8_t,const uint8_t);
void write4(FILE*,const uint8_t,const uint8_t,const uint8_t,const uint8_t);
static inline int32_t minimumI(const int32_t,const int32_t);
static inline int32_t maximumI(const int32_t,const int32_t);
static inline int32_t minimumI(const int32_t,const int32_t,const int32_t,const int32_t);

inline NTYP minimumD(const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP);
inline NTYP minimumD(const NTYP,const NTYP,const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP,const NTYP,const NTYP);

inline double minimumdouble(const double,const double);
inline double maximumdouble(const double,const double);
inline double minimumdouble(const double,const double,const double,const double);
inline double maximumdouble(const double,const double,const double,const double);

// defines used as small expressions

// AA,BB,CC,DD are VARIABLES, not expressions
// AA..BB, CC..DD may contain zero in the interior
// faster cases: AA<=BB<=0, BB>=AA>=0 and CC,DD accordingly
// then IA multiplicaiton can be performed with 2 mul-op
// instead of all 4
#define IAMUL_SIGN(ERGMIN,ERGMAX,AA,BB,CC,DD,FUNC) \
{\
	if (((AA)>=0.0)&&((AA)<=(BB))) {\
		if (((CC)>=0.0)&&((CC)<=(DD))) {\
			ERGMIN=(AA)*(CC);\
			ERGMAX=(BB)*(DD);\
		} else if (((CC)<=(DD))&&((DD)<=0.0)) {	\
			ERGMIN=(BB)*(CC);\
			ERGMAX=(AA)*(DD);\
		} else {\
			FUNC(ERGMIN,ERGMAX,(AA)*(CC),(AA)*(DD),(BB)*(CC),(BB)*(DD));\
		}\
	} \
	else if (((AA)<=(BB))&&((BB)<=0.0)) {\
		if ((0.0<=(CC))&&((CC)<=(DD))) {\
			ERGMIN=(AA)*(DD);\
			ERGMAX=(BB)*(CC);\
		} else if (((CC)<=(DD))&&((DD)<=0.0)) {\
			ERGMIN=(BB)*(DD);\
			ERGMAX=(AA)*(CC);\
		} else {\
			FUNC(ERGMIN,ERGMAX,(AA)*(CC),(AA)*(DD),(BB)*(CC),(BB)*(DD));\
		}\
	} else {\
		FUNC(ERGMIN,ERGMAX,(AA)*(CC),(AA)*(DD),(BB)*(CC),(BB)*(DD));\
	}\
}

// known: AA,BB,CC,DD are VARIABLES, not expressions
#define IAMULFPASIGNALL4(ERGMIN,ERGMAX,AA,BB,CC,DD,FUNC,ARRAY) \
{\
	FPA_mul_ZAB(ARRAY[0],&AA,&CC);\
	FPA_mul_ZAB(ARRAY[1],&AA,&DD);\
	FPA_mul_ZAB(ARRAY[2],&BB,&CC);\
	FPA_mul_ZAB(ARRAY[3],&BB,&DD);\
	FUNC(ERGMIN,ERGMAX,ARRAY);\
}

#define IAMUL_FPA_SIGN(ERGMIN,ERGMAX,AA,BB,CC,DD,FUNC,ARRAY) \
{\
	int vglAABB=FPA_vgl(&AA,&BB);\
	int vglCCDD=FPA_vgl(&CC,&DD);\
	if (((AA).vorz>=0)&&(vglAABB<=0)) {\
		if (((CC).vorz>=0)&&(vglCCDD<=0)) {\
			FPA_mul_ZAB(ARRAY[0],&AA,&CC);\
			FPA_mul_ZAB(ARRAY[1],&BB,&DD);\
			ERGMIN=&ARRAY[0];\
			ERGMAX=&ARRAY[1];\
		} else if ((vglCCDD<=0)&&((DD).vorz<=0)) {	\
			FPA_mul_ZAB(ARRAY[0],&BB,&CC);\
			FPA_mul_ZAB(ARRAY[1],&AA,&DD);\
			ERGMIN=&ARRAY[0];\
			ERGMAX=&ARRAY[1];\
		} else {\
			IAMULFPASIGNALL4(ERGMIN,ERGMAX,AA,BB,CC,DD,FUNC,ARRAY)\
		}\
	} \
	else if ((vglAABB<=0)&&((BB).vorz<=0)) {\
		if (((CC).vorz>=0)&&(vglCCDD<=0)) {\
			FPA_mul_ZAB(ARRAY[0],&AA,&DD);\
			FPA_mul_ZAB(ARRAY[1],&BB,&CC);\
			ERGMIN=&ARRAY[0];\
			ERGMAX=&ARRAY[1];\
		} else if ((vglCCDD<=0)&&((DD).vorz<=0)) {\
			FPA_mul_ZAB(ARRAY[0],&BB,&DD);\
			FPA_mul_ZAB(ARRAY[1],&AA,&CC);\
			ERGMIN=&ARRAY[0];\
			ERGMAX=&ARRAY[1];\
		} else {\
			IAMULFPASIGNALL4(ERGMIN,ERGMAX,AA,BB,CC,DD,FUNC,ARRAY)\
		}\
	} else {\
		IAMULFPASIGNALL4(ERGMIN,ERGMAX,AA,BB,CC,DD,FUNC,ARRAY)\
	}\
}

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

#define SETPCSCR(XX,YY,WX0,WX1,WY0,WY1) \
{\
	const int32_t xpos=(XX)-data5->memgrau[YY].g0;\
	data5->pcscr[YY][xpos].x0=WX0;\
	data5->pcscr[YY][xpos].x1=WX1;\
	data5->pcscr[YY][xpos].y0=WY0;\
	data5->pcscr[YY][xpos].y1=WY1;\
}

#define GETPCSCR(XX,YY,SCR) \
{\
	const int32_t xpos=(XX)-data5->memgrau[YY].g0;\
	SCR.x0=data5->pcscr[YY][xpos].x0;\
	SCR.x1=data5->pcscr[YY][xpos].x1;\
	if (SCR.x1 == -2) {\
		LOGMSG3("Error/getpcscr %i,%i\n",XX,YY);\
		exit(99);\
	}\
	SCR.y0=data5->pcscr[YY][xpos].y0;\
	SCR.y1=data5->pcscr[YY][xpos].y1;\
}

#define SETTOVISIT(XX,YY) \
{\
	data5->revcgYX[\
		(YY)*REVCGmaxnumber+\
		(XX)\
	].tovisit=1;\
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

#define SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(WW,CLEARMASKE,SETF) \
	((WW & CLEARMASKE) | SETF)


// routines

char* upper(char* s) {
	if (!s) return NULL;
	
	for(int32_t i=0;i<(int32_t)strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

void ausgabePlaneRect(FILE* f,PlaneRect& A) {
	#ifdef _FPA
	char tt[2048];
	fprintf(f,"(%s..",A.x0.str(tt));
	fprintf(f,"%s,",A.x1.str(tt));
	fprintf(f,"%s..",A.y0.str(tt));
	fprintf(f,"%s)",A.y1.str(tt));
	#else
	fprintf(f,"(%.20lg..%.20lg,%.20lg..%.20lg)",
		(double)A.x0,(double)A.x1,(double)A.y0,(double)A.y1);
	#endif
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

inline int32_t minimumI(const int32_t a,const int32_t b,const int32_t c,const int32_t d) {
	int32_t m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	return m;
}

// static reverse cell graph

RevCGBlock::RevCGBlock() {
	howmany=0;
	memused=0;
	containsgray=0;
}

void RevCGBlock::addParent(const int32_t ax,const int32_t ay) {
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

void Data5::saveRaw(const char* afn) {
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
			if ((w & 0b11) == SQUARE_BLACK) {
				//printf("read interior x=%i\n",((i<<4)+bit));
				return 1;
			}
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
	encgrayx0=encgrayy0=SCREENWIDTH-1;
	encgrayx1=encgrayy1=0;
	DDBYTE *eine=new DDBYTE[SCREENWIDTH];
	int64_t memused=0;
	
	printf("reading stored data ");

	if ((int32_t)savedlen == SCREENWIDTH) {
		// keep resolution
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
				graudensity[y]=0;
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
				int64_t ctrgrau=0;
				for(int32_t dx=0;dx<laenge;dx++) {
					if (
						(zeilen[y][dx] != SQUARE_WHITE_16_CONSECUTIVE) &&
						(zeilen[y][dx] != SQUARE_BLACK_16_CONSECUTIVE)
					) ctrgrau++;
				}
				// estimate for how many gray cells are in this row
				// needed for guided precomputing bbx's
				graudensity[y]=(int)(100*(double)ctrgrau/(double)laenge);
				memgrau[y].mem0=start;
				memgrau[y].mem1=start+laenge-1;
				memgrau[y].g0=memgrau[y].mem0 << 4;
				memgrau[y].g1=((memgrau[y].mem1+1) << 4)-1;
				if (y < encgrayy0) encgrayy0=y;
				if (y > encgrayy1) encgrayy1=y;
				if (memgrau[y].g0 < encgrayx0) encgrayx0=memgrau[y].g0;
				if (memgrau[y].g1 > encgrayx1) encgrayx1=memgrau[y].g1;
			}
		} // for y
	} else if ( (int32_t)savedlen == (SCREENWIDTH >> 1) ) {
		int32_t start,laenge;
		for(int32_t yread=0;yread<(SCREENWIDTH-1);yread+=2) {
			fread(&start,1,sizeof(start),f);
			fread(&laenge,1,sizeof(laenge),f);
			if (laenge<=0) {
				zeilen[yread]=zeilen[yread+1]=NULL;
				graudensity[yread]=0;
				graudensity[yread+1]=0;
				memgrau[yread].g0=memgrau[yread+1].g0=SCREENWIDTH;
				memgrau[yread].g1=memgrau[yread+1].g1=0;
				memgrau[yread].mem0=memgrau[yread+1].mem0=SCREENWIDTH >> 4;
				memgrau[yread].mem1=memgrau[yread+1].mem1=0;
			} else {
				// one cell into a 2x2 grid - refinement process
				int32_t readlaenge=laenge;
				laenge <<= 1;
				start <<= 1;
				zeilen[yread]=data5->datamgr->getMemory(laenge);
				zeilen[yread+1]=data5->datamgr->getMemory(laenge);
				memused += (2*laenge*sizeof(DDBYTE));
				if (!zeilen[yread+1]) {
					LOGMSG("Memory error. ReadRaw\n");
					exit(99);
				}
				fread(eine,readlaenge,sizeof(int32_t),f);
				if (interiorpresent<=0) {
					interiorpresent=interiorinrow(eine,readlaenge);
				}
				
				int64_t ctrgrau=0;
				for(int32_t dx=0;dx<readlaenge;dx++) {
					if (
						(eine[dx] != SQUARE_WHITE_16_CONSECUTIVE) &&
						(eine[dx] != SQUARE_BLACK_16_CONSECUTIVE)
					) ctrgrau++;
				}
				// estimate for how many gray cells are in this row
				int grd=(int)(100*(double)ctrgrau/(double)readlaenge);
				graudensity[yread]=grd;
				graudensity[yread+1]=grd;

				DDBYTE mem=0;
				for(int32_t k=0;k<readlaenge;k++) {
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
	
	printf("\n  %I64d GB cell memory allocated\n",1+(memused >> 30));
	
	// outward rounded enclosement of all gray
	planegrayx0=encgrayx0*scaleRangePerPixel + COMPLETE0;
	planegrayy0=encgrayy0*scaleRangePerPixel + COMPLETE0;
	planegrayx1=(encgrayx1+16)*scaleRangePerPixel + COMPLETE0;
	planegrayy1=(encgrayy1+16)*scaleRangePerPixel + COMPLETE0;
	
	if (_RESETPOTW>0) {
		// reset all potw to gray so a new propagate_XXX pass can be started fresh
		printf("resetting potw to gray ... ");
		color_changeS32(
			SQUARE_GRAY_POTENTIALLY_WHITE,
			SQUARE_GRAY,
			SQUARE_GRAYPOTW_16_CONSECUTIVE,
			SQUARE_GRAY_16_CONSECUTIVE
		);
		printf("\n");
	}
	
	return 1;
}

void Data5::saveBitmap4_twd(const char* afn,const int atwdexp) {
	// saves a trustworthily downsized version of the image: 16-fold. 
	// image format is: 8 bit Bitmap
	// 2^TWDEXP x 2^TWDEXP gives one final pixel
	
	int32_t _TWDEXPONENT=atwdexp;

	if (atwdexp<0) {
		// adjust exponent, so final image is at most 2^16 x 2^16
		_TWDEXPONENT=0;
		while ( (SCREENWIDTH >> _TWDEXPONENT) > 65536) _TWDEXPONENT++;
	}

	int32_t TWDSTEP=(1 << _TWDEXPONENT);
	char tmp[1024];
	RGB4 pal[256];
	int32_t bytes_per_row = SCREENWIDTH >> _TWDEXPONENT;
	uint8_t* rgbz=new uint8_t[bytes_per_row];
	uint32_t off
		=	14 // size of file header
		+	40 // size of bitmap header
		+	256*4; // palette entries
	uint32_t filelen
		=	off
		+	(bytes_per_row*bytes_per_row);
	
	for(int32_t i=0;i<256;i++) pal[i].R=pal[i].G=pal[i].B=pal[i].alpha=63;
	pal[SQUARE_GRAY].R=127;
	pal[SQUARE_GRAY].G=127;
	pal[SQUARE_GRAY].B=127;
	pal[SQUARE_GRAY_POTENTIALLY_WHITE].R=255;
	pal[SQUARE_GRAY_POTENTIALLY_WHITE].G=0;
	pal[SQUARE_GRAY_POTENTIALLY_WHITE].B=0;
	pal[SQUARE_BLACK].R=0;
	pal[SQUARE_BLACK].G=0;
	pal[SQUARE_BLACK].B=0;
	pal[SQUARE_WHITE].R=255;
	pal[SQUARE_WHITE].G=255;
	pal[SQUARE_WHITE].B=255;

	sprintf(tmp,"%s_2_%i-fold.bmp",afn,_TWDEXPONENT);
	FILE *fbmp=fopen(tmp,"wb");
	if (!fbmp) {
		LOGMSG("Abort error. saveBitmap4_twd.\n");
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
	
	for(int32_t y=0;y<SCREENWIDTH;y+=TWDSTEP) {
		int32_t xzeile=-1;
		for(int32_t x=0;x<SCREENWIDTH;x+=TWDSTEP) {
			xzeile++;
			
			int32_t finalf=-1;
			for(int32_t dy=0;dy<TWDSTEP;dy++) {
				for(int32_t dx=0;dx<TWDSTEP;dx++) {
					int32_t f;
					GET_SINGLE_CELLCOLOR_XY(x+dx,y+dy,f);
					if (f == SQUARE_GRAY_POTENTIALLY_WHITE) {
						finalf=SQUARE_GRAY;
						break;
					}
					
					if (finalf<0) finalf=f;
					else if (finalf != f) {
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

Data5::Data5() {
	printf("initialising main object ...\n");
	
	vgridYX=NULL;
	memgrau=new Gray_in_row[SCREENWIDTH];
	zeilen=new uint32_t*[SCREENWIDTH];
	revcgYX=new RevCGBlock[REVCGmaxnumber*REVCGmaxnumber];
	datamgr=new ArrayDDByteManager;
	graudensity=new uint8_t[SCREENWIDTH];
	for(int32_t i=0;i<SCREENWIDTH;i++) graudensity[i]=100;
	pcscr=NULL; 
	pcscrmgr=NULL;
}

Data5::~Data5() {
	delete datamgr;
	if (vgridYX) delete[] vgridYX;
	if (revcgYX) delete[] revcgYX;
	delete[] memgrau;
	delete[] zeilen;
	if (pcscr) delete[] pcscr;
	if (pcscrmgr) delete pcscrmgr;
}

// one pixel transforming into a 2x2 grid, gray or gray-potentially-white will both
// be set to gray. 
// works here on a 32bit-integer, hence 16 consecutive bits directly
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
	} // qbit
}

static inline int32_t scrcoord_as_lowerleft(const NTYP& a) {
	// calculating the screen coordinte of the pixel that contains the coordinate
	// if the coordinate lies on an edge/corner (and belongs to more than one pixel)
	// the pixel where it lies on the left,bottom edge/corner is returned
	// check for out-of-COMPLETE before flooring values
	// to remain in the representable int32_t range
	// truncates to the screen, so if it is outside is
	// checked in the calling routine

	int32_t w;
	
	#define NORMALFLOOR \
	{\
		if (a <= COMPLETE0) return 0;\
		if (a >= COMPLETE1) return (SCREENWIDTH-1);\
		w=(int)floor( (a - COMPLETE0) * scalePixelPerRange );\
		if (w >= SCREENWIDTH) return (SCREENWIDTH-1);\
		return w;\
	}
	
	#ifdef _QUADMATH
	if (a <= COMPLETE0) return 0;
	if (a >= COMPLETE1) return (SCREENWIDTH-1);
	// special floor function: floorq
	w=(int)floorq( (a - COMPLETE0) * scalePixelPerRange );
	if (w >= SCREENWIDTH) return (SCREENWIDTH-1);
	return w;
	#endif

	#ifdef _F161
	if (a <= COMPLETE0) return 0;
	if (a >= COMPLETE1) return (SCREENWIDTH-1);
	// special floor function
	w=(int)specialFloor161( (a - COMPLETE0) * scalePixelPerRange );
	if (w >= SCREENWIDTH) return (SCREENWIDTH-1);
	return w;
	#endif
	
	#ifdef _FPA
	if (FPA_vgl(a,COMPLETE0) <= 0) return 0;
	if (FPA_vgl(a,COMPLETE1) >= 0) return (SCREENWIDTH-1);
	FPA b=(a-COMPLETE0);
	// scalePixelPerRange is a power of 2
	// so bit-shifting can be used as multiplication
	//b=b*scalePixelPerRange;
	b.shiftLeft(scalePixelPerRangeExponent);
	int64_t wfl=floorFPA(b);

	if (
		(wfl < -(int64_t)UINT32MAX) ||
		(wfl > (int64_t)UINT32MAX)
	) {
		printf("Implementation error. FPA. out-of-range floor scrcoord %I64d\n",wfl);
		exit(99);
	}
	if (wfl >= SCREENWIDTH) return (SCREENWIDTH-1);
	
	return (int32_t)wfl;
	#endif
	
	#ifdef _DOUBLE
	NORMALFLOOR
	#endif
	
	#ifdef _LONGDOUBLE
	NORMALFLOOR
	#endif

	#ifdef _F107
	NORMALFLOOR
	#endif
	
	LOGMSG("Implementation error. scrcoord_as_lower_left.\n");
	exit(99);
}

static inline int32_t scrcoord_as_lowerleft_double(const double& a) {
	int32_t w;
	
	if (a <= COMPLETE0_double) return 0;
	if (a >= COMPLETE1_double) return (SCREENWIDTH-1);
	w=(int)floor( (a - COMPLETE0_double) * scalePixelPerRange_double );
	if (w >= SCREENWIDTH) return (SCREENWIDTH-1);
	return w;
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

inline NTYP maximumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}
#endif

inline double minimumdouble(const double a,const double b,const double c,const double d) {
	double m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	return m;
}

inline double minimumdouble(const double a,const double b) {
	if (a < b) return a;
	return b;
}

inline double maximumdouble(const double a,const double b) {
	if (a > b) return a;
	return b;
}

inline double maximumdouble(const double a,const double b,const double c,const double d) {
	double m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}

// function definitions
void precomputeYdep_z2c(
	PlaneRect& A,Helper *ahy
) {
	#define Z2CPRECOMPUTEY(MINMAX2,MINMAX4) \
	{\
		NTYP y02=A.y0*A.y0;\
		NTYP y12=A.y1*A.y1;\
		NTYP mi2,ma2;\
		MINMAX2(mi2,ma2,y02,y12);\
		ahy->val[_HELPER_Z2C_Ydep_C0re_minus_ma2]=seedC0re-ma2;\
		ahy->val[_HELPER_Z2C_Ydep_C1re_minus_mi2]=seedC1re-mi2;\
		return;\
	}\
	
	#ifdef _DOUBLE
	Z2CPRECOMPUTEY(minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z2CPRECOMPUTEY(minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z2CPRECOMPUTEY(minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z2CPRECOMPUTEY(minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z2CPRECOMPUTEY(minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z2CPRECOMPUTEY(minimaxFPAAB,minimaxFPAABCD)
	#endif
}

void precomputeXdep_z2c(
	PlaneRect& A,Helper *ahx
) {
	#define Z2CPRECOMPUTEX(MINMAX2,MINMAX4) \
	{\
		NTYP x02=A.x0*A.x0;\
		NTYP x12=A.x1*A.x1;\
		NTYP mi1,ma1;\
		MINMAX2(mi1,ma1,x02,x12);\
		ahx->val[_HELPER_Z2C_Xdep_mi1]=mi1;\
		ahx->val[_HELPER_Z2C_Xdep_ma1]=ma1;\
		return;\
	}\
	
	#ifdef _DOUBLE
	Z2CPRECOMPUTEX(minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z2CPRECOMPUTEX(minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z2CPRECOMPUTEX(minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z2CPRECOMPUTEX(minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z2CPRECOMPUTEX(minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z2CPRECOMPUTEX(minimaxFPAAB,minimaxFPAABCD)
	#endif
}

void getBoundingBoxfA_z2c(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re+minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1);
	fA.x1=seedC1re+maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1);
	fA.y0=seedC0im+2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
	fA.y1=seedC1im+2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
}

void getBoundingBoxfA_z2c_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper *ahx,Helper* ahy
) {
	ctrbbxfa++;

	#define Z2CHELPER(MINMAX2,MINMAX4) \
	{\
		NTYP mi3,ma3;\
		IAMUL_SIGN(mi3,ma3,A.x0,A.x1,A.y0,A.y1,MINMAX4)\
		fA.x0=ahy->val[_HELPER_Z2C_Ydep_C0re_minus_ma2]\
		+\
		ahx->val[_HELPER_Z2C_Xdep_mi1];\
		fA.x1=ahy->val[_HELPER_Z2C_Ydep_C1re_minus_mi2]\
		+\
		ahx->val[_HELPER_Z2C_Xdep_ma1];\
		fA.y0=seedC0im+2*mi3;\
		fA.y1=seedC1im+2*ma3;\
		\
		return;\
	}

	#ifdef _DOUBLE
	Z2CHELPER(minimaxdAB,minimaxdABCD);
	#endif

	#ifdef _LONGDOUBLE
	Z2CHELPER(minimaxldAB,minimaxldABCD);
	#endif

	#ifdef _F107
	Z2CHELPER(minimaxF107AB,minimaxF107ABCD);
	#endif

	#ifdef _F161
	Z2CHELPER(minimaxF161AB,minimaxF161ABCD);
	#endif

	#ifdef _QUADMATH
	Z2CHELPER(minimaxQDAB,minimaxQDABCD);
	#endif
	
	#ifdef _FPA
	PFPA mi3,ma3;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi3,ma3,
		A.x0,A.x1,A.y0,A.y1,
		minimaxFPAABCD,array1
	)
	
	FPA_add_ZAB(fA.x0,
		&ahy->val[_HELPER_Z2C_Ydep_C0re_minus_ma2],
		&ahx->val[_HELPER_Z2C_Xdep_mi1]
	);
	FPA_add_ZAB(fA.x1,
		&ahy->val[_HELPER_Z2C_Ydep_C1re_minus_mi2],
		&ahx->val[_HELPER_Z2C_Xdep_ma1]
	);
	
	// multiply by 2
	FPA dmi3,dma3;
	FPA_add_ZAB(dmi3,mi3,mi3);
	FPA_add_ZAB(dma3,ma3,ma3);
	FPA_add_ZAB(fA.y0,&seedC0im,&dmi3);
	FPA_add_ZAB(fA.y1,&seedC1im,&dma3);
	
	return;
	#endif
}

void getBoundingBoxfA_2itz2c(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=minimumD(seedC0re*seedC0re,seedC1re*seedC1re)+2*minimumD(seedC0re*minimumD(A.x0*A.x0,A.x1*A.x1),seedC0re*maximumD(A.x0*A.x0,A.x1*A.x1),seedC1re*minimumD(A.x0*A.x0,A.x1*A.x1),seedC1re*maximumD(A.x0*A.x0,A.x1*A.x1))-(2*maximumD(seedC0re*minimumD(A.y0*A.y0,A.y1*A.y1),seedC0re*maximumD(A.y0*A.y0,A.y1*A.y1),seedC1re*minimumD(A.y0*A.y0,A.y1*A.y1),seedC1re*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC0re+minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))-(4*maximumD(seedC0im*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0im*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)))+minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-maximumD(seedC0im*seedC0im,seedC1im*seedC1im);
	fA.x1=maximumD(seedC0re*seedC0re,seedC1re*seedC1re)+2*maximumD(seedC0re*minimumD(A.x0*A.x0,A.x1*A.x1),seedC0re*maximumD(A.x0*A.x0,A.x1*A.x1),seedC1re*minimumD(A.x0*A.x0,A.x1*A.x1),seedC1re*maximumD(A.x0*A.x0,A.x1*A.x1))-(2*minimumD(seedC0re*minimumD(A.y0*A.y0,A.y1*A.y1),seedC0re*maximumD(A.y0*A.y0,A.y1*A.y1),seedC1re*minimumD(A.y0*A.y0,A.y1*A.y1),seedC1re*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC1re+maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))-(4*minimumD(seedC0im*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0im*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)))+maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-minimumD(seedC0im*seedC0im,seedC1im*seedC1im);
	fA.y0=4*minimumD(seedC0re*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0re*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1))+2*minimumD(seedC0im*seedC0re,seedC0im*seedC1re,seedC1im*seedC0re,seedC1im*seedC1re)+4*minimumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)+2*minimumD(seedC0im*minimumD(A.x0*A.x0,A.x1*A.x1),seedC0im*maximumD(A.x0*A.x0,A.x1*A.x1),seedC1im*minimumD(A.x0*A.x0,A.x1*A.x1),seedC1im*maximumD(A.x0*A.x0,A.x1*A.x1))-(4*maximumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))-(2*maximumD(seedC0im*minimumD(A.y0*A.y0,A.y1*A.y1),seedC0im*maximumD(A.y0*A.y0,A.y1*A.y1),seedC1im*minimumD(A.y0*A.y0,A.y1*A.y1),seedC1im*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC0im;
	fA.y1=4*maximumD(seedC0re*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0re*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1))+2*maximumD(seedC0im*seedC0re,seedC0im*seedC1re,seedC1im*seedC0re,seedC1im*seedC1re)+4*maximumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)+2*maximumD(seedC0im*minimumD(A.x0*A.x0,A.x1*A.x1),seedC0im*maximumD(A.x0*A.x0,A.x1*A.x1),seedC1im*minimumD(A.x0*A.x0,A.x1*A.x1),seedC1im*maximumD(A.x0*A.x0,A.x1*A.x1))-(4*minimumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))-(2*minimumD(seedC0im*minimumD(A.y0*A.y0,A.y1*A.y1),seedC0im*maximumD(A.y0*A.y0,A.y1*A.y1),seedC1im*minimumD(A.y0*A.y0,A.y1*A.y1),seedC1im*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC1im;
}

#define _2ITZ2CHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP mi1,ma1;\
	IAMUL_SIGN(mi1,ma1,A.x0,A.x1,\
		ahy->val[_HELPER_2ITZ2C_Ydep_y03],\
		ahy->val[_HELPER_2ITZ2C_Ydep_y13],\
		MINMAX4);\
	NUMTYP mi2,ma2;\
	IAMUL_SIGN(mi2,ma2,\
		ahx->val[_HELPER_2ITZ2C_Xdep_x03],\
		ahx->val[_HELPER_2ITZ2C_Xdep_x13],\
		A.y0,A.y1,MINMAX4);\
	NUMTYP mi7,ma7;\
	IAMUL_SIGN(mi7,ma7,\
		ahx->val[_HELPER_2ITZ2C_Xdep_6mi3],\
		ahx->val[_HELPER_2ITZ2C_Xdep_6ma3],\
		ahy->val[_HELPER_2ITZ2C_Ydep_mi4],\
		ahy->val[_HELPER_2ITZ2C_Ydep_ma4],\
		MINMAX4);\
	NUMTYP mi11,ma11;\
	IAMUL_SIGN(mi11,ma11,A.x0,A.x1,A.y0,A.y1,MINMAX4);\
	NUMTYP mi12,ma12;\
	IAMUL_SIGN(mi12,ma12,C0IM,C1IM,mi11,ma11,MINMAX4);\
	NUMTYP mi14,ma14;\
	IAMUL_SIGN(mi14,ma14,C0RE,C1RE,mi11,ma11,MINMAX4);\
	\
	fA.x0=\
		ahx->val[_HELPER_2ITZ2C_Xdep_mi8_plus_2mi9_plus_C0RE_plus_mi5_minus_ma13]\
		-ahy->val[_HELPER_2ITZ2C_Ydep_2ma10]\
		+ahy->val[_HELPER_2ITZ2C_Ydep_mi6]\
		-ma7\
		-4*ma12;\
	fA.x1=\
		ahx->val[_HELPER_2ITZ2C_Xdep_ma8_plus_2ma9_plus_C1RE_plus_ma5_minus_mi13]\
		-ahy->val[_HELPER_2ITZ2C_Ydep_2mi10]\
		+ahy->val[_HELPER_2ITZ2C_Ydep_ma6]\
		-mi7\
		-4*mi12;\
	\
	fA.y0=\
		(((((((4*mi14\
		+ahx->val[_HELPER_2ITZ2C_Xdep_2mi15])\
		+(4*mi2))\
		+ahx->val[_HELPER_2ITZ2C_Xdep_2mi16])\
		-(4*ma1))\
		-ahy->val[_HELPER_2ITZ2C_Ydep_2ma17]))\
		+C0IM);\
	fA.y1=\
		(((((((4*ma14)\
		+ahx->val[_HELPER_2ITZ2C_Xdep_2ma15])\
		+(4*ma2))\
		+ahx->val[_HELPER_2ITZ2C_Xdep_2ma16])\
		-(4*mi1))\
		-ahy->val[_HELPER_2ITZ2C_Ydep_2mi17])\
		+C1IM);\
	\
	return;\
}

#define _2ITZ2CHELPER(NUMTYP,MINMAX2,MINMAX4) \
{\
	_2ITZ2CHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	);\
}

void getBoundingBoxfA_2itz2c_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=minimumdouble(seedC0re_double*seedC0re_double,seedC1re_double*seedC1re_double)+2*minimumdouble(seedC0re_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC0re_double*maximumdouble(A.x0*A.x0,A.x1*A.x1),seedC1re_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC1re_double*maximumdouble(A.x0*A.x0,A.x1*A.x1))-(2*maximumdouble(seedC0re_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC0re_double*maximumdouble(A.y0*A.y0,A.y1*A.y1),seedC1re_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC1re_double*maximumdouble(A.y0*A.y0,A.y1*A.y1)))+seedC0re_double+minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1)))-(4*maximumdouble(seedC0im_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0im_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)))+minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-maximumdouble(seedC0im_double*seedC0im_double,seedC1im_double*seedC1im_double);
	fA.x1=maximumdouble(seedC0re_double*seedC0re_double,seedC1re_double*seedC1re_double)+2*maximumdouble(seedC0re_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC0re_double*maximumdouble(A.x0*A.x0,A.x1*A.x1),seedC1re_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC1re_double*maximumdouble(A.x0*A.x0,A.x1*A.x1))-(2*minimumdouble(seedC0re_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC0re_double*maximumdouble(A.y0*A.y0,A.y1*A.y1),seedC1re_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC1re_double*maximumdouble(A.y0*A.y0,A.y1*A.y1)))+seedC1re_double+maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1)))-(4*minimumdouble(seedC0im_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0im_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1im_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)))+maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-minimumdouble(seedC0im_double*seedC0im_double,seedC1im_double*seedC1im_double);
	fA.y0=4*minimumdouble(seedC0re_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0re_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1))+2*minimumdouble(seedC0im_double*seedC0re_double,seedC0im_double*seedC1re_double,seedC1im_double*seedC0re_double,seedC1im_double*seedC1re_double)+4*minimumdouble((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)+2*minimumdouble(seedC0im_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC0im_double*maximumdouble(A.x0*A.x0,A.x1*A.x1),seedC1im_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC1im_double*maximumdouble(A.x0*A.x0,A.x1*A.x1))-(4*maximumdouble(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))-(2*maximumdouble(seedC0im_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC0im_double*maximumdouble(A.y0*A.y0,A.y1*A.y1),seedC1im_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC1im_double*maximumdouble(A.y0*A.y0,A.y1*A.y1)))+seedC0im_double;
	fA.y1=4*maximumdouble(seedC0re_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC0re_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re_double*minimumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),seedC1re_double*maximumdouble(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1))+2*maximumdouble(seedC0im_double*seedC0re_double,seedC0im_double*seedC1re_double,seedC1im_double*seedC0re_double,seedC1im_double*seedC1re_double)+4*maximumdouble((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)+2*maximumdouble(seedC0im_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC0im_double*maximumdouble(A.x0*A.x0,A.x1*A.x1),seedC1im_double*minimumdouble(A.x0*A.x0,A.x1*A.x1),seedC1im_double*maximumdouble(A.x0*A.x0,A.x1*A.x1))-(4*minimumdouble(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))-(2*minimumdouble(seedC0im_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC0im_double*maximumdouble(A.y0*A.y0,A.y1*A.y1),seedC1im_double*minimumdouble(A.y0*A.y0,A.y1*A.y1),seedC1im_double*maximumdouble(A.y0*A.y0,A.y1*A.y1)))+seedC1im_double;
}

void getBoundingBoxfA_2itz2c_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double *ahx,Helper_double* ahy
) {
	_2ITZ2CHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	);
}

void getBoundingBoxfA_2itz2c_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper *ahx,Helper* ahy
) {
	ctrbbxfa++;
	
	#ifdef _DOUBLE
	_2ITZ2CHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	_2ITZ2CHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	_2ITZ2CHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	_2ITZ2CHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	_2ITZ2CHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	PFPA mi1,ma1;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi1,ma1,A.x0,A.x1,
		ahy->val[_HELPER_2ITZ2C_Ydep_y03],
		ahy->val[_HELPER_2ITZ2C_Ydep_y13],
		minimaxFPAABCD,array1);
	PFPA mi2,ma2;
	FPA array2[4];
	IAMUL_FPA_SIGN(mi2,ma2,
		ahx->val[_HELPER_2ITZ2C_Xdep_x03],
		ahx->val[_HELPER_2ITZ2C_Xdep_x13],
		A.y0,A.y1,
		minimaxFPAABCD,array2);
	PFPA mi7,ma7;
	FPA array3[4];
	IAMUL_FPA_SIGN(mi7,ma7,
		ahx->val[_HELPER_2ITZ2C_Xdep_6mi3],
		ahx->val[_HELPER_2ITZ2C_Xdep_6ma3],
		ahy->val[_HELPER_2ITZ2C_Ydep_mi4],
		ahy->val[_HELPER_2ITZ2C_Ydep_ma4],
		minimaxFPAABCD,array3);
	PFPA mi11,ma11;
	FPA array4[4];
	IAMUL_FPA_SIGN(mi11,ma11,A.x0,A.x1,A.y0,A.y1,minimaxFPAABCD,array4);
	PFPA mi12,ma12;
	FPA array5[4];
	IAMUL_FPA_SIGN(mi12,ma12,seedC0im,seedC1im,*mi11,*ma11,minimaxFPAABCD,array5);
	PFPA mi14,ma14;
	FPA array6[4];
	IAMUL_FPA_SIGN(mi14,ma14,seedC0re,seedC1re,*mi11,*ma11,minimaxFPAABCD,array6);
	
	mi12->shiftLeft(2);
	ma12->shiftLeft(2);
	ma1->shiftLeft(2);
	mi1->shiftLeft(2);
	ma2->shiftLeft(2);
	mi2->shiftLeft(2);
	ma14->shiftLeft(2);
	mi14->shiftLeft(2);

	fA.x0=
		ahx->val[_HELPER_2ITZ2C_Xdep_mi8_plus_2mi9_plus_C0RE_plus_mi5_minus_ma13]
		-ahy->val[_HELPER_2ITZ2C_Ydep_2ma10]
		+ahy->val[_HELPER_2ITZ2C_Ydep_mi6]
		-(*ma7)
		-(*ma12);
		
	fA.x1=
		ahx->val[_HELPER_2ITZ2C_Xdep_ma8_plus_2ma9_plus_C1RE_plus_ma5_minus_mi13]
		-ahy->val[_HELPER_2ITZ2C_Ydep_2mi10]
		+ahy->val[_HELPER_2ITZ2C_Ydep_ma6]
		-(*mi7)
		-(*mi12);
	
	fA.y0=
		(*mi14)
		+ahx->val[_HELPER_2ITZ2C_Xdep_2mi15]
		+(*mi2)
		+ahx->val[_HELPER_2ITZ2C_Xdep_2mi16]
		-(*ma1)
		-ahy->val[_HELPER_2ITZ2C_Ydep_2ma17]
		+seedC0im;
	fA.y1=
		(*ma14)
		+ahx->val[_HELPER_2ITZ2C_Xdep_2ma15]
		+(*ma2)
		+ahx->val[_HELPER_2ITZ2C_Xdep_2ma16]
		-(*mi1)
		-ahy->val[_HELPER_2ITZ2C_Ydep_2mi17]
		+seedC1im;

	return;
	#endif
}

#define _2ITZ2CPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=A.x0*A.x0;\
	NUMTYP x03=x02*A.x0;\
	NUMTYP x04=x02*x02;\
	NUMTYP x12=A.x1*A.x1;\
	NUMTYP x13=x12*A.x1;\
	NUMTYP x14=x12*x12;\
	NUMTYP mi3,ma3;\
	ahx->val[_HELPER_2ITZ2C_Xdep_x03]=x03;\
	ahx->val[_HELPER_2ITZ2C_Xdep_x13]=x13;\
	MINMAX2(mi3,ma3,x02,x12);\
	ahx->val[_HELPER_2ITZ2C_Xdep_6mi3]=6*mi3;\
	ahx->val[_HELPER_2ITZ2C_Xdep_6ma3]=6*ma3;\
	NUMTYP mi5,ma5;\
	MINMAX2(mi5,ma5,x04,x14);\
	ahx->val[_HELPER_2ITZ2C_Xdep_mi5]=mi5;\
	ahx->val[_HELPER_2ITZ2C_Xdep_ma5]=ma5;\
	NUMTYP mi8,ma8;\
	MINMAX2(mi8,ma8,C0RE*C0RE,C1RE*C1RE);\
	ahx->val[_HELPER_2ITZ2C_Xdep_mi8]=mi8;\
	ahx->val[_HELPER_2ITZ2C_Xdep_ma8]=ma8;\
	NUMTYP mi9,ma9;\
	MINMAX4(mi9,ma9,C0RE*mi3,C0RE*ma3,C1RE*mi3,C1RE*ma3);\
	ahx->val[_HELPER_2ITZ2C_Xdep_2mi9]=2*mi9;\
	ahx->val[_HELPER_2ITZ2C_Xdep_2ma9]=2*ma9;\
	NUMTYP mi13,ma13;\
	MINMAX2(mi13,ma13,C0IM*C0IM,C1IM*C1IM);\
	ahx->val[_HELPER_2ITZ2C_Xdep_mi13]=mi13;\
	ahx->val[_HELPER_2ITZ2C_Xdep_ma13]=ma13;\
	NUMTYP mi15,ma15;\
	IAMUL_SIGN(mi15,ma15,C0IM,C1IM,C0RE,C1RE,MINMAX4);\
	ahx->val[_HELPER_2ITZ2C_Xdep_2mi15]=2*mi15;\
	ahx->val[_HELPER_2ITZ2C_Xdep_2ma15]=2*ma15;\
	NUMTYP mi16,ma16;\
	IAMUL_SIGN(mi16,ma16,C0IM,C1IM,mi3,ma3,MINMAX4);\
	ahx->val[_HELPER_2ITZ2C_Xdep_2mi16]=2*mi16;\
	ahx->val[_HELPER_2ITZ2C_Xdep_2ma16]=2*ma16;\
	ahx->val[_HELPER_2ITZ2C_Xdep_mi8_plus_2mi9_plus_C0RE_plus_mi5_minus_ma13]=\
		((((ahx->val[_HELPER_2ITZ2C_Xdep_mi8]\
		+ahx->val[_HELPER_2ITZ2C_Xdep_2mi9])\
		+C0RE)\
		+ahx->val[_HELPER_2ITZ2C_Xdep_mi5])\
		-ahx->val[_HELPER_2ITZ2C_Xdep_ma13]);\
	\
	ahx->val[_HELPER_2ITZ2C_Xdep_ma8_plus_2ma9_plus_C1RE_plus_ma5_minus_mi13]=\
		((((ahx->val[_HELPER_2ITZ2C_Xdep_ma8]\
		+ahx->val[_HELPER_2ITZ2C_Xdep_2ma9])\
		+C1RE)\
		+ahx->val[_HELPER_2ITZ2C_Xdep_ma5])\
		-ahx->val[_HELPER_2ITZ2C_Xdep_mi13]);\
	\
	return;\
}

#define _2ITZ2CPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	_2ITZ2CPRECOMPUTEX_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	);\
}



void precomputeXdep_2itz2c_double(
	PlaneRect_double& A,Helper_double *ahx
) {
	_2ITZ2CPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeXdep_2itz2c(
	PlaneRect& A,Helper *ahx
) {
	#ifdef _DOUBLE
	_2ITZ2CPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	_2ITZ2CPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	_2ITZ2CPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	_2ITZ2CPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	_2ITZ2CPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	_2ITZ2CPRECOMPUTEX(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

#define _2ITZ2CPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=A.y0*A.y0;\
	NUMTYP y03=y02*A.y0;\
	NUMTYP y04=y02*y02;\
	NUMTYP y12=A.y1*A.y1;\
	NUMTYP y13=y12*A.y1;\
	NUMTYP y14=y12*y12;\
	ahy->val[_HELPER_2ITZ2C_Ydep_y03]=y03;\
	ahy->val[_HELPER_2ITZ2C_Ydep_y13]=y13;\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,y02,y12);\
	ahy->val[_HELPER_2ITZ2C_Ydep_mi4]=mi4;\
	ahy->val[_HELPER_2ITZ2C_Ydep_ma4]=ma4;\
	NUMTYP mi6,ma6;\
	MINMAX2(mi6,ma6,y04,y14);\
	ahy->val[_HELPER_2ITZ2C_Ydep_mi6]=mi6;\
	ahy->val[_HELPER_2ITZ2C_Ydep_ma6]=ma6;\
	NUMTYP mi10,ma10;\
	IAMUL_SIGN(mi10,ma10,C0RE,C1RE,mi4,ma4,MINMAX4);\
	ahy->val[_HELPER_2ITZ2C_Ydep_2mi10]=2*mi10;\
	ahy->val[_HELPER_2ITZ2C_Ydep_2ma10]=2*ma10;\
	NUMTYP mi17,ma17;\
	IAMUL_SIGN(mi17,ma17,C0IM,C1IM,mi4,ma4,MINMAX4);\
	ahy->val[_HELPER_2ITZ2C_Ydep_2mi17]=2*mi17;\
	ahy->val[_HELPER_2ITZ2C_Ydep_2ma17]=2*ma17;\
	return;\
}\

#define _2ITZ2CPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	_2ITZ2CPRECOMPUTEY_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeYdep_2itz2c_double(
	PlaneRect_double& A,Helper_double *ahy
) {
	_2ITZ2CPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeYdep_2itz2c(
	PlaneRect& A,Helper *ahy
) {
	#ifdef _DOUBLE
	_2ITZ2CPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	_2ITZ2CPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	_2ITZ2CPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	_2ITZ2CPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	_2ITZ2CPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	_2ITZ2CPRECOMPUTEY(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

void getBoundingBoxfA_z3azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*(minimumD(A.x0*A.x0,A.x1*A.x1)-(3*maximumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre),A.x0*(maximumD(A.x0*A.x0,A.x1*A.x1)-(3*minimumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre),A.x1*(minimumD(A.x0*A.x0,A.x1*A.x1)-(3*maximumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre),A.x1*(maximumD(A.x0*A.x0,A.x1*A.x1)-(3*minimumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre));
	fA.x1=seedC1re-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*(minimumD(A.x0*A.x0,A.x1*A.x1)-(3*maximumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre),A.x0*(maximumD(A.x0*A.x0,A.x1*A.x1)-(3*minimumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre),A.x1*(minimumD(A.x0*A.x0,A.x1*A.x1)-(3*maximumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre),A.x1*(maximumD(A.x0*A.x0,A.x1*A.x1)-(3*minimumD(A.y0*A.y0,A.y1*A.y1))+FAKTORAre));
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(A.y0*(3*minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre),A.y0*(3*maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre),A.y1*(3*minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre),A.y1*(3*maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre));
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(A.y0*(3*minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre),A.y0*(3*maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre),A.y1*(3*minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre),A.y1*(3*maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+FAKTORAre));
}

void getBoundingBoxfA_z3azc_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=seedC0re_double-maximumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+minimumdouble(A.x0*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(3*maximumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double),A.x0*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(3*minimumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double),A.x1*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(3*maximumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double),A.x1*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(3*minimumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double));
	fA.x1=seedC1re_double-minimumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+maximumdouble(A.x0*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(3*maximumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double),A.x0*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(3*minimumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double),A.x1*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(3*maximumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double),A.x1*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(3*minimumdouble(A.y0*A.y0,A.y1*A.y1))+FAKTORAre_double));
	fA.y0=seedC0im_double+minimumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+minimumdouble(A.y0*(3*minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double),A.y0*(3*maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double),A.y1*(3*minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double),A.y1*(3*maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double));
	fA.y1=seedC1im_double+maximumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+maximumdouble(A.y0*(3*minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double),A.y0*(3*maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double),A.y1*(3*minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double),A.y1*(3*maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)+FAKTORAre_double));
}

#define Z3AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP mi5,ma5;\
	NUMTYP tmp1=(ahx->val[_HELPER_Z3AZC_Xdep_mi1]+ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_ma2]);\
	NUMTYP tmp2=(ahx->val[_HELPER_Z3AZC_Xdep_ma1]+ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_mi2]);\
	IAMUL_SIGN(mi5,ma5,A.x0,A.x1,tmp1,tmp2,MINMAX4)\
	NUMTYP mi6,ma6;\
	NUMTYP tmp3=(ahx->val[_HELPER_Z3AZC_Xdep_3_mul_mi1]+ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_ma2]);\
	NUMTYP tmp4=(ahx->val[_HELPER_Z3AZC_Xdep_3_mul_ma1]+ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_mi2]);\
	IAMUL_SIGN(mi6,ma6,A.y0,A.y1,tmp3,tmp4,MINMAX4);\
	\
	fA.x0=(mi5+ahy->val[_HELPER_Z3AZC_Ydep_seedC0re_minus_ma3]);\
	fA.x1=(ma5+ahy->val[_HELPER_Z3AZC_Ydep_seedC1re_minus_mi3]);\
	\
	fA.y0=(mi6+ahx->val[_HELPER_Z3AZC_Xdep_seedC0im_plus_mi4]);\
	fA.y1=(ma6+ahx->val[_HELPER_Z3AZC_Xdep_seedC1im_plus_ma4]);\
	\
	return;\
}

#define Z3AZCHELPER(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z3AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	);\
}

void getBoundingBoxfA_z3azc_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double *ahx,Helper_double* ahy
) {
	Z3AZCHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	);
}

void getBoundingBoxfA_z3azc_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper *ahx,Helper* ahy
) {
	ctrbbxfa++;
	
	#ifdef _DOUBLE
	Z3AZCHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z3AZCHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z3AZCHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z3AZCHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z3AZCHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	PFPA mi5,ma5;
	FPA array1[4];
	FPA tmp1,tmp2;
	FPA_add_ZAB(tmp1,
		&ahx->val[_HELPER_Z3AZC_Xdep_mi1],
		&ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_ma2]
	);
	FPA_add_ZAB(tmp2,
		&ahx->val[_HELPER_Z3AZC_Xdep_ma1],
		&ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_mi2]
	);
	IAMUL_FPA_SIGN(mi5,ma5,A.x0,A.x1,tmp1,tmp2,minimaxFPAABCD,array1)
	PFPA mi6,ma6;
	FPA tmp3,tmp4;
	FPA_add_ZAB(tmp3,
		&ahx->val[_HELPER_Z3AZC_Xdep_3_mul_mi1],
		&ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_ma2]
	);
	FPA_add_ZAB(tmp4,
		&ahx->val[_HELPER_Z3AZC_Xdep_3_mul_ma1],
		&ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_mi2]
	);
	FPA array2[4];
	IAMUL_FPA_SIGN(mi6,ma6,A.y0,A.y1,tmp3,tmp4,minimaxFPAABCD,array2);
	
	FPA_add_ZAB(fA.x0,
		mi5,
		&ahy->val[_HELPER_Z3AZC_Ydep_seedC0re_minus_ma3]
	);
	FPA_add_ZAB(fA.x1,
		ma5,
		&ahy->val[_HELPER_Z3AZC_Ydep_seedC1re_minus_mi3]
	);
	
	FPA_add_ZAB(fA.y0,
		mi6,
		&ahx->val[_HELPER_Z3AZC_Xdep_seedC0im_plus_mi4]
	);
	FPA_add_ZAB(fA.y1,
		ma6,
		&ahx->val[_HELPER_Z3AZC_Xdep_seedC1im_plus_ma4]
	);
	
	return;
	#endif

}

#define Z3AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=A.x0*A.x0;\
	ahx->val[_HELPER_Z3AZC_Xdep_x02]=x02;\
	NUMTYP x12=A.x1*A.x1;\
	ahx->val[_HELPER_Z3AZC_Xdep_x12]=x12;\
	NUMTYP mi1,ma1;\
	MINMAX2(mi1,ma1,x02,x12);\
	ahx->val[_HELPER_Z3AZC_Xdep_mi1]=mi1;\
	ahx->val[_HELPER_Z3AZC_Xdep_ma1]=ma1;\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,AIM*A.x0,AIM*A.x1);\
	ahx->val[_HELPER_Z3AZC_Xdep_mi4]=mi4;\
	ahx->val[_HELPER_Z3AZC_Xdep_ma4]=ma4;\
	ahx->val[_HELPER_Z3AZC_Xdep_3_mul_mi1]=3*mi1;\
	ahx->val[_HELPER_Z3AZC_Xdep_3_mul_ma1]=3*ma1;\
	ahx->val[_HELPER_Z3AZC_Xdep_seedC0im_plus_mi4]=C0IM+mi4;\
	ahx->val[_HELPER_Z3AZC_Xdep_seedC1im_plus_ma4]=C1IM+ma4;\
	return;\
}

#define Z3AZCPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z3AZCPRECOMPUTEX_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	);\
}



void precomputeXdep_z3azc_double(
	PlaneRect_double& A,Helper_double *ahx
) {
	Z3AZCPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeXdep_z3azc(
	PlaneRect& A,Helper *ahx
) {
	#ifdef _DOUBLE
	Z3AZCPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z3AZCPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z3AZCPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z3AZCPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z3AZCPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z3AZCPRECOMPUTEX(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

#define Z3AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=A.y0*A.y0;\
	ahy->val[_HELPER_Z3AZC_Ydep_y02]=y02;\
	NUMTYP y12=A.y1*A.y1;\
	ahy->val[_HELPER_Z3AZC_Ydep_y12]=y12;\
	NUMTYP mi2,ma2;\
	MINMAX2(mi2,ma2,y02,y12);\
	ahy->val[_HELPER_Z3AZC_Ydep_mi2]=mi2;\
	ahy->val[_HELPER_Z3AZC_Ydep_ma2]=ma2;\
	NUMTYP mi3,ma3;\
	MINMAX2(mi3,ma3,AIM*A.y0,AIM*A.y1);\
	ahy->val[_HELPER_Z3AZC_Ydep_mi3]=mi3;\
	ahy->val[_HELPER_Z3AZC_Ydep_ma3]=ma3;\
	ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_ma2]=ARE-3*ma2;\
	ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_3_mul_mi2]=ARE-3*mi2;\
	ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_ma2]=ARE-ma2;\
	ahy->val[_HELPER_Z3AZC_Ydep_FAKTORAre_minus_mi2]=ARE-mi2;\
	ahy->val[_HELPER_Z3AZC_Ydep_seedC0re_minus_ma3]=C0RE-ma3;\
	ahy->val[_HELPER_Z3AZC_Ydep_seedC1re_minus_mi3]=C1RE-mi3;\
	return;\
}\

#define Z3AZCPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z3AZCPRECOMPUTEY_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeYdep_z3azc_double(
	PlaneRect_double& A,Helper_double *ahy
) {
	Z3AZCPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeYdep_z3azc(
	PlaneRect& A,Helper *ahy
) {
	#ifdef _DOUBLE
	Z3AZCPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z3AZCPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z3AZCPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z3AZCPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z3AZCPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z3AZCPRECOMPUTEY(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

void getBoundingBoxfA_z4azc_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=seedC0re_double+minimumdouble(FAKTORAre_double*A.x0,FAKTORAre_double*A.x1)+minimumdouble((A.y0*A.y0*A.y0-FAKTORAim_double)*A.y0,(A.y0*A.y0*A.y0-FAKTORAim_double)*A.y1,(A.y1*A.y1*A.y1-FAKTORAim_double)*A.y0,(A.y1*A.y1*A.y1-FAKTORAim_double)*A.y1)+minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(6*maximumdouble(A.y0*A.y0,A.y1*A.y1))),minimumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(6*minimumdouble(A.y0*A.y0,A.y1*A.y1))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(6*maximumdouble(A.y0*A.y0,A.y1*A.y1))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(6*minimumdouble(A.y0*A.y0,A.y1*A.y1))));
	fA.x1=seedC1re_double+maximumdouble(FAKTORAre_double*A.x0,FAKTORAre_double*A.x1)+maximumdouble((A.y0*A.y0*A.y0-FAKTORAim_double)*A.y0,(A.y0*A.y0*A.y0-FAKTORAim_double)*A.y1,(A.y1*A.y1*A.y1-FAKTORAim_double)*A.y0,(A.y1*A.y1*A.y1-FAKTORAim_double)*A.y1)+maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(6*maximumdouble(A.y0*A.y0,A.y1*A.y1))),minimumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(6*minimumdouble(A.y0*A.y0,A.y1*A.y1))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(6*maximumdouble(A.y0*A.y0,A.y1*A.y1))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(6*minimumdouble(A.y0*A.y0,A.y1*A.y1))));
	fA.y0=seedC0im_double+minimumdouble(FAKTORAre_double*A.y0,FAKTORAre_double*A.y1)+minimumdouble(A.x0*(FAKTORAim_double+minimumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))),A.x0*(FAKTORAim_double+maximumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim_double+minimumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim_double+maximumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))));
	fA.y1=seedC1im_double+maximumdouble(FAKTORAre_double*A.y0,FAKTORAre_double*A.y1)+maximumdouble(A.x0*(FAKTORAim_double+minimumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))),A.x0*(FAKTORAim_double+maximumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim_double+minimumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim_double+maximumdouble((4*A.y0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-maximumdouble(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-minimumdouble(A.y0*A.y0,A.y1*A.y1)))));
}

void getBoundingBoxfA_z4azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	// standard un.optimized calcuation
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+minimumD((A.y0*A.y0*A.y0-FAKTORAim)*A.y0,(A.y0*A.y0*A.y0-FAKTORAim)*A.y1,(A.y1*A.y1*A.y1-FAKTORAim)*A.y0,(A.y1*A.y1*A.y1-FAKTORAim)*A.y1)+minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(6*maximumD(A.y0*A.y0,A.y1*A.y1))),minimumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(6*minimumD(A.y0*A.y0,A.y1*A.y1))),maximumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(6*maximumD(A.y0*A.y0,A.y1*A.y1))),maximumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(6*minimumD(A.y0*A.y0,A.y1*A.y1))));
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+maximumD((A.y0*A.y0*A.y0-FAKTORAim)*A.y0,(A.y0*A.y0*A.y0-FAKTORAim)*A.y1,(A.y1*A.y1*A.y1-FAKTORAim)*A.y0,(A.y1*A.y1*A.y1-FAKTORAim)*A.y1)+maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(6*maximumD(A.y0*A.y0,A.y1*A.y1))),minimumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(6*minimumD(A.y0*A.y0,A.y1*A.y1))),maximumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(6*maximumD(A.y0*A.y0,A.y1*A.y1))),maximumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(6*minimumD(A.y0*A.y0,A.y1*A.y1))));
	fA.y0=seedC0im+minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(A.x0*(FAKTORAim+minimumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))),A.x0*(FAKTORAim+maximumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim+minimumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim+maximumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))));
	fA.y1=seedC1im+maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(A.x0*(FAKTORAim+minimumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))),A.x0*(FAKTORAim+maximumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim+minimumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))),A.x1*(FAKTORAim+maximumD((4*A.y0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)),(4*A.y1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)))));
}

#define Z4AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP t2a=(ahx->val[_HELPER_Z4AZC_Xdep_mi1]-ahy->val[_HELPER_Z4AZC_Ydep_6ma2]);\
	NUMTYP t2b=(ahx->val[_HELPER_Z4AZC_Xdep_ma1]-ahy->val[_HELPER_Z4AZC_Ydep_6mi2]);\
	NUMTYP mi6,ma6;\
	IAMUL_SIGN(mi6,ma6,\
		ahx->val[_HELPER_Z4AZC_Xdep_mi1],\
		ahx->val[_HELPER_Z4AZC_Xdep_ma1],\
		t2a,t2b,MINMAX4)\
	NUMTYP t3a=(ahx->val[_HELPER_Z4AZC_Xdep_mi1]-ahy->val[_HELPER_Z4AZC_Ydep_ma2]);\
	NUMTYP t3b=(ahx->val[_HELPER_Z4AZC_Xdep_ma1]-ahy->val[_HELPER_Z4AZC_Ydep_mi2]);\
	NUMTYP mi7,ma7;\
	IAMUL_SIGN(mi7,ma7,\
		ahy->val[_HELPER_Z4AZC_Ydep_4y0],\
		ahy->val[_HELPER_Z4AZC_Ydep_4y1],\
		t3a,t3b,MINMAX4)\
	NUMTYP t4a=(AIM+mi7);\
	NUMTYP t4b=(AIM+ma7);\
	NUMTYP mi8,ma8;\
	IAMUL_SIGN(mi8,ma8,A.x0,A.x1,t4a,t4b,MINMAX4)\
	\
	fA.x0=\
	(\
	ahx->val[_HELPER_Z4AZC_Xdep_C0re_plus_mi3]\
	+ahy->val[_HELPER_Z4AZC_Ydep_mi5]\
	)\
	+mi6;\
	fA.x1=\
	(\
	ahx->val[_HELPER_Z4AZC_Xdep_C1re_plus_ma3]\
	+ahy->val[_HELPER_Z4AZC_Ydep_ma5]\
	)\
	+ma6;\
	\
	fA.y0=(ahy->val[_HELPER_Z4AZC_Ydep_C0im_plus_mi4]+mi8);\
	fA.y1=(ahy->val[_HELPER_Z4AZC_Ydep_C1im_plus_ma4]+ma8);\
	\
	return;\
}

#define Z4AZCHELPER(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z4AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,FAKTORAre,FAKTORAim,seedC0re,seedC1re,seedC0im,seedC1im)\
}

void getBoundingBoxfA_z4azc_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double* ahx,Helper_double* ahy
 ) {
	Z4AZCHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	);
}

#define Z4AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=(A.y0*A.y0);\
	NUMTYP y03=(y02*A.y0);\
	NUMTYP y12=(A.y1*A.y1);\
	NUMTYP y13=(y12*A.y1);\
	NUMTYP mi2,ma2;\
	MINMAX2(mi2,ma2,y02,y12);\
	ahy->val[_HELPER_Z4AZC_Ydep_mi2]=mi2;\
	ahy->val[_HELPER_Z4AZC_Ydep_ma2]=ma2;\
	ahy->val[_HELPER_Z4AZC_Ydep_6mi2]=(6*mi2);\
	ahy->val[_HELPER_Z4AZC_Ydep_6ma2]=(6*ma2);\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,ARE*A.y0,ARE*A.y1);\
	ahy->val[_HELPER_Z4AZC_Ydep_mi4]=mi4;\
	ahy->val[_HELPER_Z4AZC_Ydep_ma4]=ma4;\
	NUMTYP t1a=(y03-AIM);\
	NUMTYP t1b=(y13-AIM);\
	NUMTYP mi5,ma5;\
	IAMUL_SIGN(mi5,ma5,t1a,t1b,A.y0,A.y1,MINMAX4)\
	ahy->val[_HELPER_Z4AZC_Ydep_mi5]=mi5;\
	ahy->val[_HELPER_Z4AZC_Ydep_ma5]=ma5;\
	ahy->val[_HELPER_Z4AZC_Ydep_4y0]=(4*A.y0);\
	ahy->val[_HELPER_Z4AZC_Ydep_4y1]=(4*A.y1);\
	ahy->val[_HELPER_Z4AZC_Ydep_C0im_plus_mi4]=(C0IM+mi4);\
	ahy->val[_HELPER_Z4AZC_Ydep_C1im_plus_ma4]=(C1IM+ma4);\
	return;\
}\

#define Z4AZCPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z4AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,FAKTORAre,FAKTORAim,seedC0re,seedC1re,seedC0im,seedC1im);\
}

void precomputeYdep_z4azc_double(
	PlaneRect_double& A,Helper_double *ahy
) {
	Z4AZCPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeYdep_z4azc(
	PlaneRect& A,Helper *ahy
) {
	#ifdef _DOUBLE
	Z4AZCPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z4AZCPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z4AZCPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z4AZCPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z4AZCPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z4AZCPRECOMPUTEY(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif

}

#define Z4AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=(A.x0*A.x0);\
	NUMTYP x12=(A.x1*A.x1);\
	NUMTYP mi1,ma1;\
	MINMAX2(mi1,ma1,x02,x12);\
	ahx->val[_HELPER_Z4AZC_Xdep_mi1]=mi1;\
	ahx->val[_HELPER_Z4AZC_Xdep_ma1]=ma1;\
	NUMTYP mi3,ma3;\
	MINMAX2(mi3,ma3,(ARE*A.x0),(ARE*A.x1));\
	ahx->val[_HELPER_Z4AZC_Xdep_mi3]=mi3;\
	ahx->val[_HELPER_Z4AZC_Xdep_ma3]=ma3;\
	ahx->val[_HELPER_Z4AZC_Xdep_C0re_plus_mi3]=(C0RE+mi3);\
	ahx->val[_HELPER_Z4AZC_Xdep_C1re_plus_ma3]=(C1RE+ma3);\
	return;\
}

#define Z4AZCPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z4AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im);\
}

void precomputeXdep_z4azc_double(
	PlaneRect_double& A,Helper_double *ahx
) {
	Z4AZCPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeXdep_z4azc(
	PlaneRect& A,Helper *ahx
) {
	#ifdef _DOUBLE
	Z4AZCPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z4AZCPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z4AZCPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z4AZCPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z4AZCPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z4AZCPRECOMPUTEX(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

void getBoundingBoxfA_z4azc_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper* ahx,Helper* ahy
) {
	ctrbbxfa++;
	
	#ifdef _DOUBLE
	Z4AZCHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z4AZCHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _F107
	Z4AZCHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z4AZCHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif

	#ifdef _QUADMATH
	Z4AZCHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _FPA
	// optomized version without operator overloading
	FPA t2a,t2b;
	FPA_sub_ZAB(t2a,
		&ahx->val[_HELPER_Z4AZC_Xdep_mi1],
		&ahy->val[_HELPER_Z4AZC_Ydep_6ma2]
	);
	FPA_sub_ZAB(t2b,
		&ahx->val[_HELPER_Z4AZC_Xdep_ma1],
		&ahy->val[_HELPER_Z4AZC_Ydep_6mi2]
	);

	PFPA mi6,ma6;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi6,ma6,
		ahx->val[_HELPER_Z4AZC_Xdep_mi1],
		ahx->val[_HELPER_Z4AZC_Xdep_ma1],
		t2a,t2b,minimaxFPAABCD,
		array1)
	FPA t3a,t3b;
	FPA_sub_ZAB(t3a,
		&ahx->val[_HELPER_Z4AZC_Xdep_mi1],
		&ahy->val[_HELPER_Z4AZC_Ydep_ma2]
	);
	FPA_sub_ZAB(t3b,
		&ahx->val[_HELPER_Z4AZC_Xdep_ma1],
		&ahy->val[_HELPER_Z4AZC_Ydep_mi2]
	);
	PFPA mi7,ma7;
	FPA array2[4];
	IAMUL_FPA_SIGN(mi7,ma7,
		ahy->val[_HELPER_Z4AZC_Ydep_4y0],
		ahy->val[_HELPER_Z4AZC_Ydep_4y1],
		t3a,t3b,minimaxFPAABCD,array2)
		
	FPA t4a,t4b;
	FPA_add_ZAB(t4a,&FAKTORAim,mi7);
	FPA_add_ZAB(t4b,&FAKTORAim,ma7);
	
	PFPA mi8,ma8;
	FPA array3[4];
	IAMUL_FPA_SIGN(mi8,ma8,
		A.x0,A.x1,t4a,t4b,
		minimaxFPAABCD,array3)
	
	FPA e1;
	FPA_add_ZAB(e1,
		&ahx->val[_HELPER_Z4AZC_Xdep_C0re_plus_mi3],
		&ahy->val[_HELPER_Z4AZC_Ydep_mi5]
	);
	FPA_add_ZAB(fA.x0,&e1,mi6);

	FPA_add_ZAB(e1,
		&ahx->val[_HELPER_Z4AZC_Xdep_C1re_plus_ma3],
		&ahy->val[_HELPER_Z4AZC_Ydep_ma5]
	);
	FPA_add_ZAB(fA.x1,&e1,ma6);
	
	FPA_add_ZAB(fA.y0,
		&ahy->val[_HELPER_Z4AZC_Ydep_C0im_plus_mi4],
		mi8);
	FPA_add_ZAB(fA.y1,
		&ahy->val[_HELPER_Z4AZC_Ydep_C1im_plus_ma4],
		ma8);
	
	return;
	#endif
}

void getBoundingBoxfA_z5azc_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=seedC0re_double-maximumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+minimumdouble(A.x0*(5*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double),A.x0*(5*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double),A.x1*(5*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double),A.x1*(5*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double))+minimumdouble((A.x0*A.x0*A.x0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))),(A.x0*A.x0*A.x0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))));
	fA.x1=seedC1re_double-minimumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+maximumdouble(A.x0*(5*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double),A.x0*(5*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double),A.x1*(5*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double),A.x1*(5*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double))+maximumdouble((A.x0*A.x0*A.x0)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))),(A.x0*A.x0*A.x0)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(minimumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(maximumdouble(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))));
	fA.y0=seedC0im_double+minimumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+minimumdouble(A.y0*(5*minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y0*(5*maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y1*(5*minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y1*(5*maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double))+minimumdouble((A.y0*A.y0*A.y0)*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),(A.y0*A.y0*A.y0)*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))));
	fA.y1=seedC1im_double+maximumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+maximumdouble(A.y0*(5*minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y0*(5*maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y1*(5*minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y1*(5*maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double))+maximumdouble((A.y0*A.y0*A.y0)*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),(A.y0*A.y0*A.y0)*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))));
}

void getBoundingBoxfA_z5azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*(5*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre),A.x0*(5*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre),A.x1*(5*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre),A.x1*(5*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre))+minimumD((A.x0*A.x0*A.x0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))),(A.x0*A.x0*A.x0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))));
	fA.x1=seedC1re-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*(5*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre),A.x0*(5*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre),A.x1*(5*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre),A.x1*(5*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+FAKTORAre))+maximumD((A.x0*A.x0*A.x0)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))),(A.x0*A.x0*A.x0)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(minimumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))),(A.x1*A.x1*A.x1)*(maximumD(A.x0*A.x0,A.x1*A.x1)-(2*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))));
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(A.y0*(5*minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y0*(5*maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y1*(5*minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y1*(5*maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre))+minimumD((A.y0*A.y0*A.y0)*(minimumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),(A.y0*A.y0*A.y0)*(maximumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(minimumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(maximumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))));
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(A.y0*(5*minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y0*(5*maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y1*(5*minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y1*(5*maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+FAKTORAre))+maximumD((A.y0*A.y0*A.y0)*(minimumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),(A.y0*A.y0*A.y0)*(maximumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(minimumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),(A.y1*A.y1*A.y1)*(maximumD(A.y0*A.y0,A.y1*A.y1)-(2*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))));
}

#define Z5AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP mi7,ma7;\
	IAMUL_SIGN(mi7,ma7,\
		A.x0,A.x1,\
		ahy->val[_HELPER_Z5AZC_Ydep_tmp1],\
		ahy->val[_HELPER_Z5AZC_Ydep_tmp2],\
		MINMAX4);\
	NUMTYP mi8,ma8;\
	NUMTYP tmp3=(ahx->val[_HELPER_Z5AZC_Xdep_mi1]-ahy->val[_HELPER_Z5AZC_Ydep_10_mul_ma2]);\
	NUMTYP tmp4=(ahx->val[_HELPER_Z5AZC_Xdep_ma1]-ahy->val[_HELPER_Z5AZC_Ydep_10_mul_mi2]);\
	IAMUL_SIGN(mi8,ma8,\
		tmp3,tmp4,\
		ahx->val[_HELPER_Z5AZC_Xdep_x03],\
		ahx->val[_HELPER_Z5AZC_Xdep_x13],\
		MINMAX4);\
	NUMTYP mi9,ma9;\
	IAMUL_SIGN(mi9,ma9,\
		A.y0,A.y1,\
		ahx->val[_HELPER_Z5AZC_Xdep_tmp5],\
		ahx->val[_HELPER_Z5AZC_Xdep_tmp6],\
		MINMAX4);\
	NUMTYP mi10,ma10;\
	NUMTYP tmp7=(ahy->val[_HELPER_Z5AZC_Ydep_mi2]-ahx->val[_HELPER_Z5AZC_Xdep_10_mul_ma1]);\
	NUMTYP tmp8=(ahy->val[_HELPER_Z5AZC_Ydep_ma2]-ahx->val[_HELPER_Z5AZC_Xdep_10_mul_mi1]);\
	IAMUL_SIGN(mi10,ma10,\
		tmp7,tmp8,\
		ahy->val[_HELPER_Z5AZC_Ydep_y03],\
		ahy->val[_HELPER_Z5AZC_Ydep_y13],\
		MINMAX4);\
	\
	fA.x0=\
		(ahy->val[_HELPER_Z5AZC_Ydep_C0re_minus_ma5]+mi7) \
		+mi8;\
	fA.x1=\
		(ahy->val[_HELPER_Z5AZC_Ydep_C1re_minus_mi5]+ma7)\
		+ma8;\
	\
	fA.y0=\
		(ahx->val[_HELPER_Z5AZC_Xdep_C0im_plus_mi6]+mi9)\
		+mi10;\
	fA.y1=\
		(ahx->val[_HELPER_Z5AZC_Xdep_C1im_plus_ma6]+ma9)\
		+ma10;\
	\
	return;\
}

#define Z5AZCHELPER(NUMLIST,MINMAX2,MINMAX4) \
{\
	Z5AZCHELPER_VARLIST(\
		NUMLIST,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void getBoundingBoxfA_z5azc_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double* ahx,Helper_double* ahy
) {
	Z5AZCHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void getBoundingBoxfA_z5azc_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper* ahx,Helper* ahy
) {
	ctrbbxfa++;
	
	#ifdef _DOUBLE
	Z5AZCHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z5AZCHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _QUADMATH
	Z5AZCHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _F107
	Z5AZCHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z5AZCHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif

	#ifdef _FPA
	PFPA mi7,ma7;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi7,ma7,
		A.x0,A.x1,
		ahy->val[_HELPER_Z5AZC_Ydep_tmp1],
		ahy->val[_HELPER_Z5AZC_Ydep_tmp2],
		minimaxFPAABCD,array1);
	PFPA mi8,ma8;
	FPA array2[4];
	FPA tmp3=ahx->val[_HELPER_Z5AZC_Xdep_mi1]-ahy->val[_HELPER_Z5AZC_Ydep_10_mul_ma2];
	FPA tmp4=ahx->val[_HELPER_Z5AZC_Xdep_ma1]-ahy->val[_HELPER_Z5AZC_Ydep_10_mul_mi2];
	IAMUL_FPA_SIGN(mi8,ma8,
		tmp3,tmp4,
		ahx->val[_HELPER_Z5AZC_Xdep_x03],
		ahx->val[_HELPER_Z5AZC_Xdep_x13],
		minimaxFPAABCD,array2);
	PFPA mi9,ma9;
	FPA array3[4];
	IAMUL_FPA_SIGN(mi9,ma9,
		A.y0,A.y1,
		ahx->val[_HELPER_Z5AZC_Xdep_tmp5],
		ahx->val[_HELPER_Z5AZC_Xdep_tmp6],
		minimaxFPAABCD,array3);
	PFPA mi10,ma10;
	FPA tmp7,tmp8;
	FPA_sub_ZAB(tmp7,
		&ahy->val[_HELPER_Z5AZC_Ydep_mi2],
		&ahx->val[_HELPER_Z5AZC_Xdep_10_mul_ma1]
	);
	FPA_sub_ZAB(tmp8,
		&ahy->val[_HELPER_Z5AZC_Ydep_ma2],
		&ahx->val[_HELPER_Z5AZC_Xdep_10_mul_mi1]
	);
	FPA array4[4];
	IAMUL_FPA_SIGN(mi10,ma10,
		tmp7,tmp8,
		ahy->val[_HELPER_Z5AZC_Ydep_y03],
		ahy->val[_HELPER_Z5AZC_Ydep_y13],
		minimaxFPAABCD,array4);
	
	FPA e1,e2;
	FPA_add_ZAB(e1,
		&ahy->val[_HELPER_Z5AZC_Ydep_C0re_minus_ma5],
		mi7);
	FPA_add_ZAB(fA.x0,&e1,mi8);
	FPA_add_ZAB(e2,
		&ahy->val[_HELPER_Z5AZC_Ydep_C1re_minus_mi5],
		ma7);
	FPA_add_ZAB(fA.x1,&e2,ma8);
		
	FPA_add_ZAB(e1,
		&ahx->val[_HELPER_Z5AZC_Xdep_C0im_plus_mi6],
		mi9);
	FPA_add_ZAB(fA.y0,&e1,mi10);
	FPA_add_ZAB(e2,
		&ahx->val[_HELPER_Z5AZC_Xdep_C1im_plus_ma6],
		ma9);
	FPA_add_ZAB(fA.y1,&e2,ma10);
	
	return;
	#endif
}

#define Z5AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=A.y0*A.y0;\
	NUMTYP y03=y02*A.y0;\
	ahy->val[_HELPER_Z5AZC_Ydep_y03]=y03;\
	NUMTYP y04=y02*y02;\
	NUMTYP y12=A.y1*A.y1;\
	NUMTYP y13=y12*A.y1;\
	ahy->val[_HELPER_Z5AZC_Ydep_y13]=y13;\
	NUMTYP y14=y12*y12;\
	NUMTYP mi2,ma2;\
	MINMAX2(mi2,ma2,y02,y12);\
	ahy->val[_HELPER_Z5AZC_Ydep_mi2]=mi2;\
	ahy->val[_HELPER_Z5AZC_Ydep_ma2]=ma2;\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,y04,y14);\
	ahy->val[_HELPER_Z5AZC_Ydep_mi4]=mi4;\
	ahy->val[_HELPER_Z5AZC_Ydep_ma4]=ma4;\
	NUMTYP tmp1=(5*mi4)+ARE;\
	NUMTYP tmp2=(5*ma4)+ARE;\
	ahy->val[_HELPER_Z5AZC_Ydep_tmp1]=tmp1;\
	ahy->val[_HELPER_Z5AZC_Ydep_tmp2]=tmp2;\
	NUMTYP mi5,ma5;\
	MINMAX2(mi5,ma5,AIM*A.y0,AIM*A.y1);\
	ahy->val[_HELPER_Z5AZC_Ydep_mi5]=mi5;\
	ahy->val[_HELPER_Z5AZC_Ydep_ma5]=ma5;\
	ahy->val[_HELPER_Z5AZC_Ydep_10_mul_ma2]=10*ma2;\
	ahy->val[_HELPER_Z5AZC_Ydep_10_mul_mi2]=10*mi2;\
	ahy->val[_HELPER_Z5AZC_Ydep_C0re_minus_ma5]=C0RE-ma5;\
	ahy->val[_HELPER_Z5AZC_Ydep_C1re_minus_mi5]=C1RE-mi5;\
	return;\
}

#define Z5AZCPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z5AZCPRECOMPUTEY_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeYdep_z5azc_double(
	PlaneRect_double& A,Helper_double *ahy
) {
	Z5AZCPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeYdep_z5azc(
	PlaneRect& A,Helper *ahy
) {
	#ifdef _DOUBLE
	Z5AZCPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z5AZCPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z5AZCPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z5AZCPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z5AZCPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z5AZCPRECOMPUTEY(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

#define Z5AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=A.x0*A.x0;\
	NUMTYP x03=x02*A.x0;\
	NUMTYP x04=x02*x02;\
	ahx->val[_HELPER_Z5AZC_Xdep_x03]=x03;\
	NUMTYP x12=A.x1*A.x1;\
	NUMTYP x13=x12*A.x1;\
	ahx->val[_HELPER_Z5AZC_Xdep_x13]=x13;\
	NUMTYP x14=x12*x12;\
	NUMTYP mi1,ma1;\
	MINMAX2(mi1,ma1,x02,x12);\
	ahx->val[_HELPER_Z5AZC_Xdep_mi1]=mi1;\
	ahx->val[_HELPER_Z5AZC_Xdep_ma1]=ma1;\
	NUMTYP mi3,ma3;\
	MINMAX2(mi3,ma3,x04,x14);\
	ahx->val[_HELPER_Z5AZC_Xdep_mi3]=mi3;\
	ahx->val[_HELPER_Z5AZC_Xdep_ma3]=ma3;\
	NUMTYP mi6,ma6;\
	MINMAX2(mi6,ma6,AIM*A.x0,AIM*A.x1);\
	ahx->val[_HELPER_Z5AZC_Xdep_mi6]=mi6;\
	ahx->val[_HELPER_Z5AZC_Xdep_ma6]=ma6;\
	NUMTYP tmp5=(5*mi3)+ARE;\
	NUMTYP tmp6=(5*ma3)+ARE;\
	ahx->val[_HELPER_Z5AZC_Xdep_tmp5]=tmp5;\
	ahx->val[_HELPER_Z5AZC_Xdep_tmp6]=tmp6;\
	ahx->val[_HELPER_Z5AZC_Xdep_10_mul_mi1]=10*mi1;\
	ahx->val[_HELPER_Z5AZC_Xdep_10_mul_ma1]=10*ma1;\
	ahx->val[_HELPER_Z5AZC_Xdep_C0im_plus_mi6]=C0IM+mi6;\
	ahx->val[_HELPER_Z5AZC_Xdep_C1im_plus_ma6]=C1IM+ma6;\
	return;\
}

#define Z5AZCPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z5AZCPRECOMPUTEX_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeXdep_z5azc_double(
	PlaneRect_double& A,Helper_double *ahx
) {
	Z5AZCPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeXdep_z5azc(
	PlaneRect& A,Helper *ahx
) {
	#ifdef _DOUBLE
	Z5AZCPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _F161
	Z5AZCPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z5AZCPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z5AZCPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _QUADMATH
	Z5AZCPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	#ifdef _FPA
	Z5AZCPRECOMPUTEX(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif
}

void getBoundingBoxfA_z6azc_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=seedC0re_double+minimumdouble(FAKTORAre_double*A.x0,FAKTORAre_double*A.x1)-maximumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)-maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))));
	fA.x1=seedC1re_double+maximumdouble(FAKTORAre_double*A.x0,FAKTORAre_double*A.x1)-minimumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)-minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(A.x0*A.x0,A.x1*A.x1)*(maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumdouble((3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(minimumdouble(A.y0*A.y0,A.y1*A.y1)-maximumdouble(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumdouble(A.y0*A.y0,A.y1*A.y1)))*(maximumdouble(A.y0*A.y0,A.y1*A.y1)-minimumdouble(A.x0*A.x0,A.x1*A.x1)))));
	fA.y0=seedC0im_double+minimumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+minimumdouble(A.y0*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre_double),A.y0*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y1*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre_double),A.y1*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double))+minimumdouble(minimumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))));
	fA.y1=seedC1im_double+maximumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+maximumdouble(A.y0*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre_double),A.y0*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double),A.y1*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre_double),A.y1*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre_double))+maximumdouble(minimumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumdouble(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumdouble(A.x0*A.x0,A.x1*A.x1)))));
}

void getBoundingBoxfA_z6azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))));
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(A.x0*A.x0,A.x1*A.x1)*(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+minimumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(A.x0*A.x0,A.x1*A.x1)*(maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)+maximumD((3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*minimumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(minimumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1)),(3*(5*maximumD(A.y0*A.y0,A.y1*A.y1)))*(maximumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1)))));
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(A.y0*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre),A.y0*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y1*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre),A.y1*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre))+minimumD(minimumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))));
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(A.y0*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre),A.y0*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre),A.y1*(6*(A.x0*A.x0*A.x0*A.x0*A.x0)+FAKTORAre),A.y1*(6*(A.x1*A.x1*A.x1*A.x1*A.x1)+FAKTORAre))+maximumD(minimumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*minimumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1)*(6*maximumD(A.y0*A.y0,A.y1*A.y1)-(4*(5*minimumD(A.x0*A.x0,A.x1*A.x1)))));
}

#define Z6AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM)\
{\
	NUMTYP mi8,ma8;\
	NUMTYP tmp21=(ahy->val[_HELPER_Z6AZC_Ydep_mi2]-ahx->val[_HELPER_Z6AZC_Xdep_ma1]);\
	NUMTYP tmp22=(ahy->val[_HELPER_Z6AZC_Ydep_ma2]-ahx->val[_HELPER_Z6AZC_Xdep_mi1]);\
	IAMUL_SIGN(mi8,ma8,\
		tmp21,tmp22,\
		ahy->val[_HELPER_Z6AZC_Ydep_15_mul_mi2],\
		ahy->val[_HELPER_Z6AZC_Ydep_15_mul_ma2],\
		MINMAX4);\
	NUMTYP mi9,ma9;\
	NUMTYP tmp23=(ahx->val[_HELPER_Z6AZC_Xdep_mi3]+mi8);\
	NUMTYP tmp24=(ahx->val[_HELPER_Z6AZC_Xdep_ma3]+ma8);\
	IAMUL_SIGN(mi9,ma9,\
		tmp23,tmp24,\
		ahx->val[_HELPER_Z6AZC_Xdep_mi1],\
		ahx->val[_HELPER_Z6AZC_Xdep_ma1],\
		MINMAX4);\
	NUMTYP mi10,ma10;\
	IAMUL_SIGN(mi10,ma10,\
		A.y0,A.y1,\
		ahx->val[_HELPER_Z6AZC_Xdep_tmp1],\
		ahx->val[_HELPER_Z6AZC_Xdep_tmp2],\
		MINMAX4);\
	NUMTYP mi11,ma11;\
	IAMUL_SIGN(mi11,ma11,\
		A.x0,A.x1,\
		ahy->val[_HELPER_Z6AZC_Ydep_y03],\
		ahy->val[_HELPER_Z6AZC_Ydep_y13],\
		MINMAX4);\
	NUMTYP mi12,ma12;\
	NUMTYP tmp3=(ahy->val[_HELPER_Z6AZC_Ydep_6_mul_mi2]-ahx->val[_HELPER_Z6AZC_Xdep_20_mul_ma1]);\
	NUMTYP tmp4=(ahy->val[_HELPER_Z6AZC_Ydep_6_mul_ma2]-ahx->val[_HELPER_Z6AZC_Xdep_20_mul_mi1]);\
	IAMUL_SIGN(mi12,ma12,mi11,ma11,tmp3,tmp4,MINMAX4);\
	\
	fA.x0=\
		(ahx->val[_HELPER_Z6AZC_Xdep_mi4]\
		+mi9\
		)\
		+ahy->val[_HELPER_Z6AZC_Ydep_C0re_minus_ma5_minus_ma7];\
	fA.x1=\
		(ahx->val[_HELPER_Z6AZC_Xdep_ma4]\
		+ma9\
		)\
		+ahy->val[_HELPER_Z6AZC_Ydep_C1re_minus_mi5_minus_mi7];\
	fA.y0=\
		(ahx->val[_HELPER_Z6AZC_Xdep_C0im_plus_mi6]\
		+mi10\
		)\
		+mi12;\
	fA.y1=\
		(ahx->val[_HELPER_Z6AZC_Xdep_C1im_plus_ma6]\
		+ma10\
		)\
		+ma12;\
	\
	return;\
}

#define Z6AZCHELPER(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z6AZCHELPER_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void getBoundingBoxfA_z6azc_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double* ahx,Helper_double* ahy
) {
	Z6AZCHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void getBoundingBoxfA_z6azc_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper* ahx,Helper* ahy
) {
	ctrbbxfa++;
	
	#ifdef _DOUBLE
	Z6AZCHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z6AZCHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _QUADMATH
	Z6AZCHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _F107
	Z6AZCHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z6AZCHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
	
	#ifdef _FPA
	FPA tmp21,tmp22;
	FPA_sub_ZAB(tmp21,
		&ahy->val[_HELPER_Z6AZC_Ydep_mi2],
		&ahx->val[_HELPER_Z6AZC_Xdep_ma1]
	);
	FPA_sub_ZAB(tmp22,
		&ahy->val[_HELPER_Z6AZC_Ydep_ma2],
		&ahx->val[_HELPER_Z6AZC_Xdep_mi1]
	);
	PFPA mi8,ma8;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi8,ma8,
		tmp21,tmp22,
		ahy->val[_HELPER_Z6AZC_Ydep_15_mul_mi2],
		ahy->val[_HELPER_Z6AZC_Ydep_15_mul_ma2],
		minimaxFPAABCD,array1);
	FPA tmp23,tmp24;
	FPA_add_ZAB(tmp23,
		&ahx->val[_HELPER_Z6AZC_Xdep_mi3],
		mi8
	);
	FPA_add_ZAB(tmp24,
		&ahx->val[_HELPER_Z6AZC_Xdep_ma3],
		ma8
	);
	PFPA mi9,ma9;
	FPA array2[4];
	IAMUL_FPA_SIGN(mi9,ma9,
		tmp23,tmp24,
		ahx->val[_HELPER_Z6AZC_Xdep_mi1],
		ahx->val[_HELPER_Z6AZC_Xdep_ma1],
		minimaxFPAABCD,array2);
	PFPA mi10,ma10;
	FPA array3[4];
	IAMUL_FPA_SIGN(mi10,ma10,
		A.y0,A.y1,
		ahx->val[_HELPER_Z6AZC_Xdep_tmp1],
		ahx->val[_HELPER_Z6AZC_Xdep_tmp2],
		minimaxFPAABCD,array3);
	PFPA mi11,ma11;
	FPA array4[4];
	IAMUL_FPA_SIGN(mi11,ma11,
		A.x0,A.x1,
		ahy->val[_HELPER_Z6AZC_Ydep_y03],
		ahy->val[_HELPER_Z6AZC_Ydep_y13],
		minimaxFPAABCD,array4);
	FPA tmp3,tmp4;
	FPA_sub_ZAB(tmp3,
		&ahy->val[_HELPER_Z6AZC_Ydep_6_mul_mi2],
		&ahx->val[_HELPER_Z6AZC_Xdep_20_mul_ma1]
	);
	FPA_sub_ZAB(tmp4,
		&ahy->val[_HELPER_Z6AZC_Ydep_6_mul_ma2],
		&ahx->val[_HELPER_Z6AZC_Xdep_20_mul_mi1]
	);
	PFPA mi12,ma12;
	FPA array5[4];
	IAMUL_FPA_SIGN(mi12,ma12,
		*mi11,*ma11,tmp3,tmp4,
		minimaxFPAABCD,array5);
	
	FPA e1,e2;
	FPA_add_ZAB(e1,
		&ahx->val[_HELPER_Z6AZC_Xdep_mi4],
		mi9
	);
	FPA_add_ZAB(fA.x0,
		&e1,
		&ahy->val[_HELPER_Z6AZC_Ydep_C0re_minus_ma5_minus_ma7]
	);

	FPA_add_ZAB(e2,
		&ahx->val[_HELPER_Z6AZC_Xdep_ma4],
		ma9
	);
	FPA_add_ZAB(fA.x1,
		&e2,
		&ahy->val[_HELPER_Z6AZC_Ydep_C1re_minus_mi5_minus_mi7]
	);
	
	FPA_add_ZAB(e1,
		&ahx->val[_HELPER_Z6AZC_Xdep_C0im_plus_mi6],
		mi10
	);
	FPA_add_ZAB(fA.y0,&e1,mi12);

	FPA_add_ZAB(e2,
		&ahx->val[_HELPER_Z6AZC_Xdep_C1im_plus_ma6],
		ma10
	);
	FPA_add_ZAB(fA.y1,&e2,ma12);
		
	return;
	#endif

}

#define Z6AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=A.y0*A.y0;\
	NUMTYP y03=y02*A.y0;\
	ahy->val[_HELPER_Z6AZC_Ydep_y03]=y03;\
	NUMTYP y12=A.y1*A.y1;\
	NUMTYP y13=y12*A.y1;\
	NUMTYP y04=y02*y02;\
	NUMTYP y06=y04*y02;\
	NUMTYP y14=y12*y12;\
	NUMTYP y16=y14*y12;\
	ahy->val[_HELPER_Z6AZC_Ydep_y13]=y13;\
	NUMTYP mi5,ma5;\
	MINMAX2(mi5,ma5,AIM*A.y0,AIM*A.y1);\
	ahy->val[_HELPER_Z6AZC_Ydep_mi5]=mi5;\
	ahy->val[_HELPER_Z6AZC_Ydep_ma5]=ma5;\
	NUMTYP mi2,ma2;\
	MINMAX2(mi2,ma2,y02,y12);\
	ahy->val[_HELPER_Z6AZC_Ydep_mi2]=mi2;\
	ahy->val[_HELPER_Z6AZC_Ydep_ma2]=ma2;\
	NUMTYP mi7,ma7;\
	MINMAX2(mi7,ma7,y06,y16);\
	ahy->val[_HELPER_Z6AZC_Ydep_mi7]=mi7;\
	ahy->val[_HELPER_Z6AZC_Ydep_ma7]=ma7;\
	ahy->val[_HELPER_Z6AZC_Ydep_15_mul_mi2]=15*mi2;\
	ahy->val[_HELPER_Z6AZC_Ydep_15_mul_ma2]=15*ma2;\
	ahy->val[_HELPER_Z6AZC_Ydep_6_mul_mi2]=6*mi2;\
	ahy->val[_HELPER_Z6AZC_Ydep_6_mul_ma2]=6*ma2;\
	ahy->val[_HELPER_Z6AZC_Ydep_C0re_minus_ma5_minus_ma7]=(C0RE-ma5)-ma7;\
	ahy->val[_HELPER_Z6AZC_Ydep_C1re_minus_mi5_minus_mi7]=(C1RE-mi5)-mi7;\
	return;\
}

#define Z6AZCPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z6AZCPRECOMPUTEY_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeYdep_z6azc_double(
	PlaneRect_double& A,Helper_double* ahy
) {
	Z6AZCPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeYdep_z6azc(
	PlaneRect& A,Helper* ahy
) {
	#ifdef _DOUBLE
	Z6AZCPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _FPA
	Z6AZCPRECOMPUTEY(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z6AZCPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _QUADMATH
	Z6AZCPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _F107
	Z6AZCPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z6AZCPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
}

#define Z6AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=A.x0*A.x0;\
	NUMTYP x12=A.x1*A.x1;\
	NUMTYP mi1,ma1;\
	MINMAX2(mi1,ma1,x02,x12);\
	ahx->val[_HELPER_Z6AZC_Xdep_mi1]=mi1;\
	ahx->val[_HELPER_Z6AZC_Xdep_ma1]=ma1;\
	NUMTYP x03=x02*A.x0;\
	NUMTYP x13=x12*A.x1;\
	NUMTYP x05=x02*x03;\
	NUMTYP x15=x12*x13;\
	NUMTYP x04=x02*x02;\
	NUMTYP x14=x12*x12;\
	NUMTYP mi3,ma3;\
	MINMAX2(mi3,ma3,x04,x14);\
	ahx->val[_HELPER_Z6AZC_Xdep_mi3]=mi3;\
	ahx->val[_HELPER_Z6AZC_Xdep_ma3]=ma3;\
	NUMTYP mi6,ma6;\
	MINMAX2(mi6,ma6,AIM*A.x0,AIM*A.x1);\
	ahx->val[_HELPER_Z6AZC_Xdep_mi6]=mi6;\
	ahx->val[_HELPER_Z6AZC_Xdep_ma6]=ma6;\
	NUMTYP tmp1=(6*x05)+ARE;\
	NUMTYP tmp2=(6*x15)+ARE;\
	ahx->val[_HELPER_Z6AZC_Xdep_tmp1]=tmp1;\
	ahx->val[_HELPER_Z6AZC_Xdep_tmp2]=tmp2;\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,ARE*A.x0,ARE*A.x1);\
	ahx->val[_HELPER_Z6AZC_Xdep_mi4]=mi4;\
	ahx->val[_HELPER_Z6AZC_Xdep_ma4]=ma4;\
	ahx->val[_HELPER_Z6AZC_Xdep_20_mul_ma1]=20*ma1;\
	ahx->val[_HELPER_Z6AZC_Xdep_20_mul_mi1]=20*mi1;\
	ahx->val[_HELPER_Z6AZC_Xdep_C0im_plus_mi6]=C0IM+mi6;\
	ahx->val[_HELPER_Z6AZC_Xdep_C1im_plus_ma6]=C1IM+ma6;\
	return;\
}

#define Z6AZCPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z6AZCPRECOMPUTEX_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeXdep_z6azc_double(
	PlaneRect_double& A,Helper_double* ahx
) {
	Z6AZCPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeXdep_z6azc(
	PlaneRect& A,Helper* ahx
) {
	#ifdef _DOUBLE
	Z6AZCPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _FPA
	Z6AZCPRECOMPUTEX(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z6AZCPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _QUADMATH
	Z6AZCPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _F107
	Z6AZCPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z6AZCPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
}

void getBoundingBoxfA_z7azc_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=seedC0re_double-maximumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+minimumdouble(A.x0*(minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.x0*(maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))));
	fA.x1=seedC1re_double-minimumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+maximumdouble(A.x0*(minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.x0*(maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))));
	fA.y0=seedC0im_double+minimumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+minimumdouble(A.y0*(7*minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.y0*(7*maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))));
	fA.y1=seedC1im_double+maximumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+maximumdouble(A.y0*(7*minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.y0*(7*maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+minimumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre_double+maximumdouble(minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumdouble(A.x0*A.x0,A.x1*A.x1)))),maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumdouble(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumdouble(A.x0*A.x0,A.x1*A.x1)))))));
}

void getBoundingBoxfA_z7azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*(minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.x0*(maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))));
	fA.x1=seedC1re-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*(minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.x0*(maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.x1*(maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(7*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(5*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))));
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(A.y0*(7*minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.y0*(7*maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))));
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(A.y0*(7*minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.y0*(7*maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+minimumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))),A.y1*(7*maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)+FAKTORAre+maximumD(minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*minimumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*maximumD(A.x0*A.x0,A.x1*A.x1)))),maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))*(3*(7*maximumD(A.y0*A.y0,A.y1*A.y1))-(5*(7*minimumD(A.x0*A.x0,A.x1*A.x1)))))));
}

#define Z7AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=A.x0*A.x0;\
	NUMTYP x04=x02*x02;\
	NUMTYP x06=x04*x02;\
	NUMTYP x12=A.x1*A.x1;\
	NUMTYP x14=x12*x12;\
	NUMTYP x16=x14*x12;\
	NUMTYP mi1,ma1;\
	MINMAX2(mi1,ma1,x02,x12);\
	ahx->val[_HELPER_Z7AZC_Xdep_mi1]=mi1;\
	ahx->val[_HELPER_Z7AZC_Xdep_ma1]=ma1;\
	ahx->val[_HELPER_Z7AZC_Xdep_21mi1]=21*mi1;\
	ahx->val[_HELPER_Z7AZC_Xdep_21ma1]=21*ma1;\
	ahx->val[_HELPER_Z7AZC_Xdep_35mi1]=35*mi1;\
	ahx->val[_HELPER_Z7AZC_Xdep_35ma1]=35*ma1;\
	NUMTYP mi6,ma6;\
	MINMAX2(mi6,ma6,x06,x16);\
	ahx->val[_HELPER_Z7AZC_Xdep_mi6]=mi6;\
	ahx->val[_HELPER_Z7AZC_Xdep_ma6]=ma6;\
	ahx->val[_HELPER_Z7AZC_Xdep_7mi6]=7*mi6;\
	ahx->val[_HELPER_Z7AZC_Xdep_7ma6]=7*ma6;\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,AIM*A.x0,AIM*A.x1);\
	ahx->val[_HELPER_Z7AZC_Xdep_C0im_plus_mi4]=C0IM+mi4;\
	ahx->val[_HELPER_Z7AZC_Xdep_C1im_plus_ma4]=C1IM+ma4;\
	ahx->val[_HELPER_Z7AZC_Xdep_Are_plus_mi6]=ARE+mi6;\
	ahx->val[_HELPER_Z7AZC_Xdep_Are_plus_ma6]=ARE+ma6;\
	\
	return;\
}

#define Z7AZCPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z7AZCPRECOMPUTEX_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeXdep_z7azc_double(
	PlaneRect_double& A,Helper_double *ahx
) {
	Z7AZCPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	);
}

void precomputeXdep_z7azc(
	PlaneRect& A,Helper *ahx
) {
	#ifdef _F107
	Z7AZCPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _DOUBLE
	Z7AZCPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _QUADMATH
	Z7AZCPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _FPA
	Z7AZCPRECOMPUTEX(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z7AZCPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _F161
	Z7AZCPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
}

#define Z7AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=A.y0*A.y0;\
	NUMTYP y04=y02*y02;\
	NUMTYP y06=y04*y02;\
	NUMTYP y12=A.y1*A.y1;\
	NUMTYP y14=y12*y12;\
	NUMTYP y16=y14*y12;\
	NUMTYP mi2,ma2;\
	MINMAX2(mi2,ma2,y02,y12);\
	ahy->val[_HELPER_Z7AZC_Ydep_mi2]=mi2;\
	ahy->val[_HELPER_Z7AZC_Ydep_ma2]=ma2;\
	ahy->val[_HELPER_Z7AZC_Ydep_21mi2]=21*mi2;\
	ahy->val[_HELPER_Z7AZC_Ydep_21ma2]=21*ma2;\
	ahy->val[_HELPER_Z7AZC_Ydep_35mi2]=35*mi2;\
	ahy->val[_HELPER_Z7AZC_Ydep_35ma2]=35*ma2;\
	NUMTYP mi7,ma7;\
	MINMAX2(mi7,ma7,y06,y16);\
	ahy->val[_HELPER_Z7AZC_Ydep_mi7]=mi7;\
	ahy->val[_HELPER_Z7AZC_Ydep_ma7]=ma7;\
	ahy->val[_HELPER_Z7AZC_Ydep_7mi7]=7*mi7;\
	ahy->val[_HELPER_Z7AZC_Ydep_7ma7]=7*ma7;\
	NUMTYP mi3,ma3;\
	MINMAX2(mi3,ma3,AIM*A.y0,AIM*A.y1);\
	ahy->val[_HELPER_Z7AZC_Ydep_C0re_minus_ma3]=C0RE-ma3;\
	ahy->val[_HELPER_Z7AZC_Ydep_C1re_minus_mi3]=C1RE-mi3;\
	ahy->val[_HELPER_Z7AZC_Ydep_Are_minus_mi7]=ARE-mi7;\
	ahy->val[_HELPER_Z7AZC_Ydep_Are_minus_ma7]=ARE-ma7;\
	\
	return;\
}

#define Z7AZCPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z7AZCPRECOMPUTEY_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeYdep_z7azc_double(
	PlaneRect_double& A,Helper_double *ahy
) {
	Z7AZCPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void precomputeYdep_z7azc(
	PlaneRect& A,Helper *ahy
) {
	#ifdef _F107
	Z7AZCPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _DOUBLE
	Z7AZCPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _QUADMATH
	Z7AZCPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _FPA
	Z7AZCPRECOMPUTEY(NTYP,minimaxFPAAB,minimaxFPAABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z7AZCPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _F161
	Z7AZCPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
}

#define Z7AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP mi5,ma5;\
	IAMUL_SIGN(mi5,ma5,\
		ahx->val[_HELPER_Z7AZC_Xdep_mi1],\
		ahx->val[_HELPER_Z7AZC_Xdep_ma1],\
		ahy->val[_HELPER_Z7AZC_Ydep_mi2],\
		ahy->val[_HELPER_Z7AZC_Ydep_ma2],\
		MINMAX4);\
	NUMTYP mi8,ma8;\
	NUMTYP tmp1=(ahy->val[_HELPER_Z7AZC_Ydep_35mi2]-ahx->val[_HELPER_Z7AZC_Xdep_21ma1]);\
	NUMTYP tmp2=(ahy->val[_HELPER_Z7AZC_Ydep_35ma2]-ahx->val[_HELPER_Z7AZC_Xdep_21mi1]);\
	IAMUL_SIGN(mi8,ma8,mi5,ma5,tmp1,tmp2,MINMAX4)\
	NUMTYP mi9,ma9;\
	NUMTYP tmp3=\
		(ahx->val[_HELPER_Z7AZC_Xdep_Are_plus_mi6]\
		-ahy->val[_HELPER_Z7AZC_Ydep_7ma7]\
		)\
		+mi8;\
	NUMTYP tmp4=\
		(ahx->val[_HELPER_Z7AZC_Xdep_Are_plus_ma6]\
		-ahy->val[_HELPER_Z7AZC_Ydep_7mi7]\
		)\
		+ma8;\
	IAMUL_SIGN(mi9,ma9,A.x0,A.x1,tmp3,tmp4,MINMAX4)\
	NUMTYP mi10,ma10;\
	NUMTYP tmp5=(ahy->val[_HELPER_Z7AZC_Ydep_21mi2]-ahx->val[_HELPER_Z7AZC_Xdep_35ma1]);\
	NUMTYP tmp6=(ahy->val[_HELPER_Z7AZC_Ydep_21ma2]-ahx->val[_HELPER_Z7AZC_Xdep_35mi1]);\
	IAMUL_SIGN(mi10,ma10,mi5,ma5,tmp5,tmp6,MINMAX4)\
	NUMTYP mi11,ma11;\
	NUMTYP tmp7=\
		(ahx->val[_HELPER_Z7AZC_Xdep_7mi6]\
		+ahy->val[_HELPER_Z7AZC_Ydep_Are_minus_ma7]\
		)\
		+mi10;\
	NUMTYP tmp8=\
		(ahy->val[_HELPER_Z7AZC_Ydep_Are_minus_mi7]\
		+ahx->val[_HELPER_Z7AZC_Xdep_7ma6]\
		)\
		+ma10;\
	IAMUL_SIGN(mi11,ma11,A.y0,A.y1,tmp7,tmp8,MINMAX4)\
	\
	fA.x0=ahy->val[_HELPER_Z7AZC_Ydep_C0re_minus_ma3]+mi9;\
	fA.x1=ahy->val[_HELPER_Z7AZC_Ydep_C1re_minus_mi3]+ma9;\
	\
	fA.y0=ahx->val[_HELPER_Z7AZC_Xdep_C0im_plus_mi4]+mi11;\
	fA.y1=ahx->val[_HELPER_Z7AZC_Xdep_C1im_plus_ma4]+ma11;\
	\
	return;\
}

#define Z7AZCHELPER(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z7AZCHELPER_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void getBoundingBoxfA_z7azc_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double* ahx,Helper_double* ahy
) {
	Z7AZCHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void getBoundingBoxfA_z7azc_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper* ahx,Helper* ahy
) {
	ctrbbxfa++;
	
	#ifdef _F107
	Z7AZCHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _DOUBLE
	Z7AZCHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _QUADMATH
	Z7AZCHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _FPA
	PFPA mi5,ma5;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi5,ma5,
		ahx->val[_HELPER_Z7AZC_Xdep_mi1],
		ahx->val[_HELPER_Z7AZC_Xdep_ma1],
		ahy->val[_HELPER_Z7AZC_Ydep_mi2],
		ahy->val[_HELPER_Z7AZC_Ydep_ma2],
		minimaxFPAABCD,array1);
	FPA tmp1,tmp2;
	FPA_sub_ZAB(tmp1,
		&ahy->val[_HELPER_Z7AZC_Ydep_35mi2],
		&ahx->val[_HELPER_Z7AZC_Xdep_21ma1]
	);
	FPA_sub_ZAB(tmp2,
		&ahy->val[_HELPER_Z7AZC_Ydep_35ma2],
		&ahx->val[_HELPER_Z7AZC_Xdep_21mi1]
	);
	FPA array2[4];
	PFPA mi8,ma8;
	IAMUL_FPA_SIGN(mi8,ma8,
		*mi5,*ma5,tmp1,tmp2,
		minimaxFPAABCD,array2)
	FPA e1,e2,tmp3,tmp4;
	FPA_sub_ZAB(e1,
		&ahx->val[_HELPER_Z7AZC_Xdep_Are_plus_mi6],
		&ahy->val[_HELPER_Z7AZC_Ydep_7ma7]
	);
	FPA_add_ZAB(tmp3,&e1,mi8);
	FPA_sub_ZAB(e2,
		&ahx->val[_HELPER_Z7AZC_Xdep_Are_plus_ma6],
		&ahy->val[_HELPER_Z7AZC_Ydep_7mi7]
	);
	FPA_add_ZAB(tmp4,&e2,ma8);
	PFPA mi9,ma9;
	FPA array3[4];
	IAMUL_FPA_SIGN(mi9,ma9,
		A.x0,A.x1,tmp3,tmp4,
		minimaxFPAABCD,array3)
	FPA tmp5,tmp6;
	FPA_sub_ZAB(tmp5,
		&ahy->val[_HELPER_Z7AZC_Ydep_21mi2],
		&ahx->val[_HELPER_Z7AZC_Xdep_35ma1]
	);
	FPA_sub_ZAB(tmp6,
		&ahy->val[_HELPER_Z7AZC_Ydep_21ma2],
		&ahx->val[_HELPER_Z7AZC_Xdep_35mi1]
	);
	PFPA mi10,ma10;
	FPA array4[4];
	IAMUL_FPA_SIGN(mi10,ma10,
		*mi5,*ma5,tmp5,tmp6,
		minimaxFPAABCD,array4)
	FPA tmp7,tmp8;
	FPA_add_ZAB(e1,
		&ahx->val[_HELPER_Z7AZC_Xdep_7mi6],
		&ahy->val[_HELPER_Z7AZC_Ydep_Are_minus_ma7]
	);
	FPA_add_ZAB(tmp7,&e1,mi10);
	FPA_add_ZAB(e2,
		&ahy->val[_HELPER_Z7AZC_Ydep_Are_minus_mi7],
		&ahx->val[_HELPER_Z7AZC_Xdep_7ma6]
	);
	FPA_add_ZAB(tmp8,&e2,ma10);
	PFPA mi11,ma11;
	FPA array5[4];
	IAMUL_FPA_SIGN(mi11,ma11,
		A.y0,A.y1,tmp7,tmp8,
		minimaxFPAABCD,array5)
	
	FPA_add_ZAB(fA.x0,
		&ahy->val[_HELPER_Z7AZC_Ydep_C0re_minus_ma3],
		mi9
	);
	FPA_add_ZAB(fA.x1,
		&ahy->val[_HELPER_Z7AZC_Ydep_C1re_minus_mi3],
		ma9
	);
	
	FPA_add_ZAB(fA.y0,
		&ahx->val[_HELPER_Z7AZC_Xdep_C0im_plus_mi4],
		mi11
	);
	FPA_add_ZAB(fA.y1,
		&ahx->val[_HELPER_Z7AZC_Xdep_C1im_plus_ma4],
		ma11
	);
	
	return;
	#endif

	#ifdef _LONGDOUBLE
	Z7AZCHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _F161
	Z7AZCHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif
}

void getBoundingBoxfA_z8azc_double_oh(PlaneRect_double& A,PlaneRect_double& fA) {
	fA.x0=seedC0re_double+minimumdouble(FAKTORAre_double*A.x0,FAKTORAre_double*A.x1)-maximumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(4*(7*maximumdouble(minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))))+2*(5*(7*minimumdouble(minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))))-(4*(7*maximumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))))+minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.x1=seedC1re_double+maximumdouble(FAKTORAre_double*A.x0,FAKTORAre_double*A.x1)-minimumdouble(FAKTORAim_double*A.y0,FAKTORAim_double*A.y1)+maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(4*(7*minimumdouble(minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),minimumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0,A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0,A.y1*A.y1))))+2*(5*(7*maximumdouble(minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))))-(4*(7*minimumdouble(minimumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),minimumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*minimumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumdouble(A.x0*A.x0,A.x1*A.x1)*maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))))+maximumdouble(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.y0=minimumdouble(FAKTORAre_double*A.y0,FAKTORAre_double*A.y1)+minimumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+8*minimumdouble((A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(7*(8*maximumdouble((A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+7*(8*minimumdouble((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1*A.y1*A.y1)))-(8*maximumdouble(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)))+seedC0im_double;
	fA.y1=maximumdouble(FAKTORAre_double*A.y0,FAKTORAre_double*A.y1)+maximumdouble(FAKTORAim_double*A.x0,FAKTORAim_double*A.x1)+8*maximumdouble((A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(7*(8*minimumdouble((A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+7*(8*maximumdouble((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1*A.y1*A.y1)))-(8*minimumdouble(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)))+seedC1im_double;
}

void getBoundingBoxfA_z8azc(PlaneRect& A,PlaneRect& fA) {
	ctrbbxfa++;
	
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(4*(7*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+2*(5*(7*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))))-(4*(7*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))))+minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(4*(7*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+2*(5*(7*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))))-(4*(7*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1))))+maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+8*minimumD((A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(7*(8*maximumD((A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+7*(8*minimumD((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1*A.y1*A.y1)))-(8*maximumD(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)))+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+8*maximumD((A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(7*(8*minimumD((A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+7*(8*maximumD((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1*A.y1*A.y1)))-(8*minimumD(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1*A.y1*A.y1)))+seedC1im;
}

#define Z8AZCHELPER_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP mi9,ma9;\
	IAMUL_SIGN(mi9,ma9,\
		ahx->val[_HELPER_Z8AZC_Xdep_mi370],\
		ahx->val[_HELPER_Z8AZC_Xdep_ma370],\
		ahy->val[_HELPER_Z8AZC_Ydep_mi4],\
		ahy->val[_HELPER_Z8AZC_Ydep_ma4],\
		MINMAX4)\
	NUMTYP mi12,ma12;\
	IAMUL_SIGN(mi12,ma12,\
		ahx->val[_HELPER_Z8AZC_Xdep_mi11],\
		ahx->val[_HELPER_Z8AZC_Xdep_ma11],\
		ahy->val[_HELPER_Z8AZC_Ydep_mi228],\
		ahy->val[_HELPER_Z8AZC_Ydep_ma228],\
		MINMAX4)\
	NUMTYP mi14,ma14;\
	IAMUL_SIGN(mi14,ma14,\
		ahx->val[_HELPER_Z8AZC_Xdep_mi128],\
		ahx->val[_HELPER_Z8AZC_Xdep_ma128],\
		ahy->val[_HELPER_Z8AZC_Ydep_mi13],\
		ahy->val[_HELPER_Z8AZC_Ydep_ma13],\
		MINMAX4)\
	NUMTYP mi17,ma17;\
	IAMUL_SIGN(mi17,ma17,\
		ahx->val[_HELPER_Z8AZC_Xdep_t1a],\
		ahx->val[_HELPER_Z8AZC_Xdep_t1b],\
		A.y0,A.y1,MINMAX4)\
	NUMTYP mi18,ma18;\
	IAMUL_SIGN(mi18,ma18,\
		ahx->val[_HELPER_Z8AZC_Xdep_t2a],\
		ahx->val[_HELPER_Z8AZC_Xdep_t2b],\
		ahy->val[_HELPER_Z8AZC_Ydep_y03],\
		ahy->val[_HELPER_Z8AZC_Ydep_y13],\
		MINMAX4)\
	NUMTYP mi19,ma19;\
	IAMUL_SIGN(mi19,ma19,\
		ahx->val[_HELPER_Z8AZC_Xdep_t3a],\
		ahx->val[_HELPER_Z8AZC_Xdep_t3b],\
		ahy->val[_HELPER_Z8AZC_Ydep_y05],\
		ahy->val[_HELPER_Z8AZC_Ydep_y15],\
		MINMAX4);\
	NUMTYP mi20,ma20;\
	IAMUL_SIGN(mi20,ma20,\
		ahy->val[_HELPER_Z8AZC_Ydep_t4a],\
		ahy->val[_HELPER_Z8AZC_Ydep_t4b],\
		A.x0,A.x1,MINMAX4)\
	\
	fA.x0=\
		((ahx->val[_HELPER_Z8AZC_Xdep_C0re_plus_mi5_plus_mi10]\
		-ma14)\
		+(mi9-ma12))\
		+ahy->val[_HELPER_Z8AZC_Ydep_mi15_minus_ma6];\
	fA.x1=\
		((ahx->val[_HELPER_Z8AZC_Xdep_C1re_plus_ma5_plus_ma10]\
		-mi12)\
		+(ma9-mi14))\
		+ahy->val[_HELPER_Z8AZC_Ydep_ma15_minus_mi6];\
	\
	fA.y0=\
		(ahy->val[_HELPER_Z8AZC_Ydep_mi7]\
		+ahx->val[_HELPER_Z8AZC_Xdep_C0im_plus_mi8]\
		)\
		+((mi17-ma18)+(mi19-ma20));\
	fA.y1=\
		(ahy->val[_HELPER_Z8AZC_Ydep_ma7]\
		+ahx->val[_HELPER_Z8AZC_Xdep_C1im_plus_ma8])\
		+((ma17-mi18)+(ma19-mi20));\
	\
	return;\
}

#define Z8AZCHELPER(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z8AZCHELPER_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void getBoundingBoxfA_z8azc_helper(
	PlaneRect& A,PlaneRect& fA,
	Helper* ahx,Helper* ahy
) {
	ctrbbxfa++;

	#ifdef _DOUBLE
	Z8AZCHELPER(NTYP,minimaxdAB,minimaxdABCD)
	#endif
	
	#ifdef _LONGDOUBLE
	Z8AZCHELPER(NTYP,minimaxldAB,minimaxldABCD)
	#endif
	
	#ifdef _F107
	Z8AZCHELPER(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif
	
	#ifdef _F161
	Z8AZCHELPER(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif

	#ifdef _QUADMATH
	Z8AZCHELPER(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif

	#ifdef _FPA
	// optimierte Aufrufroutine
	PFPA mi9,ma9;
	FPA array1[4];
	IAMUL_FPA_SIGN(mi9,ma9,
		ahx->val[_HELPER_Z8AZC_Xdep_mi370],
		ahx->val[_HELPER_Z8AZC_Xdep_ma370],
		ahy->val[_HELPER_Z8AZC_Ydep_mi4],
		ahy->val[_HELPER_Z8AZC_Ydep_ma4],
		minimaxFPAABCD,array1);
	PFPA mi12,ma12;
	FPA array2[4];
	IAMUL_FPA_SIGN(mi12,ma12,
		ahx->val[_HELPER_Z8AZC_Xdep_mi11],
		ahx->val[_HELPER_Z8AZC_Xdep_ma11],
		ahy->val[_HELPER_Z8AZC_Ydep_mi228],
		ahy->val[_HELPER_Z8AZC_Ydep_ma228],
		minimaxFPAABCD,array2);
	PFPA mi14,ma14;
	FPA array3[4];
	IAMUL_FPA_SIGN(mi14,ma14,
		ahx->val[_HELPER_Z8AZC_Xdep_mi128],
		ahx->val[_HELPER_Z8AZC_Xdep_ma128],
		ahy->val[_HELPER_Z8AZC_Ydep_mi13],
		ahy->val[_HELPER_Z8AZC_Ydep_ma13],
		minimaxFPAABCD,array3);
	PFPA mi17,ma17;
	FPA array4[4];
	IAMUL_FPA_SIGN(mi17,ma17,
		ahx->val[_HELPER_Z8AZC_Xdep_t1a],
		ahx->val[_HELPER_Z8AZC_Xdep_t1b],
		A.y0,A.y1,minimaxFPAABCD,array4);
	PFPA mi18,ma18;
	FPA array5[4];
	IAMUL_FPA_SIGN(mi18,ma18,
		ahx->val[_HELPER_Z8AZC_Xdep_t2a],
		ahx->val[_HELPER_Z8AZC_Xdep_t2b],
		ahy->val[_HELPER_Z8AZC_Ydep_y03],
		ahy->val[_HELPER_Z8AZC_Ydep_y13],
		minimaxFPAABCD,array5);
	PFPA mi19,ma19;
	FPA array6[4];
	IAMUL_FPA_SIGN(mi19,ma19,
		ahx->val[_HELPER_Z8AZC_Xdep_t3a],
		ahx->val[_HELPER_Z8AZC_Xdep_t3b],
		ahy->val[_HELPER_Z8AZC_Ydep_y05],
		ahy->val[_HELPER_Z8AZC_Ydep_y15],
		minimaxFPAABCD,array6);
	PFPA mi20,ma20;
	FPA array7[4];
	IAMUL_FPA_SIGN(mi20,ma20,
		ahy->val[_HELPER_Z8AZC_Ydep_t4a],
		ahy->val[_HELPER_Z8AZC_Ydep_t4b],
		A.x0,A.x1,minimaxFPAABCD,array7);
	
	FPA e1,e2,e3;
	FPA_sub_ZAB(e1,
		&ahx->val[_HELPER_Z8AZC_Xdep_C0re_plus_mi5_plus_mi10],
		ma14);
	FPA_sub_ZAB(e2,mi9,ma12);
	FPA_add_ZAB(e3,e1,e2);
	FPA_add_ZAB(fA.x0,&e3,&ahy->val[_HELPER_Z8AZC_Ydep_mi15_minus_ma6]);

	FPA f1,f2,f3;
	FPA_sub_ZAB(f1,
		&ahx->val[_HELPER_Z8AZC_Xdep_C1re_plus_ma5_plus_ma10],
		mi14);
	FPA_sub_ZAB(f2,ma9,mi12);
	FPA_add_ZAB(f3,f1,f2);
	FPA_add_ZAB(fA.x1,&f3,&ahy->val[_HELPER_Z8AZC_Ydep_ma15_minus_mi6]);
	
	FPA_add_ZAB(e1,
		&ahy->val[_HELPER_Z8AZC_Ydep_mi7],
		&ahx->val[_HELPER_Z8AZC_Xdep_C0im_plus_mi8]);
	FPA_sub_ZAB(e2,mi17,ma18);
	FPA_sub_ZAB(e3,mi19,ma20);
	FPA e4;
	FPA_add_ZAB(e4,&e2,&e3);
	FPA_add_ZAB(fA.y0,&e1,&e4);

	FPA_add_ZAB(f1,
		&ahy->val[_HELPER_Z8AZC_Ydep_ma7],
		&ahx->val[_HELPER_Z8AZC_Xdep_C1im_plus_ma8]);
	FPA_sub_ZAB(f2,ma17,mi18);
	FPA_sub_ZAB(f3,ma19,mi20);
	FPA f4;
	FPA_add_ZAB(f4,&f2,&f3);
	FPA_add_ZAB(fA.y1,&f1,&f4);
	
	return;
	#endif
}

#define Z8AZCPRECOMPUTEY_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP y02=A.y0*A.y0;\
	NUMTYP y03=y02*A.y0;\
	NUMTYP y04=y02*y02;\
	NUMTYP y05=y02*y03;\
	NUMTYP y06=y04*y02;\
	NUMTYP y12=A.y1*A.y1;\
	NUMTYP y13=y12*A.y1;\
	NUMTYP y14=y12*y12;\
	NUMTYP y15=y12*y13;\
	NUMTYP y16=y14*y12;\
	NUMTYP tmi2,tma2;\
	MINMAX2(tmi2,tma2,y02,y12);\
	NUMTYP mi228,ma228;\
	mi228=28*tmi2;\
	ma228=28*tma2;\
	NUMTYP mi4,ma4;\
	MINMAX2(mi4,ma4,y04,y14);\
	NUMTYP mi6,ma6;\
	MINMAX2(mi6,ma6,AIM*A.y0,AIM*A.y1);\
	NUMTYP mi7,ma7;\
	MINMAX2(mi7,ma7,ARE*A.y0,ARE*A.y1);\
	NUMTYP mi13,ma13;\
	MINMAX2(mi13,ma13,y06,y16);\
	NUMTYP mi15,ma15;\
	MINMAX2(mi15,ma15,y04*y04,y14*y14);\
	NUMTYP t4a=8*y04*y03;\
	NUMTYP t4b=8*y14*y13;\
	ahy->val[_HELPER_Z8AZC_Ydep_mi4]=mi4;\
	ahy->val[_HELPER_Z8AZC_Ydep_ma4]=ma4;\
	ahy->val[_HELPER_Z8AZC_Ydep_mi228]=mi228;\
	ahy->val[_HELPER_Z8AZC_Ydep_ma228]=ma228;\
	ahy->val[_HELPER_Z8AZC_Ydep_mi13]=mi13;\
	ahy->val[_HELPER_Z8AZC_Ydep_ma13]=ma13;\
	ahy->val[_HELPER_Z8AZC_Ydep_y03]=y03;\
	ahy->val[_HELPER_Z8AZC_Ydep_y13]=y13;\
	ahy->val[_HELPER_Z8AZC_Ydep_y05]=y05;\
	ahy->val[_HELPER_Z8AZC_Ydep_y15]=y15;\
	ahy->val[_HELPER_Z8AZC_Ydep_t4a]=t4a;\
	ahy->val[_HELPER_Z8AZC_Ydep_t4b]=t4b;\
	ahy->val[_HELPER_Z8AZC_Ydep_mi7]=mi7;\
	ahy->val[_HELPER_Z8AZC_Ydep_ma7]=ma7;\
	ahy->val[_HELPER_Z8AZC_Ydep_ma15_minus_mi6]=ma15-mi6;\
	ahy->val[_HELPER_Z8AZC_Ydep_mi15_minus_ma6]=mi15-ma6;\
	return;\
}

#define Z8AZCPRECOMPUTEY(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z8AZCPRECOMPUTEY_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeYdep_z8azc(
	PlaneRect& A,Helper *ahy
) {
	#ifdef _FPA
	Z8AZCPRECOMPUTEY(FPA,minimaxFPAAB,minimaxFPAABCD)
	#endif

	#ifdef _F107
	Z8AZCPRECOMPUTEY(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z8AZCPRECOMPUTEY(NTYP,minimaxF161AB,minimaxF107ABCD)
	#endif

	#ifdef _DOUBLE
	Z8AZCPRECOMPUTEY(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z8AZCPRECOMPUTEY(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _QUADMATH
	Z8AZCPRECOMPUTEY(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	LOGMSG("Missing implementation. z8azc precompute Ydep helper\n");
	exit(99);
}

#define Z8AZCPRECOMPUTEX_VARLIST(NUMTYP,MINMAX2,MINMAX4,ARE,AIM,C0RE,C1RE,C0IM,C1IM) \
{\
	NUMTYP x02=A.x0*A.x0;\
	NUMTYP x03=x02*A.x0;\
	NUMTYP x04=x02*x02;\
	NUMTYP x05=x02*x03;\
	NUMTYP x06=x04*x02;\
	NUMTYP x12=A.x1*A.x1;\
	NUMTYP x13=x12*A.x1;\
	NUMTYP x14=x12*x12;\
	NUMTYP x15=x12*x13;\
	NUMTYP x16=x14*x12;\
	NUMTYP tmi1,tma1;\
	MINMAX2(tmi1,tma1,x02,x12);\
	NUMTYP mi128,ma128;\
	mi128=28*tmi1;\
	ma128=28*tma1;\
	ahx->val[_HELPER_Z8AZC_Xdep_mi128]=mi128;\
	ahx->val[_HELPER_Z8AZC_Xdep_ma128]=ma128;\
	NUMTYP tmi3,tma3;\
	MINMAX2(tmi3,tma3,x04,x14);\
	ahx->val[_HELPER_Z8AZC_Xdep_mi370]=70*tmi3;\
	ahx->val[_HELPER_Z8AZC_Xdep_ma370]=70*tma3;\
	NUMTYP mi5,ma5;\
	MINMAX2(mi5,ma5,ARE*A.x0,ARE*A.x1);\
	NUMTYP mi8,ma8;\
	MINMAX2(mi8,ma8,AIM*A.x0,AIM*A.x1);\
	NUMTYP mi10,ma10;\
	MINMAX2(mi10,ma10,x04*x04,x14*x14);\
	NUMTYP mi11,ma11;\
	MINMAX2(mi11,ma11,x06,x16);\
	ahx->val[_HELPER_Z8AZC_Xdep_mi11]=mi11;\
	ahx->val[_HELPER_Z8AZC_Xdep_ma11]=ma11;\
	NUMTYP t1a=8*(x04*x03);\
	NUMTYP t1b=8*(x14*x13);\
	NUMTYP t2a=56*x05;\
	NUMTYP t2b=56*x15;\
	NUMTYP t3a=56*x03;\
	NUMTYP t3b=56*x13;\
	ahx->val[_HELPER_Z8AZC_Xdep_t1a]=t1a;\
	ahx->val[_HELPER_Z8AZC_Xdep_t1b]=t1b;\
	ahx->val[_HELPER_Z8AZC_Xdep_t2a]=t2a;\
	ahx->val[_HELPER_Z8AZC_Xdep_t2b]=t2b;\
	ahx->val[_HELPER_Z8AZC_Xdep_t3a]=t3a;\
	ahx->val[_HELPER_Z8AZC_Xdep_t3b]=t3b;\
	ahx->val[_HELPER_Z8AZC_Xdep_C0im_plus_mi8]=C0IM+mi8;\
	ahx->val[_HELPER_Z8AZC_Xdep_C1im_plus_ma8]=C1IM+ma8;\
	ahx->val[_HELPER_Z8AZC_Xdep_C0re_plus_mi5_plus_mi10]=(C0RE+mi5)+mi10;\
	ahx->val[_HELPER_Z8AZC_Xdep_C1re_plus_ma5_plus_ma10]=(C1RE+ma5)+ma10;\
	return;\
}

#define Z8AZCPRECOMPUTEX(NUMTYP,MINMAX2,MINMAX4) \
{\
	Z8AZCPRECOMPUTEX_VARLIST(\
		NUMTYP,MINMAX2,MINMAX4,\
		FAKTORAre,FAKTORAim,\
		seedC0re,seedC1re,\
		seedC0im,seedC1im\
	)\
}

void precomputeXdep_z8azc(
	PlaneRect& A,Helper *ahx
) {
	#ifdef _FPA
	Z8AZCPRECOMPUTEX(FPA,minimaxFPAAB,minimaxFPAABCD)
	#endif

	#ifdef _F107
	Z8AZCPRECOMPUTEX(NTYP,minimaxF107AB,minimaxF107ABCD)
	#endif

	#ifdef _F161
	Z8AZCPRECOMPUTEX(NTYP,minimaxF161AB,minimaxF161ABCD)
	#endif

	#ifdef _DOUBLE
	Z8AZCPRECOMPUTEX(NTYP,minimaxdAB,minimaxdABCD)
	#endif

	#ifdef _LONGDOUBLE
	Z8AZCPRECOMPUTEX(NTYP,minimaxldAB,minimaxldABCD)
	#endif

	#ifdef _QUADMATH
	Z8AZCPRECOMPUTEX(NTYP,minimaxQDAB,minimaxQDABCD)
	#endif
	
	LOGMSG("Missing implementation. z8azc precompute Xdep helper\n");
	exit(99);
}

void precomputeXdep_z8azc_double(
	PlaneRect_double& A,Helper_double *ahx
) {
	Z8AZCPRECOMPUTEX_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	);
}

void precomputeYdep_z8azc_double(
	PlaneRect_double& A,Helper_double *ahy
) {
	Z8AZCPRECOMPUTEY_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
}

void getBoundingBoxfA_z8azc_double(
	PlaneRect_double& A,PlaneRect_double& fA,
	Helper_double* ahx,Helper_double* ahy
) {
	Z8AZCHELPER_VARLIST(
		double,minimaxdAB,minimaxdABCD,
		FAKTORAre_double,FAKTORAim_double,
		seedC0re_double,seedC1re_double,
		seedC0im_double,seedC1im_double
	)
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
		
		data5->revcgYX[i].containsgray=0;
	}
	
	// check which revcg vertices may contains gray
	int32_t parentx,parenty;
	for(int32_t y=0;y<SCREENWIDTH;y+=REVCGBLOCKWIDTH) {
		parenty=(y >> REVCGBITS);
		int32_t poffsety=(int64_t)parenty*REVCGmaxnumber;

		for(int32_t x=0;x<SCREENWIDTH;x+=REVCGBLOCKWIDTH) {
			parentx=(x >> REVCGBITS);
			
			// do the coordinates of that revcg vertex
			// lie in memgrau.g0..g1 of row y
			int hasgray=0;
			
			for(int y2=y;y2<(y+REVCGBLOCKWIDTH);y2++) {
				if (!data5->zeilen[y2]) continue;
			
				// revcg block starting at x must OVERLAP with
				// gray in row, but not necessarily be contained
				// fully within
				// [x..x+REVCGBLOCKWIDTH-1] geschniten [g0..g1]
				int32_t xe=x+REVCGBLOCKWIDTH-1;
				
				if (
					(xe < data5->memgrau[y2].g0) ||
					(x > data5->memgrau[y2].g1)
				) {
					// not overlapping
				} 
				else {
					hasgray=1;
					break;
				}
			}
			
			if (hasgray>0) {
				data5->revcgYX[poffsety+parentx].containsgray=1;
			} else {
				data5->revcgYX[poffsety+parentx].containsgray=0;
			}
		}
	} // y
	
	PlaneRect A,bbxfA;
	const NTYP DD=REVCGBLOCKWIDTH*scaleRangePerPixel;
	
	for(int32_t dl=1;dl<=2;dl++) {
		A.y1=COMPLETE0;
		// first pass: calculuate how many are nedded per square
		// 2nd pass: build cell graph as array allocating exact amount of memory
		if (dl==1) printf("\ncounting parents ...");
		else printf("\nsetting parents to squares ... ");
		
		for(int32_t y=0;y<SCREENWIDTH;y+=REVCGBLOCKWIDTH) {
			parenty=(y >> REVCGBITS);
			#ifdef _FPA
			A.y0=A.y1;
			FPA_add_ZAB(A.y1,A.y0,DD);
			#else
			A.y0=y*scaleRangePerPixel + COMPLETE0;
			A.y1=A.y0 + DD;
			#endif
		
			A.x1=COMPLETE0;
			for(int32_t x=0;x<SCREENWIDTH;x+=REVCGBLOCKWIDTH) {
				parentx=(x >> REVCGBITS);
				#ifdef _FPA
				A.x0=A.x1;
				FPA_add_ZAB(A.x1,A.x1,DD);
				#else
				A.x0=x*scaleRangePerPixel + COMPLETE0;
				A.x1=A.x0+DD;
				#endif
			
				// no use of helper object here
				// as size of A16 is NOT scaleRangePerPixel
				getBoundingBoxfA(A,bbxfA);
			
				if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA) > 0) {
					continue;
				}

				ScreenRect scr;
				// truncated, outside part is not relevant
				// for rev. cel graph
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
					for(int32_t by=scr.y0;by<=scr.y1;by++) {
						int64_t yoffset=(int64_t)by*REVCGmaxnumber;
						for(int32_t bx=scr.x0;bx<=scr.x1;bx++) {
							if (data5->revcgYX[yoffset+bx].containsgray>0) {
								data5->revcgYX[yoffset+bx].memused++;
							}
						}
					} // by
				} else {
					// setting parents to square
					for(int32_t by=scr.y0;by<=scr.y1;by++) {
						int64_t yoffset=(int64_t)by*REVCGmaxnumber;
						for(int32_t bx=scr.x0;bx<=scr.x1;bx++) {
							// if the revcg vertex by,bx contains gray
							// => a parent needs to be set
							if (data5->revcgYX[yoffset+bx].containsgray>0) {
								data5->revcgYX[yoffset+bx].addParent(parentx,parenty);
							}
						}
					} // by
				}
			} // x
		} // y
	} // passes
}

void compute(void) {
	if (
		(_PROPAGATEPOTW>0) ||
		(_PROPAGATEDEF>0) 
	) {
		construct_static_reverse_cellgraph();
	}
	
	// propagating definite colors
	if (_PROPAGATEDEF>0) {
		propagate_definite();
	} 
		
	// propagate potentially white
	// takes most of time - especially if no black
	// will finally emerge
	if (_PROPAGATEPOTW>0) {
		if (_PROPAGATEDEF>0) {
			printf("saving raw data ... ");
			data5->saveRaw("_temp");
			printf("done\n");
		}
		propagate_potw();
		printf("\nsearching for interior cells ... ");
		int32_t res=color_changeS32(
			SQUARE_GRAY,
			SQUARE_BLACK,
			SQUARE_GRAY_16_CONSECUTIVE,
			SQUARE_BLACK_16_CONSECUTIVE
		);
		if (res>0) interiorpresent=1;
	} else {
		printf("\nskipping interior coloring (potw not propagated)\n");
	}
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
	
	int32_t noch0=(SCREENWIDTH >> 4) >> 2;
	int32_t noch=1;
	const NTYP DD16SCALE=16.0*scaleRangePerPixel;

	A16.y1=COMPLETE0; 
	for(int32_t y16=0;y16<SCREENWIDTH;y16+=16) {
		A16.y0=A16.y1;
		#ifdef _FPA
		FPA_add_ZAB(A16.y1,A16.y0,DD16SCALE);
		#else
		A16.y1=A16.y1+DD16SCALE;
		#endif
		if ( (--noch)==0) {
			printf("%i ",SCREENWIDTH-y16);
			noch=noch0;
		}
		int32_t gray0=SCREENWIDTH-1,gray1=0;

		A16.x1=COMPLETE0; 
		for(int32_t x16=0;x16<SCREENWIDTH;x16+=16) {
			A16.x0=A16.x1;
			#ifdef _FPA
			FPA_add_ZAB(A16.x1,A16.x0,DD16SCALE);
			#else
			A16.x1=A16.x1+DD16SCALE;
			#endif
					
			// no use of helper object here
			// as size of A16 is NOT scaleRangePerPixel
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
			if (gray1 >= gray0) {
				data5->graudensity[yy]=(uint8_t)(100.0*(double)(gray1-gray0+1)/(double)SCREENWIDTH);
			}  else {
				data5->graudensity[yy]=0;
			}

			// whole row was allocated outside in calling routine
			// => mem spans the entire pixel coordinate range
			data5->memgrau[yy].mem0=0;
			data5->memgrau[yy].mem1=(SCREENWIDTH >> 4)-1;
		} // yy
	} // y16
	
	// adjusting image gray enclosement
	// one 16-block left/riht/up/down as buffer
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
	
	// check the border of the data: no gray cell
	// must touch the image's boundary
	// if so => larger RANGE value is necessary
	
	int8_t touches=0;
	
	#define CHECK(XX0,XX1,YY0,YY1) \
	{\
		for(int32_t y=(YY0);y<=(YY1);y++) {\
			for(int32_t x=(XX0);x<=(XX1);x++) {\
				int32_t f;\
				GET_SINGLE_CELLCOLOR_XY(x,y,f)\
				if (f != SQUARE_WHITE) {\
					touches=1;\
					break;\
				}\
			} \
			if (touches>0) break;\
		} \
	}
	
	if (touches<=0) {
		CHECK(0,0,0,SCREENWIDTH-1)
	}
	if (touches<=0) {
		CHECK(SCREENWIDTH-1,SCREENWIDTH-1,0,SCREENWIDTH-1)
	}
	if (touches<=0) {
		CHECK(0,SCREENWIDTH-1,0,0)
	}
	if (touches<=0) {
		CHECK(0,SCREENWIDTH-1,SCREENWIDTH-1,SCREENWIDTH-1)
	}
	
	if (touches>0) {
		LOGMSG("\n\nGray region touches image border which is not possible in the current implementation.\nIncreased RANGE value is recommended.\n");
		exit(99);
	}
}

void propagate_definite(void) {
	PlaneRect A,bbxfA;
	ScreenRect scr;
	int32_t noch0=((SCREENWIDTH >> REVCGBITS) >> 1);
	int32_t noch=6;
	int8_t changed=1;
	int32_t lastsavetime=0;
	int64_t checkclockat=ctrbbxfa+checkclockatbbxcount0;
	
	// mark every tile as to be visited
	for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
		data5->revcgYX[i].tovisit=1;
		data5->revcgYX[i].containsgray=1;
	}
	
	changed=1;
	
	while (changed>0) {
		changed=0;
		printf("\npropagating definite color ... ");

		for(int32_t y256=0,YBLOCK=0;y256<SCREENWIDTH;y256+=REVCGBLOCKWIDTH,YBLOCK++) {
			int32_t yrevoffset=YBLOCK*REVCGmaxnumber;
			if ( (--noch) <= 0) {
				printf("%i ",SCREENWIDTH-y256);
				noch=noch0;
			}
			
			if (ctrbbxfa > checkclockat) {
				checkclockat += checkclockatbbxadd;
				int t2=clock();
				if ((t2-lastsavetime) > CLOCKHOURSTOSAVE) {
					printf("saving raw data ... ");
					data5->saveRaw("_temp");
					printf("done\n");
					lastsavetime=t2;
				}
			}
		
			// outside gray enclosement => jump out
			if ( (y256+REVCGBLOCKWIDTH) < encgrayy0) continue;
			if (y256 > encgrayy1) break;
	
			int32_t lastBX=-1; // last vertex set to VISITED-AGAIN
			int32_t GLOBALBY=y256 >> REVCGBITS;
			int32_t GLOBALBYOFFSET=GLOBALBY*REVCGmaxnumber;

			for(int32_t x256=0,XBLOCK=0;x256<SCREENWIDTH;x256+=REVCGBLOCKWIDTH,XBLOCK++) {
				if ( (data5->revcgYX[yrevoffset+XBLOCK].tovisit<=0) ) {
					continue;
				}
			
				// visit, but is there still gray
				if (data5->revcgYX[yrevoffset+XBLOCK].containsgray<=0) {
					continue;
				}
				
				// block has now been checked
				data5->revcgYX[yrevoffset+XBLOCK].tovisit=0;
	
				const int32_t Y256ENDE=y256+REVCGBLOCKWIDTH;
			
				// ATTN: even if bbxprecomputed
				// A.y HAS TO BE calculated
				// as it is INCREMENTED in every yloop
				#ifdef _FPA
				FPA_mul_ZAuvlong(A.y1,scaleRangePerPixel,y256);
				FPA_add_ZAB(A.y1,A.y1,COMPLETE0);
				#else
				A.y1=y256*scaleRangePerPixel + COMPLETE0;
				#endif
				
				int8_t blockhasgray=0;
				
				for(int32_t y=y256;y<Y256ENDE;y++) {
					A.y0=A.y1;
					#ifdef _FPA
					FPA_add_ZAB(A.y1,A.y0,scaleRangePerPixel);
					#else
					A.y1=A.y0+scaleRangePerPixel;
					#endif

					int8_t bbxprecomputedrow=0;
					if (data5->pcscr) {
						if (data5->pcscr[y]) bbxprecomputedrow=1;
					}

					const int32_t xanf=data5->memgrau[y].g0;
					const int32_t xende=data5->memgrau[y].g1;
					if (xende < xanf) continue;
			
					// gray in that row lies outside of area 
					// of this revcg vertex ?
					if ( 
						(x256 > xende) ||
						( (x256+REVCGBLOCKWIDTH) < xanf) 
					) continue;
				
					Helper* helperY=helperYdep->getHelper(y);
					int32_t wmem=-1 + (x256 >> 4);
		
					for(int32_t x=x256;x<(x256+REVCGBLOCKWIDTH);x+=16) {
						wmem++;
						uint32_t w;
						GETDATA5BYMEM_MY(wmem,y,w)
			
						// no gray square in this 32-bit integer
						if (
							(w == SQUARE_WHITE_16_CONSECUTIVE) ||
							(w == SQUARE_BLACK_16_CONSECUTIVE) ||
							(w == SQUARE_GRAYPOTW_16_CONSECUTIVE) 
						) continue; 
			
						uint32_t wneu=w;
						int32_t w_changed=0;
						if (bbxprecomputedrow<=0) {
							#ifdef _FPA
							FPA tmp;
							FPA_mul_ZAuvlong(tmp,scaleRangePerPixel,x);
							FPA_add_ZAB(A.x1,tmp,COMPLETE0);
							#else
							A.x1=x*scaleRangePerPixel + COMPLETE0;
							#endif
						} 
						
						for(int32_t wbith=0;wbith<16;wbith++) {
							if (bbxprecomputedrow<=0) {
								A.x0=A.x1;
								#ifdef _FPA
								FPA_add_ZAB(A.x1,A.x0,scaleRangePerPixel);
								#else
								A.x1=A.x0+scaleRangePerPixel;
								#endif
							}
						
							uint32_t globalf=w & 0b11;
							w >>= 2;
				
							if (globalf != SQUARE_GRAY) continue;
							
							blockhasgray=1;
							int8_t hits_white=0;
							int8_t hits_black=0;
						
							if (bbxprecomputedrow>0) {
								// bbx already computed befprehand
								// bounding box in special exterior
								GETPCSCR(x+wbith,y,scr)
							
								if (scr.x0<0) {
									// partially or fully outside GRAY ENCLOSMEENT
									if (scr.x1<0) {
										// fully
										wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
										w_changed=1;
										continue;
									} else {
										// partially
										hits_white=1;
										scr.x0 = (-scr.x0)-1; // make it valid agfain
									}
								}
							} else {
								getBoundingBoxfA_helper(
									A,bbxfA,
									helperXdep->getHelper(x+wbith),
									helperY
								);
		
								// bounding box in special exterior
								if ((SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(bbxfA))>0) {
									if (bbxprecomputedrow>0) {
										// scr ist der gespeicherte Wert
										if (scr.x1>=0) {
											LOGMSG("Implementation error def1\n");
											exit(99);
										}
									}
									wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
									w_changed=1;
									continue;
								}
					
								if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
									// overlaps with white region
									hits_white=1;
									if (bbxprecomputedrow>0) {
										// scr ist der gespeicherte Wert
										if (scr.x0>=0) {
											LOGMSG("Implementation error. def2\n");
											exit(99);
										}
									}
								}
					
								scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
								scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
								scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
								scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
							} // newly computed screenrect done
						
							for(int32_t ty=scr.y0;ty<=scr.y1;ty++) {
								for(int32_t tx=scr.x0;tx<=scr.x1;tx++) {
									int32_t f;
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
							SETDATA5BYMEM_MY(wmem,y,wneu)
							int32_t BX=x >> REVCGBITS;
							if (lastBX != BX) {
								int32_t offset=GLOBALBYOFFSET+BX;
								for(int32_t i=0;i<data5->revcgYX[offset].howmany;i++) {
									SETTOVISIT(
										data5->revcgYX[offset].parent[i].BX,
										data5->revcgYX[offset].parent[i].BY
									)
								}
								lastBX=BX;
							}
							changed=1;
						}
					} // x
				} // y
				
				if (blockhasgray<=0) {
					data5->revcgYX[yrevoffset+XBLOCK].containsgray=0;
				}
			} // X256
		} // Y256
	} // while
}

void propagate_potw(void) {
	PlaneRect A,bbxfA;
	ScreenRect scr;

	int8_t changed=1;
	int32_t lastsavetime=0;
	int64_t checkclockat=checkclockatbbxcount0;
	int32_t noch0=(SCREENWIDTH >> 4) >> 2;
	int32_t noch=1;
	
	// mark every tile as to be visited
	for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
		data5->revcgYX[i].tovisit=1;
		data5->revcgYX[i].containsgray=1;
	}

	changed=1;
	
	while (changed>0) {
		changed=0;
		printf("\npropagating potentially white ... ");
	
		for(int32_t y256=0,YBLOCK=0;y256<SCREENWIDTH;y256+=REVCGBLOCKWIDTH,YBLOCK++) {
			int32_t yrevoffset=YBLOCK*REVCGmaxnumber;
			if ( (--noch) <= 0) {
				printf("%i ",SCREENWIDTH-y256);
				noch=noch0;
			}
			if (ctrbbxfa > checkclockat) {
				checkclockat += checkclockatbbxadd;
				int32_t t2=clock();
				if ((t2-lastsavetime) > CLOCKHOURSTOSAVE) {
					printf("saving raw data ... ");
					data5->saveRaw("_temp");
					printf("done\n");
					lastsavetime=t2;
				}
			}

			// outside gray enclosement => jump out
			if ( (y256+REVCGBLOCKWIDTH) < encgrayy0) continue;
			if (y256 > encgrayy1) break;
	
			int32_t lastBX=-1;
			int32_t GLOBALBY=y256 >> REVCGBITS;
			int32_t GLOBALBYOFFSET=GLOBALBY*REVCGmaxnumber;

			for(int32_t x256=0,XBLOCK=0;x256<SCREENWIDTH;x256+=REVCGBLOCKWIDTH,XBLOCK++) {
				if ( (data5->revcgYX[yrevoffset+XBLOCK].tovisit<=0) ) {
					continue;
				}
				
				if (data5->revcgYX[yrevoffset+XBLOCK].containsgray<=0) {
					continue;
				}
		
				// block has now been checked
				data5->revcgYX[yrevoffset+XBLOCK].tovisit=0;
	
				const int32_t Y256ENDE=y256+REVCGBLOCKWIDTH;
				#ifdef _FPA
				FPA_mul_ZAuvlong(A.y1,scaleRangePerPixel,y256);
				FPA_add_ZAB(A.y1,A.y1,COMPLETE0);
				#else
				A.y1=y256*scaleRangePerPixel + COMPLETE0;
				#endif
				
				int8_t blockhasgray=0;
				
				for(int32_t y=y256;y<Y256ENDE;y++) {
					A.y0=A.y1;
					#ifdef _FPA
					FPA_add_ZAB(A.y1,A.y0,scaleRangePerPixel);
					#else
					A.y1=A.y0+scaleRangePerPixel;
					#endif
					
					int8_t bbxprecomputedrow=0;
					if (data5->pcscr) {
						if (data5->pcscr[y]) bbxprecomputedrow=1;
					}
					const int32_t xanf=data5->memgrau[y].g0;
					const int32_t xende=data5->memgrau[y].g1;
					if (xende < xanf) continue;
		
					// gray in that row lies outside of area to be checked
					if ( 
						(x256 > xende) ||
						( (x256+REVCGBLOCKWIDTH) < xanf) 
					) continue;
				
					Helper *helperY=helperYdep->getHelper(y);
					int32_t wmem=-1 + (x256 >> 4);
		
					for(int32_t x=x256;x<(x256+REVCGBLOCKWIDTH);x+=16) {
						wmem++;
						uint32_t w;
						GETDATA5BYMEM_MY(wmem,y,w)
			
						// no gray square in this 32-bit integer
						if (
							(w == SQUARE_WHITE_16_CONSECUTIVE) ||
							(w == SQUARE_BLACK_16_CONSECUTIVE) ||
							(w == SQUARE_GRAYPOTW_16_CONSECUTIVE) 
						) continue; 
			
						uint32_t wneu=w;
						int32_t w_changed=0;
						if (bbxprecomputedrow<=0) {
							#ifdef _FPA
							FPA_mul_ZAuvlong(A.x1,scaleRangePerPixel,x);
							FPA_add_ZAB(A.x1,A.x1,COMPLETE0);
							#else
							A.x1=x*scaleRangePerPixel + COMPLETE0;
							#endif
						}
					
						for(int32_t wbith=0;wbith<16;wbith++) {
							// continuous adding here even if
							// pixel is gray (continue)
							if (bbxprecomputedrow<=0) {
								A.x0=A.x1;
								#ifdef _FPA
								FPA_add_ZAB(A.x1,A.x0,scaleRangePerPixel);
								#else
								A.x1=A.x0+scaleRangePerPixel;
								#endif
							} 
						
							uint32_t globalf=w & 0b11;
							w >>= 2;
				
							// current pixel is already gray potw
							if (globalf != SQUARE_GRAY) continue;

							blockhasgray=1;
							int32_t pathtowhite=0;
						
							if (bbxprecomputedrow>0) {
								GETPCSCR(x+wbith,y,scr)
								if (scr.x0<0) {
									if (scr.x1<0) {
										// fully outside GRAY ENCLOSMEENT
										wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
										w_changed=1;
										continue;
									} else {
										// partially outside
										pathtowhite=1;
										scr.x0 = (-scr.x0)-1; // valid
										// scr is entirely in screen 
									}
								} 
							} else {
								getBoundingBoxfA_helper(
									A,bbxfA,
									helperXdep->getHelper(x+wbith),
									helperY
								);
			
								// bounding box in special exterior
								if ((SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(bbxfA))>0) {
									if (bbxprecomputedrow>0) {
										// scr ist der gespeicherte Wert
										if (scr.x1>=0) {
											LOGMSG("Implementation error. potw1\n");
											exit(99);
										}
									}
									wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
									w_changed=1;
									continue;
								}
				
								if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
									// at least overlaps with white region
									// if completely outisde => would've been set to white in _def
									pathtowhite=1;
									if (bbxprecomputedrow>0) {
										// scr ist der gespeicherte Wert
										if (scr.x0>=0) {
											LOGMSG("Implementation error. potw/2\n");
											exit(99);
										}
									}
								} 
								
								scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
								scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
								scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
								scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
							} // newly computed bbx
						
							if (pathtowhite <= 0) {
								for(int32_t ty=scr.y0;ty<=scr.y1;ty++) {
									for(int32_t tx=scr.x0;tx<=scr.x1;tx++) {
										int32_t f;
										GET_SINGLE_CELLCOLOR_XY(tx,ty,f);
										if (
											(f==SQUARE_WHITE) ||
											(f==SQUARE_GRAY_POTENTIALLY_WHITE)
										) {
											pathtowhite=1;
											break;
										}
									
										if (pathtowhite>0) break;
									} // tx
							
									if (pathtowhite>0) break;
								} // ty
							}
				
							if (pathtowhite>0) {
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_GRAYPOTW[wbith]);
								w_changed=1;
							}
						} // wbith
			
						if (w_changed>0) {
							SETDATA5BYMEM_MY(wmem,y,wneu)
							int32_t BX=x >> REVCGBITS;
							if (lastBX != BX) {
								int32_t offset=GLOBALBYOFFSET+BX;
								for(int32_t i=0;i<data5->revcgYX[offset].howmany;i++) {
									SETTOVISIT(
										data5->revcgYX[offset].parent[i].BX,
										data5->revcgYX[offset].parent[i].BY
									)
								}
								lastBX=BX;
							} 
							changed=1;
						}
					} // x
				} // y
				
				if (blockhasgray<=0) {
					data5->revcgYX[yrevoffset+XBLOCK].containsgray=0;
				}
			} // X256
		} // Y256
	} // while
}

int color_changeS32(
	const DDBYTE source1,const DDBYTE target1,
	const DDBYTE source16,const DDBYTE target16
) {
	int32_t res=0;
	int32_t noch0=SCREENWIDTH >> 3;
	int32_t noch=1;
	
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
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

// struct RefList

RefPoint* RefList::getRefPtr(const int32_t ax) {
	// list is sorted increasingly, so binary search is possible
	int32_t left=0,right=anz-1;
	
	while ( (left<(right-5)) ) {
		int32_t m=(int32_t)(((int64_t)right + left) >> 1);
		if (m<left) m=left;
		if (m>right) m=right;
		
		if (points[m].x == ax) return &points[m];
		
		if (points[m].x < ax) {
			left=m+1;
		} else {
			right=m-1;
		}
	}
	
	for(int i=left;i<=right;i++) {
		if (points[i].x == ax) {
			return &points[i];
		}
	}
	
	return NULL;
}

void RefList::addXB(const int32_t ax,const int32_t ablobid) {
	if (anz >= memused) {
		LOGMSG("Implementation error RefList. Too many XB values.\n");
		exit(99);
	}
	points[anz].x=ax;
	points[anz].blobid=ablobid;
	anz++;
}

RefList::RefList() {
	anz=0;
	points=NULL;
	memused=0;
}

// ParentManager
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
	#else
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
	#else
	sprintf(erg,"A_%.20lg_%.20lg",
		(double)FAKTORAre,(double)FAKTORAim);
	#endif
	
	return erg;
}

char* seedCstr225(char* erg) {
	#ifdef _FPA
	char t1[1024],t2[1024],t3[1024],t4[1024];
	sprintf(erg,"c_ia_%s_%s_x_%s_%s",
		seedC0re.str(t1),
		seedC1re.str(t2),
		seedC0im.str(t3),
		seedC1im.str(t4)
	);
	#else
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
	sprintf(erg,"A_%s_%s",
		FAKTORAre.str(t1),
		FAKTORAim.str(t2)
	);
	#else
	sprintf(erg,"A_%I64d_%I64d",
		(int64_t)floor(DENOM225*(double)FAKTORAre),
		(int64_t)floor(DENOM225*(double)FAKTORAim)
	);
	#endif
		
	return erg;
}

void periodicity_m3(char* afn) {
	printf("initializing Fatou component search ... ");
	// change every graypotw to gray
	// graypotw gis later used as flag "visited"
	color_changeS32(
		SQUARE_GRAY_POTENTIALLY_WHITE,
		SQUARE_GRAY,
		SQUARE_GRAYPOTW_16_CONSECUTIVE,
		SQUARE_GRAY_16_CONSECUTIVE
	);
	
	uint32_t SQUARE_VISITED=SQUARE_GRAY_POTENTIALLY_WHITE;
	
	// build hashed list of reference points
	// going once over the image, checking black
	// If there's no black direct south neighbour
	// => found new reference point
	// the list is then automatically increasingly sorted
	
	RefPointArray *refpoints=new RefPointArray;
	
	int32_t noch0=SCREENWIDTH >> 3;
	int32_t noch=1;
	RefPointManager *rpmgr=new RefPointManager;

	printf("\nsearching for reference points ... ");
	int64_t refpointsctr=0;
	// interior is and remains BLACK colored throughout this step
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		refpoints->listY[y].anz=0;
		refpoints->listY[y].points=NULL;
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
			noch=noch0;
		}
		
		if (!data5->zeilen) continue;
		if (data5->memgrau[y].g1 < data5->memgrau[y].g0) continue;
	
		// counting in a first pass
		int32_t localctr=0;
		for(int x=data5->memgrau[y].g0;x<=data5->memgrau[y].g1;x++) {
			int32_t f;
			GET_SINGLE_CELLCOLOR_XY(x,y,f)
			if (f != SQUARE_BLACK) continue;
			
			if (y==0) {
				LOGMSG("Error. m3/1\n");
				exit(99);
			}
			
			int32_t fsouth;
			GET_SINGLE_CELLCOLOR_XY(x,y-1,fsouth);
			if (fsouth == SQUARE_BLACK) continue;
			
			// new reference point found
			localctr++;
		} // x
		
		// allocate
		if (localctr <= 0) continue;
		refpointsctr += localctr;

		refpoints->listY[y].points=rpmgr->getMemory(localctr);
		refpoints->listY[y].anz=0;
		refpoints->listY[y].memused=localctr;
		
		// second pass: add ref points
		for(int x=data5->memgrau[y].g0;x<=data5->memgrau[y].g1;x++) {
			int32_t f;
			GET_SINGLE_CELLCOLOR_XY(x,y,f)
			if (f != SQUARE_BLACK) continue;
			
			if (y==0) {
				LOGMSG("Error/2 m3.\n");
				exit(99);
			}
			
			int32_t fsouth;
			GET_SINGLE_CELLCOLOR_XY(x,y-1,fsouth);
			if (fsouth == SQUARE_BLACK) continue;
			
			// new reference point found (x,y)
			refpoints->addRefPoint(x,y,0);
		} // x
	} // y
	
	printf("\n  %I64d reference points identified\n",refpointsctr);
	// go over all ref points and determine their
	// blobid number
	
	printf("searching for Fatou components ... ");
	int32_t nextblobid=1; // not zero
	int32_t MAXBLOBID=(UINT32MAX >> 1);
	
	noch=1;
	StreakArray *streaks=new StreakArray;
	
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
			noch=noch0;
		}
		
		if (
			(refpoints->listY[y].anz==0) ||
			(!refpoints->listY[y].points)
		) continue;
		
		int32_t rpy=y;
		for(int32_t rp=0;rp<refpoints->listY[y].anz;rp++) {
			// already identified
			if (refpoints->listY[y].points[rp].blobid>0) continue;
			
			int32_t rpx=refpoints->listY[y].points[rp].x;
			
			// floodFill blob
			refpoints->listY[y].points[rp].blobid=nextblobid;
			streaks->fastEmpty();
			
			int32_t xx0,xx1;

			#define SEARCHSTREAK_AND_MARK_REFP(BLOBID,REFX,REFY,ERGX0,ERGX1) \
			{\
				ERGX0=REFX;\
				while (ERGX0>0) {\
					int32_t f;\
					GET_SINGLE_CELLCOLOR_XY(ERGX0,REFY,f);\
					if (\
						(ERGX0 != (REFX)) &&\
						(f != SQUARE_BLACK) \
					) {\
						ERGX0++;\
						break;\
					}\
					SET_SINGLE_CELLCOLOR_XYFF(ERGX0,REFY,SQUARE_VISITED);\
					if (REFY>0) {\
						GET_SINGLE_CELLCOLOR_XY(ERGX0,REFY-1,f)\
						if (f == SQUARE_GRAY) {\
							RefPoint *p=refpoints->getRefPtr(ERGX0,REFY);\
							if (p) {\
								if ( (p->blobid>0) && (p->blobid != (BLOBID)) ) {\
									LOGMSG3("Error/3. blobid inconsistency: has %i, set to %i\n",\
										p->blobid,BLOBID);\
									exit(99);\
								} else {\
									p->blobid=BLOBID;\
								}\
							} else {\
								LOGMSG3("Error/M3. refpoint ptr null %i,%i\n",\
									ERGX0,REFY);\
								exit(99);\
							}\
						}\
					} else {\
						LOGMSG("Error/periodM3. Not able to analyze reference points in bottom row.\n");\
							exit(99);\
					}\
					ERGX0--;\
				}\
				\
				ERGX1=REFX;\
				while (ERGX1<SCREENWIDTH) {\
					int f;\
					GET_SINGLE_CELLCOLOR_XY(ERGX1,REFY,f);\
					if (\
						(ERGX1 != (REFX) ) &&\
						(f != SQUARE_BLACK) \
					) {\
						ERGX1--;\
						break;\
					}\
					SET_SINGLE_CELLCOLOR_XYFF(ERGX1,REFY,SQUARE_VISITED);\
					if (REFY>0) {\
						GET_SINGLE_CELLCOLOR_XY(ERGX1,REFY-1,f)\
						if (f == SQUARE_GRAY) {\
							RefPoint *p=refpoints->getRefPtr(ERGX1,REFY);\
							if (p) {\
								if ( (p->blobid>0) && (p->blobid != (BLOBID)) ) {\
									LOGMSG3("Error/3. blobid inconsistency: has %i, set to %i\n",\
										p->blobid,BLOBID);\
									exit(99);\
								} else {\
									p->blobid=BLOBID;\
								}\
							} else {\
								LOGMSG3("Error/M3. refpoint ptr null %i,%i\n",\
									ERGX1,REFY);\
								exit(99);\
							}\
						}\
					} else {\
						LOGMSG("Error/periodM3. Not able to analyze reference points in bottom row.\n");\
							exit(99);\
					}\
					ERGX1++;\
				}\
			}

			// rpx,rpy might already be VISITED => no following here
			SEARCHSTREAK_AND_MARK_REFP(nextblobid,rpx,rpy,xx0,xx1)
			
			streaks->pushStreak(xx0,xx1,rpy);
			
			Streak str;
			while (1) {
				streaks->popStreak(str);
				if (str.y<0) break; // stack empty
				
				// use leftmost point
				int actx=str.x0;
				int acty=str.y;
				str.x0++;
				if (str.x0 <= str.x1) {
					// pushback streak shortend by one at the left
					streaks->pushStreak(str.x0,str.x1,str.y);
				}

				// go one row up: if pixel black
				// => new streak pushed
				
				int32_t newx0,newx1;
				if (acty>0) {
					int32_t f;
					GET_SINGLE_CELLCOLOR_XY(actx,acty-1,f)
					if (f == SQUARE_BLACK) {
						SEARCHSTREAK_AND_MARK_REFP(
							nextblobid,
							actx,acty-1,
							newx0,newx1
						);
						streaks->pushStreak(newx0,newx1,acty-1);
					}
				} else {
					LOGMSG("Periodicity/m3. error. streak at y=0 not possible\n");
					exit(99);
				}
				
				// go one row down: if pixel black
				if (acty < (SCREENWIDTH-1) ) {
					int32_t f;
					GET_SINGLE_CELLCOLOR_XY(actx,acty+1,f)
					if (f == SQUARE_BLACK) {
						SEARCHSTREAK_AND_MARK_REFP(
							nextblobid,
							actx,acty+1,
							newx0,newx1
						);
						streaks->pushStreak(newx0,newx1,acty+1);
					}
				} else {
					LOGMSG("Periodicity/m3. error. streak at topmost row not possible\n");
					exit(99);
				}
				// => new streak pushed
			} // while streaks
			
			nextblobid++;
			if (nextblobid > (MAXBLOBID-8)) {
				LOGMSG("Error. Too many blobs.\n");
				exit(99);
			}
		} // rp
	} // y
	
	printf("\n  %i Fatou components found\n",nextblobid);
	
	// now go over all blob ids and construct their orbit
	// it always ends in a cycle
	// if it is new => add
	printf("searching cycles ... ");

	int32_t *oneorbit=new int32_t[M3MAXORBITLEN];
	int32_t orbitlen=0;
	int8_t *blobvisited=new int8_t[nextblobid];
	for(int32_t i=0;i<nextblobid;i++) blobvisited[i]=0;
	
	#define FREEMEM \
	{\
		delete[] oneorbit;\
		delete[] cycles;\
		delete[] blobvisited;\
		delete refpoints;\
		delete rpmgr;\
		delete perpoimgr;\
	}

	noch=1;
	noch0=SCREENWIDTH >> 3;
	uint8_t CYCLECOLOROFFSET=16;
	Palette4 periodpal;
	periodpal.setPaletteRGB(SQUARE_BLACK,0,0,0);
	periodpal.setPaletteRGB(SQUARE_WHITE,255,255,255);
	periodpal.setPaletteRGB(SQUARE_GRAY,127,127,127);
	periodpal.setPaletteRGB(SQUARE_GRAY_POTENTIALLY_WHITE,255,0,0);
	// some shuffled version of a heat-map
	double d=0.0;
	double dst=0.19;
	for(int32_t i=CYCLECOLOROFFSET;i<=255;i++) {
		int32_t r,g,b;
		basinpal.getColor(d,r,g,b);
		periodpal.setPaletteRGB(i,r,g,b);
		d += dst;
		while (d >= 1.0) d -= 1.0;
	}

	#define SETRGB(IDX,RR,GG,BB) \
	{\
		periodpal.setPaletteRGB(IDX,RR,GG,BB);\
	}
	SETRGB(0+CYCLECOLOROFFSET,0,255,255)
	SETRGB(1+CYCLECOLOROFFSET,255,0,255)
	SETRGB(2+CYCLECOLOROFFSET,255,0,0)
	SETRGB(3+CYCLECOLOROFFSET,0,255,0)
	SETRGB(4+CYCLECOLOROFFSET,255,255,0)
	SETRGB(5+CYCLECOLOROFFSET,193,193,255)
	SETRGB(6+CYCLECOLOROFFSET,193,63,255)
	
	int32_t anzcycles=0;
	CycleM3 *cycles=new CycleM3[M3MAXCYCLES];
	IntManager *perpoimgr=new IntManager;
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",y);
			noch=noch0;
		}

		for(int32_t rp=0;rp<refpoints->listY[y].anz;rp++) {
			if (refpoints->listY[y].points[rp].blobid <= 0) {
				LOGMSG("Error. Period/M3 at item4: blobid not determined.\n");
				exit(99);
			}
			
			RefPoint* pcurrblob=&refpoints->listY[y].points[rp];
			
			if (blobvisited[pcurrblob->blobid] > 0) continue;
			
			// new orbit
			orbitlen=1;
			int32_t currentx=pcurrblob->x;;
			int32_t currenty=y;
			oneorbit[0]=pcurrblob->blobid;
			blobvisited[pcurrblob->blobid]=1;
			
			int32_t orbit0=-1,orbit1=-1;
			PlaneRect A,bbxfA;
			int8_t skipit=0;
			while (1) {
				// take orbit end's blobid ref point
				A.x0=currentx*scaleRangePerPixel+ COMPLETE0;
				A.x1=A.x0+scaleRangePerPixel;
				A.y0=currenty*scaleRangePerPixel+ COMPLETE0;
				A.y1=A.y0+scaleRangePerPixel;
				
				getBoundingBoxfA_helper(
					A,bbxfA,
					helperXdep->getHelper(currentx),
					helperYdep->getHelper(currenty)
				);
				
				// is completely in interior, so no check necessary
				ScreenRect ascr;
				ascr.x0=scrcoord_as_lowerleft(bbxfA.x0);
				ascr.x1=scrcoord_as_lowerleft(bbxfA.x1);
				ascr.y0=scrcoord_as_lowerleft(bbxfA.y0);
				ascr.y1=scrcoord_as_lowerleft(bbxfA.y1);

				int32_t rx=ascr.x0;
				int32_t ry=ascr.y0;

				// bbx is in one fatou component
				// so just take the lower left corner and find
				// its blobid as the next orbit component
				
				// go down until ref point is found
				while (ry>0) {
					int32_t f;
					GET_SINGLE_CELLCOLOR_XY(rx,ry,f);
					if (f == SQUARE_GRAY) {
						ry++; // ref point coordinates found
						break;
					}
					ry--;
				} // while
				
				RefPoint *ptarget=refpoints->getRefPtr(rx,ry);
				if (
					(!ptarget) ||
					(ptarget->blobid <= 0)
				) {
					LOGMSG3("Error. Period/M3. Target ref point %i,%i not found\n",rx,ry);
					exit(99);
				}
				
				currentx=rx;
				currenty=ry;
				oneorbit[orbitlen]=ptarget->blobid;
				orbitlen++;
				
				if (orbitlen >= (M3MAXORBITLEN-8)) {
					LOGMSG("Orbit too long. Periodicity skipped.\n");
					FREEMEM
					return;
				}
				
				// is a cycle just reached
				// i.e [orbitlen-1] already present earlier
				// can only be a cycle if the last blob
				// has already been visited
				
				if (blobvisited[ptarget->blobid]<=0) {
					// mark now as visited
					blobvisited[ptarget->blobid]=1;				
					continue;
				}
				
				// is it in a cycle
				for(int o=(orbitlen-2);o>=0;o--) {
					if (oneorbit[o] == oneorbit[orbitlen-1]) {
						// cycle found
						orbit0=o;
						orbit1=orbitlen-1;
						break;
					}
				} // o
				
				if (orbit0>=0) break;
				
				// no cycle in current orbit, but
				// that means the last blob has been visited in
				// a previous orbit, which already
				// has been followed to its cycle
				// so this orbit can be skipped as it
				// does not land in a new cycle
				
				skipit=1;
				break;
			} // while
			
			if (skipit>0) continue; // next rp
			
			if (orbit0<0) {
				LOGMSG("Error. Period/M3 no orbit found\n");
				exit(99);
			}
			
			// [o0].blobid==[o1].blobid, d.h. laenge ist diff
			// new cycle
			int32_t found=-1;
			
			for(int32_t cyc=0;cyc<anzcycles;cyc++) {
				for(int32_t k=0;k<cycles[cyc].len;k++) {
					if (cycles[cyc].perblobs[k] == oneorbit[orbitlen-1]) {
						found=cyc;
						break;
					}
				}
				if (found>=0) break;
			} // cyc
			
			if (found<0) {
				printf(" \n  cycle len=%i found\n",orbit1-orbit0);
				cycles[anzcycles].len=orbit1-orbit0;
				cycles[anzcycles].color=CYCLECOLOROFFSET+anzcycles;
				cycles[anzcycles].perblobs=perpoimgr->getMemory(cycles[anzcycles].len);
				for(int32_t o=orbit0;o<orbit1;o++) {
					cycles[anzcycles].perblobs[o-orbit0]=oneorbit[o];
				}
				anzcycles++;
				if (anzcycles >= (M3MAXCYCLES-8)) {
					LOGMSG("Too many cycles detect. Periodicity check skipped.\n");
					FREEMEM
					return;
				}
			}
		} // rp
	} // y
	
	printf("\n%i cycles detected\n",anzcycles);
	for(int32_t cyc=0;cyc<anzcycles;cyc++) {
		LOGMSG6("  cycle #%i len=%i immediate RGB(%i,%i,%i)\n",
			cyc,cycles[cyc].len,
			periodpal.rgbs[cycles[cyc].color].R,
			periodpal.rgbs[cycles[cyc].color].G,
			periodpal.rgbs[cycles[cyc].color].B);
	}
	
	// construct image
	int32_t TWD=0;
	// adjust exponentn so final image is at most 2^16 x 2^16
	while ( (SCREENWIDTH >> TWD) > 65536) TWD++;
	
	char tmp[1024];
	int32_t TWDSTEP=1 << TWD;
	int32_t bytes_per_row = SCREENWIDTH >> TWD;
	uint8_t* rgbz=new uint8_t[bytes_per_row];
	uint8_t* lastrgbz=new uint8_t[bytes_per_row];
	for(int32_t i=0;i<bytes_per_row;i++) {
		rgbz[i]=lastrgbz[i]=SQUARE_GRAY;
	}
	uint32_t off
		=	14 // size of file header
		+	40 // size of bitmap header
		+	256*4; // palette entries
	uint32_t filelen
		=	off
		+	(bytes_per_row*bytes_per_row);
	
	sprintf(tmp,"%s_period_2_%i-fold.bmp",afn,TWD);
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
	for(int32_t i=0;i<256;i++) write4(fbmp,periodpal.rgbs[i].B,periodpal.rgbs[i].G,periodpal.rgbs[i].R,periodpal.rgbs[i].alpha);
	
	int32_t lastblob=-1;
	noch=1;
	noch0=(SCREENWIDTH >> 3) >> TWD;
	printf("saving image ... ");
	
	#define GETREFPOINT(XSTART,YSTART,ERGPTR) \
	{\
		int32_t rx=XSTART;\
		int32_t ry=YSTART;\
		\
		while (ry>0) {\
			int32_t f;\
			GET_SINGLE_CELLCOLOR_XY(rx,ry,f)\
			if (f==SQUARE_GRAY) {\
				ry++;\
				break;\
			}\
			ry--;\
		}\
		\
		if (ry <= 0) {\
			LOGMSG("PeriodM3. save image, reference point below bottom.\n");\
			exit(99);\
		}\
		ERGPTR=refpoints->getRefPtr(rx,ry);\
	}
	
	// POTW color here means VISITED
	for(int32_t y=0;y<SCREENWIDTH;y+=TWDSTEP) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
			noch=noch0;
		}
		int32_t xwrite=-1;
		lastblob=-1;
		for(int32_t x=0;x<SCREENWIDTH;x+=TWDSTEP) {	
			xwrite++;
			
			int32_t finalf=-1;
			for(int32_t dy=0;dy<TWDSTEP;dy++) {
				for(int32_t dx=0;dx<TWDSTEP;dx++) {
					int32_t ftwd;
					GET_SINGLE_CELLCOLOR_XY(x+dx,y+dy,ftwd);
					if (ftwd==SQUARE_GRAY) {
						finalf=SQUARE_GRAY;
						break;
					}
					if (finalf<0) finalf=ftwd;
					else if (finalf != ftwd) {
						finalf=SQUARE_GRAY;
						break;
					}
				}
				if (finalf==SQUARE_GRAY) break;
			}
			
			if (finalf != SQUARE_VISITED) {
				lastblob=-1;
			} else {
				finalf=-1;

				if (y>0) {
					if (lastrgbz[xwrite] == SQUARE_BLACK) {
						finalf=SQUARE_BLACK; // connected
					} else if (lastrgbz[xwrite] >= CYCLECOLOROFFSET) {
						finalf=lastrgbz[xwrite]; // connected
					}
				}
				
				if (finalf<0) {
					if (lastblob<0) {
						// finding ONE reference point and its
						// blob is sufficient
						
						RefPoint* ptr;
						GETREFPOINT(x,y,ptr)
	
						if (!ptr) {
							LOGMSG("PeriodM3. save image, reference point not found as pointer.\n");
							exit(99);
						}
						
						lastblob=ptr->blobid;
					} // interior pixel in image
					
					// find cycle of that blob
					for(int32_t cyc=0;cyc<anzcycles;cyc++) {
						for(int32_t k=0;k<cycles[cyc].len;k++) {
							if (lastblob == cycles[cyc].perblobs[k]) {
								finalf=cycles[cyc].color;
								break;
							}
						}
					}
					if (finalf<0) {
						// non-immediate => black
						finalf=SQUARE_BLACK;
					}
				} // finalf to dtermine
			}
			
			rgbz[xwrite]=finalf;
		} // x
		fwrite(rgbz,sizeof(BYTE),bytes_per_row,fbmp);
		// swap rgbz pointers
		PBYTE tmpp=rgbz;
		rgbz=lastrgbz;
		lastrgbz=tmpp;
	} // y
	
	fclose(fbmp);
	FREEMEM
}

// another routine, memory-inefficient but can detect
// attraction basins and periodic points
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

	printf("converting data to periodicity structure ... ");
	int64_t memused=0;
	int32_t noch0=SCREENWIDTH >> 3;
	int32_t noch=1;
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
			noch=noch0;
		}
		if (data5->memgrau[y].g1 < data5->memgrau[y].g0) {
			dbY[y]=NULL;
		} else {
			int32_t m0=(data5->memgrau[y].g0 >> 4);
			int32_t m1=( (data5->memgrau[y].g1 >> 4) + 1);
			data5->memgrau[y].g0 = (m0 << 4);
			data5->memgrau[y].g1 = ( ((m1+1) << 4) - 1);
			
			int32_t xlen=data5->memgrau[y].g1-data5->memgrau[y].g0+1;
			memused += (xlen*sizeof(DBYTE));

			if (xlen>0) {
				dbY[y]=mgr->getMemory(xlen);
				if (!dbY[y]) {
					LOGMSG("Error/periodicity.\n");
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
	} // y
	printf("\nperiodicity memory used %I64d GB\n",1+(memused >> 30));
	
	FatouComponent *oneorbit=new FatouComponent[MAXFATOUCOMPONENTS];
	int32_t anzfatouinorbit=0;
	ibfcomponents=new FatouComponent[MAXFATOUCOMPONENTS];
	anzibf=0;
	cycles=new Cycle[MAXCYCLES];
	anzcycles=0;
	
	const int32_t MINTEMPCOLOR=256;
	int32_t orbitfcnbr=MINTEMPCOLOR; // >= 256
	int32_t cyclesetnbrimmediate=FATOUCOMPONENTCOLOROFFSET;
	int32_t cyclesetnbrattraction=FATOUCOMPONENTCOLOROFFSET+1;
	DBYTE blobaktiv=FATOUCOMPONENTCOLOROFFSET-1;
	
	noch0=SCREENWIDTH >> 2;
	noch=1;
	const int64_t MAXINLISTE=( (int64_t)1 << 26); // a 1 GB Speicher
	Int2 *liste=new Int2[MAXINLISTE];
	int64_t anzliste=0;
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
	
	SETCYCLE(0,0,255,255)
	SETCYCLE(1,255,0,255)
	SETCYCLE(2,255,0,0)
	SETCYCLE(3,0,255,0)
	SETCYCLE(4,255,255,0)
	SETCYCLE(5,193,193,255)
	SETCYCLE(6,193,63,255)
	
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
		for(int32_t xx=(BX-1);xx>=0;xx--) {\
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
					
					getBoundingBoxfA_helper(
						A,bbxfA,
						helperXdep->getHelper(x),
						helperYdep->getHelper(y)
					);
					
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
	int64_t totalsize=(int64_t)SCREENWIDTH * SCREENWIDTH;
	int64_t fourgb=( (int64_t )1 << 32 );
	while (totalsize > fourgb) {
		DSEXPONENT++;
		totalsize >>= 1;
	}
	int32_t TWDSTEP=(1 << DSEXPONENT); // n x n pixels are trustworthily downscaled
	int32_t dslen=(SCREENWIDTH >> DSEXPONENT);

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
			=	14 // FILEHeader
			+	40 // Bitmapheader
			+	256*4; // ColorPalette
	unsigned int filelen
			=	off
			+	(ybytes*dslen);

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
	
	delete[] liste;
	liste=NULL;
	anzliste=0;

	anzliste=0;
	ListeDFS orbit;
	ListeFIFO possibleper; 
	int32_t ppx,ppy;
	ScreenRect* ppscr=new ScreenRect[MAXPERIODICPOINTS];
	int32_t rausausfc=0;
	
	for(int32_t cyc=0;cyc<anzcycles;cyc++) {
		printf("\ncycle #%i periodic start ",cyc);
		int32_t PRELEN=-1 + cycles[cyc].len;
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
		int32_t noch0=(ibfcomponents[fc].scrc.y1 - ibfcomponents[fc].scrc.y0) >> 3;
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
					
					getBoundingBoxfA_helper(
						A,bbxfA,
						helperXdep->getHelper(ox),
						helperYdep->getHelper(oy)
					);
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
						
			getBoundingBoxfA_helper(
				A,bbxfA,
				helperXdep->getHelper(wx),
				helperYdep->getHelper(wy)
			);

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
			#else
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
	
	// BLACK = periodic point
	if (TWDSTEP>1) {
		sprintf(tt,"%s_periodic_points_twd_%i_fold.bmp",afn,DSEXPONENT);
	} else {
		sprintf(tt,"%s_periodic_points.bmp",afn);
	}
	
	if (DSEXPONENT>0) printf("  trustworthily downscaled 2^%i-fold\n",DSEXPONENT);

	fbmp=fopen(tt,"wb");
	BMPHEADER1(fbmp)
	
	for(int32_t y=0;y<SCREENWIDTH;y+=TWDSTEP) {
		int32_t setx=-1;
		for(int32_t x=0;x<SCREENWIDTH;x+=TWDSTEP) {
			setx++;
			int32_t f=-1;
			int raus=0;
			
			for(int32_t dy=0;dy<TWDSTEP;dy++) {
				for(int32_t dx=0;dx<TWDSTEP;dx++) {
					int32_t tmpf;
					GETDBY(x+dx,y+dy,tmpf);
					if (tmpf == SQUARE_GRAY_POTENTIALLY_WHITE) tmpf=SQUARE_GRAY;
					
					if ( (tmpf & PP_POSSIBLEPER) != 0) {
						// a periodic point => keep it
						f=SQUARE_BLACK;
						raus=1;
						break;
					}
					
					if (f<0) f=tmpf;
					else if (f != tmpf) {
						f=SQUARE_GRAY;
						// no break here, as BLACK (i.e. meaning here = PEIORIDCPOINT containing possibly => takes precedence
					}
				}
				if (raus>0) break;				
			}
			if ( (f<0) || (f>=256) ) {
				LOGMSG2("Periodicity point. Farbfehler %i\n",f);
				exit(99);
			}
			rgbz[setx]=f;
		}
		fwrite(rgbz,sizeof(char),dslen,fbmp);
	}
		
	fclose(fbmp);
	
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

// struct IntManager
IntManager::IntManager() {
	current=NULL;
	allocatedIdx=0;
	freeFromIdx=-1;
	double d=CHUNKSIZE; 
	d /= sizeof(int32_t);
	allocatePerBlockIdx=(int32_t)floor(d);
	anzptr=0;
}

int32_t* IntManager::getMemory(const int32_t aanz) {
	if (
		(!current) ||
		((freeFromIdx + aanz) >= allocatedIdx)
	) {
		if (anzptr >= MAXPTR) {
			LOGMSG("Error. Memory-IntManager.\n");
			exit(99);
		}
		ptr[anzptr]=current=new int32_t[allocatePerBlockIdx];
		anzptr++;
		if (!current) {
			LOGMSG("Error/2. Memory-IntManager.\n");
			exit(99);
		}
		freeFromIdx=0;
		allocatedIdx=allocatePerBlockIdx;
	}
	
	int32_t *p=&current[freeFromIdx];
	freeFromIdx += aanz;
	
	return p;
}

IntManager::~IntManager() {
	for(int32_t i=0;i<anzptr;i++) {
		if (ptr[i]) delete[] ptr[i];
	}
	anzptr=0;
}

// struct ByteManager
ByteManager::ByteManager() {
	current=NULL;
	allocatedIdx=0;
	freeFromIdx=-1;
	double d=CHUNKSIZE; 
	d /= sizeof(BYTE);
	allocatePerBlockIdx=(int32_t)floor(d);
	anzptr=0;
}

BYTE* ByteManager::getMemory(const int32_t aanz) {
	if (
		(!current) ||
		((freeFromIdx + aanz) >= allocatedIdx)
	) {
		if (anzptr >= MAXPTR) {
			LOGMSG("Error. Memory-ByteManager.\n");
			exit(99);
		}
		ptr[anzptr]=current=new BYTE[allocatePerBlockIdx];
		anzptr++;
		if (!current) {
			LOGMSG("Error/2. Memory-ByteManager.\n");
			exit(99);
		}
		freeFromIdx=0;
		allocatedIdx=allocatePerBlockIdx;
	}
	
	BYTE *p=&current[freeFromIdx];
	freeFromIdx += aanz;
	
	return p;
}

ByteManager::~ByteManager() {
	for(int32_t i=0;i<anzptr;i++) {
		if (ptr[i]) delete[] ptr[i];
	}
	anzptr=0;
}

// struct DByteMgr
DByteMgr::DByteMgr() {
	current=NULL;
	allocatedIdx=0;
	freeFromIdx=-1;
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

int32_t getfuncidx(const char* s) {
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
	if (strstr(p,rl)) {
		delete[] test;
		return 1;
	}

	delete[] test;
	return 0;
}

int32_t testA(void) {
	#ifdef _FPA
	FPA two,minustwo;
	two.set_vlong(2);
	minustwo.set_vlong(-2);
	if (
		(FAKTORAre > two) ||
		(FAKTORAre < minustwo) ||
		(FAKTORAim > two) ||
		(FAKTORAim < minustwo)
	) return 0;
	#else
	if (
		( fabs( (double)FAKTORAre ) > 2.0 ) ||
		( fabs( (double)FAKTORAim ) > 2.0 )
	) return 0;
	#endif
	
	return 1;
}

void setfunc_and_bitprecision(const int afunc,char* afn) {
	char tmp2[1024],tmp3[1024];
	
	int8_t bitprecision=1; 
	
	// |Creal|<=2, |Cimag|<=2
	#ifdef _FPA
	FPA two,minustwo;
	two.set_vlong(2);
	minustwo.set_vlong(-2);
	if (
		(seedC0re < minustwo) ||
		(seedC0re > two) ||
		(seedC1re < minustwo) ||
		(seedC1re > two) ||
		(seedC0im < minustwo) ||
		(seedC0im > two) ||
		(seedC1im < minustwo) ||
		(seedC1im > two)
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
			
			precompute_helperYdep=precomputeYdep_z3azc;
			precompute_helperXdep=precomputeXdep_z3azc;
			getBoundingBoxfA_helper=getBoundingBoxfA_z3azc_helper;
			
			precompute_helperYdep_double=precomputeYdep_z3azc_double;
			precompute_helperXdep_double=precomputeXdep_z3azc_double;
			getBoundingBoxfA_double=getBoundingBoxfA_z3azc_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_z3azc_double_oh;

			if (
				(_HELPER_Z3AZC_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z3AZC_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z3azc.\n");
				exit(99);
			}
			
			if (bitprecision>0) if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,D,LD,F1,QD,F6,FP,;R2L12,A,D,LD,F1,QD,F6,FP,;R2L13,A,D,LD,F1,QD,F6,FP,;R2L14,A,D,LD,F1,QD,F6,FP,;R2L15,A,D,LD,F1,QD,F6,FP,;R2L16,A,D,LD,F1,QD,F6,FP,;R2L17,A,D,LD,F1,QD,F6,FP,;R2L18,A,LD,F1,QD,F6,FP,;R2L19,A,LD,F1,QD,F6,FP,;R2L20,A,LD,F1,QD,F6,FP,;R2L21,A,F1,QD,F6,FP,;R2L22,A,F1,QD,F6,FP,;R2L23,A,F1,QD,F6,FP,;R2L24,A,F1,QD,F6,FP,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,D,LD,F1,QD,F6,FP,;R4L12,A,D,LD,F1,QD,F6,FP,;R4L13,A,D,LD,F1,QD,F6,FP,;R4L14,A,D,LD,F1,QD,F6,FP,;R4L15,A,D,LD,F1,QD,F6,FP,;R4L16,A,D,LD,F1,QD,F6,FP,;R4L17,A,D,LD,F1,QD,F6,FP,;R4L18,A,LD,F1,QD,F6,FP,;R4L19,A,LD,F1,QD,F6,FP,;R4L20,A,LD,F1,QD,F6,FP,;R4L21,A,F1,QD,F6,FP,;R4L22,A,F1,QD,F6,FP,;R4L23,A,F1,QD,F6,FP,;R4L24,A,F1,QD,F6,FP,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,D,LD,F1,QD,F6,FP,;R8L12,A,D,LD,F1,QD,F6,FP,;R8L13,A,D,LD,F1,QD,F6,FP,;R8L14,A,D,LD,F1,QD,F6,FP,;R8L15,A,D,LD,F1,QD,F6,FP,;R8L16,A,D,LD,F1,QD,F6,FP,;R8L17,A,D,LD,F1,QD,F6,FP,;R8L18,A,LD,F1,QD,F6,FP,;R8L19,A,LD,F1,QD,F6,FP,;R8L20,A,LD,F1,QD,F6,FP,;R8L21,A,F1,QD,F6,FP,;R8L22,A,F1,QD,F6,FP,;R8L23,A,F1,QD,F6,FP,;R8L24,A,F1,QD,F6,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,D,LD,F1,QD,F6,FP,;R2L12,A,D,LD,F1,QD,F6,FP,;R2L13,A,D,LD,F1,QD,F6,FP,;R2L14,A,D,LD,F1,QD,F6,FP,;R2L15,A,D,LD,F1,QD,F6,FP,;R2L16,A,D,LD,F1,QD,F6,FP,;R2L17,A,D,LD,F1,QD,F6,FP,;R2L18,A,LD,F1,QD,F6,FP,;R2L19,A,LD,F1,QD,F6,FP,;R2L20,A,LD,F1,QD,F6,FP,;R2L21,A,F1,QD,F6,FP,;R2L22,A,F1,QD,F6,FP,;R2L23,A,F1,QD,F6,FP,;R2L24,A,F1,QD,F6,FP,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,D,LD,F1,QD,F6,FP,;R4L12,A,D,LD,F1,QD,F6,FP,;R4L13,A,D,LD,F1,QD,F6,FP,;R4L14,A,D,LD,F1,QD,F6,FP,;R4L15,A,D,LD,F1,QD,F6,FP,;R4L16,A,D,LD,F1,QD,F6,FP,;R4L17,A,D,LD,F1,QD,F6,FP,;R4L18,A,LD,F1,QD,F6,FP,;R4L19,A,LD,F1,QD,F6,FP,;R4L20,A,LD,F1,QD,F6,FP,;R4L21,A,F1,QD,F6,FP,;R4L22,A,F1,QD,F6,FP,;R4L23,A,F1,QD,F6,FP,;R4L24,A,F1,QD,F6,FP,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,D,LD,F1,QD,F6,FP,;R8L12,A,D,LD,F1,QD,F6,FP,;R8L13,A,D,LD,F1,QD,F6,FP,;R8L14,A,D,LD,F1,QD,F6,FP,;R8L15,A,D,LD,F1,QD,F6,FP,;R8L16,A,D,LD,F1,QD,F6,FP,;R8L17,A,D,LD,F1,QD,F6,FP,;R8L18,A,LD,F1,QD,F6,FP,;R8L19,A,LD,F1,QD,F6,FP,;R8L20,A,LD,F1,QD,F6,FP,;R8L21,A,F1,QD,F6,FP,;R8L22,A,F1,QD,F6,FP,;R8L23,A,F1,QD,F6,FP,;R8L24,A,F1,QD,F6,FP,;"
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

			precompute_helperYdep=precomputeYdep_z4azc;
			precompute_helperXdep=precomputeXdep_z4azc;
			getBoundingBoxfA_helper=getBoundingBoxfA_z4azc_helper;

			precompute_helperYdep_double=precomputeYdep_z4azc_double;
			precompute_helperXdep_double=precomputeXdep_z4azc_double;
			getBoundingBoxfA_double=getBoundingBoxfA_z4azc_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_z4azc_double_oh;

			if (
				(_HELPER_Z4AZC_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z4AZC_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z4azc.\n");
				exit(99);
			}

			if (bitprecision>0) if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,D,LD,F1,QD,F6,FP,;R2L12,A,D,LD,F1,QD,F6,FP,;R2L13,A,D,LD,F1,QD,F6,FP,;R2L14,A,LD,F1,QD,F6,FP,;R2L15,A,LD,F1,QD,F6,FP,;R2L16,A,F1,QD,F6,FP,;R2L17,A,F1,QD,F6,FP,;R2L18,A,F1,QD,F6,FP,;R2L19,A,F1,QD,F6,FP,;R2L20,A,F1,QD,F6,FP,;R2L21,A,F1,QD,F6,FP,;R2L22,A,F1,QD,F6,FP,;R2L23,A,F1,QD,F6,FP,;R2L24,A,F1,QD,F6,FP,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,D,LD,F1,QD,F6,FP,;R4L12,A,D,LD,F1,QD,F6,FP,;R4L13,A,D,LD,F1,QD,F6,FP,;R4L14,A,LD,F1,QD,F6,FP,;R4L15,A,LD,F1,QD,F6,FP,;R4L16,A,F1,QD,F6,FP,;R4L17,A,F1,QD,F6,FP,;R4L18,A,F1,QD,F6,FP,;R4L19,A,F1,QD,F6,FP,;R4L20,A,F1,QD,F6,FP,;R4L21,A,F1,QD,F6,FP,;R4L22,A,F1,QD,F6,FP,;R4L23,A,F1,QD,F6,FP,;R4L24,A,F1,QD,F6,FP,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,D,LD,F1,QD,F6,FP,;R8L12,A,D,LD,F1,QD,F6,FP,;R8L13,A,D,LD,F1,QD,F6,FP,;R8L14,A,LD,F1,QD,F6,FP,;R8L15,A,LD,F1,QD,F6,FP,;R8L16,A,F1,QD,F6,FP,;R8L17,A,F1,QD,F6,FP,;R8L18,A,F1,QD,F6,FP,;R8L19,A,F1,QD,F6,FP,;R8L20,A,F1,QD,F6,FP,;R8L21,A,F1,QD,F6,FP,;R8L22,A,F1,QD,F6,FP,;R8L23,A,F1,QD,F6,FP,;R8L24,A,F1,QD,F6,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,D,LD,F1,QD,F6,FP,;R2L12,A,D,LD,F1,QD,F6,FP,;R2L13,A,D,LD,F1,QD,F6,FP,;R2L14,A,LD,F1,QD,F6,FP,;R2L15,A,LD,F1,QD,F6,FP,;R2L16,A,F1,QD,F6,FP,;R2L17,A,F1,QD,F6,FP,;R2L18,A,F1,QD,F6,FP,;R2L19,A,F1,QD,F6,FP,;R2L20,A,F1,QD,F6,FP,;R2L21,A,F1,QD,F6,FP,;R2L22,A,F1,QD,F6,FP,;R2L23,A,F1,QD,F6,FP,;R2L24,A,F1,QD,F6,FP,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,D,LD,F1,QD,F6,FP,;R4L12,A,D,LD,F1,QD,F6,FP,;R4L13,A,D,LD,F1,QD,F6,FP,;R4L14,A,LD,F1,QD,F6,FP,;R4L15,A,LD,F1,QD,F6,FP,;R4L16,A,F1,QD,F6,FP,;R4L17,A,F1,QD,F6,FP,;R4L18,A,F1,QD,F6,FP,;R4L19,A,F1,QD,F6,FP,;R4L20,A,F1,QD,F6,FP,;R4L21,A,F1,QD,F6,FP,;R4L22,A,F1,QD,F6,FP,;R4L23,A,F1,QD,F6,FP,;R4L24,A,F1,QD,F6,FP,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,D,LD,F1,QD,F6,FP,;R8L12,A,D,LD,F1,QD,F6,FP,;R8L13,A,D,LD,F1,QD,F6,FP,;R8L14,A,LD,F1,QD,F6,FP,;R8L15,A,LD,F1,QD,F6,FP,;R8L16,A,F1,QD,F6,FP,;R8L17,A,F1,QD,F6,FP,;R8L18,A,F1,QD,F6,FP,;R8L19,A,F1,QD,F6,FP,;R8L20,A,F1,QD,F6,FP,;R8L21,A,F1,QD,F6,FP,;R8L22,A,F1,QD,F6,FP,;R8L23,A,F1,QD,F6,FP,;R8L24,A,F1,QD,F6,FP,;"
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
			
			precompute_helperYdep=precomputeYdep_z5azc;
			precompute_helperXdep=precomputeXdep_z5azc;
			getBoundingBoxfA_helper=getBoundingBoxfA_z5azc_helper;
			
			precompute_helperYdep_double=precomputeYdep_z5azc_double;
			precompute_helperXdep_double=precomputeXdep_z5azc_double;
			getBoundingBoxfA_double=getBoundingBoxfA_z5azc_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_z5azc_double_oh;

			if (
				(_HELPER_Z5AZC_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z5AZC_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z5azc.\n");
				exit(99);
			}

			if (bitprecision>0) if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,LD,F1,QD,F6,FP,;R2L12,A,LD,F1,QD,F6,FP,;R2L13,A,F1,QD,F6,FP,;R2L14,A,F1,QD,F6,FP,;R2L15,A,F1,QD,F6,FP,;R2L16,A,F1,QD,F6,FP,;R2L17,A,F1,QD,F6,FP,;R2L18,A,F1,QD,F6,FP,;R2L19,A,F1,QD,F6,FP,;R2L20,A,F1,QD,F6,FP,;R2L21,A,F1,QD,F6,FP,;R2L22,A,QD,F6,;R2L23,A,QD,F6,;R2L24,A,QD,F6,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,LD,F1,QD,F6,FP,;R4L12,A,LD,F1,QD,F6,FP,;R4L13,A,F1,QD,F6,FP,;R4L14,A,F1,QD,F6,FP,;R4L15,A,F1,QD,F6,FP,;R4L16,A,F1,QD,F6,FP,;R4L17,A,F1,QD,F6,FP,;R4L18,A,F1,QD,F6,FP,;R4L19,A,F1,QD,F6,FP,;R4L20,A,F1,QD,F6,FP,;R4L21,A,F1,QD,F6,FP,;R4L22,A,QD,F6,FP,;R4L23,A,QD,F6,;R4L24,A,QD,F6,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,LD,F1,QD,F6,FP,;R8L12,A,LD,F1,QD,F6,FP,;R8L13,A,F1,QD,F6,FP,;R8L14,A,F1,QD,F6,FP,;R8L15,A,F1,QD,F6,FP,;R8L16,A,F1,QD,F6,FP,;R8L17,A,F1,QD,F6,FP,;R8L18,A,F1,QD,F6,FP,;R8L19,A,F1,QD,F6,FP,;R8L20,A,F1,QD,F6,FP,;R8L21,A,F1,QD,F6,FP,;R8L22,A,QD,F6,FP,;R8L23,A,QD,F6,FP,;R8L24,A,QD,F6,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,LD,F1,QD,F6,FP,;R2L12,A,LD,F1,QD,F6,FP,;R2L13,A,F1,QD,F6,FP,;R2L14,A,F1,QD,F6,FP,;R2L15,A,F1,QD,F6,FP,;R2L16,A,F1,QD,F6,FP,;R2L17,A,F1,QD,F6,FP,;R2L18,A,F1,QD,F6,FP,;R2L19,A,F1,QD,F6,FP,;R2L20,A,F1,QD,F6,FP,;R2L21,A,F1,QD,F6,FP,;R2L22,A,QD,F6,;R2L23,A,QD,F6,;R2L24,A,QD,F6,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,LD,F1,QD,F6,FP,;R4L12,A,LD,F1,QD,F6,FP,;R4L13,A,F1,QD,F6,FP,;R4L14,A,F1,QD,F6,FP,;R4L15,A,F1,QD,F6,FP,;R4L16,A,F1,QD,F6,FP,;R4L17,A,F1,QD,F6,FP,;R4L18,A,F1,QD,F6,FP,;R4L19,A,F1,QD,F6,FP,;R4L20,A,F1,QD,F6,FP,;R4L21,A,F1,QD,F6,FP,;R4L22,A,QD,F6,FP,;R4L23,A,QD,F6,;R4L24,A,QD,F6,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,LD,F1,QD,F6,FP,;R8L12,A,LD,F1,QD,F6,FP,;R8L13,A,F1,QD,F6,FP,;R8L14,A,F1,QD,F6,FP,;R8L15,A,F1,QD,F6,FP,;R8L16,A,F1,QD,F6,FP,;R8L17,A,F1,QD,F6,FP,;R8L18,A,F1,QD,F6,FP,;R8L19,A,F1,QD,F6,FP,;R8L20,A,F1,QD,F6,FP,;R8L21,A,F1,QD,F6,FP,;R8L22,A,QD,F6,FP,;R8L23,A,QD,F6,FP,;R8L24,A,QD,F6,;"
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
			
			precompute_helperYdep=precomputeYdep_z6azc;
			precompute_helperXdep=precomputeXdep_z6azc;
			getBoundingBoxfA_helper=getBoundingBoxfA_z6azc_helper;

			precompute_helperYdep_double=precomputeYdep_z6azc_double;
			precompute_helperXdep_double=precomputeXdep_z6azc_double;
			getBoundingBoxfA_double=getBoundingBoxfA_z6azc_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_z6azc_double_oh;

			if (
				(_HELPER_Z6AZC_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z6AZC_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z6azc.\n");
				exit(99);
			}

			if (bitprecision>0) if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,LD,F1,QD,F6,FP,;R2L10,A,LD,F1,QD,F6,FP,;R2L11,A,F1,QD,F6,FP,;R2L12,A,F1,QD,F6,FP,;R2L13,A,F1,QD,F6,FP,;R2L14,A,F1,QD,F6,FP,;R2L15,A,F1,QD,F6,FP,;R2L16,A,F1,QD,F6,FP,;R2L17,A,F1,QD,F6,FP,;R2L18,A,QD,F6,;R2L19,A,QD,F6,;R2L20,A,QD,F6,;R2L21,A,QD,F6,;R2L22,A,F6,;R2L23,A,F6,;R2L24,A,F6,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,LD,F1,QD,F6,FP,;R4L10,A,LD,F1,QD,F6,FP,;R4L11,A,F1,QD,F6,FP,;R4L12,A,F1,QD,F6,FP,;R4L13,A,F1,QD,F6,FP,;R4L14,A,F1,QD,F6,FP,;R4L15,A,F1,QD,F6,FP,;R4L16,A,F1,QD,F6,FP,;R4L17,A,F1,QD,F6,FP,;R4L18,A,QD,F6,FP,;R4L19,A,QD,F6,;R4L20,A,QD,F6,;R4L21,A,QD,F6,;R4L22,A,F6,;R4L23,A,F6,;R4L24,A,F6,;R8L8,A,LD,F1,QD,F6,FP,;R8L9,A,LD,F1,QD,F6,FP,;R8L10,A,LD,F1,QD,F6,FP,;R8L11,A,F1,QD,F6,FP,;R8L12,A,F1,QD,F6,FP,;R8L13,A,F1,QD,F6,FP,;R8L14,A,F1,QD,F6,FP,;R8L15,A,F1,QD,F6,FP,;R8L16,A,F1,QD,F6,FP,;R8L17,A,F1,QD,F6,FP,;R8L18,A,QD,F6,FP,;R8L19,A,QD,F6,FP,;R8L20,A,QD,F6,;R8L21,A,QD,F6,;R8L22,A,F6,;R8L23,A,F6,;R8L24,A,F6,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,LD,F1,QD,F6,FP,;R2L10,A,LD,F1,QD,F6,FP,;R2L11,A,F1,QD,F6,FP,;R2L12,A,F1,QD,F6,FP,;R2L13,A,F1,QD,F6,FP,;R2L14,A,F1,QD,F6,FP,;R2L15,A,F1,QD,F6,FP,;R2L16,A,F1,QD,F6,FP,;R2L17,A,F1,QD,F6,FP,;R2L18,A,QD,F6,;R2L19,A,QD,F6,;R2L20,A,QD,F6,;R2L21,A,QD,F6,;R2L22,A,F6,;R2L23,A,F6,;R2L24,A,F6,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,LD,F1,QD,F6,FP,;R4L10,A,LD,F1,QD,F6,FP,;R4L11,A,F1,QD,F6,FP,;R4L12,A,F1,QD,F6,FP,;R4L13,A,F1,QD,F6,FP,;R4L14,A,F1,QD,F6,FP,;R4L15,A,F1,QD,F6,FP,;R4L16,A,F1,QD,F6,FP,;R4L17,A,F1,QD,F6,FP,;R4L18,A,QD,F6,FP,;R4L19,A,QD,F6,;R4L20,A,QD,F6,;R4L21,A,QD,F6,;R4L22,A,F6,;R4L23,A,F6,;R4L24,A,F6,;R8L8,A,LD,F1,QD,F6,FP,;R8L9,A,LD,F1,QD,F6,FP,;R8L10,A,LD,F1,QD,F6,FP,;R8L11,A,F1,QD,F6,FP,;R8L12,A,F1,QD,F6,FP,;R8L13,A,F1,QD,F6,FP,;R8L14,A,F1,QD,F6,FP,;R8L15,A,F1,QD,F6,FP,;R8L16,A,F1,QD,F6,FP,;R8L17,A,F1,QD,F6,FP,;R8L18,A,QD,F6,FP,;R8L19,A,QD,F6,FP,;R8L20,A,QD,F6,;R8L21,A,QD,F6,;R8L22,A,F6,;R8L23,A,F6,;R8L24,A,F6,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
		case FUNC_2ITZ2C: {
			checkclockatbbxadd >>= 5;
			sprintf(afn,"_L%02i_%s_2itz2c_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2));
			getBoundingBoxfA=getBoundingBoxfA_2itz2c;
			
			precompute_helperYdep=precomputeYdep_2itz2c;
			precompute_helperXdep=precomputeXdep_2itz2c;
			getBoundingBoxfA_helper=getBoundingBoxfA_2itz2c_helper;

			precompute_helperYdep_double=precomputeYdep_2itz2c_double;
			precompute_helperXdep_double=precomputeXdep_2itz2c_double;
			getBoundingBoxfA_double=getBoundingBoxfA_2itz2c_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_2itz2c_double_oh;

			if (
				(_HELPER_2ITZ2C_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_2ITZ2C_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices 2itz2c.\n");
				exit(99);
			}

			break;
		}
		case FUNC_Z7AZC: {
			checkclockatbbxadd >>= 6;
			checkclockatbbxcount0 >>= 1;
			sprintf(afn,"_L%02i_%s_z7azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z7azc;
			
			precompute_helperYdep=precomputeYdep_z7azc;
			precompute_helperXdep=precomputeXdep_z7azc;
			getBoundingBoxfA_helper=getBoundingBoxfA_z7azc_helper;

			precompute_helperYdep_double=precomputeYdep_z7azc_double;
			precompute_helperXdep_double=precomputeXdep_z7azc_double;
			getBoundingBoxfA_double=getBoundingBoxfA_z7azc_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_z7azc_double_oh;

			if (
				(_HELPER_Z7AZC_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z7AZC_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z7azc.\n");
				exit(99);
			}

			if (bitprecision>0) if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,LD,F1,QD,F6,FP,;R2L9,A,F1,QD,F6,FP,;R2L10,A,F1,QD,F6,FP,;R2L11,A,F1,QD,F6,FP,;R2L12,A,F1,QD,F6,FP,;R2L13,A,F1,QD,F6,FP,;R2L14,A,F1,QD,F6,FP,;R2L15,A,F1,QD,F6,FP,;R2L16,A,QD,F6,;R2L17,A,QD,F6,;R2L18,A,QD,F6,;R2L19,A,F6,;R2L20,A,F6,;R2L21,A,F6,;R2L22,A,F6,;R2L23,A,;R2L24,A,;R4L8,A,LD,F1,QD,F6,FP,;R4L9,A,F1,QD,F6,FP,;R4L10,A,F1,QD,F6,FP,;R4L11,A,F1,QD,F6,FP,;R4L12,A,F1,QD,F6,FP,;R4L13,A,F1,QD,F6,FP,;R4L14,A,F1,QD,F6,FP,;R4L15,A,F1,QD,F6,FP,;R4L16,A,QD,F6,FP,;R4L17,A,QD,F6,;R4L18,A,QD,F6,;R4L19,A,F6,;R4L20,A,F6,;R4L21,A,F6,;R4L22,A,F6,;R4L23,A,;R4L24,A,;R8L8,A,LD,F1,QD,F6,FP,;R8L9,A,F1,QD,F6,FP,;R8L10,A,F1,QD,F6,FP,;R8L11,A,F1,QD,F6,FP,;R8L12,A,F1,QD,F6,FP,;R8L13,A,F1,QD,F6,FP,;R8L14,A,F1,QD,F6,FP,;R8L15,A,F1,QD,F6,FP,;R8L16,A,QD,F6,FP,;R8L17,A,QD,F6,FP,;R8L18,A,QD,F6,;R8L19,A,F6,;R8L20,A,F6,;R8L21,A,F6,;R8L22,A,F6,;R8L23,A,;R8L24,A,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,LD,F1,QD,F6,FP,;R2L9,A,F1,QD,F6,FP,;R2L10,A,F1,QD,F6,FP,;R2L11,A,F1,QD,F6,FP,;R2L12,A,F1,QD,F6,FP,;R2L13,A,F1,QD,F6,FP,;R2L14,A,F1,QD,F6,FP,;R2L15,A,F1,QD,F6,FP,;R2L16,A,QD,F6,;R2L17,A,QD,F6,;R2L18,A,QD,F6,;R2L19,A,F6,;R2L20,A,F6,;R2L21,A,F6,;R2L22,A,F6,;R2L23,A,;R2L24,A,;R4L8,A,LD,F1,QD,F6,FP,;R4L9,A,F1,QD,F6,FP,;R4L10,A,F1,QD,F6,FP,;R4L11,A,F1,QD,F6,FP,;R4L12,A,F1,QD,F6,FP,;R4L13,A,F1,QD,F6,FP,;R4L14,A,F1,QD,F6,FP,;R4L15,A,F1,QD,F6,FP,;R4L16,A,QD,F6,FP,;R4L17,A,QD,F6,;R4L18,A,QD,F6,;R4L19,A,F6,;R4L20,A,F6,;R4L21,A,F6,;R4L22,A,F6,;R4L23,A,;R4L24,A,;R8L8,A,LD,F1,QD,F6,FP,;R8L9,A,F1,QD,F6,FP,;R8L10,A,F1,QD,F6,FP,;R8L11,A,F1,QD,F6,FP,;R8L12,A,F1,QD,F6,FP,;R8L13,A,F1,QD,F6,FP,;R8L14,A,F1,QD,F6,FP,;R8L15,A,F1,QD,F6,FP,;R8L16,A,QD,F6,FP,;R8L17,A,QD,F6,FP,;R8L18,A,QD,F6,;R8L19,A,F6,;R8L20,A,F6,;R8L21,A,F6,;R8L22,A,F6,;R8L23,A,;R8L24,A,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
		case FUNC_Z8AZC: {
			checkclockatbbxadd >>= 7;
			checkclockatbbxcount0 >>= 1;
			sprintf(afn,"_L%02i_%s_z8azc_%s_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
			getBoundingBoxfA=getBoundingBoxfA_z8azc;
			
			precompute_helperYdep=precomputeYdep_z8azc;
			precompute_helperXdep=precomputeXdep_z8azc;
			getBoundingBoxfA_helper=getBoundingBoxfA_z8azc_helper;

			precompute_helperYdep_double=precomputeYdep_z8azc_double;
			precompute_helperXdep_double=precomputeXdep_z8azc_double;
			getBoundingBoxfA_double=getBoundingBoxfA_z8azc_double;
			getBoundingBoxfA_double_oh=getBoundingBoxfA_z8azc_double_oh;

			if (
				(_HELPER_Z8AZC_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z8AZC_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z7azc.\n");
				exit(99);
			}

			if (bitprecision>0) if (testA() <= 0) bitprecision=0;
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
		default: {
			getBoundingBoxfA=getBoundingBoxfA_z2c;
			sprintf(afn,"_L%02i_%s_z2c_%s.bmp",
				REFINEMENTLEVEL,
				NNTYPSTR,seedCstr(tmp2));
			precompute_helperYdep=precomputeYdep_z2c;
			precompute_helperXdep=precomputeXdep_z2c;
			getBoundingBoxfA_helper=getBoundingBoxfA_z2c_helper;
			// no _double functions here
			if (
				(_HELPER_Z2C_Xdep_ANZ >= MAXHELPERVALUES) ||
				(_HELPER_Z2C_Ydep_ANZ >= MAXHELPERVALUES)
			) {
				LOGMSG("Implementation error. Too many helper value indices z2c.\n");
				exit(99);
			}
			
			// real part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,D,LD,F1,QD,F6,FP,;R2L12,A,D,LD,F1,QD,F6,FP,;R2L13,A,D,LD,F1,QD,F6,FP,;R2L14,A,D,LD,F1,QD,F6,FP,;R2L15,A,D,LD,F1,QD,F6,FP,;R2L16,A,D,LD,F1,QD,F6,FP,;R2L17,A,D,LD,F1,QD,F6,FP,;R2L18,A,D,LD,F1,QD,F6,FP,;R2L19,A,D,LD,F1,QD,F6,FP,;R2L20,A,D,LD,F1,QD,F6,FP,;R2L21,A,D,LD,F1,QD,F6,FP,;R2L22,A,D,LD,F1,QD,F6,FP,;R2L23,A,D,LD,F1,QD,F6,FP,;R2L24,A,D,LD,F1,QD,F6,FP,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,D,LD,F1,QD,F6,FP,;R4L12,A,D,LD,F1,QD,F6,FP,;R4L13,A,D,LD,F1,QD,F6,FP,;R4L14,A,D,LD,F1,QD,F6,FP,;R4L15,A,D,LD,F1,QD,F6,FP,;R4L16,A,D,LD,F1,QD,F6,FP,;R4L17,A,D,LD,F1,QD,F6,FP,;R4L18,A,D,LD,F1,QD,F6,FP,;R4L19,A,D,LD,F1,QD,F6,FP,;R4L20,A,D,LD,F1,QD,F6,FP,;R4L21,A,D,LD,F1,QD,F6,FP,;R4L22,A,D,LD,F1,QD,F6,FP,;R4L23,A,D,LD,F1,QD,F6,FP,;R4L24,A,D,LD,F1,QD,F6,FP,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,D,LD,F1,QD,F6,FP,;R8L12,A,D,LD,F1,QD,F6,FP,;R8L13,A,D,LD,F1,QD,F6,FP,;R8L14,A,D,LD,F1,QD,F6,FP,;R8L15,A,D,LD,F1,QD,F6,FP,;R8L16,A,D,LD,F1,QD,F6,FP,;R8L17,A,D,LD,F1,QD,F6,FP,;R8L18,A,D,LD,F1,QD,F6,FP,;R8L19,A,D,LD,F1,QD,F6,FP,;R8L20,A,D,LD,F1,QD,F6,FP,;R8L21,A,D,LD,F1,QD,F6,FP,;R8L22,A,D,LD,F1,QD,F6,FP,;R8L23,A,D,LD,F1,QD,F6,FP,;R8L24,A,D,LD,F1,QD,F6,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;
			
			// imaginary part
			if (bitprecision>0) if (bitsSufficient(
				";R2L8,A,D,LD,F1,QD,F6,FP,;R2L9,A,D,LD,F1,QD,F6,FP,;R2L10,A,D,LD,F1,QD,F6,FP,;R2L11,A,D,LD,F1,QD,F6,FP,;R2L12,A,D,LD,F1,QD,F6,FP,;R2L13,A,D,LD,F1,QD,F6,FP,;R2L14,A,D,LD,F1,QD,F6,FP,;R2L15,A,D,LD,F1,QD,F6,FP,;R2L16,A,D,LD,F1,QD,F6,FP,;R2L17,A,D,LD,F1,QD,F6,FP,;R2L18,A,D,LD,F1,QD,F6,FP,;R2L19,A,D,LD,F1,QD,F6,FP,;R2L20,A,D,LD,F1,QD,F6,FP,;R2L21,A,D,LD,F1,QD,F6,FP,;R2L22,A,D,LD,F1,QD,F6,FP,;R2L23,A,D,LD,F1,QD,F6,FP,;R2L24,A,D,LD,F1,QD,F6,FP,;R4L8,A,D,LD,F1,QD,F6,FP,;R4L9,A,D,LD,F1,QD,F6,FP,;R4L10,A,D,LD,F1,QD,F6,FP,;R4L11,A,D,LD,F1,QD,F6,FP,;R4L12,A,D,LD,F1,QD,F6,FP,;R4L13,A,D,LD,F1,QD,F6,FP,;R4L14,A,D,LD,F1,QD,F6,FP,;R4L15,A,D,LD,F1,QD,F6,FP,;R4L16,A,D,LD,F1,QD,F6,FP,;R4L17,A,D,LD,F1,QD,F6,FP,;R4L18,A,D,LD,F1,QD,F6,FP,;R4L19,A,D,LD,F1,QD,F6,FP,;R4L20,A,D,LD,F1,QD,F6,FP,;R4L21,A,D,LD,F1,QD,F6,FP,;R4L22,A,D,LD,F1,QD,F6,FP,;R4L23,A,D,LD,F1,QD,F6,FP,;R4L24,A,D,LD,F1,QD,F6,FP,;R8L8,A,D,LD,F1,QD,F6,FP,;R8L9,A,D,LD,F1,QD,F6,FP,;R8L10,A,D,LD,F1,QD,F6,FP,;R8L11,A,D,LD,F1,QD,F6,FP,;R8L12,A,D,LD,F1,QD,F6,FP,;R8L13,A,D,LD,F1,QD,F6,FP,;R8L14,A,D,LD,F1,QD,F6,FP,;R8L15,A,D,LD,F1,QD,F6,FP,;R8L16,A,D,LD,F1,QD,F6,FP,;R8L17,A,D,LD,F1,QD,F6,FP,;R8L18,A,D,LD,F1,QD,F6,FP,;R8L19,A,D,LD,F1,QD,F6,FP,;R8L20,A,D,LD,F1,QD,F6,FP,;R8L21,A,D,LD,F1,QD,F6,FP,;R8L22,A,D,LD,F1,QD,F6,FP,;R8L23,A,D,LD,F1,QD,F6,FP,;R8L24,A,D,LD,F1,QD,F6,FP,;"
				,RANGE1,REFINEMENTLEVEL,NTS
			) <= 0) bitprecision=0;

			break;
		}
	} // switch
}

// struct ListeFIFO
ListeFIFO::ListeFIFO() {
	next_readpos=next_writepos=0;
	werte=NULL;
	allokiereAmStueck=(CHUNKSIZE / sizeof(Int2));
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
		werte=new Int2[allokiereAmStueck];
		if (!werte) {
			LOGMSG("ListeFIFO: No memory\n");
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
		werte=new DFSPunkt[allokiereAmStueck];
		if (!werte) {
			LOGMSG("ListeDFS: No memory\n");
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

// RefPointArray
void RefPointArray::addRefPoint(const int32_t ax,const int32_t ay,const int32_t ablobid) {
	if (!listY[ay].points) {
		LOGMSG2("Implementation error. RefPointArray::addRefPoint not allocated at y=%i\n",ay);
		exit(99);
	}
	listY[ay].addXB(ax,ablobid);
}

RefPoint* RefPointArray::getRefPtr(const int32_t ax,const int32_t ay) {
	if (!listY[ay].points) {
		LOGMSG2("Error. getBloId of non-refpoint-row %i\n",ay);
		exit(99);
	}
	
	return listY[ay].getRefPtr(ax);
}

RefPointArray::RefPointArray() {
	listY=new RefList[SCREENWIDTH];
}

RefPointArray::~RefPointArray() {
	delete[] listY; // content is freed by managaer
}

// struct Int2Manager
Int2Manager::Int2Manager() {
	current=NULL;
	allokierteIdx=0;
	freiAbIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(Int2);
	allokierePerBlockIdx=(int32_t)floor(d);
}

Int2Manager::~Int2Manager() {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
}

PInt2 Int2Manager::getMemory(const int aanz) {
	if (anzptr >= MAXPTR) {
		LOGMSG("Int2Manager:: Memory error.\n");
		exit(99);
	}
	if (
		(!current) ||
		((freiAbIdx + aanz + 2) >= allokierteIdx)
	) {
		ptr[anzptr]=current=new Int2[allokierePerBlockIdx];
		anzptr++;
		if (!current) {
			printf("Memory-Fehler. Int2Manager.\n");
			exit(99);
		}
		freiAbIdx=0;
		allokierteIdx=allokierePerBlockIdx;
	}
	
	PInt2 p=&current[freiAbIdx];
	freiAbIdx += aanz;
	return p;
}

// struct ScreenRectManager
ScreenRectManager::ScreenRectManager() {
	current=NULL;
	allokierteIdx=0;
	freiAbIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(ScreenRect);
	allokierePerBlockIdx=(int32_t)floor(d);
}

ScreenRectManager::~ScreenRectManager() {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
}

ScreenRect* ScreenRectManager::getMemory(const int aanz) {
	if (anzptr >= MAXPTR) {
		LOGMSG("ScreenRectManagerManager:: Memory error.\n");
		exit(99);
	}
	if (
		(!current) ||
		((freiAbIdx + aanz + 2) >= allokierteIdx)
	) {
		ptr[anzptr]=current=new ScreenRect[allokierePerBlockIdx];
		anzptr++;
		if (!current) {
			printf("Memory-Fehler. ScreenRectManager.\n");
			exit(99);
		}
		freiAbIdx=0;
		allokierteIdx=allokierePerBlockIdx;
	}
	
	ScreenRect* p=&current[freiAbIdx];
	freiAbIdx += aanz;
	return p;
}

// struct RefPointManager
RefPointManager::RefPointManager() {
	current=NULL;
	allokierteIdx=0;
	freiAbIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(RefPoint);
	allokierePerBlockIdx=(int32_t)floor(d);
}

RefPointManager::~RefPointManager() {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
}

RefPoint* RefPointManager::getMemory(const int aanz) {
	if (anzptr >= MAXPTR) {
		LOGMSG("RefPointManager:: Memory error.\n");
		exit(99);
	}
	if (
		(!current) ||
		((freiAbIdx + aanz + 2) >= allokierteIdx)
	) {
		ptr[anzptr]=current=new RefPoint[allokierePerBlockIdx];
		anzptr++;
		if (!current) {
			printf("Memory-Fehler. RefPointManager.\n");
			exit(99);
		}
		freiAbIdx=0;
		allokierteIdx=allokierePerBlockIdx;
	}
	
	PRefPoint p=&current[freiAbIdx];
	freiAbIdx += aanz;
	return p;
}

// struct HelperManager
HelperManager::HelperManager() {
	current=NULL;
	allokierteIdx=0;
	freiAbIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(Helper);
	allokierePerBlockIdx=(int32_t)floor(d);
}

HelperManager::~HelperManager() {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
}

PHelper HelperManager::getMemory(const int aanz) {
	if (anzptr >= MAXPTR) {
		LOGMSG("HelperManager:: Memory error.\n");
		exit(99);
	}
	if (
		(!current) ||
		((freiAbIdx + aanz + 2) >= allokierteIdx)
	) {
		ptr[anzptr]=current=new Helper[allokierePerBlockIdx];
		anzptr++;
		if (!current) {
			printf("Memory-error. HelperManager.\n");
			exit(99);
		}
		freiAbIdx=0;
		allokierteIdx=allokierePerBlockIdx;
	}
	
	PHelper p=&current[freiAbIdx];
	freiAbIdx += aanz;
	
	return p;
}

// struct Helper_doublemanager
Helper_doubleManager::Helper_doubleManager() {
	current=NULL;
	allokierteIdx=0;
	freiAbIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(Helper_double);
	allokierePerBlockIdx=(int32_t)floor(d);
}

Helper_doubleManager::~Helper_doubleManager() {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
}

PHelper_double Helper_doubleManager::getMemory(const int aanz) {
	if (anzptr >= MAXPTR) {
		LOGMSG("HelperManager_double:: Memory error.\n");
		exit(99);
	}
	if (
		(!current) ||
		((freiAbIdx + aanz + 2) >= allokierteIdx)
	) {
		ptr[anzptr]=current=new Helper_double[allokierePerBlockIdx];
		anzptr++;
		if (!current) {
			printf("Memory-error. HelperManager_double.\n");
			exit(99);
		}
		freiAbIdx=0;
		allokierteIdx=allokierePerBlockIdx;
	}
	
	PHelper_double p=&current[freiAbIdx];
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
inline bool operator!=(const FPA a,const FPA b) {
	return (FPA_vgl(a,b) != 0);
}

inline bool operator>(const FPA a,const FPA b) {
	int vgl=FPA_vgl(a,b);
	if (vgl > 0) return 1;
	return 0;
}

inline bool operator>=(const FPA a,const FPA b) {
	return (FPA_vgl(a,b)>=0);
}

inline bool operator<=(const FPA a,const FPA b) {
	return (FPA_vgl(a,b)<=0);
}

inline bool operator==(const FPA a,const FPA b) {
	return (FPA_vgl(a,b)==0);
}

inline bool operator<(const FPA a,const FPA b) {
	int vgl=FPA_vgl(a,b);
	if (vgl < 0) return 1;
	return 0;
}

inline FPA operator+(const FPA a,const FPA b) {
    FPA ret;
    FPA_add_ZAB(ret,a,b);

    return ret;
}

inline FPA operator-(const FPA a,const FPA b) {
    FPA ret;
    FPA_sub_ZAB(ret,a,b);

    return ret;
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
		erg.vorz=1;
	} else {
		// |a| < |b|
		FPA_sub_abs_ovgl_ZAB(erg,b,a);
		erg.vorz=-1;
	}
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

void minimaxFPAAB(FPA& mi,FPA& ma,const FPA term1,const FPA term2) {
	int vgl=FPA_vgl(term1,term2);
	if (vgl>=0) { mi=term2; ma=term1; }
	else { mi=term1; ma=term2; }
}

void minimaxFPAAB(PFPA& mi,PFPA& ma,PFPA term1,PFPA term2) {
	int vgl=FPA_vgl(term1,term2);
	if (vgl>=0) { mi=term2; ma=term1; }
	else { mi=term1; ma=term2; }
}

void minimaxFPAABCD(FPA& mi,FPA& ma,const FPA term1,const FPA term2,const FPA term3,const FPA term4) {
	FPA miab,maab,micd,macd;

	minimaxFPAAB(miab,maab,term1,term2);
	minimaxFPAAB(micd,macd,term3,term4);
	if (FPA_vgl(miab,micd) <= 0) mi=miab; else mi=micd;
	if (FPA_vgl(maab,macd) >= 0) ma=maab; else ma=macd;
}

void minimaxFPAABCD(PFPA& mi,PFPA& ma,PFPA term1,PFPA term2,PFPA term3,PFPA term4) {
	PFPA miab,maab,micd,macd;

	minimaxFPAAB(miab,maab,term1,term2);
	minimaxFPAAB(micd,macd,term3,term4);
	if (FPA_vgl(miab,micd) <= 0) mi=miab; else mi=micd;
	if (FPA_vgl(maab,macd) >= 0) ma=maab; else ma=macd;
}

void minimaxFPAABCD(PFPA& mi,PFPA& ma,FPA* array4) {
	PFPA miab,maab,micd,macd;

	minimaxFPAAB(miab,maab,&array4[0],&array4[1]);
	minimaxFPAAB(micd,macd,&array4[2],&array4[3]);
	if (FPA_vgl(miab,micd) <= 0) mi=miab; else mi=micd;
	if (FPA_vgl(maab,macd) >= 0) ma=maab; else ma=macd;
}

inline void FPA_sub_abs_ZAB(FPA& erg,FPA* a,FPA* b) {
	// only absolute value
	int32_t vgl=FPA_vgl_abs(a,b);
	if (vgl==0) {
		erg.setNull();
		return;
	}
	
	if (vgl>0) {
		// |a| > |b|
		FPA_sub_abs_ovgl_ZAB(erg,a,b);
		erg.vorz=1;
	} else {
		// |a| < |b|
		FPA_sub_abs_ovgl_ZAB(erg,b,a);
		erg.vorz=-1;
	}
}

inline void FPA_sub_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if (term1.vorz == 0) {
		if (term2.vorz == 0) {
			erg.setNull();
		} else {
			erg.copyFrom(term2);
			erg.vorz*=-1;
		}
		return;
	} else if (term2.vorz == 0) {
		erg.copyFrom(term1);
		// term1 != 0
		return;
	}

	if (term1.vorz>0) {
		// a>0
		if (term2.vorz>0) {
			// a>0, b>0
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
}

inline void FPA_sub_ZAB(FPA& erg,FPA* term1,FPA* term2) {
	if (term1->vorz == 0) {
		if (term2->vorz == 0) {
			erg.setNull();
		} else {
			erg.copyFrom(term2);
			erg.vorz*=-1;
			// checkNull not necessary as term2 is not zero
		}
		return;
	} else if (term2->vorz == 0) {
		erg.copyFrom(term1);
		return;
	}

	if (term1->vorz>0) {
		// a>=0
		if (term2->vorz>0) {
			// a>=0, b>0
			FPA_sub_abs_ZAB(erg,term1,term2);
		} else {
			// a>=0, b<0: a+|b|
			FPA_add_abs_ZAB(erg,term1,term2);
		}
	} else {
		// a<0
		if (term2->vorz>0) {
			// a<0,b>=0:  -|a|-|b| = -(|a|+|b|)
			FPA_add_abs_ZAB(erg,term1,term2);
			erg.vorz=-1;
		} else {
			// a<0,b<0: -|a| - -|b| = |b|-|a|
			FPA_sub_abs_ZAB(erg,term2,term1);
		}
	}
}

inline void FPA_sub2tpos_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	// b is assumed positive (sign ignored=

	if (term1.vorz == 0) {
		if (term2.vorz == 0) {
			erg.setNull();
			return;
		} else {
			erg.copyFrom(term2);
			erg.vorz=-1;
			return;
		}
	} else if (term2.vorz == 0) {
		erg.copyFrom(term1);
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
}

inline void FPA_sub2tpos_ZAB(FPA& erg,FPA* term1,FPA* term2) {
	// b is assumed positive (sign ignored=

	if (term1->vorz == 0) {
		if (term2->vorz == 0) {
			erg.setNull();
		} else {
			erg.copyFrom(term2);
			erg.vorz=-1;
		}
		return;
	} else if (term2->vorz == 0) {
		erg.copyFrom(term1);
		return;
	}

	if (term1->vorz>0) {
		// a>=0
		FPA_sub_abs_ZAB(erg,term1,term2);
	} else {
		// a<0
		FPA_add_abs_ZAB(erg,term1,term2);
		erg.vorz=-1;
	}
}

inline void FPA_add_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if (term1.vorz==0) {
		if (term2.vorz==0) {
			erg.setNull();
			return;
		} else {
			erg.copyFrom(term2);
			return;
		}
	} else 
	if (term2.vorz==0) {
		erg.copyFrom(term1);
		// is not zero
		return;
	}
	
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
}

inline void FPA_add_ZAB(FPA& erg,FPA* term1,FPA* term2) {
	if (term1->vorz==0) {
		if (term2->vorz==0) {
			erg.setNull();
			return;
		} else {
			erg.copyFrom(term2);
			return;
		}
	} else 
	if (term2->vorz==0) {
		erg.copyFrom(term1);
		return;
	}
	
	// check sign
	if (term1->vorz == term2->vorz) {
		FPA_add_abs_ZAB(erg,term1,term2);
		erg.vorz=term1->vorz;
		return; 
	} else {
		if ( (term1->vorz > 0) && (term2->vorz < 0) ) {
			// erg=a-|b|
			FPA_sub2tpos_ZAB(erg,term1,term2);
		} else {
			// a < 0, b > 0
			FPA_sub2tpos_ZAB(erg,term2,term1);
		}
	}
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

inline void FPA_mul_ZAB(FPA& erg,const FPA term1,const FPA term2) {
	if ( (term1.vorz==0) || (term2.vorz==0) ) {
		erg.setNull();
		return;
	}
	
	// result cannot be zero
	
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
	uint64_t tmp;
	erg.setNull();
	w.a=w.b=w.c=w.d=0;
	w.vorz=1; // not 0, because later digits are set into

	// R^6
	FPA testerg;
	testerg.a=testerg.b=testerg.c=testerg.d=0;
	testerg.vorz=1;
	
	if ((term1.d!=0)&&(term2.d!=0)) {
		//d*h*R^6 => R^3
		tmp=(uint64_t)term1.d*(uint64_t)term2.d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	w.d=0;
	if ((term1.c!=0)&&(term2.d!=0)) {
		//c*h*R^5 + 
		tmp=(uint64_t)term1.c*(uint64_t)term2.d;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	if ((term1.d!=0)&&(term2.c!=0)) {
		//d*g*R^5
		tmp=(uint64_t)term1.d*(uint64_t)term2.c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	w.c=0;
	if ((term1.b!=0)&&(term2.d!=0)) {
		// b*h*R^4 
		tmp=(uint64_t)term1.b*(uint64_t)term2.d;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if ((term1.c!=0)&&(term2.c!=0)) {
		// + c*g*R^4 
		tmp=(uint64_t)term1.c*(uint64_t)term2.c;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if ((term1.d!=0)&&(term2.b!=0)) {
		// + d*f*R^4
		tmp=(uint64_t)term1.d*(uint64_t)term2.b;
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
		tmp=(uint64_t)term1.d*(uint64_t)term2.a;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	
	// a*h*R^3 + 
	if ( (term1.a!=0) && (term2.d!=0) ) {
		tmp=(uint64_t)term1.a*(uint64_t)term2.d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// b*g*R^3 + 
	if ( (term1.b!=0) && (term2.c!=0) ) {
		tmp=(uint64_t)term1.b*(uint64_t)term2.c;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// c*f*R^3 + 
	if ( (term1.c!=0) && (term2.b!=0) ) {
		tmp=(uint64_t)term1.c*(uint64_t)term2.b;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R^2
	w.d=0;
	// a*g*R^2 + 
	if ( (term1.a!=0) && (term2.c!=0) ) {
		tmp=(uint64_t)term1.a*(uint64_t)term2.c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// b*f*R^2 + 
	if ( (term1.b!=0) && (term2.b!=0) ) {
		tmp=(uint64_t)term1.b*(uint64_t)term2.b;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	// c*e*R^2
	if ( (term1.c!=0) && (term2.a!=0) ) {
		tmp=(uint64_t)term1.c*(uint64_t)term2.a;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R
	w.c=0;
	// a*f*R 
	if ( (term1.a!=0) && (term2.b!=0) ) {
		tmp=(uint64_t)term1.a*(uint64_t)term2.b;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	// + b*e*R
	if ( (term1.b!=0) && (term2.a!=0) ) {
		tmp=(uint64_t)term1.b*(uint64_t)term2.a;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
		
	// integer-Part
	w.b=0;
	// a*e
	if ( (term1.a!=0) && (term2.a!=0) ) {
		tmp=(uint64_t)term1.a*(uint64_t)term2.a;
		w.a=tmp & MAXDDBYTE;
		if ( (tmp >> 32) != 0) {
			LOGMSG("Overflow FPA mul.\n");
			exit(99);
		}
		FPA_add_ZAB(erg,erg,w);
	}
				
	if (term1.vorz==term2.vorz) erg.vorz=1;
	else erg.vorz=-1;
}

inline void FPA_mul_ZAuvlong(FPA& erg,const FPA term,const uint64_t intmul) {
	if ( (intmul==0) || (term.vorz==0) ) {
		erg.setNull();
		return;
	}
	
	// result cannot be zero
	
	uint64_t w=(uint64_t)term.d * (uint64_t)intmul;
	erg.d=(w & MAXDDBYTE);
	w=(uint64_t)term.c * (uint64_t)intmul + (w >> 32);
	erg.c=(w & MAXDDBYTE);
	w=(uint64_t)term.b * (uint64_t)intmul + (w >> 32);
	erg.b=(w & MAXDDBYTE);
	w=(uint64_t)term.a * (uint64_t)intmul + (w >> 32);
	
	if (w > MAXDDBYTE) {
		LOGMSG("mul_zap out of range\n");
		exit(99);
	}
	erg.a=(w & MAXDDBYTE);
	if (intmul<0) erg.vorz=-term.vorz;
	else erg.vorz=term.vorz;
}

inline int64_t floorFPA(const FPA term) {
	if (term.vorz==0) return 0;

	if (term.vorz>0) {
		return (int64_t)term.a;
	}

	// negative
	// add one
	if (
		(term.b != 0) ||
		(term.c != 0) ||
		(term.d != 0)
	) return (-(int64_t)term.a-(int64_t)1);

	return -(int64_t)term.a;
}

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

inline void FPA::squareTo(FPA& erg) {
	// works faster than simple multipliucation
	if (vorz==0) {
		erg.setNull();
		return;
	}
	
	// result cannot be zero
	
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
	uint64_t tmp;
	erg.setNull();
	w.a=w.b=w.c=w.d=0;
	w.vorz=1; // not 0, as values are set into w later

	// R^6
	FPA testerg;
	testerg.a=testerg.b=testerg.c=testerg.d=0;
	testerg.vorz=1;
	
	if (d!=0) {
		// d^2R^6 => R^3 locally here in testerg
		tmp=(uint64_t)d*(uint64_t)d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	w.d=0;
	if ((c!=0)&&(d!=0)) {
		//2cdR^5: c < 2^32, d <2^32
		tmp=(uint64_t)c*(uint64_t)d;
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
		tmp=(uint64_t)b*(uint64_t)d;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		w.shiftLeft(1);
		FPA_add_ZAB(testerg,testerg,w);
		// only a and b can be set
	}

	if (c!=0) {
		// + c^2*R^4 
		tmp=(uint64_t)c*(uint64_t)c;
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
		tmp=2*(uint64_t)d*(uint64_t)a;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	
	if ( (b!=0) && (c!=0) ) {
		// + 2bc*R^3
		tmp=2*(uint64_t)b*(uint64_t)c;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R^2
	w.d=0;
	if ( (a!=0) && (c!=0) ) {
		//2acR^2 + 
		tmp=2*(uint64_t)a*(uint64_t)c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	if (b!=0) {
		// b^2*R^2
		tmp=(uint64_t)b*(uint64_t)b;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R
	w.c=0;
	if ( (a!=0) && (b!=0) ) {
		// 2abR
		tmp=2*(uint64_t)a*(uint64_t)b;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
		
	// integer-Part
	w.b=0;
	if (a!=0) {
		//a^2
		tmp=(uint64_t)a*(uint64_t)a;
		w.a=tmp & MAXDDBYTE;
		if ( (tmp >> 32) != 0) {
			FEHLER("a overflow");
		}
		FPA_add_ZAB(erg,erg,w);
	}
				
	erg.vorz=1; 
}

inline void FPA_mul_ZAB(FPA& erg,FPA* term1,FPA* term2) {
	if ( (term1->vorz==0) || (term2->vorz==0) ) {
		erg.setNull();
		return;
	}
	
	// result cannot be zero
	
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
	uint64_t tmp;
	erg.setNull();
	w.a=w.b=w.c=w.d=0;
	w.vorz=1; // not 0, because later digits are set into

	// R^6
	FPA testerg;
	testerg.a=testerg.b=testerg.c=testerg.d=0;
	testerg.vorz=1;
	
	if ((term1->d!=0)&&(term2->d!=0)) {
		//d*h*R^6 => R^3
		tmp=(uint64_t)term1->d*(uint64_t)term2->d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	w.d=0;
	if ((term1->c!=0)&&(term2->d!=0)) {
		//c*h*R^5 + 
		tmp=(uint64_t)term1->c*(uint64_t)term2->d;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}
	
	if ((term1->d!=0)&&(term2->c!=0)) {
		//d*g*R^5
		tmp=(uint64_t)term1->d*(uint64_t)term2->c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	w.c=0;
	if ((term1->b!=0)&&(term2->d!=0)) {
		// b*h*R^4 
		tmp=(uint64_t)term1->b*(uint64_t)term2->d;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if ((term1->c!=0)&&(term2->c!=0)) {
		// + c*g*R^4 
		tmp=(uint64_t)term1->c*(uint64_t)term2->c;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(testerg,testerg,w);
	}

	if ((term1->d!=0)&&(term2->b!=0)) {
		// + d*f*R^4
		tmp=(uint64_t)term1->d*(uint64_t)term2->b;
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
	if ( (term1->d!=0) && (term2->a!=0) ) {
		tmp=(uint64_t)term1->d*(uint64_t)term2->a;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	
	// a*h*R^3 + 
	if ( (term1->a!=0) && (term2->d!=0) ) {
		tmp=(uint64_t)term1->a*(uint64_t)term2->d;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// b*g*R^3 + 
	if ( (term1->b!=0) && (term2->c!=0) ) {
		tmp=(uint64_t)term1->b*(uint64_t)term2->c;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// c*f*R^3 + 
	if ( (term1->c!=0) && (term2->b!=0) ) {
		tmp=(uint64_t)term1->c*(uint64_t)term2->b;
		w.d=tmp & MAXDDBYTE;
		w.c=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R^2
	w.d=0;
	// a*g*R^2 + 
	if ( (term1->a!=0) && (term2->c!=0) ) {
		tmp=(uint64_t)term1->a*(uint64_t)term2->c;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// b*f*R^2 + 
	if ( (term1->b!=0) && (term2->b!=0) ) {
		tmp=(uint64_t)term1->b*(uint64_t)term2->b;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	// c*e*R^2
	if ( (term1->c!=0) && (term2->a!=0) ) {
		tmp=(uint64_t)term1->c*(uint64_t)term2->a;
		w.c=tmp & MAXDDBYTE;
		w.b=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}

	// R
	w.c=0;
	// a*f*R 
	if ( (term1->a!=0) && (term2->b!=0) ) {
		tmp=(uint64_t)term1->a*(uint64_t)term2->b;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
	// + b*e*R
	if ( (term1->b!=0) && (term2->a!=0) ) {
		tmp=(uint64_t)term1->b*(uint64_t)term2->a;
		w.b=tmp & MAXDDBYTE;
		w.a=(tmp >> 32);
		FPA_add_ZAB(erg,erg,w);
	}
		
	// integer-Part
	w.b=0;
	// a*e
	if ( (term1->a!=0) && (term2->a!=0) ) {
		tmp=(uint64_t)term1->a*(uint64_t)term2->a;
		w.a=tmp & MAXDDBYTE;
		if ( (tmp >> 32) != 0) {
			LOGMSG("Overflow FPA mul.\n");
			exit(99);
		}
		FPA_add_ZAB(erg,erg,w);
	}
				
	if (term1->vorz==term2->vorz) erg.vorz=1;
	else erg.vorz=-1;
}

FPA::FPA(const int64_t a) {
	set_vlong(a);
}

FPA::FPA(const double a) {
	set_double32(a);
}

void FPA::shiftLeft(const int32_t bits) {
	// especially for multiplication in scrcoord_as_lowerleft
	// as scalePixelPerRange is a power of 2
	// sign might be affected to zero
	if (bits==0) return;
	if ((bits<0)||(bits>30)) {
		LOGMSG2("Error. shiftLeft. %i not implemented\n",bits);
	}
	uint64_t tmpd=( (uint64_t)d << bits);
	d=tmpd & MAXDDBYTE;
	uint64_t tmpc=(((uint64_t)c << bits) | (tmpd >> 32));
	c=tmpc & MAXDDBYTE;
	uint64_t tmpb=(((uint64_t)b << bits) | (tmpc >> 32));
	b=tmpb & MAXDDBYTE;
	uint64_t tmpa=(((uint64_t)a << bits) | (tmpb >> 32));
	a=tmpa & MAXDDBYTE;

	if ((tmpa >> 32) != 0) {
		// overflow
		LOGMSG("Implementation error. Shift-FPA Left overflow.\n");
		exit(99);
	}
	
	checkNull();
}

void FPA::shiftRight(const int32_t bits) {
	// especially for division in scrcoord_as_lowerleft
	// as scalePixelPerRange is a power of 2
	// sign might be affected to zero
	
	if (bits==0) return; // nothing
	if ((bits<0) || (bits > 16)) {
		LOGMSG2("Error. shiftRight. Bit %i not implemented.\n",bits);
		exit(99);
	}
	
	uint32_t tmpa=(a & power2[bits]) << (32-bits);
	a >>= bits;
	uint32_t tmpb=(b & power2[bits]) << (32-bits);
	b = (b >> bits) | tmpa;
	uint32_t tmpc=(c & power2[bits]) << (32-bits);
	c = (c >> bits) | tmpb;
	uint32_t tmpd=(d & power2[bits]) << (32-bits);
	d = (d >> bits) | tmpc;
	
	if (tmpd) {
		// underflow
		LOGMSG2("Error. Underflow FPA: %i\n",tmpd);
		char tt[2048];
		LOGMSG2("Result: %s\n",str(tt));
		exit(99);
	}

	checkNull();
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

inline void FPA::copyFrom(FPA* w) {
	vorz=w->vorz;
	a=w->a;
	b=w->b;
	c=w->c;
	d=w->d;
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
	
	int64_t vl=(int64_t)floor(w);
	if (
		(vl > (int64_t)MAXDDBYTE) ||
		(vl < (-(int64_t)MAXDDBYTE) )
	) {
		printf("setdouble32 out of range. FAP %I64d\n",vl);
		exit(99);
	}
	
	a=(vl & MAXDDBYTE);
	w -= vl;
	w *= ZWEIHOCH32;
	
	vl=(int64_t)floor(w);
	b=(vl & MAXDDBYTE);
	w -= vl;
	w *= ZWEIHOCH32;

	vl=(int64_t)floor(w);
	c=(vl & MAXDDBYTE);

	d=0; // sonst sind's zuviele Nachkommastellen
	
	// d ist ja 0
	if (
		(b==0) &&
		(a==0) &&
		(c==0)
	) vorz=0;
}

void FPA::set_vlong(const int64_t avalue) {
	uint64_t wert;
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
		erg.vorz=1; // zero not possible
		uint64_t sumdd=(uint64_t)term1.d+(uint64_t)term2.d;
		uint64_t sumcc=(uint64_t)term1.c+(uint64_t)term2.c;
		uint64_t sumbb=(uint64_t)term1.b+(uint64_t)term2.b;
		uint64_t sumaa=(uint64_t)term1.a+(uint64_t)term2.a;
	
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
		erg.checkNull();
	}
}

inline void FPA_add_abs_ZAB(FPA& erg,FPA* term1,FPA* term2) {
	if (term1->vorz==0) {
		if (term2->vorz==0) erg.setNull();
		else erg.copyFrom(*term2);
	}
	else if (term2->vorz==0) erg.copyFrom(*term1);
	else {
		erg.vorz=1; 
		uint64_t sumdd=(uint64_t)term1->d+(uint64_t)term2->d;
		uint64_t sumcc=(uint64_t)term1->c+(uint64_t)term2->c;
		uint64_t sumbb=(uint64_t)term1->b+(uint64_t)term2->b;
		uint64_t sumaa=(uint64_t)term1->a+(uint64_t)term2->a;
	
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
		erg.checkNull();
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
	
	int64_t w=(int64_t)term1.d-(int64_t)term2.d;
	
	if (w < 0) {
		w += (int64_t)ZWEIHOCH32;
		erg.d=(w & MAXDDBYTE);	
		// carry over
		w=(int64_t)term1.c-((int64_t)term2.c+(int64_t)1);
	} else {
		erg.d=(w & MAXDDBYTE);	
		w=(int64_t)term1.c-(int64_t)term2.c;
	}

	if (w < 0) {
		w += (int64_t)ZWEIHOCH32;
		erg.c=(w & MAXDDBYTE);	
		// 1 Uebertrag
		w=(int64_t)term1.b-((int64_t)term2.b+(int64_t)1);
	} else {
		erg.c=(w & MAXDDBYTE);	
		w=(int64_t)term1.b-(int64_t)term2.b;
	}

	if (w < 0) {
		w += (int64_t)ZWEIHOCH32;
		erg.b=(w & MAXDDBYTE);	
		// carry over
		w=(int64_t)term1.a-((int64_t)term2.a+(int64_t)1);
	} else {
		erg.b=(w & MAXDDBYTE);	
		w=(int64_t)term1.a-(int64_t)term2.a;
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

inline void FPA_sub_abs_ovgl_ZAB(FPA& erg,FPA* term1,FPA* term2) {
	// only absolute part, |a| >= |b|
	if (term1->vorz==0) {
		if (term2->vorz==0) {
			erg.setNull();
			return;
		}
		else {
			// Fehler
			LOGMSG("sub_abs_vgl out of range\n");
			exit(99);
		}
	} else if (term2->vorz==0) {
		erg.copyFrom(term1);
		return;
	}
	
	int64_t w=(int64_t)term1->d-(int64_t)term2->d;
	
	if (w < 0) {
		w += (int64_t)ZWEIHOCH32;
		erg.d=(w & MAXDDBYTE);	
		// carry over
		w=(int64_t)term1->c-((int64_t)term2->c+(int64_t)1);
	} else {
		erg.d=(w & MAXDDBYTE);	
		w=(int64_t)term1->c-(int64_t)term2->c;
	}

	if (w < 0) {
		w += (int64_t)ZWEIHOCH32;
		erg.c=(w & MAXDDBYTE);	
		// 1 Uebertrag
		w=(int64_t)term1->b-((int64_t)term2->b+(int64_t)1);
	} else {
		erg.c=(w & MAXDDBYTE);	
		w=(int64_t)term1->b-(int64_t)term2->b;
	}

	if (w < 0) {
		w += (int64_t)ZWEIHOCH32;
		erg.b=(w & MAXDDBYTE);	
		// carry over
		w=(int64_t)term1->a-((int64_t)term2->a+(int64_t)1);
	} else {
		erg.b=(w & MAXDDBYTE);	
		w=(int64_t)term1->a-(int64_t)term2->a;
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

FPA::FPA() {
	// uninitialized
}
#endif

#ifdef _F107
f107_o operator*(int a,const f107_o b) {
	f107_o w=a;
	f107_o erg=w*b;
	
	return erg;
}

inline void	minimaxF107AB(f107_o& mi,f107_o& ma,const f107_o a,const f107_o b) {
	if (a < b) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

inline void	minimaxF107ABCD(f107_o& mi,f107_o& ma,
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

inline void	minimaxF107AB(Pf107& mi,Pf107& ma,f107_o* a,f107_o* b) {
	if ( (*a) < (*b) ) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

inline void	minimaxF107ABCD(Pf107& mi,Pf107& ma,
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

inline void minimaxdABCD(double& mi,double& ma,
	const double a,const double b,
	const double c,const double d
) {
	double miab,micd,maab,macd;
	minimaxdAB(miab,maab,a,b);
	minimaxdAB(micd,macd,c,d);
	
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

inline void minimaxdAB(double& mi,double& ma,const double a,const double b) {
	if (a < b) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

#ifdef _LONGDOUBLE
inline void minimaxldABCD(long double& mi,long double& ma,
	const long double a,const long double b,
	const long double c,const long double d
) {
	long double miab,micd,maab,macd;
	minimaxldAB(miab,maab,a,b);
	minimaxldAB(micd,macd,c,d);
	
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

inline void minimaxldAB(
	long double& mi,
	long double& ma,
	const long double a,const long double b) {
	if (a < b) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}
#endif

#ifdef _QUADMATH
inline void minimaxQDABCD(__float128& mi,__float128& ma,
	const __float128 a,const __float128 b,
	const __float128 c,const __float128 d
) {
	__float128 miab,micd,maab,macd;
	minimaxQDAB(miab,maab,a,b);
	minimaxQDAB(micd,macd,c,d);
	
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

inline void minimaxQDAB(
	__float128& mi,
	__float128& ma,
	const __float128 a,const __float128 b) {
	if (a < b) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}
#endif

#ifdef _F161
f161_o operator*(int a,const f161_o b) {
	f161_o w=a;
	f161_o erg=w*b;
	
	return erg;
}

inline void minimaxF161AB(f161_o& mi,f161_o& ma,const f161_o a,const f161_o b) {
	if (a < b) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

inline void minimaxF161ABCD(f161_o& mi,f161_o& ma,
	const f161_o a,const f161_o b,
	const f161_o c,const f161_o d
) {
	f161_o miab,micd,maab,macd;
	minimaxF161AB(miab,maab,a,b);
	minimaxF161AB(micd,macd,c,d);
	
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

inline void minimaxF161AB(Pf161& mi,Pf161& ma,f161_o* a,f161_o* b) {
	if ( (*a) < (*b) ) {
		mi=a; ma=b;
	} else {
		mi=b; ma=a;
	}
}

inline void minimaxF161ABCD(Pf161& mi,Pf161& ma,
	f161_o* a,f161_o* b,
	f161_o* c,f161_o* d
) {
	Pf161 miab,micd,maab,macd;
	minimaxF161AB(miab,maab,a,b);
	minimaxF161AB(micd,macd,c,d);
	
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

int32_t makePowerOf2(const double avalue) {
	return (1 << ( (int32_t)ceil(log(ceil(avalue))/log(2.0)) ) );
}

#ifdef _F161
double specialFloor161(const f161_o a) {
	// non-overlapping parts in f161_o
	// a is definitely positive AND a <= 2^32
	// i.e. if a=(b_h,b_m,b_l)
	// and |b_m| < 2^-53*|b_h|
	// and |b_l| < 2^-53*|b_m|
	// highest (order-wise) non-zero value floored
	if (a.c0 != 0) return floor(a.c0);
	if (a.c1 != 0) return floor(a.c1);
	
	return floor(a.c2);
}
#endif

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

#define CALCSHRINKAGE \
{\
	shrinkage=0;\
	int32_t touchesborder=0;\
		\
	for(int32_t yk=0;yk<split;yk++) {\
		Aklein.y0=A.y0 + yk*scale2;\
		Aklein.y1=Aklein.y0 + scale2;\
		\
		for(int32_t xk=0;xk<split;xk++) {\
			Aklein.x0=A.x0 + xk*scale2;\
			Aklein.x1=Aklein.x0 + scale2;\
			\
			getBoundingBoxfA(Aklein,bbxfAklein);\
			\
			if (bbxfAklein.x0 < tight.x0) {\
				tight.x0=bbxfAklein.x0;\
				if (tight.x0 <= bbxfA.x0) {\
					touchesborder=1;\
					break;\
				}\
			}\
			if (bbxfAklein.x1 > tight.x1) {\
				tight.x1=bbxfAklein.x1;\
				if (tight.x1 >= bbxfA.x1) {\
					touchesborder=1;\
					break;\
				}\
			}\
			if (bbxfAklein.y0 < tight.y0) {\
				tight.y0=bbxfAklein.y0;\
				if (tight.y0 <= bbxfA.y0) {\
					touchesborder=1;\
					break;\
				}\
			}\
			if (bbxfAklein.y1 > tight.y1) {\
				tight.y1=bbxfAklein.y1;\
				if (tight.y1 >= bbxfA.y1) {\
					touchesborder=1;\
					break;\
				}\
			}\
		}\
		\
		if (touchesborder>0) break;\
	}\
	\
	if (touchesborder<=0) {\
		ScreenRect scr;\
		scr.x0=scrcoord_as_lowerleft(bbxfA.x0);\
		scr.x1=scrcoord_as_lowerleft(bbxfA.x1);\
		scr.y0=scrcoord_as_lowerleft(bbxfA.y0);\
		scr.y1=scrcoord_as_lowerleft(bbxfA.y1);\
		\
		ScreenRect delta;\
		delta.x0=scrcoord_as_lowerleft(tight.x0) - scr.x0;\
		delta.x1=scr.x1 - scrcoord_as_lowerleft(tight.x1);\
		delta.y0=scrcoord_as_lowerleft(tight.y0) - scr.y0;\
		delta.y1=scr.y1 - scrcoord_as_lowerleft(tight.y1);\
		\
		shrinkage=minimumI(\
			delta.x0,delta.x1,delta.y0,delta.y1\
		);\
	}\
}

// struct HelperAccess
HelperAccess::HelperAccess() {
	helperblocks=NULL;
	blockanz=0;
}

HelperAccess::~HelperAccess() {
	if (helperblocks) delete[] helperblocks;
}
	
void HelperAccess::initMemory(void) {
	blockanz=1 + (SCREENWIDTH >> HELPERPERBLOCKBITS);
	helperblocks=new PHelper[blockanz];
	if (!helperblocks) {
		LOGMSG("Memory error. HelperAccess::initMemory\n");
		exit(99);
	}
	
	for(int32_t i=0;i<blockanz;i++) {
		helperblocks[i]=helpermgr->getMemory(MAXHELPERPERBLOCK);
		// no initialization
	} // i
}

inline PHelper HelperAccess::getHelper(const int32_t acoord) {
	return 
		&helperblocks
		[acoord >> HELPERPERBLOCKBITS]
		[acoord & HELPERPERBLOCKMODULO];
}

void HelperAccess::precompute(const int32_t adir) {
	// uses precompute
	// precompute_helperYdep)(const int32_t,PlaneRect&,Helper*) = NULL;
	// precompute_helperXdep)(const int32_t,PlaneRect&,Helper*) = NULL;
	
	PlaneRect A;
	PHelper helper;
	
	for(int32_t idx=0;idx<SCREENWIDTH;idx++) {
		helper=getHelper(idx);
		if (adir==DIRECTIONX) {
			A.x0=idx*scaleRangePerPixel + COMPLETE0;
			A.x1=A.x0+scaleRangePerPixel;
			precompute_helperXdep(A,helper);
		} else if (adir==DIRECTIONY) {
			A.y0=idx*scaleRangePerPixel + COMPLETE0;
			A.y1=A.y0+scaleRangePerPixel;
			precompute_helperYdep(A,helper);
		}
	} // idx
}

// helperacess_double
HelperAccess_double::HelperAccess_double() {
	helperblocks=NULL;
	blockanz=0;
}

HelperAccess_double::~HelperAccess_double() {
	if (helperblocks) delete[] helperblocks;
}
	
void HelperAccess_double::initMemory(void) {
	blockanz=1 + (SCREENWIDTH >> HELPERPERBLOCKBITS);
	helperblocks=new PHelper_double[blockanz];
	if (!helperblocks) {
		LOGMSG("Memory error. HelperAccess_double::initMemory\n");
		exit(99);
	}
	
	for(int32_t i=0;i<blockanz;i++) {
		helperblocks[i]=helper_doublemgr->getMemory(MAXHELPERPERBLOCK);
		// no initialization
	} // i
}

inline PHelper_double HelperAccess_double::getHelper(const int32_t acoord) {
	return 
		&helperblocks
		[acoord >> HELPERPERBLOCKBITS]
		[acoord & HELPERPERBLOCKMODULO];
}

void HelperAccess_double::precompute(const int32_t adir) {
	PlaneRect_double A;
	PHelper_double helper;
	
	for(int32_t idx=0;idx<SCREENWIDTH;idx++) {
		helper=getHelper(idx);
		if (adir==DIRECTIONX) {
			A.x0=idx*scaleRangePerPixel_double + COMPLETE0_double;
			A.x1=A.x0+scaleRangePerPixel_double;
			precompute_helperXdep_double(A,helper);
		} else if (adir==DIRECTIONY) {
			A.y0=idx*scaleRangePerPixel_double + COMPLETE0_double;
			A.y1=A.y0+scaleRangePerPixel_double;
			precompute_helperYdep_double(A,helper);
		}
	} // idx
}

// main
int32_t fastdtcheck_double(void) {
	// the current NTYP provides by assumption enough
	// precision to store all calculation (intermediate) results during
	// comutation of PlaneRect
	// screenRect with NTYP are compared with 
	// screenRects calculated with double only
	// only gray cells are checked
	
	if (
		(!getBoundingBoxfA_double) ||
		(!getBoundingBoxfA_double_oh)
	) {
		LOGMSG("fastdtcheck: function pointer not defined.\n");
		exit(99);
	}
	
	PlaneRect A,bbxfA;
	PlaneRect_double Ad,bbxfAd;
	ScreenRect scr,scrd;
	
	int32_t noch0=SCREENWIDTH >> 3;
	if (noch0 >= 8192) noch0=8192;
	int32_t noch=1;
	
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
			noch=noch0;
		}
		
		if (!data5->zeilen) continue;
		if (data5->memgrau[y].g1 < data5->memgrau[y].g0) continue;
		
		A.y0=y*scaleRangePerPixel + COMPLETE0;
		A.y1=A.y0+scaleRangePerPixel;
		
		Ad.y0=y*scaleRangePerPixel_double + COMPLETE0_double;
		Ad.y1=Ad.y0+scaleRangePerPixel_double;

		Helper* helperY=helperYdep->getHelper(y);
		Helper_double* helperY_double=helperYdep_double->getHelper(y);

		int32_t xlast=-1;
		
		for(int32_t x=data5->memgrau[y].g0;x<=data5->memgrau[y].g1;x++) {
			int32_t f;
			GET_SINGLE_CELLCOLOR_XY(x,y,f);
			
			if (
				(f != SQUARE_GRAY) &&
				(f != SQUARE_GRAY_POTENTIALLY_WHITE)
			) continue;
			
			if (x==(xlast+1)) {
				A.x0=A.x1;
				Ad.x0=Ad.x1;
				Ad.x1=Ad.x0+scaleRangePerPixel_double;
				#ifdef _FPA
				FPA_add_ZAB(A.x1,&A.x0,&scaleRangePerPixel);
				#else
				A.x1=A.x0+scaleRangePerPixel;
				#endif
			} else {
				#ifdef _FPA
				FPA tmp;
				FPA_mul_ZAuvlong(tmp,scaleRangePerPixel,x);
				FPA_add_ZAB(A.x0,&tmp,&COMPLETE0);
				FPA_add_ZAB(A.x1,&A.x0,&scaleRangePerPixel);
				#else
				A.x0=x*scaleRangePerPixel + COMPLETE0;
				A.x1=A.x0+scaleRangePerPixel;
				#endif
			
				Ad.x0=x*scaleRangePerPixel_double + COMPLETE0_double;
				Ad.x1=Ad.x0+scaleRangePerPixel_double;
			}
			xlast=x;
			
			// POINTHERE
			getBoundingBoxfA_helper(
				A,bbxfA,
				helperXdep->getHelper(x),
				helperY);
			getBoundingBoxfA_double(
				Ad,bbxfAd,
				helperXdep_double->getHelper(x),
				helperY_double);

			const int32_t PLACE_ENTIRELYOUTSIDERANGE=1;
			const int32_t PLACE_INTERSECTSTRUELYOUTSIDERANGE=2;
			const int32_t PLACE_COMPLETELYINRANGE=3;
			
			int32_t place=-1;
			int32_t placed=-1;
			
			if (
				(bbxfA.x0 > COMPLETE1) ||
				(bbxfA.x1 < COMPLETE1) ||
				(bbxfA.y0 > COMPLETE1) ||
				(bbxfA.y1 < COMPLETE1)
			) place=PLACE_ENTIRELYOUTSIDERANGE;
			else 
			if (
				(COMPLETE0 <= bbxfA.x0) &&
				(bbxfA.x1 <= COMPLETE1) &&
				(COMPLETE0 <= bbxfA.y0) &&
				(bbxfA.y1 <= COMPLETE1)
			) place=PLACE_COMPLETELYINRANGE;
			else place=PLACE_INTERSECTSTRUELYOUTSIDERANGE;
			
			if (
				(bbxfAd.x0 > COMPLETE1_double) ||
				(bbxfAd.x1 < COMPLETE1_double) ||
				(bbxfAd.y0 > COMPLETE1_double) ||
				(bbxfAd.y1 < COMPLETE1_double)
			) placed=PLACE_ENTIRELYOUTSIDERANGE;
			else
			if (
				(COMPLETE0_double <= bbxfAd.x0) &&
				(bbxfAd.x1 <= COMPLETE1_double) &&
				(COMPLETE0_double <= bbxfAd.y0) &&
				(bbxfAd.y1 <= COMPLETE1_double)
			) placed=PLACE_COMPLETELYINRANGE;
			else placed=PLACE_INTERSECTSTRUELYOUTSIDERANGE;

			if (
				(place<0) ||
				(placed<0) ||
				(place != placed)
			) {
				LOGMSG("\nINVALID placing fastdtchk\n");
				return 0;
			}
			
			if (place==PLACE_ENTIRELYOUTSIDERANGE) {
				continue;
			}
			
			scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
			scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
			scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
			scr.y1=scrcoord_as_lowerleft(bbxfA.y1);

			scrd.x0=scrcoord_as_lowerleft_double(bbxfAd.x0);
			scrd.x1=scrcoord_as_lowerleft_double(bbxfAd.x1);
			scrd.y0=scrcoord_as_lowerleft_double(bbxfAd.y0);
			scrd.y1=scrcoord_as_lowerleft_double(bbxfAd.y1);
			
			if (
				(scr.x0 != scrd.x0) ||
				(scr.x1 != scrd.x1) ||
				(scr.y0 != scrd.y0) ||
				(scr.y1 != scrd.y1)
			) {
				LOGMSG("INVALID screenRects\n");
				return 0;
			}
		} // x
	} // y
	
	// analyze reverse cell graph
	printf("\n  validating reverse cell graph coordinates ... ");
	const NTYP DD=REVCGBLOCKWIDTH*scaleRangePerPixel;
	const double DDd=REVCGBLOCKWIDTH*scaleRangePerPixel_double;

	noch=1;
	noch0=(SCREENWIDTH >> REVCGBITS) >> 3;
	for(int32_t y=0;y<SCREENWIDTH;y+=REVCGBLOCKWIDTH) {
		if ((--noch)<=0) {
			printf("%i ",y);
			noch=noch0;
		}

		A.y0=y*scaleRangePerPixel + COMPLETE0;
		A.y1=A.y0+DD;
		Ad.y0=y*scaleRangePerPixel_double + COMPLETE0_double;
		Ad.y1=Ad.y0+DDd;

		for(int32_t x=0;x<SCREENWIDTH;x+=REVCGBLOCKWIDTH) {
			int8_t hasgray=0;
			for(int32_t y2=y;y2<(y+REVCGBLOCKWIDTH);y2++) {
				if (!data5->zeilen[y2]) continue;
			
				int32_t xe=x+REVCGBLOCKWIDTH-1;
				
				if (
					(xe < data5->memgrau[y2].g0) ||
					(x > data5->memgrau[y2].g1)
				) {
					// not overlapping
				} 
				else {
					hasgray=1;
					break;
				}
			}
			
			if (hasgray<=0) continue;

			A.x0=x*scaleRangePerPixel + COMPLETE0;
			A.x1=A.x0+DD;
			Ad.x0=x*scaleRangePerPixel_double + COMPLETE0_double;
			Ad.x1=Ad.x0+DDd;

			// no use of helper object here
			// as size of A16 is NOT scaleRangePerPixel
			getBoundingBoxfA(A,bbxfA);
			getBoundingBoxfA_double_oh(Ad,bbxfAd);
			
			if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA) > 0) {
				if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfAd) > 0) {
					// correct
					continue;
				} else {
					LOGMSG("INVALID. Reverse cell graph not correct in double numbertype\n");
					return 0;
				}
			}
			
			// now only compare part in COMPLETE

			ScreenRect scr,scrd;
			// scroord_as_lowerleft trims to SCREENWIDTH
			
			scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
			scr.x0 >>= REVCGBITS;
			scrd.x0=scrcoord_as_lowerleft_double(bbxfAd.x0);
			scrd.x0 >>= REVCGBITS;
			if (scr.x0 != scrd.x0) return 0; // not usable

			scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
			scr.x1 >>= REVCGBITS;
			scrd.x1=scrcoord_as_lowerleft_double(bbxfAd.x1);
			scrd.x1 >>= REVCGBITS;
			if (scr.x1 != scrd.x1) return 0; // not usable

			scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
			scr.y0 >>= REVCGBITS;
			scrd.y0=scrcoord_as_lowerleft_double(bbxfAd.y0);
			scrd.y0 >>= REVCGBITS;
			if (scr.y0 != scrd.y0) return 0; // not usable

			scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
			scr.y1 >>= REVCGBITS;
			scrd.y1=scrcoord_as_lowerleft_double(bbxfAd.y1);
			scrd.y1 >>= REVCGBITS;
			if (scr.y1 != scrd.y1) return 0; // not usable
		}
	} // y
	
	return 1;
}

// streakArray
void StreakArray::fastEmpty(void) {
	writeidx=0;
	writenextpos=0;
	for(int32_t i=0;i<STREAK_MAXBLOCKS;i++) anz[i]=0;
}

void StreakArray::pushStreak(const int32_t ax0,const int32_t ax1,const int32_t ay) {
	if (writenextpos >= STREAK_PERBLOCK) {
		writeidx++;
		writenextpos=0;
	}

	if (writeidx >= STREAK_MAXBLOCKS) {
		LOGMSG("Memory-error. pushStreak full\n");
		exit(99);
	}
	if (!ptr[writeidx]) {
		ptr[writeidx]=new Streak[STREAK_PERBLOCK];
		writenextpos=0;
	}
	
	ptr[writeidx][writenextpos].x0=ax0;
	ptr[writeidx][writenextpos].x1=ax1;
	ptr[writeidx][writenextpos].y=ay;;
	writenextpos++;
	anz[writeidx]=writenextpos;
}

void StreakArray::popStreak(Streak& str) {
	if (!ptr[writeidx]) {
		LOGMSG2("Implementation error. StreakArray.pop %i\n",writeidx);
		exit(99);
	}
	if (writenextpos==0) {
		writeidx--;
		if (writeidx<0) {
			writeidx=0;
			writenextpos=0;
			str.y=-1;
			return;
		}
		anz[writeidx]--;
		writenextpos=anz[writeidx];
		str.x0=ptr[writeidx][writenextpos].x0;
		str.x1=ptr[writeidx][writenextpos].x1;
		str.y=ptr[writeidx][writenextpos].y;
		return;
	}
	
	anz[writeidx]--;
	writenextpos--;
	str.x0=ptr[writeidx][writenextpos].x0;
	str.x1=ptr[writeidx][writenextpos].x1;
	str.y=ptr[writeidx][writenextpos].y;
}

StreakArray::StreakArray() {
	writeidx=0;
	writenextpos=0;
	for(int32_t i=0;i<STREAK_MAXBLOCKS;i++) {
		ptr[i]=NULL;
		anz[i]=0;
	}
}

StreakArray::~StreakArray() {
	for(int32_t i=0;i<STREAK_MAXBLOCKS;i++) {
		delete[] ptr[i];
	}
}

int8_t validateInterior(void) {
	// check if every black pixel lands in only black pixels
	PlaneRect A,bbxfA;
	
	int32_t noch=1;
	int32_t noch0=SCREENWIDTH >> 3;
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-y);
			noch=noch0;
		}
		if (!data5->zeilen[y]) continue; // all white
		
		A.y0=y*scaleRangePerPixel + COMPLETE0;
		A.y1=A.y0+scaleRangePerPixel;
		
		Helper* helperY=helperYdep->getHelper(y);

		for(int32_t x=0;x<SCREENWIDTH;x++) {
			int32_t f;
			GET_SINGLE_CELLCOLOR_XY(x,y,f);
			
			if (f != SQUARE_BLACK) continue;
			
			A.x0=x*scaleRangePerPixel + COMPLETE0;
			A.x1=A.x0+scaleRangePerPixel;
			
			getBoundingBoxfA_helper(
				A,bbxfA,
				helperXdep->getHelper(x),
				helperY);
			
			if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
				return 0; 
			}
			
			// bbx completely in screen
			ScreenRect scr;
			scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
			scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
			scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
			scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
			
			for(int32_t ty=scr.y0;ty<=scr.y1;ty++) {
				for(int32_t tx=scr.x0;tx<=scr.x1;tx++) {
					int32_t flokal;
					GET_SINGLE_CELLCOLOR_XY(tx,ty,flokal);
					if (flokal != SQUARE_BLACK) {
						return 0;
					}
				} // tx
			} // ty
		} // x
	} // y
	
	return 1;
}

void Data5::precomputeScreenRect(void) {
	// precomputing ScreenRects up to a specific memory consumption
	int64_t maxmemory=_PRECOMPUTEBBXMEMORYGB * ( (int64_t)1 << 30 );
	// in bytes here
	int64_t memoryused=0;
	
	// from the center of the image upwards and
	// downwards one row repeatedly
	
	PlaneRect A,bbxfA;
	pcscr=new PScreenRect[SCREENWIDTH];
	pcscrmgr=new ScreenRectManager;
	int32_t noch=1;
	int32_t noch0=SCREENWIDTH >> 2;

	for(int32_t i=0;i<SCREENWIDTH;i++) pcscr[i]=NULL;
	printf("precomputing some screenRects ... ");
	
	for(int32_t threshold=75;threshold>=0;threshold-=25) {
		printf("%i ",threshold);
		for(int32_t y=0;y<SCREENWIDTH;y++) {
			if ((--noch)<=0) {
				printf(".");
				noch=noch0;
			}

			if (pcscr[y]) continue;
			if (!zeilen[y]) continue; // empty 
			if (memgrau[y].g1 < memgrau[y].g0) continue;
			
			// only specific gray density cells
			if (graudensity[y] < threshold) continue;
		
			int64_t touse=((int64_t)memgrau[y].g1-memgrau[y].g0+1);
			memoryused += (touse*sizeof(ScreenRect));
			if (memoryused > maxmemory) break;
			pcscr[y]=pcscrmgr->getMemory(touse);
			
			if (!pcscr[y]) {
				LOGMSG("precompute memory exhausted\n");
				memoryused=maxmemory; // exit
				threshold=0;
				break;
			}
		
			A.y0=y*scaleRangePerPixel + COMPLETE0;
			A.y1=A.y0+scaleRangePerPixel;
			Helper* helperY=helperYdep->getHelper(y);
		
			for(int32_t x=memgrau[y].g0;x<=memgrau[y].g1;x++) {
				int32_t f;
				GET_SINGLE_CELLCOLOR_XY(x,y,f);
				if (f != SQUARE_GRAY) {
					SETPCSCR(x,y,-2,-2,-2,-2);
					continue;
				}

				A.x0=x*scaleRangePerPixel + COMPLETE0;
				A.x1=A.x0+scaleRangePerPixel;
				
				getBoundingBoxfA_helper(
					A,bbxfA,
					helperXdep->getHelper(x),
					helperY
				);
				
				if ((SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(bbxfA))>0) {
					SETPCSCR(x,y,-1,-1,0,0);
				} else {
					ScreenRect scr;
					// scr is trimmed to the screen: [0..SCREENWIDT-1]
					scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
					scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
					scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
					scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
					
					// but does it partially overlap with the outside
					if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
						scr.x0 = -(scr.x0+1); // strictly negative
					} 
					SETPCSCR(
						x,y,
						scr.x0,scr.x1,scr.y0,scr.y1
					);
				} // not totally outside
			}
		} // y
		
		if (memoryused >= maxmemory) {
			break;
		}
	} // while

	printf(" used %I64d GB\n",
		(memoryused >> 30)+1);
}

int32_t main(int32_t argc,char** argv) {
	int32_t c0=clock();
	
	flog=fopen("juliatsacoredyn.log.txt","at");
	fprintf(flog,"\n-----------------\n");
	printf("juliatsacoredyn\n");

	#ifdef _USESPLITGRIDBBX
	LOGMSG2("_USESPLITGRIDBBX %i\n",_SPLITBBX);
	#endif
	
	parentmgr=new ParentManager;
	int32_t cmd=CMD_CALC;
	getBoundingBoxfA=getBoundingBoxfA_z2c;
	_FUNC=FUNC_Z2C;
	RANGE0=-2.0;
	RANGE1=2.0;
	seedC0re=seedC1re=-1.0;
	seedC0im=seedC1im=0.0;
	FAKTORAre=FAKTORAim=0.0;
	REVCGBITS=4;
	SCREENWIDTH=(1 << 10);
	_RESETPOTW=0;
	
	for(int32_t i=1;i<argc;i++) {
		upper(argv[i]);
		if (strstr(argv[i],"FUNC=")==argv[i]) {
			_FUNC=getfuncidx(&argv[i][5]);
		} else
		if (strstr(argv[i],"CMD=")==argv[i]) {
			if (strstr(&argv[i][4],"PERIOD")==&argv[i][4]) {
				cmd=CMD_PERIOD;
				if (strstr(argv[i],",PP")) _PERIODICPOINTS=1;
				if (strstr(argv[i],",M2")) _PERIODICITYMETHOD=2;
				if (strstr(argv[i],",M3")) _PERIODICITYMETHOD=3;
			} else
			if (strstr(&argv[i][4],"CONVERT")==&argv[i][4]) {
				convert_raw_into_newstructure();
				fclose(flog);
				return 0;
			} else
			if (strstr(&argv[i][4],"FASTDTCHK")==&argv[i][4]) {
				cmd=CMD_FASTDTCHECK;
			} 
		} else 
		if (strstr(argv[i],"PRECOMPUTE=")) {
			int a;
			if (sscanf(&argv[i][11],"%i",&a) == 1) {
				if (a<1) a=1;
				_PRECOMPUTEBBXMEMORYGB=a;
			}
		} else
		if (strstr(argv[i],"C=")==argv[i]) {
			double r0,r1,i0,i1; // not NTYP
			// command line parameters are always considered double no matter the datatype used
			if (sscanf(&argv[i][2],"%lf,%lf,%lf,%lf",&r0,&r1,&i0,&i1) == 4) {
				double w1=floor(r0*DENOM225)/DENOM225;
				double w2=floor(r1*DENOM225)/DENOM225;
				if (w1 > w2) {
					seedC0re=w2;
					seedC1re=w1;
					seedC0re_double=w2;
					seedC1re_double=w1;
				} else {
					seedC0re=w1;
					seedC1re=w2;
					seedC0re_double=w1;
					seedC1re_double=w2;
				}
				w1=floor(i0*DENOM225)/DENOM225;
				w2=floor(i1*DENOM225)/DENOM225;
				if (w1 > w2) {
					seedC0im=w2;
					seedC1im=w1;
					seedC0im_double=w2;
					seedC1im_double=w1;
				} else {
					seedC0im=w1;
					seedC1im=w2;
					seedC0im_double=w1;
					seedC1im_double=w2;
				}
			} else
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				double ra=floor(r0*DENOM225)/DENOM225;;
				double ib=floor(i0*DENOM225)/DENOM225;;
				seedC0re=seedC1re=ra;
				seedC0im=seedC1im=ib;
				seedC0re_double=seedC1re_double=ra;
				seedC0im_double=seedC1im_double=ib;
			}
		} else if (strstr(argv[i],"CD=")==argv[i]) {
			int r0,r1,i0,i1; 
			// command line parameters are always considered double no matter the datatype used
			if (sscanf(&argv[i][3],"%i,%i,%i,%i",&r0,&r1,&i0,&i1) == 4) {
				double w1=(double)r0/DENOM225;
				double w2=(double)r1/DENOM225;
				if (w1 > w2) {
					seedC0re=w2;
					seedC1re=w1;
					seedC0re_double=w2;
					seedC1re_double=w1;
				} else {
					seedC0re=w1;
					seedC1re=w2;
					seedC0re_double=w1;
					seedC1re_double=w2;
				}
				w1=(double)i0/DENOM225;
				w2=(double)i1/DENOM225;
				if (w1 > w2) {
					seedC0im=w2;
					seedC1im=w1;
					seedC0im_double=w2;
					seedC1im_double=w1;
				} else {
					seedC0im=w1;
					seedC1im=w2;
					seedC0im_double=w1;
					seedC1im_double=w2;
				}
			} else
			if (sscanf(&argv[i][3],"%i,%i",&r0,&i0) == 2) {
				double ra=(double)r0/DENOM225;
				double ib=(double)i0/DENOM225;;
				seedC0re=seedC1re=ra;
				seedC0im=seedC1im=ib;
				seedC0re_double=seedC1re_double=ra;
				seedC0im_double=seedC1im_double=ib;
			}
		} else if (strstr(argv[i],"PROP=")==argv[i]) {
			_PROPAGATEDEF=0;
			_PROPAGATEPOTW=0;

			if (strstr(argv[i],"DEF")) _PROPAGATEDEF=1;
			if (strstr(argv[i],"POTW")) _PROPAGATEPOTW=1;
		} else if (strstr(argv[i],"A=")==argv[i]) {
			double r0,i0;
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				r0=floor(r0*DENOM225)/DENOM225;
				i0=floor(i0*DENOM225)/DENOM225;

				#ifdef _FPA
				FAKTORAre.set_double32(r0);
				FAKTORAim.set_double32(i0);
				#else
				FAKTORAre=r0;
				FAKTORAim=i0;
				#endif

				FAKTORAre_double=r0;
				FAKTORAim_double=i0;
			}
		} else if (strstr(argv[i],"AD=")==argv[i]) {
			int a,b;
			if (sscanf(&argv[i][3],"%i,%i",&a,&b) == 2) {
				double rw=(double)a/DENOM225;
				double iw=(double)b/DENOM225;

				FAKTORAre=rw;
				FAKTORAim=iw;
				FAKTORAre_double=rw;
				FAKTORAim_double=iw;
			}
		} else
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
				if (a<0) {
					// a is a negative exponent
					double w=1.0 / (double)((int64_t)1 << (-a) );
					RANGE0=-w;
					RANGE1=w;
				} else {
					RANGE1=makePowerOf2(a);
					if (RANGE1 != a) LOGMSG2("range adjusted to next-bigger power of 2: %.2lg\n",RANGE1);
					RANGE0=-RANGE1;
				}
				
				printf("RANGE: %lg\n",RANGE1);
			}
		}
	} // i
	
	COMPLETE0_double=RANGE0;
	COMPLETE1_double=RANGE1;
	#ifdef _FPA
	COMPLETE0.set_vlong(RANGE0);
	COMPLETE1.set_vlong(RANGE1);
	#else
	COMPLETE0=RANGE0;
	COMPLETE1=RANGE1;
	#endif
	
	if (SCREENWIDTH < 256) SCREENWIDTH=256;
	if (REVCGBITS < 4) REVCGBITS=4;
	
	REFINEMENTLEVEL=(int)ceil(log(SCREENWIDTH)/log(2.0));
	
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
	scaleRangePerPixel_double=w;

	w=(double)SCREENWIDTH / (RANGE1-RANGE0);
	scalePixelPerRange=w;
	scalePixelPerRange_double=w;
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

		// squares whose bounding box lies completely in the special exterior
		find_special_exterior_hitting_squares();
	} else {
		// add/substract 16 for safetly
		plane.x0=(encgrayx0-16) * scaleRangePerPixel + COMPLETE0;
		plane.x1=(encgrayx1+16) * scaleRangePerPixel + COMPLETE0;
		plane.y0=(encgrayy0-16) * scaleRangePerPixel + COMPLETE0;
		plane.y1=(encgrayy1+16) * scaleRangePerPixel + COMPLETE0;
	}
	
	#ifdef _FPA
	LOGMSG5("  roughly %.20lg..%.20lg x %.5lg..%.5lg used\n",
		plane.x0.convert_to_double(),
		plane.x1.convert_to_double(),
		plane.y0.convert_to_double(),
		plane.y1.convert_to_double());
	#else
	LOGMSG5("  roughly %.20lg..%.20lg x %.20lg..%.20lg used\n",
		(double)plane.x0,(double)plane.x1,(double)plane.y0,(double)plane.y1);
	#endif
	
	if (
		(encgrayx0 > (SCREENWIDTH >> 2) ) &&
		(encgrayx1 < (3*(SCREENWIDTH >> 2) ) ) &&
		(encgrayy0 > (SCREENWIDTH >> 2) ) &&
		(encgrayy1 < (3*(SCREENWIDTH >> 2) ) ) 
	) {
		LOGMSG5("  gray in pixel region [%i..%i] x [%i..%i]\n",
			encgrayx0,encgrayx1,encgrayy0,encgrayy1);
		LOGMSG("  range could be adjusted (half is enough)\n");
	}
	
	char tmp[4096];
	fprintf(flog,"%s * 2^-%i\n",seedCstr225(tmp),BASEDENOMINATOR);
	fprintf(flog,"(if needed): %s * 2^-%i\n",FAKTORAstr225(tmp),BASEDENOMINATOR);
	
	#define CLOCK1 \
	{\
		int32_t c1=clock();\
		LOGMSG2("\nduration %.0lf sec\n",(double)(c1-c0)/CLOCKS_PER_SEC);\
	}
	
	if (!getBoundingBoxfA_helper) {
		LOGMSG("Error. No helper bbx function defined.\n");
		exit(99);
	}
	helpermgr=new HelperManager;
	helperYdep=new HelperAccess;
	helperXdep=new HelperAccess;
	helperYdep->initMemory();
	helperXdep->initMemory();
	printf("precomputing sub-expressions ... Y ");
	helperYdep->precompute(DIRECTIONY);
	printf("X\n");
	helperXdep->precompute(DIRECTIONX);

	if (cmd==CMD_FASTDTCHECK) {
		// compute respective _double variants
		helper_doublemgr=new Helper_doubleManager;
		helperYdep_double=new HelperAccess_double;
		helperXdep_double=new HelperAccess_double;
		helperYdep_double->initMemory();
		helperXdep_double->initMemory();
		printf("precomputing number type double sub-expressions ... Y ");
		helperYdep_double->precompute(DIRECTIONY);
		printf("X ");
		helperXdep_double->precompute(DIRECTIONX);
		printf("\n");
	}
	
	if (cmd==CMD_FASTDTCHECK) {
		printf("checking if double can be used instead of %s ... ",NNTYPSTR);
		char tt[2048];
		if (fastdtcheck_double() > 0) {
			LOGMSG("\n  PASSED: At current formula/parameters/level, number type 'double' results in correct screenRects\n");
			sprintf(tt,"__L%i_fastdtcheck_PASSED",REFINEMENTLEVEL);
			FILE *f=fopen(tt,"wt");
			fprintf(f,"\n  Checking DOUBLE against sufficient %s type\n",NNTYPSTR);
			fprintf(f,"\n  PASSED: At current formula/parameters/level number type 'double' results in correct screenRects\n");
			fclose(f);
		} else {
			LOGMSG("\n  FAILED. Using double is discouraged due to rounding errors.\n");
			sprintf(tt,"__L%i_fastdtcheck_FAILED",REFINEMENTLEVEL);
			FILE *f=fopen(tt,"wt");
			fprintf(f,"\n  Checking DOUBLE against sufficient %s type\n",NNTYPSTR);
			fprintf(f,"\n  FAILED. Using double is discouraged due to rounding errors.\n");
			fclose(f);
		}
		CLOCK1
		
		// dirty exit
		return 0;
	}

	// if paramter provided => bbx can be precomputed uop
	// to a certain memory consumption
	if (_PRECOMPUTEBBXMEMORYGB>0) {
		data5->precomputeScreenRect();
	}
	
	// //////////////////////////////////////
	compute(); 
	// //////////////////////////////////////

	if (interiorpresent>0) {
		LOGMSG("\nINTERIOR present\n");
	}
	
	// storing raw data
	printf("saving raw data ... ");
	data5->saveRaw(fn);
	printf("done\n");

	// storing a trustworthily downscaled image
	if (SAVEIMAGE>0) {
		if (SCREENWIDTH > 65536) {
			printf("\nsaving trustworthily downscaled image ... ");
		} else {
			printf("\nsaving image ... ");
		}
		data5->saveBitmap4_twd(fn,-1);
	}
	
	//followallgray(fn);

	// free memory to get enough to allocate for the periodicty check
	printf("freeing non-image memory ...\n");
	freeRevCGMem();

	if (data5->pcscr) {
		delete[] data5->pcscr;
		data5->pcscr=NULL;
	}
	if (data5->pcscrmgr) {
		// frees the allocated screenrects
		delete data5->pcscrmgr;
		data5->pcscrmgr=NULL;
	}

	// data is now computed or loaded
	if (cmd==CMD_PERIOD) {
		if (interiorpresent>0) {
			if (_PERIODICITYMETHOD==3) {
				periodicity_m3(fn);
				// data is no longer valid as potw has
				// a different meaning there
			} else {
				periodicity(fn);
			}
		} else {
			LOGMSG("No interior present. Periodicity check skipped.\n");
		}
	} 
	
	delete data5;
	
	CLOCK1
	
	LOGMSG2("%I64d bounding boxes calculated\n",ctrbbxfa);
	
	delete helpermgr;
	delete helperYdep;
	delete helperXdep;
	if (helper_doublemgr) delete helper_doublemgr;
	if (helperYdep_double) delete helperYdep_double;
	if (helperXdep_double) delete helperXdep_double;
	
	fclose(flog);
	
    return 0;
}

