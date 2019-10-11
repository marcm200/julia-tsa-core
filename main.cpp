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
#include "quadmath.h"


// used floating type
// comment out or in what is needed

#define _DOUBLE
//#define _LONGDOUBLE
//#define _QUADMATH

#ifdef _QUADMATH
typedef __float128 NTYP;
const char NNTYPSTR[]="f128_";
const int _BITSNTYP=113;
#endif

#ifdef  _LONGDOUBLE
typedef long double NTYP;
const char NNTYPSTR[]="ld_";
const int _BITSNTYP=63;
#endif

#ifdef _DOUBLE
typedef double NTYP;
const char NNTYPSTR[]="";
const int _BITSNTYP=53;
#endif

typedef uint16_t DBYTE;
typedef DBYTE *PDBYTE;

const int32_t BLOBOFFSET=16;
const int32_t MAXBLOBS=(1 << 15);

const int32_t DBYTEMAX = (1 << 16);
const uint8_t FIXED2FLAG_INUSE=1;
const uint8_t FIXED2FLAG_TODELETE=2;
const uint8_t FIXED2FLAG_SETTOACTIVECOLOR=4;
const uint8_t FIXED2FLAG_TMPCOLORDA=8;

const int32_t FIXED2_ACTIVETOVISIT=16;
const int32_t COLORSTART=24;


// const definitions

const int32_t MAXPARENT=1024;

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

enum { CMD_CALC=1,CMD_PERIOD,CMD_BASIN };
enum { FUNC_Z2C=1,FUNC_Z2AZC,FUNC_Z3AZC,FUNC_Z4AZC,FUNC_Z5AZC,FUNC_Z6AZC };

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
const int64_t DENOM225=( (int64_t)1 << 25 );

// structs

struct SCCColor {
	DBYTE alt,neu;
	uint8_t flag;
	uint8_t palidx;
};

struct SCCData {
	PDBYTE* scc;
	int32_t grayencx0,grayencx1,grayxlen;
	int32_t grayency0,grayency1,grayylen;
	SCCColor* scccolor;
};

// palette-entry for bitmaps
struct RGB4 {
	uint8_t R,G,B,alpha;
};

// a rectangle in the complex plane - used for a square and its bounding box
struct PlaneRect {
	NTYP x0,x1,y0,y1;
};

// for the reverse cell graph: position of one tile
struct Parent {
	int32_t BX,BY;
};

// one reverse cell graph tile and its parents (preimages before the iteration)
struct RevCGBlock {
	int32_t howmany; 
	// flag, whether this tile has to be checked (again) for gray pixels and the bounding
	// boxes' hits after one iteration
	int32_t tovisit;
	// how many parents are needed here? First pass result
	int32_t memused;
	// parents as variable size array
	Parent* parent;
		
	RevCGBlock();
	void addParent(const int32_t,const int32_t);
};

char* upper(char* s) {
	if (!s) return NULL;
	for(int i=(strlen(s)-1);i>=0;i--) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

// pixel-coordinate rectangle
struct ScreenRect {
	int32_t x0,x1,y0,y1;
};

struct Gray_in_row {
	int32_t g0,g1;
};

struct ColorPalette {
	int32_t anz;
	RGB4* rgbs;

	ColorPalette();
	virtual ~ColorPalette();
	void setlen(const int32_t);
	void setInterval(const double,const double,const int32_t,const int32_t,const int32_t,const int32_t,const int32_t,const int32_t);
	void getColor(const double ,int32_t&,int32_t&,int32_t&);
};

struct Blob {
	ScreenRect rect;
	int32_t mx,my;
	int32_t zyklusnr;
	DBYTE color; 
	int32_t targetblob; 
	uint8_t visited;
};

struct Charmap {
	int64_t xlen,ylen; // würde zwar int reichen, aber da man xlen*y macht, könnte man ausserhalb des Speichers laden
	int64_t memused;
	uint8_t *cmp;
	RGB4 palette[256];
		
	Charmap();
	virtual ~Charmap();
		
	void fill(const uint8_t);
	void clearPalette(void);
	void setlenxy(const int32_t,const int32_t);
	void save(const char*);
	void setPaletteRGB(const int32_t,const uint8_t,const uint8_t,const uint8_t);
	void setP(const int32_t,const int32_t,const uint8_t);
	uint8_t getP(const int32_t,const int32_t);
	void line(const uint8_t,const int32_t,const int32_t,const int32_t,const int32_t,const uint8_t);
};

// main object
struct Data5 {
	uint32_t ** pixelrow;
	Gray_in_row* gray;
	RevCGBlock* revcgYX;
		
	Data5();
	virtual ~Data5();
		
	void saveBitmap4(const char*);
	void saveBitmap4_trustworthily_downscaled_16fold(const char*);
	void saveRaw(const char*);
	int32_t readRawBlowUp(const char*);
};

struct ParentManager {
	Parent* lastallocated;
	int memused;
	int freefrom;
	
	ParentManager();
	// dirty cleanup
	Parent* getParentSpace(const int);
};


// globals

// stores the 
ColorPalette basinpal;
int32_t bitpower=2;
int32_t MAXCOLOREXPECTED=1;
void (*getBoundingBoxfA)(PlaneRect&,PlaneRect&) = NULL;
int _FUNC;
ParentManager parentmgr;
int32_t REFINEMENTLEVEL=0;
int32_t _BUDDHADEPTH=0;
NTYP seedC0re,seedC1re,seedC0im,seedC1im; 
NTYP FAKTORAre,FAKTORAim;
Data5 *data5;
NTYP scaleRangePerPixel,scalePixelPerRange;
int64_t countsquares_white,countsquares_gray;
int64_t countsquares_black,countsquares_graypotw;
int32_t SPLITAFTER;
// region in the complex plane where all the currently still gray squares reside
NTYP planegrayx0,planegrayx1;
NTYP planegrayy0,planegrayy1;
// region in screen coordinates
int32_t encgrayx0,encgrayx1;
int32_t encgrayy0,encgrayy1;
// width of screen in pixel, must be a power of 2
int32_t SCREENWIDTH;
// number of 32bit-integers a screen row holds, equal to (SCREENWIDTH >> 4)
int32_t MEMWIDTH;
int32_t RANGE0=-2,RANGE1=2;
// complte holds the same values as range
NTYP COMPLETE0,COMPLETE1;
// for the low-resolution reverse cell graph working on (usually) 64x64 pixel squares or bigger
int32_t REVCGBITS,REVCGBLOCKWIDTH;
int32_t REVCGmaxnumber,REVCGmaxnumberQ;


// forward declarations

// main function to compute the Julia set
void compute(void);
// constructing the low-level reverse cell-graph
void construct_static_reverse_cellgraph(void);
// squares that directly hit the special exterior outside COMPLETE
void find_special_exterior_hitting_squares(void);
// white and potentially-white cells wiul
void propagate_white(void);
void coloring_interior_black(void);
inline int32_t scrcoord_as_lowerleft(const NTYP);
void copy_pixel_to_2x2grid(const uint32_t,uint32_t*);

void write2(FILE*,const uint8_t,const uint8_t);
void write4(FILE*,const uint8_t,const uint8_t,const uint8_t,const uint8_t);
inline NTYP minimumD(const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP);
inline NTYP minimumD(const NTYP,const NTYP,const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP,const NTYP,const NTYP);


// defines used as small expressions

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

#define SQUARE_LIES_ENTIRELY_IN_COMPLETE(BBX) \
	(\
		(BBX.x1 < COMPLETE1) &&\
		(BBX.x0 > COMPLETE0) &&\
		(BBX.y1 < COMPLETE1) &&\
		(BBX.y0 > COMPLETE0)\
	)

#define SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(BBX) \
	(\
		(BBX.x1 < planegrayx1) &&\
		(BBX.x0 > planegrayx0) &&\
		(BBX.y1 < planegrayy1) &&\
		(BBX.y0 > planegrayy0)\
	)

#define GET_SINGLE_PIXELCOLOR_FROM_4BYTEINTEGER(WW,BPOS) \
	( ( (WW) >> (BPOS) ) & (uint32_t)3 )
	
#define SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(WW,CLEARMASKE,SETF) \
	((WW & CLEARMASKE) | SETF)


// struct definitions

// RevCGBlock
// static reverse cell graph at the given SCREENWIDTH and REVCGBLOCKWIDTH

RevCGBlock::RevCGBlock() {
	howmany=0;
}

void RevCGBlock::addParent(const int ax,const int ay) {
	if (
		(ax < 0) ||
		(ay < 0) ||
		(ax >= REVCGmaxnumber) ||
		(ay >= REVCGmaxnumber) 
	) return; 
	
	if (parent == NULL) {
		// Speicher allokieren
		// geht über einen speziellen ArrayManager
		// die Anzahl wurde vorab bereits ja bestimmt
		parent=parentmgr.getParentSpace(memused);
		howmany=0;
		if (!parent) {
			printf("Memory failure for parent.\n");
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
void Data5::saveRaw(const char* afn) {
	char fn[1024];
	
	// header file
	sprintf(fn,"%s.raw_header",afn);
	FILE *f=fopen(fn,"wb");
	uint32_t w=SCREENWIDTH;
	fwrite(&w,sizeof(w),1,f);
	fwrite(gray,sizeof(Gray_in_row),SCREENWIDTH,f);
	fwrite(&encgrayx0,sizeof(encgrayx0),1,f);
	fwrite(&encgrayx1,sizeof(encgrayx1),1,f);
	fwrite(&encgrayy0,sizeof(encgrayy0),1,f);
	fwrite(&encgrayy1,sizeof(encgrayy1),1,f);
	fclose(f);
	
	// saving the actual pixel data in several files of at most 2 GB size
	int32_t fctr=0;
	uint32_t slen=SPLITAFTER*MEMWIDTH;
	uint32_t spla=SPLITAFTER;

	for(int32_t z=0;z<SCREENWIDTH;z+=SPLITAFTER) {
		fctr++;
		sprintf(fn,"%s.raw_%04i",afn,fctr);
		f=fopen(fn,"wb");
		fwrite(&spla,sizeof(spla),1,f);
		fwrite(pixelrow[z],sizeof(uint32_t),slen,f);
		fclose(f);
	}
	
	fclose(f);
}

int32_t Data5::readRawBlowUp(
	const char* afn) {
	#define ENLARGEL(WW)\
	if (add>=2) {\
		WW <<= 1;\
		WW -= 16;\
		if (WW < 0) WW=0;\
		WW=( (WW >> 4) << 4);\
	}
	#define ENLARGER(WW)\
	if (add>=2) {\
		WW <<= 1;\
		WW += 16;\
		if (WW >= SCREENWIDTH) WW=SCREENWIDTH-16;\
		WW=( (WW >> 4) << 4);\
	}

	char fn[1024];
	sprintf(fn,"%s.raw_header",afn);
	FILE *f=fopen(fn,"rb");
	if (!f) {
		// if file does not exist =>
		// no error, then computation has to start from scratch
		return 0;
	}
	
	printf("reading stored data ...\n");
	
	int32_t qlen,add;

	fread(&qlen,sizeof(qlen),1,f);
	if (qlen == (SCREENWIDTH >> 1)) add=2; 
	else if (qlen == SCREENWIDTH) add=1;
	else {
		fprintf(stderr,"Error. data5::readRawBlowUp\n");
		exit(99);
	}
	
	Gray_in_row g;

	for(int32_t i=0;i<SCREENWIDTH;i+=add) {
		fread(&g,sizeof(Gray_in_row),1,f);
		if (add>=2) {
			ENLARGEL(g.g0)
			ENLARGER(g.g1);
			gray[i].g0=gray[i+1].g0=g.g0;
			gray[i].g1=gray[i+1].g1=g.g1;
		} else {
			gray[i].g0=g.g0;
			gray[i].g1=g.g1;
		}
	}
	
	fread(&encgrayx0,sizeof(encgrayx0),1,f);
	fread(&encgrayx1,sizeof(encgrayx1),1,f);
	fread(&encgrayy0,sizeof(encgrayy0),1,f);
	fread(&encgrayy1,sizeof(encgrayy1),1,f);
	fclose(f);

	// if image is to be doubled in width, gray-coordinates have to be adjusted
	ENLARGEL(encgrayx0)
	ENLARGEL(encgrayy0)
	ENLARGER(encgrayx1)
	ENLARGER(encgrayy1)
	
	planegrayx0=encgrayx0*scaleRangePerPixel + COMPLETE0;
	planegrayy0=encgrayy0*scaleRangePerPixel + COMPLETE0;
	planegrayx1=(encgrayx1+16)*scaleRangePerPixel + COMPLETE0;
	planegrayy1=(encgrayy1+16)*scaleRangePerPixel + COMPLETE0;
	
	uint32_t qpixelrowlen=qlen >> 4;
	uint32_t *qpixelrow=new uint32_t[qpixelrowlen];
	
	int32_t fctr=0;
	int32_t rows_read=0;
	uint32_t zidx=0;
	while (1) {
		fctr++;
		sprintf(fn,"%s.raw_%04i",afn,fctr);
		f=fopen(fn,"rb");
		// when there is no more file existent, exit loop
		if (!f) break;
		uint32_t so_many_rows_in_file;
		fread(&so_many_rows_in_file,sizeof(so_many_rows_in_file),1,f);
		for(uint32_t i=0;i<so_many_rows_in_file;i++) {
			fread(qpixelrow,sizeof(uint32_t),qpixelrowlen,f);
		
			uint32_t ziel[2];
			int32_t mem=0;
			for(uint32_t k=0;k<qpixelrowlen;k++) {
				uint32_t quelle=qpixelrow[k];
				if (add>=2) {
					// any gray-potentially-white will
					// be reduced to gray
					copy_pixel_to_2x2grid(quelle,ziel);
					pixelrow[zidx][mem]=ziel[0];
					pixelrow[zidx][mem+1]=ziel[1];
					pixelrow[zidx+1][mem]=ziel[0];
					pixelrow[zidx+1][mem+1]=ziel[1];
					mem+=2;
				} else {
					// keep current juding including potentially white
					// to continue computation
					pixelrow[zidx][mem]=quelle;
					mem++;
				}
			
			} // setting pixels of one row
			
			zidx+=add;
			rows_read++;
		} // reading all rows in file

		fclose(f);
	} // while next raw split file
	
	if (rows_read != qlen) {
		fprintf(stderr,"Error. readBlowUp2\n");
		exit(99);
	}
	
	delete[] qpixelrow;

	return 1;
}

void Data5::saveBitmap4_trustworthily_downscaled_16fold(const char* afn) {
	// saves a trustworthily downsized version of the image: 16-fold. 
	// image format is: 8 bit Bitmap

	int32_t scr16 = SCREENWIDTH >> 4;
	int32_t bytes_per_row=scr16; 
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
	for(int i=0;i<256;i++) pal[i].R=pal[i].G=pal[i].B=pal[i].alpha=0;
	pal[SQUARE_GRAY].R=127;
	pal[SQUARE_GRAY].G=127;
	pal[SQUARE_GRAY].B=127;
	pal[SQUARE_BLACK].R=0;
	pal[SQUARE_BLACK].G=0;
	pal[SQUARE_BLACK].B=0;
	pal[SQUARE_WHITE].R=255;
	pal[SQUARE_WHITE].G=255;
	pal[SQUARE_WHITE].B=255;

	sprintf(tmp,"%s_twd16.bmp",afn);
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
	fwrite(&scr16,sizeof(scr16),1,fbmp);
	fwrite(&scr16,sizeof(scr16),1,fbmp);
	write2(fbmp,1,0);
	write2(fbmp,8,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,19,11,0,0);
	write4(fbmp,19,11,0,0);
	write4(fbmp,0,1,0,0);
	write4(fbmp,0,0,0,0);

	for(int32_t i=0;i<256;i++) write4(fbmp,pal[i].B,pal[i].G,pal[i].R,pal[i].alpha);
			
	// 16x16 pixels into one: 16 rows and each with one 32bit-integer
	for(int32_t y=0;y<SCREENWIDTH;y+=16) {
		for(int32_t x=0;x<scr16;x++) {
			uint32_t f=pixelrow[y][x];
			// 16 consecutive pixels starting at screen x*16,y
			if (
				(f==SQUARE_GRAY_16_CONSECUTIVE) ||
				(f==SQUARE_GRAYPOTW_16_CONSECUTIVE) ||
				(f==SQUARE_WHITE_16_CONSECUTIVE) ||
				(f==SQUARE_BLACK_16_CONSECUTIVE) 
			) {
				// first row is uniformly colored
				// are the other 15 rows the same color
				for(int32_t y2=1;y2<16;y2++) {
					if (pixelrow[y+y2][x] != f) {
						// no: resulting color is gray
						f=SQUARE_GRAY_16_CONSECUTIVE;
						break;
					}
				}
				switch (f) {
					case SQUARE_BLACK_16_CONSECUTIVE: rgbz[x]=SQUARE_BLACK; break;
					case SQUARE_WHITE_16_CONSECUTIVE: rgbz[x]=SQUARE_WHITE; break;
					default: rgbz[x]=SQUARE_GRAY; break;
				};
				
			} else {
				// if the first pixel row is not of one single color
				// the 16x16 grid cannot be either and
				// hence must be colored gray when downscaled
				rgbz[x]=SQUARE_GRAY;
			}
		} // x
		
		fwrite(rgbz,scr16,sizeof(uint8_t),fbmp);
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
	for(int i=0;i<16;i++) pal[i].R=pal[i].G=pal[i].B=pal[i].alpha=0;
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
					uint32_t w=data5->pixelrow[y][mem];
					
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

	printf("%I64d interior, %I64d exterior, %I64d (%.0lf%%) gray squares\n",
		countsquares_black,
		countsquares_white,
		countsquares_gray+countsquares_graypotw,
		100.0*(double)(countsquares_gray+countsquares_graypotw)/(countsquares_black+countsquares_graypotw+countsquares_gray+countsquares_white));

	if (countsquares_graypotw>0) {
		printf("(of which %I64d are gray but potentially white)\n",countsquares_graypotw);
	}
}

Data5::Data5() {
	printf("initialising main object ...\n");
	
	gray=new Gray_in_row[SCREENWIDTH];
	pixelrow=new uint32_t*[SCREENWIDTH];
	// at most 2 GB continuous 32bit integers
	int64_t gesamt=(1 << (31-2)); 
	gesamt /= (SCREENWIDTH >> 4);
	SPLITAFTER = gesamt; 
	if (SPLITAFTER > SCREENWIDTH) SPLITAFTER=SCREENWIDTH;

	// allokate enough memory to store the necessary amount of 32 bit ints
	// in as few continuous blocks as possible (avoiding allocating memory for each row)
	int32_t still_todo=SCREENWIDTH; 
	int32_t y=0;
	
	while (still_todo>0) {
		int64_t allok=SPLITAFTER;
		if (allok > still_todo) allok=still_todo;
		
		uint32_t *p=new uint32_t[allok*MEMWIDTH];
		if (!p) {
			fprintf(stderr,"Memory error. data5 consstructor\n");
			exit(99);
		}
		
		// setting the points to the array
		// portion that belongs to the row
		for(int32_t i=0;i<SPLITAFTER;i++) {
			pixelrow[y]=&p[i*MEMWIDTH];
			y++;
		}
		
		still_todo -= allok;
	} // while
	
	revcgYX=new RevCGBlock[REVCGmaxnumber*REVCGmaxnumber];
}

Data5::~Data5() {
	for(int32_t i=0;i<SCREENWIDTH;i+=SPLITAFTER) {
		if (pixelrow[i]) delete[] pixelrow[i];
	}
	delete[] revcgYX;
	delete[] gray;
	delete[] pixelrow;
}


// non-struct function

// functions

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

inline int32_t scrcoord_as_lowerleft(const NTYP a) {
	// calculating the screen coordinte of the pixel that contains the coordinate
	// if the coordinate lies on an edge/corner (and belongs to more than one pixel)
	// the pixel where it lies on the left,bottom edge/corner is returned
	#ifdef _QUADMATH
	return (int)floorq( (a - COMPLETE0) * scalePixelPerRange );
	#else
	return (int)floor( (a - COMPLETE0) * scalePixelPerRange );
	#endif
}

// maximum of 4 doubles. often used by the boundingbox function
inline NTYP maximumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}

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

// IMPORTANT FUNCTION
// - function that performs the iteration of the quadratic case: z := z*z + seedC
// - function can be replaced by any other boundiong box generating routine to accommodate
// for other types of iteration formulas like z^3+c, z^6+c etc.
// - CAVE: Care must be taken that the C++ double floating point type used here can handle
// the produced numbers to absolute accuracy
// - expression is semi-automatically generated and left un-optimized for ease of substitution

// iterting function z^2+c
void getBoundingBoxfA_z2c(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+seedC0re;
	fA.x1=maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+seedC1re;
	fA.y0=2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedC0im;
	fA.y1=2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedC1im;
}

// z^2+A*z+c
void getBoundingBoxfA_z2azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-maximumD(A.y0*A.y0,A.y1*A.y1);
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-minimumD(A.y0*A.y0,A.y1*A.y1);
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
}

// z^3+A*z+c
void getBoundingBoxfA_z3azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x0*A.x0*A.x0-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x1*A.x1*A.x1-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y1*A.y1*A.y1)+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y0*A.y0*A.y0)+seedC1im;
}

// z^4+A*z+c
void getBoundingBoxfA_z4azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+4*minimumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*maximumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+4*maximumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*minimumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedC1im;
}

// z^5+A*z+c
void getBoundingBoxfA_z5azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x0*A.x0*A.x0*A.x0*A.x0-(2*(5*maximumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*minimumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x1*A.x1*A.x1*A.x1*A.x1-(2*(5*minimumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*maximumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y0*A.y0*A.y0*A.y0*A.y0+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y1*A.y1*A.y1*A.y1*A.y1+seedC1im;
}

// z^6+A*z+c
void getBoundingBoxfA_z6azc(PlaneRect& A,PlaneRect& fA) {
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
	
	const NTYP  DD=REVCGBLOCKWIDTH*scaleRangePerPixel;
	int32_t parentx,parenty;
	
	for(int dl=1;dl<=2;dl++) {
		A.y1=COMPLETE0;
		// first pass: calculuate how many are nedded per square
		// 2nd pass: build cell graph as array
		if (dl==1) printf("counting parents ...\n");
		else printf("setting parents to squares ...\n");

		for(int32_t y=0;y<SCREENWIDTH;y+=REVCGBLOCKWIDTH) {
			parenty=(y >> REVCGBITS);
			A.y0=A.y1;
			A.y1 += DD;
		
			A.x1=COMPLETE0;
			for(int32_t x=0;x<SCREENWIDTH;x+=REVCGBLOCKWIDTH) {
				parentx=(x >> REVCGBITS);
				A.x0=A.x1;
				A.x1 += DD;
			
				getBoundingBoxfA(A,bbxfA);
			
				if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA) > 0) {
					continue;
				}

				ScreenRect scr;
				scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
				scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
				scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
				scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
			
				#define RANGE(WW) \
				if (WW < 0) WW=0; \
				else if (WW >= SCREENWIDTH) WW=SCREENWIDTH-1;
			
				RANGE(scr.x0)
				RANGE(scr.x1)
				RANGE(scr.y0)
				RANGE(scr.y1)
			
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
	
	// if raw data file exists: read it and if necessary blow up the pixels 2fold
	if (data5->readRawBlowUp("_in") <= 0) {
		// data5 object is - no matter what data it holds - considered uninitialised
		printf("searching for special exterior ... ");

		// at the start: everything is considered gray
		encgrayx0=encgrayy0=0;
		encgrayx1=encgrayy1=SCREENWIDTH-16;
		planegrayx0=planegrayy0=COMPLETE0;
		planegrayx1=planegrayy1=COMPLETE1;

		// squares(pixels) whose bounding box lies completely in the special exterior
		find_special_exterior_hitting_squares();
	}
	
	// mark every tile as to be visited
	for(int32_t i=0;i<REVCGmaxnumberQ;i++) {
		data5->revcgYX[i].tovisit=1;
	}
	
	// propagating white and potentially white
	propagate_white();
	
	printf("\ncoloring interior cells ...\n");
	coloring_interior_black();
}

void find_special_exterior_hitting_squares(void) {
	PlaneRect bbxfA,A16,A;
	
	encgrayx0=encgrayy0=SCREENWIDTH-16;
	encgrayx1=encgrayy1=0;
	
	// initialising gray enclosement per row
		
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		data5->gray[y].g0=SCREENWIDTH-16;
		data5->gray[y].g1=0;
	}
	
	// some step size to output progress
	int32_t still_todo0=(SCREENWIDTH >> 4) >> 2;
	int32_t still_todo=1;
	
	const NTYP  DD16SCALE=16.0*scaleRangePerPixel;
	A16.y1=COMPLETE0; 
	for(int32_t y16=0;y16<SCREENWIDTH;y16+=16) {
		A16.y0=A16.y1;
		A16.y1 += DD16SCALE;
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
			A16.x1 += DD16SCALE;
					
			getBoundingBoxfA(A16,bbxfA);

			// does the 16x16 sqzare lie completely in the special exterior
			if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA)>0) {
				// yes => all can be colored white
				int32_t memx0=(x16 >> 4);			
				for(int32_t y2=y16;y2<(y16+16);y2++) {
					data5->pixelrow[y2][memx0]=SQUARE_WHITE_16_CONSECUTIVE;
				} // y2

				continue;
			}
					
			// now the individual points have to be checked
			const int32_t yendebit=y16+16;
			A.y0=y16*scaleRangePerPixel + COMPLETE0;
			A.y1=A.y0 + scaleRangePerPixel;
			for(int32_t ybit=y16;ybit<yendebit;++ybit) {
				int32_t mem=x16 >> 4;
				uint32_t w=0;
				uint8_t new_gray_set=0;
				A.x0=x16*scaleRangePerPixel + COMPLETE0;
				A.x1=A.x0 + scaleRangePerPixel;
				for(int32_t bitpos=0;bitpos<32;bitpos+=2) {
					getBoundingBoxfA(A,bbxfA);
					
					if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA)>0) {
						w |= (SQUARE_WHITE << bitpos);
					} else {
						new_gray_set=1;
						if (SQUARE_LIES_ENTIRELY_IN_COMPLETE(bbxfA) <= 0) {
							w |= (SQUARE_GRAY_POTENTIALLY_WHITE << bitpos);
						} else {
							w |= (SQUARE_GRAY << bitpos);
						}
					} // gray
					
					A.x0=A.x1;
					A.x1 += scaleRangePerPixel;
				} // bitpos
				
				data5->pixelrow[ybit][mem]=w;
				
				if (new_gray_set>0) {
					if (ybit < encgrayy0) encgrayy0=ybit;
					if (ybit > encgrayy1) encgrayy1=ybit;
					if (x16 < encgrayx0) encgrayx0=x16;
					if (x16 > encgrayx1) encgrayx1=x16;
					if (x16 < gray0) gray0=x16;
					if (x16 > gray1) gray1=x16;
				}
				
				A.y0=A.y1;
				A.y1 += scaleRangePerPixel;
			} // ybit
		} // x16
		
		// adjusting gray enclosement per row
		for(int32_t yy=y16;yy<y16+16;yy++) {
			data5->gray[yy].g0=gray0;
			data5->gray[yy].g1=gray1;
		} // yy
		
	} // y16
	
	// adjusting image gray enclosement
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

void propagate_white(void) {
	PlaneRect A,bbxfA;
	ScreenRect scr;
	int32_t still_todo0=((SCREENWIDTH >> REVCGBITS) >> 1);
	int32_t still_todo=1;
	int32_t loops_till_saving_raw_data=10;
	
	int32_t changed=1;
	
	while (changed>0) {
		changed=0;
		if ( (--loops_till_saving_raw_data) <= 0) {
			#ifdef _QUADMATH
			loops_till_saving_raw_data=6;
			#else
			loops_till_saving_raw_data=16;
			#endif
			printf("saving raw data as temporary ... ");
			data5->saveRaw("_temp");
		}
		
		printf("\npropagating (potentially) white ... ");
		for(int32_t y256=0,YBLOCK=0;y256<SCREENWIDTH;y256+=REVCGBLOCKWIDTH,YBLOCK++) {
			int32_t yrevoffset=YBLOCK*REVCGmaxnumber;
			if ( (--still_todo) <= 0) {
				printf("%i ",SCREENWIDTH-y256);
				still_todo=still_todo0;
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
				for(int32_t y=y256;y<Y256ENDE;y++) {
					A.y0=y*scaleRangePerPixel + COMPLETE0;
					A.y1=A.y0 + scaleRangePerPixel;

					const int32_t xanf=data5->gray[y].g0;
					const int32_t xende=data5->gray[y].g1;
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
						uint32_t w=data5->pixelrow[y][wmem];
				
						// no gray square in this 32-bit integer
						if (
							(w == SQUARE_WHITE_16_CONSECUTIVE) ||
							(w == SQUARE_BLACK_16_CONSECUTIVE) 
						) continue; 
				
						uint32_t wneu=w;
						int32_t w_changed=0;
						for(int32_t wbith=0;wbith<16;++wbith) {
							uint32_t f=w & 0b11;
							w >>= 2;
					
							// current pixel is already gray
							if (
								(f != SQUARE_GRAY) &&
								(f != SQUARE_GRAY_POTENTIALLY_WHITE) 
							) continue;
					
							A.x0=(x+wbith)*scaleRangePerPixel + COMPLETE0;
							A.x1=A.x0+scaleRangePerPixel;
				
							getBoundingBoxfA(A,bbxfA);
				
							// bounding box in special exterior
							if ((SQUARE_LIES_ENTIRELY_OUTSIDE_GRAY_ENCLOSEMENT(bbxfA))>0) {
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
								w_changed=1;
								continue;
							}
					
							int32_t potw=(f == SQUARE_GRAY_POTENTIALLY_WHITE);
							int32_t hits_white=0;
							int32_t hits_nonwhite=0;
					
							if (SQUARE_LIES_ENTIRELY_IN_GRAY_ENCLOSEMENT(bbxfA) <= 0) {
								// overlaps with white region
								potw=1;
								hits_white=1;
							}
					
							scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
							scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
							scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
							scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
					
							#define RANGE2(WW) \
							if (WW < 0) { WW=0; hits_white=1; }\
							else if (WW >= SCREENWIDTH) { WW=SCREENWIDTH-1; hits_white=1; }
					
							RANGE2(scr.x0)
							RANGE2(scr.x1)
							RANGE2(scr.y0)
							RANGE2(scr.y1)

							int32_t goon=1;
							uint32_t wbbx;

							for(int32_t ty=scr.y0;((ty<=scr.y1)&&(goon>0));ty++) {
								int32_t bbxbith=scr.x0 % (1 << 4);
								int32_t bbxmem=scr.x0 >> 4;
						
								int32_t still_todobits=scr.x1-scr.x0+1;
								wbbx=data5->pixelrow[ty][bbxmem];
								// sometimes the whole 32 bit are white
								if (wbbx == SQUARE_WHITE_16_CONSECUTIVE) {
									hits_white=1;
									still_todobits -= (16-bbxbith);
									bbxbith=16; 
								} else if (wbbx == SQUARE_BLACK_16_CONSECUTIVE) {
									hits_nonwhite=1;
									still_todobits -= (16-bbxbith);
									bbxbith=16; 
								} else if (wbbx == SQUARE_GRAY_16_CONSECUTIVE) {
									hits_nonwhite=1;
									still_todobits -= (16-bbxbith);
									bbxbith=16; 
								} else if (wbbx == SQUARE_GRAYPOTW_16_CONSECUTIVE) {
									hits_nonwhite=1;
									potw=1;
									still_todobits -= (16-bbxbith);
									bbxbith=16; 
								} else {
									wbbx >>= (bbxbith << 1);
								}
								
								while (still_todobits>0) {
									if (bbxbith >= 16) {
										bbxmem++;
										bbxbith=0;
										wbbx=data5->pixelrow[ty][bbxmem];
										if (wbbx == SQUARE_WHITE_16_CONSECUTIVE) {
											hits_white=1;
											still_todobits -= (16-bbxbith);
											bbxbith=16; 
											continue;
										} else if (wbbx == SQUARE_BLACK_16_CONSECUTIVE) {
											hits_nonwhite=1;
											still_todobits -= (16-bbxbith);
											bbxbith=16; 
											continue;
										} else if (wbbx == SQUARE_GRAY_16_CONSECUTIVE) {
											hits_nonwhite=1;
											still_todobits -= (16-bbxbith);
											bbxbith=16; 
											continue;
										} else if (wbbx == SQUARE_GRAYPOTW_16_CONSECUTIVE) {
											hits_nonwhite=1;
											potw=1;
											still_todobits -= (16-bbxbith);
											bbxbith=16; 
											continue;
										}
									} 
							
									int32_t f=wbbx & 0b11;
									if (f == SQUARE_WHITE) hits_white=1;
									else if (f == SQUARE_GRAY_POTENTIALLY_WHITE) { potw=1; hits_nonwhite=1; }
									else hits_nonwhite=1;
							
									if (
										(potw > 0) && 
										(hits_white > 0) && 
										(hits_nonwhite > 0) 
									) {
										// no further information to be gained
										// checking can be surpassed
										goon=0;
										break;
									}
									still_todobits--;
									bbxbith++;
									wbbx >>= 2;
								} // while
							} // ty
					
							if ( (hits_white>0) && (hits_nonwhite<=0) ) {
								// only white pixels in the bounding box
								wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_WHITE[wbith]);
								w_changed=1;
							} else if ( (hits_white>0) || (potw>0) ) {
								if (f != SQUARE_GRAY_POTENTIALLY_WHITE) {
									// if not already marked as potentially
									// white => new information generated
									wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[wbith],ARRAY_SQUARE_GRAYPOTW[wbith]);
									w_changed=1;
								}
							}
						} // wbith
				
						if (w_changed>0) {
							data5->pixelrow[y][wmem]=wneu;
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

void coloring_interior_black(void) {
	// every gray but not potentially gray pixel
	// can be colored black as it never hits any
	// pixel that is either white or leads to white
	
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		const int32_t xanf=data5->gray[y].g0;
		const int32_t xende=data5->gray[y].g1;
		int32_t mem=-1 + (xanf >> 4);

		for(int32_t x=xanf;x<=xende;x+=16) {
			mem++; 
			uint32_t w=data5->pixelrow[y][mem];
			uint32_t wneu=w;
			uint8_t w_changed=0;
			
			for(int32_t b=0;b<16;b++) {
				if ( (w & 0b11) == SQUARE_GRAY) {
					w_changed=1;
					wneu=SET_SINGLE_PIXELCOLOR_INTO_4BYTEINTEGER(wneu,COLOR_CLEARMASK[b],ARRAY_SQUARE_BLACK[b]);
				}
				w >>= 2;
			} // b
			
			if (w_changed>0) data5->pixelrow[y][mem]=wneu;
		} // x
		
	} // y
}

ParentManager::ParentManager() {
	lastallocated=NULL;
	memused=0;
	freefrom=0;
}

Parent* ParentManager::getParentSpace(const int32_t aneeded) {
	if (
		(!lastallocated) ||
		( (memused-freefrom) < aneeded)
	) {
		// allocate new memory in 1 GB chunks
		// deleting is done in a dirty manner when
		// program ends
		int64_t all=( (int64_t)1 << 30) / sizeof(Parent);
		lastallocated=new Parent[all];
		if (!lastallocated) {
			printf("Memory error. Parentmanager.\n");
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
	sprintf(erg,"c_ia_%lg_%lg_x_%lg_%lg",
		(double)seedC0re,(double)seedC1re,(double)seedC0im,(double)seedC1im);
	return erg;
}

char* FAKTORAstr(char* erg) {
	sprintf(erg,"A_%lg_%lg",
		(double)FAKTORAre,(double)FAKTORAim);
	return erg;
}

void test_for_sufficient_precision_bits(void) {
	// using the current datatype, the current screen width
	// the chosen function and the range in the complex plane
	// Are there enough bits precision in NTYP for the
	// mantissa to support the boundingbox calculation
	PlaneRect A,bbxfA;

	double d=0.0;
	
	#ifdef _QUADMATH
	NTYP dtmp=0.0;
	#define MAXVALUE(AA,BB,CC,DD) \
	{\
		A.x0=AA;\
		A.x1=BB;\
		A.y0=CC;\
		A.y1=DD;\
		getBoundingBoxfA(A,bbxfA);\
		NTYP d2=ceilq(maximumD(\
			fabsq(bbxfA.x0),\
			fabsq(bbxfA.x1),\
			fabsq(bbxfA.y0),\
			fabsq(bbxfA.y1)));\
		if (d2 > dtmp) dtmp=d2;\
	}
	#else
	NTYP dtmp=0.0;
	#define MAXVALUE(AA,BB,CC,DD) \
	{\
		A.x0=AA;\
		A.x1=BB;\
		A.y0=CC;\
		A.y1=DD;\
		getBoundingBoxfA(A,bbxfA);\
		NTYP d2=ceil(maximumD(\
			fabs(bbxfA.x0),\
			fabs(bbxfA.x1),\
			fabs(bbxfA.y0),\
			fabs(bbxfA.y1)));\
		if (d2 > dtmp) dtmp=d2;\
	}
	#endif
	
	// Quadrant-wise 1
	// since bounding boxes are calculated
	// via interval arithmetics and exponentiation,
	// zero is not allowed to be within an interval
	// just at its borders
	MAXVALUE(0,RANGE1,0,RANGE1); // upper right quadrant
	MAXVALUE(0,RANGE1,RANGE0,0); 
	MAXVALUE(RANGE0,0,0,RANGE1); 
	MAXVALUE(RANGE0,0,RANGE0,0); 
	d=(double)dtmp;
	
	// how many bits are necessary for the max (absolute) value
	int bitsd=(int)ceil(log(d)/log(2.0));

	// how many bits are necessary for the fractional part
	int bitsfract;
	#ifdef _QUADMATH
	bitsfract=(int)ceilq(
		fabsq(
			logq( (__float128)(RANGE1-RANGE0) / (__float128)SCREENWIDTH )
			/ logq(2.0)
		) );
	#else
	bitsfract=(int)ceil(
		fabs(
			log( (NTYP)(RANGE1-RANGE0) / (NTYP)SCREENWIDTH )
			/ log(2.0)
		) );
	#endif
	
	int bitsused=bitsd + bitpower*bitsfract;
	printf("mantissa precision is ");
	// -1 as buffer
	if (bitsused >= (_BITSNTYP-1)) printf("not ");
	printf("sufficient (%i needed, %i provided. Decimal part's maximum value=%.0lf\n",
		bitsused,_BITSNTYP,d);
}

// allocates chunks of memory
void arrayManagerDbyte(
	PDBYTE* memptr,
	const int32_t aylen,
	const int32_t axlen,
	const DBYTE avalue
) {
	int64_t maxmem=( (int64_t)1 << 30); // 1-GB chunks
	
	int32_t rows_at_once=(int)floor( (double)maxmem / (sizeof(DBYTE)*axlen) );
	int32_t todoy=aylen;
	
	int32_t pos=0;
	while (todoy>0) {
		int32_t alloky;
		if (todoy >= rows_at_once) alloky=rows_at_once;
		else alloky=todoy;
		DBYTE *p=new DBYTE[alloky*axlen];
		for(int32_t i=0;i<alloky;i++) {
			if (pos < aylen) {
				memptr[pos]=&p[i*axlen];
				for(int k=0;k<axlen;k++) {
					memptr[pos][k]=avalue;
				}
				pos++;
			}
		}
		
		todoy -= alloky;
	} // while
}

#define SETSCCREL(SCCDATA,XX,YY,COL) \
SCCDATA.scc[(YY-SCCDATA.grayency0)] [(XX-SCCDATA.grayencx0)]=COL;\
	
#define GETSCCREL(SCCDATA,XX,YY,ERG) \
ERG=SCCDATA.scc[YY-SCCDATA.grayency0] [XX-SCCDATA.grayencx0];

#define GETFREECOLOR(SCCDATA,ERG) \
{\
	ERG=0;\
	for(int32_t i=COLORSTART;i<DBYTEMAX;i++) {\
		if ( (SCCDATA.scccolor[i].flag & FIXED2FLAG_INUSE) == 0) {\
			ERG=i;\
			SCCDATA.scccolor[i].flag |= FIXED2FLAG_INUSE;\
			SCCDATA.scccolor[i].neu=i;\
			break;\
		}\
	}\
}

int32_t isInterior(SCCData& ain,const int32_t ax,const int32_t ay,const int32_t bx,const int32_t by) {
	for(int32_t dy=ay;dy<=by;dy++) {
		for(int32_t dx=ax;dx<=bx;dx++) {
			DBYTE f;
			GETSCCREL(ain,dx,dy,f);
			if (f != SQUARE_BLACK) return 0;
		}
	} // dy
	return 1;
}

void generateBasins(SCCData& ain) {
	int activecolor=0;
	
	PlaneRect A,bbxfA;
	ScreenRect scr;
	
	char* settoactivecolor=new char[DBYTEMAX];
	for(int32_t i=0;i<DBYTEMAX;i++) settoactivecolor[i]=0;
	const int32_t MAXNUMBERTOVISIT=16384;
	Gray_in_row *tovisit=new Gray_in_row[MAXNUMBERTOVISIT];
	int32_t counttovisit=0;
	
	#define MERGECOLORS(SCCDATA) \
	{\
		printf("merging colors ...\n");\
		for(int32_t i=COLORSTART;i<DBYTEMAX;i++) {\
			if ((SCCDATA.scccolor[i].flag & FIXED2FLAG_INUSE) == 0) continue;\
			if (SCCDATA.scccolor[i].neu == SCCDATA.scccolor[i].alt) continue;\
			int32_t finalf=i;\
			while (SCCDATA.scccolor[finalf].alt != SCCDATA.scccolor[finalf].neu) {\
				finalf=SCCDATA.scccolor[finalf].neu;\
			}\
			SCCDATA.scccolor[i].neu=finalf;\
		}\
		int32_t a0=DBYTEMAX,a1=0;\
		for(int32_t dy=SCCDATA.grayency0;dy<=SCCDATA.grayency1;dy++) {\
			for(int dx=SCCDATA.grayencx0;dx<=SCCDATA.grayencx1;dx++) {\
				DBYTE f;\
				GETSCCREL(SCCDATA,dx,dy,f)\
				if (f < COLORSTART) continue;\
				if (SCCDATA.scccolor[f].neu == SCCDATA.scccolor[f].alt) continue;\
				SETSCCREL(SCCDATA,dx,dy,SCCDATA.scccolor[f].neu)\
				SCCDATA.scccolor[f].flag |= FIXED2FLAG_TODELETE;\
				if (f < a0) a0=f;\
				if (f > a1) a1=f;\
			}\
		}\
		for(int32_t i=a0;i<=a1;i++) {\
			if ( (SCCDATA.scccolor[i].flag & FIXED2FLAG_TODELETE) == 0) continue;\
			SCCDATA.scccolor[i].flag=0;\
			SCCDATA.scccolor[i].neu=i;\
		}\
	}
	
	int todo0=SCREENWIDTH >> 6;
	int todo=1;
	for(int32_t y=ain.grayency0;y<=ain.grayency1;y++) {
		if ((--todo)<=0) {
			printf("%i ",y);
			todo=todo0;
		}
		for(int32_t x=ain.grayencx0;x<=ain.grayencx1;x++) {
			DBYTE f;
			GETSCCREL(ain,x,y,f);
			if (f != SQUARE_BLACK) continue;
			
			// neuen Punkt gefunden
			// neue Farbe suchen
			activecolor=0;
			GETFREECOLOR(ain,activecolor);
			if (activecolor < COLORSTART) {
				// keine mehr frei
				// Faerben
				MERGECOLORS(ain)
				activecolor=0;
				GETFREECOLOR(ain,activecolor)
				if (activecolor < COLORSTART) {
					printf("Not implemented. Color table full - even after merging.\n");
					exit(99);
				}
			}
			
			// activecolor contains index to an unused color
			
			counttovisit=0;
			tovisit[counttovisit].g0=x;
			tovisit[counttovisit].g1=y;
			counttovisit=1;
					
			int32_t ers0=DBYTEMAX;
			int32_t ers1=0;
			
			while (counttovisit>0) {
				int xx=tovisit[counttovisit-1].g0;
				int yy=tovisit[counttovisit-1].g1;
				counttovisit--;
				SETSCCREL(ain,xx,yy,activecolor);
				
				A.y0=yy*scaleRangePerPixel + COMPLETE0;
				A.y1=A.y0+scaleRangePerPixel;
				A.x0=xx*scaleRangePerPixel + COMPLETE0;
				A.x1=A.x0+scaleRangePerPixel;
						
				getBoundingBoxfA(A,bbxfA);
						
				scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
				scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
				scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
				scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
						
				if (
					(scr.x0 < ain.grayencx0) ||
					(scr.x1 > ain.grayencx1) ||
					(scr.y0 < ain.grayency0) ||
					(scr.y1 > ain.grayency1) 
				) {
					printf("Implementation error Interior cells not pointing to other interiro cells.\n");
					exit(99);
				}

				#define TRIMSCRX(WW) \
				if (WW < ain.grayencx0) WW=ain.grayencx0;\
				if (WW > ain.grayencx1) WW=ain.grayencx1;
						
				#define TRIMSCRY(WW) \
				if (WW < ain.grayency0) WW=ain.grayency0;\
				if (WW > ain.grayency1) WW=ain.grayency1;

				TRIMSCRX(scr.x0)
				TRIMSCRX(scr.x1)
				TRIMSCRY(scr.y0)
				TRIMSCRY(scr.y1)
				
				int area=(scr.x1-scr.x0+1)*(scr.y1-scr.y0+1);
				if ((counttovisit+area) > (MAXNUMBERTOVISIT-128)) continue; 
						
				for(int32_t by=scr.y0;by<=scr.y1;by++) {
					for(int32_t bx=scr.x0;bx<=scr.x1;bx++) {
						DBYTE fb;
						GETSCCREL(ain,bx,by,fb);
						if (fb == activecolor) continue;
								
						if (fb == SQUARE_BLACK) {
							// new reached pixel
							tovisit[counttovisit].g0=bx;
							tovisit[counttovisit].g1=by;
							counttovisit++;
						} else if (fb >= COLORSTART) {
							// trifft andere Farbe
							// fusionieren
							settoactivecolor[fb]=1;
							if (fb < ers0) ers0=fb;
							if (fb > ers1) ers1=fb;
						} else if (
							(fb != FIXED2_ACTIVETOVISIT)
						) {
							printf("Implementation error. Wrong pixel color hit.\n");
							exit(99);
						}
					} // bx
				} // by
			} // Liste
			
			// die Fusionen ent-transitivieren im Feld scccolor
			// aber noch nicht im Bild
			// alle bereits in setzeaktiv vorhanden
			// Farben: gehe ihren Pfad nach und nimm
			// die dortigen NEU auch auf
			int ch=1;
			while (ch>0) {
				ch=0;
				for(int i=ers0;i<=ers1;i++) {
					if (settoactivecolor[i]<=0) continue;
					// folow
					int idx=i;
					while (ain.scccolor[idx].neu != ain.scccolor[idx].alt) {
						if (settoactivecolor[idx]<=0) {
							settoactivecolor[idx]=1;
							ch=1;
						}
						idx=ain.scccolor[idx].neu;
						if (idx < ers0) ers0=idx;
						if (idx > ers1) ers1=idx;
					}
					if (settoactivecolor[idx]<=0) {
						settoactivecolor[idx]=1;
						ch=1;
					}
				} // i
			} // while
			
			settoactivecolor[activecolor]=0;
			for(int32_t i=ers0;i<=ers1;i++) {
				if (settoactivecolor[i]<=0) continue;
				settoactivecolor[i]=0;
				ain.scccolor[i].neu=activecolor;
			}
		} // x
	} // y
	
	// check if merging is possible
	MERGECOLORS(ain)
	
	delete[] settoactivecolor;
	delete[] tovisit;
}

// chunk precoloring - speed up
void preColoring(SCCData& ain,const int32_t STEP) {
	const int32_t STEPM1=STEP-1;
	int32_t WALK=16;
	if (WALK > STEP) WALK=STEP;
	
	printf("\nprecoloring ... ");
	PlaneRect A,bbxfA;
	ScreenRect scr;
	NTYP RECTSIZE=STEP*scaleRangePerPixel;
	
	int32_t done=0;
	
	const DBYTE HALVCOLOR=(1 << 15);

	#define SETCOLORRECT(SCCDATA,XX0,YY0,XX1,YY1,FF)\
	{\
		for(int32_t ty=(YY0);ty<=(YY1);ty++) {\
			for(int32_t tx=(XX0);tx<=(XX1);tx++) {\
				SETSCCREL(SCCDATA,tx,ty,FF);\
			}\
		}\
	}
	
	int32_t todo0=SCREENWIDTH >> 3;
	int32_t todo=1;
	for(int32_t y=ain.grayency0;y<=(ain.grayency1-WALK-STEP);y+=WALK) {
		if (((todo-=WALK))<=0) {
			printf("%i ",y);
			todo=todo0;
		}
		if (done > 0) break;

		A.y0=y*scaleRangePerPixel + COMPLETE0;
		A.y1=A.y0+RECTSIZE;
		
		for(int32_t x=ain.grayencx0;x<=(ain.grayencx1-STEP-WALK);x+=WALK) {
			if (isInterior(ain,x,y,x+STEPM1,y+STEPM1)
				<= 0) continue;
			
			A.x0=x*scaleRangePerPixel + COMPLETE0;
			A.x1=A.x0 + RECTSIZE;
			
			getBoundingBoxfA(A,bbxfA);
			scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
			#define CONTX(XX)\
				if (\
					((XX) < ain.grayencx0) ||\
					((XX) > ain.grayencx1)\
				) continue;
			
			#define CONTY(YY)\
				if (\
					((YY) < ain.grayency0) ||\
					((YY) > ain.grayency1)\
				) continue;

			// not pre-colorable
			CONTX(scr.x0)
			scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
			CONTX(scr.x1)
			scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
			CONTY(scr.y0)
			scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
			CONTY(scr.y1)
			
			if (
				isInterior(ain,scr.x0,scr.y0,scr.x1,scr.y1)
				<= 0
			) continue;
			
			DBYTE ff;
			GETFREECOLOR(ain,ff)

			SETCOLORRECT(ain,x,y,x+STEPM1,y+STEPM1,ff);
			SETCOLORRECT(ain,scr.x0,scr.y0,scr.x1,scr.y1,ff);

			if (ff >= HALVCOLOR) {
				done=1; // leave some colors for the actual routine
				break;
			}

		} // x
	} // y
}

void generateBasinImage(SCCData& ain,const char* afn) {
	printf("creating basin image ...\n");
	RGB4 pal[256];
	char tmp[1024];
	for(int32_t i=0;i<256;++i) pal[i].R=pal[i].G=pal[i].B=pal[i].alpha=0;

	pal[SQUARE_GRAY].R=127;
	pal[SQUARE_GRAY].G=127;
	pal[SQUARE_GRAY].B=127;
	pal[SQUARE_GRAY_POTENTIALLY_WHITE].R=127;
	pal[SQUARE_GRAY_POTENTIALLY_WHITE].G=127;
	pal[SQUARE_GRAY_POTENTIALLY_WHITE].B=127;
	pal[SQUARE_BLACK].R=0;
	pal[SQUARE_BLACK].G=0;
	pal[SQUARE_BLACK].B=0;
	pal[SQUARE_WHITE].R=255;
	pal[SQUARE_WHITE].G=255;
	pal[SQUARE_WHITE].B=255;
	pal[COLORRED].R=255;
	pal[COLORRED].G=0;
	pal[COLORRED].B=0;
	pal[FIXED2_ACTIVETOVISIT].R=255;
	pal[FIXED2_ACTIVETOVISIT].R=255;
	pal[FIXED2_ACTIVETOVISIT].R=0;

	int32_t anzbasin=0;
	for(int32_t i=COLORSTART;i<DBYTEMAX;i++) {
		if ( (ain.scccolor[i].flag & FIXED2FLAG_INUSE) != 0) {
			anzbasin++;
		}
	}
	printf("\n\n%i basins of attraction found\n",anzbasin);
	
	char pref[1024];
	pref[0]=0;
	if (anzbasin > 200) {
		printf("NICHT möglich. Zuviele Basins gefunden.\n");
		printf("die meisten SIND GLEICHFARBIG.\n");
		anzbasin=200;
		strcpy(pref,"_FEHLER_ZUVIELE_BASINS_");
	}
	
	int32_t idx=0;
	int32_t maxidx=255-COLORSTART;
	for(int32_t i=COLORSTART;i<DBYTEMAX;i++) {
		if ((ain.scccolor[i].flag & FIXED2FLAG_INUSE) != 0) {
			double d=idx; d /= (anzbasin+1);
			int r,g,b;
			basinpal.getColor(d,r,g,b);
			pal[idx+COLORSTART].R=r;
			pal[idx+COLORSTART].G=g;
			pal[idx+COLORSTART].B=b;
			ain.scccolor[i].palidx=idx+COLORSTART;
			idx++;
			if (idx > maxidx) {
				printf("Too many (>240) basins. Some are not colored correctly.\n");
				idx=maxidx;
			}
		}
	}
	
	int ybytes=SCREENWIDTH; // 1 Pixel pro Byte
	uint8_t* rgbz=new uint8_t[ybytes];
	uint32_t off
	=		14 // FILEHeader
		+	40 // Bitmapheader
		+	256*4 // ColorPalette
	;
	uint32_t filelen
		=	off
		+	(ybytes*SCREENWIDTH);
	;
	
	sprintf(tmp,"%s_basins.bmp",afn);
	FILE *fbmp=fopen(tmp,"wb");
	if (!fbmp) {
		printf("Error. Bitmap-file not writeable.\n");
		exit(99);
	}
	write2(fbmp,66,77);
	fwrite(&filelen,1,sizeof(filelen),fbmp);
	write4(fbmp,0,0,0,0); 
	fwrite(&off,1,sizeof(off),fbmp);
	write4(fbmp,40,0,0,0);
	
	uint32_t a=SCREENWIDTH;
	fwrite(&a,sizeof(a),1,fbmp);
	fwrite(&a,sizeof(a),1,fbmp);
	write2(fbmp,1,0);
	write2(fbmp,8,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,19,11,0,0);
	write4(fbmp,19,11,0,0);
	write4(fbmp,0,1,0,0);
	write4(fbmp,0,0,0,0);

	// Palette
	for(int i=0;i<256;++i) write4(fbmp,pal[i].B,pal[i].G,pal[i].R,pal[i].alpha);
	
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		for(int32_t x=0;x<SCREENWIDTH;x++) {
			DBYTE f;
			if (
				(y < ain.grayency0) || 
				(y > ain.grayency1) ||
				(x < ain.grayencx0) || 
				(x > ain.grayencx1)
			) {
				f=SQUARE_WHITE;
			} else {
				GETSCCREL(ain,x,y,f)
				if (f >= COLORSTART) {
					if (ain.scccolor[f].palidx < COLORSTART) {
						printf("Implementation error. palidx too low.\n");
						exit(99);
					}
					f=ain.scccolor[f].palidx;
				} 
			}
			
			rgbz[x]=f;
		} // x
		
		fwrite(rgbz,SCREENWIDTH,sizeof(uint8_t),fbmp);
	} // y
	
	fclose(fbmp);
	
	delete[] rgbz;
}

#define GETFREIEFARBE(SCCDATA,ERG) \
{\
	ERG=0;\
	for(int i=COLORSTART;i<DBYTEMAX;i++) {\
		if ( (SCCDATA.scccolor[i].flag & FIXED2FLAG_INBENUTZUNG) == 0) {\
			ERG=i;\
			SCCDATA.scccolor[i].flag |= FIXED2FLAG_INBENUTZUNG;\
			SCCDATA.scccolor[i].neu=i;\
			break;\
		}\
	}\
}

// if interior cells are present, tries to color
// them according to basins of attraction
void basin(const char* fnbase) {
	printf("searching for basins of attraction ...\n");
	// set up data structure
	SCCData sccd;
	
	SCCColor *pt=new SCCColor[DBYTEMAX];
	for(int i=0;i<DBYTEMAX;i++) {
		pt[i].alt=pt[i].neu=i;
		pt[i].flag=0;
		pt[i].palidx=COLORRED;
	}
	sccd.scccolor=pt;

	int32_t gm0=encgrayx0 >> 4;
	sccd.grayencx0=(gm0 << 4);
	int gm1=encgrayx1 >> 4;
	if (gm1 >= MEMWIDTH) gm1=MEMWIDTH-1;
	sccd.grayencx1=( (gm1+1) << 4) - 1;
	if (sccd.grayencx1 >= SCREENWIDTH) sccd.grayencx1=SCREENWIDTH-1;
	sccd.grayency0=encgrayy0;
	sccd.grayency1=encgrayy1;
	
	sccd.grayxlen=sccd.grayencx1-sccd.grayencx0 + 1;
	sccd.grayylen=sccd.grayency1-sccd.grayency0 + 1;
		
	PDBYTE *sccp=new PDBYTE[sccd.grayylen];
	arrayManagerDbyte(sccp,sccd.grayylen,sccd.grayxlen,SQUARE_WHITE);
	sccd.scc=sccp;
	printf(".");
	
	// copy bits of data5 into DBYTE words here
	// everything else has already be filled with SQUARE_WHITE
	for(int32_t y=encgrayy0;y<=encgrayy1;y++) {
		for(int mem=gm0;mem<=gm1;mem++) {
			uint32_t w=data5->pixelrow[y][mem];
			
			int32_t xx=mem*(1 << 4);
			for(int32_t bph=0;bph<16;bph++) {
				uint8_t f=w & 0b11;
				// identically treated here
				if (f == SQUARE_GRAY_POTENTIALLY_WHITE) f=SQUARE_GRAY;
				w >>= 2;
				SETSCCREL(sccd,xx,y,f);
				xx++;
			}
			
		} // x
	} // y
	printf(".");
	
	// preColoring
	// check for whole blocks of 128x128 black pixels
	// where their bounding box falls => coloring
	preColoring(sccd,128);
	
	// searching for basins
	generateBasins(sccd);
	
	// and saving
	generateBasinImage(sccd,fnbase);
	
	delete[] pt;
	delete[] sccp;
	// dirty exit: memory reusage done by program exit
}

// if interior cells are present and the set posseses
// attracting cycle(s) the immediate basins of periodic
// points are colored
void periodicity(const char* fnbase) {
	Blob *blobs=new Blob[MAXBLOBS];
	int anzblobs=0;
	
	Charmap md;
	md.setlenxy(SCREENWIDTH,SCREENWIDTH);
	md.clearPalette();
	md.setPaletteRGB(SQUARE_BLACK,0,0,0);
	md.setPaletteRGB(SQUARE_WHITE,255,255,255);
	md.setPaletteRGB(SQUARE_GRAY,127,127,127);
	md.setPaletteRGB(SQUARE_GRAY_POTENTIALLY_WHITE,255,255,0);
	md.setPaletteRGB(COLORRED,255,0,0); 

	// some shuffled version of a heat-map
	double d=0.0;
	double dst=0.19;
	for(int i=BLOBOFFSET;i<=255;i++) {
		int r,g,b;
		basinpal.getColor(d,r,g,b);
		md.setPaletteRGB(i,r,g,b);
		d += dst;
		while (d >= 1.0) d -= 1.0;
	}
	
	DBYTE *db=new DBYTE[SCREENWIDTH*SCREENWIDTH];
	
	int32_t schwencx0=SCREENWIDTH,schwencx1=0;
	int32_t schwency0=SCREENWIDTH,schwency1=0;
	
	int64_t offset=-1;
	printf("generating cycles ... \n");
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		int32_t x=-1;
		for(int32_t m=0;m<MEMWIDTH;m++) {
			uint32_t w=data5->pixelrow[y][m];
			for(int32_t b=0;b<16;b++) {
				offset++;
				x++;
				int f=w & 0b11;
				w >>= 2;
				if (
					(f==SQUARE_GRAY) ||
					(f==SQUARE_GRAY_POTENTIALLY_WHITE)
				) db[offset]=SQUARE_GRAY;
				else
				if (f==SQUARE_BLACK) {
					db[offset]=SQUARE_BLACK;
					if (x < schwencx0) schwencx0=x;
					if (x > schwencx1) schwencx1=x;
					if (y < schwency0) schwency0=y;
					if (y > schwency1) schwency1=y;
				} else db[offset]=f;
			}
		} // x
	} // y
	
	// nun über das Bild gehen und Blobs aufbauen
	anzblobs=0;
	DBYTE blobaktiv=BLOBOFFSET-1;
	DBYTE blobfertig=blobaktiv-1;
	printf("looking for interior blobs ...\n");
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		int64_t offy=y*SCREENWIDTH;
		for(int32_t x=0;x<SCREENWIDTH;x++) {
			if (db[offy+x] != SQUARE_BLACK) continue;
			
			// new blob
			if (anzblobs > (MAXBLOBS-8)) {
				printf("Error. Too many blobs.\n");
				exit(99);
			}
			int bx0=x,bx1=x,by0=y,by1=y;
			int setcol=BLOBOFFSET+anzblobs;
			int setidx=anzblobs;
			blobs[setidx].targetblob=-1;
			blobs[setidx].color=setcol;
			blobs[setidx].visited=0;
			blobs[setidx].zyklusnr=-1;
			anzblobs++;
			
			db[offy+x]=blobaktiv;
			
			// Floodfill
			int changed=1;
			while (changed>0) {
				changed=0;
				
				int32_t ey0=by0,ey1=by1,ex0=bx0,ex1=bx1;
				for(int32_t by=ey0;by<=ey1;by++) {
					int64_t offby=by*SCREENWIDTH;
					for(int32_t bx=ex0;bx<=ex1;bx++) {
						if (db[offby+bx] != blobaktiv) continue;
						db[offby+bx]=setcol;
						changed=1;
						if (bx < bx0) bx0=bx;
						if (bx > bx1) bx1=bx;
						if (by < by0) by0=by;
						if (by > by1) by1=by;
						
						for(int32_t xx=(bx-1);xx>=0;xx--) {
							if (db[by*SCREENWIDTH+xx] == SQUARE_BLACK) {
								db[by*SCREENWIDTH+xx]=blobaktiv;
								if (xx < bx0) bx0=xx;
								if (xx > bx1) bx1=xx;
							} else break;
						}
						for(int32_t xx=(bx+1);xx<SCREENWIDTH;xx++) {
							if (db[by*SCREENWIDTH+xx] == SQUARE_BLACK) {
								db[by*SCREENWIDTH+xx]=blobaktiv;
								if (xx < bx0) bx0=xx;
								if (xx > bx1) bx1=xx;
							} else break;
						}
						for(int32_t yy=(by+1);yy<SCREENWIDTH;yy++) {
							if (db[yy*SCREENWIDTH+bx] == SQUARE_BLACK) {
								db[yy*SCREENWIDTH+bx]=blobaktiv;
								if (yy < by0) by0=yy;
								if (yy > by1) by1=yy;
							} else break;
						}
						for(int32_t yy=(by-1);yy>=0;yy--) {
							if (db[yy*SCREENWIDTH+bx] == SQUARE_BLACK) {
								db[yy*SCREENWIDTH+bx]=blobaktiv;
								if (yy < by0) by0=yy;
								if (yy > by1) by1=yy;
							} else break;
						}
					} // bx
				} // by
			} // while
			
			blobs[setidx].rect.x0=bx0;
			blobs[setidx].rect.x1=bx1;
			blobs[setidx].rect.y0=by0;
			blobs[setidx].rect.y1=by1;
			blobs[setidx].mx=((bx0+bx1) >> 1);
			blobs[setidx].my=((by0+by1) >> 1);
			
		}
	} // y
	
	int64_t targetfree=0;
	
	// finding target blobs
	PlaneRect A,bbxfA;
	for(int32_t b=0;b<anzblobs;b++) {
		int32_t zf=-1;
		for(int32_t y=blobs[b].rect.y0;y<=blobs[b].rect.y1;y++) {
			A.y0=y*scaleRangePerPixel+COMPLETE0;
			A.y1=A.y0+scaleRangePerPixel;
			for(int32_t x=blobs[b].rect.x0;x<=blobs[b].rect.x1;x++) {
				if (db[y*SCREENWIDTH+x] != blobs[b].color) continue;
				
				A.x0=x*scaleRangePerPixel+COMPLETE0;
				A.x1=A.x0+scaleRangePerPixel;
				
				getBoundingBoxfA(A,bbxfA);
				
				if (SQUARE_LIES_ENTIRELY_IN_SPECEXT(bbxfA) > 0) {
					printf("Implementation error. Interior points to special exterior.\n");
					exit(99);
				}
				
				ScreenRect scr;
				scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
				scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
				scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
				scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
				
				for(int32_t bby=scr.y0;bby<=scr.y1;bby++) {
					for(int32_t bbx=scr.x0;bbx<=scr.x1;bbx++) {
						int32_t bf=db[bby*SCREENWIDTH+bbx];
						if (bf < BLOBOFFSET) {
							printf("Implementation error. Blob 2: bf=%i\n",bf);
							exit(99);
						}
						if (zf < 0) zf=bf;
						else if (zf != bf) { 
							zf=-2; 
							printf("zf-2\n");
							break; 
						}
					}
					if (zf == -2) break;
				} // bby
				
				if (zf == -2) break;
			} // x
			if (zf == -2) break;
		} // y
		
		if (zf>0) {
			blobs[b].targetblob=zf-BLOBOFFSET; 
		} else {
			blobs[b].targetblob=-1;
			targetfree++;
		}
	} // b
	
	// look for cycles in blob->targetblob->atrgetblob...
	int32_t zyklusanz=0;
	char *visited=new char[anzblobs+16];
	
	#define COLORZYKLUS(BNR,ZIELF) \
	{\
		if (((ZIELF)<BLOBOFFSET)||((ZIELF)>=256)) {\
			printf("Error. Too many cycle colors\n");\
			exit(99);\
		}\
		int32_t qf=blobs[BNR].color;\
		for(int32_t y=blobs[BNR].rect.y0;y<=blobs[BNR].rect.y1;y++) {\
			int64_t doff=y*SCREENWIDTH;\
			for(int32_t x=blobs[BNR].rect.x0;x<=blobs[BNR].rect.x1;x++) {\
				int32_t f=db[doff+x];\
				if (db[doff+x] == qf) {\
					md.setP(x,y,ZIELF);\
				}\
			}\
		}\
	}

	printf("identifying cycles ...\n");
	offset=-1;
	// coloring everything interior in the first round
	for(int32_t y=0;y<SCREENWIDTH;y++) {
		for(int32_t x=0;x<SCREENWIDTH;x++) {
			offset++;
			int32_t f=db[offset];
			if (
				(f <= 4)
			) md.setP(x,y,f);
			else if ((f >= BLOBOFFSET) && (f< 256)) {
				md.setP(x,y,SQUARE_BLACK);
			} else {
				if (f >= 256) {
					md.setP(x,y,SQUARE_BLACK);
				} else {
					printf("Implementation error. Blob color %i,%i,F%i\n",x,y,f);
					md.setP(x,y,COLORRED);
				}
			}
		}
	}

	for(int32_t b0=0;b0<anzblobs;b0++) {
		if (blobs[b0].zyklusnr>=0) continue;
		
		for(int32_t i=0;i<anzblobs;i++) visited[i]=0;
		// is b0 in a cycle?
		// is there a targetblob path back to b0 ?
		visited[b0]=1;
		int32_t b=b0;
		int32_t inz=1;
		while (1) {
			if (blobs[b].targetblob == b0) {
				// cycle found
				break;
			}
			int zb=blobs[b].targetblob;
			if (zb < 0) {
				inz=0;
				break;
			}
			if (visited[zb] > 0) {
				inz=0;
				break;
			}
			visited[zb]=1;
			b=zb;
		} // while
		if (inz>0) {
			b=b0;
			printf("Cycle #%i found: ",zyklusanz);
			int l=1;
			while (1) {
				blobs[b].zyklusnr=zyklusanz;
				COLORZYKLUS(b,zyklusanz+BLOBOFFSET);
				if (blobs[b].targetblob == b0) break;
				b=blobs[b].targetblob;
				l++;
			}
			printf("length %i\n",l);
			zyklusanz++;
		}
	} // b0
	
	// lines to indicate diretion of the cycles
	for(int32_t b=0;b<anzblobs;b++) {
		if (blobs[b].zyklusnr<0) continue;
		if (blobs[b].targetblob<0) continue;

		int32_t zx=blobs[blobs[b].targetblob].mx;
		int32_t zy=blobs[blobs[b].targetblob].my;
		
		int32_t f=blobs[b].zyklusnr + BLOBOFFSET;
		if ((f>=BLOBOFFSET)&&(f<256)) {
			int lf=f+13;
			while (lf>255) lf-=256;
			while (lf<BLOBOFFSET) lf+=BLOBOFFSET;
			md.line(f,blobs[b].mx,blobs[b].my,zx,zy,lf);
			md.line(f,blobs[b].mx+4,blobs[b].my,zx,zy,lf);
			md.line(f,blobs[b].mx-4,blobs[b].my,zx,zy,lf);
			md.line(f,blobs[b].mx+2,blobs[b].my,zx,zy,lf);
			md.line(f,blobs[b].mx-2,blobs[b].my,zx,zy,lf);
		}	}
	
	delete[] blobs;
	delete[] db;
	delete[] visited;
	
	char tt[1024];
	sprintf(tt,"%s_periodicity.bmp",fnbase);
	md.save(tt);
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

// Charmap

void Charmap::setlenxy(const int32_t ax,const int32_t ay) {
	if ((xlen>0)&&(cmp)) delete[] cmp;
	
	xlen=ax;
	ylen=ay;
	memused=xlen*ylen;
	cmp=new uint8_t[memused];
}

Charmap::Charmap() {
	xlen=ylen=0;
	memused=0;
	cmp=NULL;
}

Charmap::~Charmap() {
	if ((xlen>0)&&(cmp)) delete[] cmp;
}

void Charmap::fill(const uint8_t swert) {
	if (!cmp) return;
	for(int64_t i=0;i<memused;i++) cmp[i]=swert;
}

void Charmap::clearPalette(void) {
	for(int i=0;i<256;i++) palette[i].R=palette[i].G=palette[i].B=0;
}

void Charmap::setPaletteRGB(const int32_t pos,const uint8_t ar,const uint8_t ag,const uint8_t ab) {
	if ((pos<0)||(pos>255)) return;
	palette[pos].R=ar;
	palette[pos].G=ag;
	palette[pos].B=ab;
}

void Charmap::setP(const int32_t ax,const int32_t ay,const uint8_t awert) {
	if (!cmp) return;
	int64_t pos=(int64_t)ay*xlen+ax;
	if ((pos<0)||(pos>=memused)) return;
	cmp[pos]=awert;
}

uint8_t Charmap::getP(const int32_t ax,const int32_t ay) {
	if (!cmp) return 0;
	int64_t pos=(int64_t)ay*xlen+ax;
	if ((pos<0)||(pos>=memused)) return 0;
	return cmp[pos];
}

void Charmap::save(const char* afn) {
	FILE *fbmp=fopen(afn,"wb");
	write2(fbmp,66,77); 
	int ybytes=ylen; 
	ybytes=(int)(4*ceil(ybytes*0.25));
	
	uint32_t off
		=		14 
			+	40 
			+	256*4
		;
	uint32_t filelen
			=	off
			+	(ybytes*xlen);
		;
			
	fwrite(&filelen,1,sizeof(filelen),fbmp);
	write4(fbmp,0,0,0,0);
	fwrite(&off,1,sizeof(off),fbmp); 
	write4(fbmp,40,0,0,0);
	
	uint32_t w= xlen;
	fwrite(&w,sizeof(w),1,fbmp);
	w = ylen;
	fwrite(&w,sizeof(w),1,fbmp);
	write2(fbmp,1,0); 
	write2(fbmp,8,0); 
	write4(fbmp,0,0,0,0); 
	write4(fbmp,0,0,0,0); 
	write4(fbmp,19,10,0,0);
	write4(fbmp,19,10,0,0);
	write4(fbmp,0,1,0,0); 
	write4(fbmp,0,0,0,0); 
	uint8_t puffer[4];
	for(int32_t i=0;i<256;i++) {
		puffer[0]=palette[i].B;
		puffer[1]=palette[i].G;
		puffer[2]=palette[i].R;
		puffer[3]=0;
		fwrite(puffer,4,sizeof(uint8_t),fbmp);
	}
	
	fwrite(cmp,memused,sizeof(uint8_t),fbmp);

	fclose(fbmp);	
}

void Charmap::line(const uint8_t keepf,const int32_t vx,const int32_t vy,const int32_t zx,const int32_t zy,const uint8_t awert) {
	int32_t dx = zx-vx;
	int32_t dy = zy-vy;

	double x=vx, y=vy;
	double ma;
	if (fabs(dx) > fabs(dy)) ma=(int32_t)ceil(fabs(dx)*3); else ma=(int32_t)ceil(fabs(dy)*3);
	double steps=1.0; steps /= ma;
	double rvx=dx,rvy=dy;
	int32_t draussenq=0;
	int32_t draussenz=0;
	
	for(double  t=0;t<=0.5;t+=steps) {
		// Q -> Z
		y=vy+t*rvy;
		x=vx+t*rvx;
		int32_t xx=(int32_t)floor(x);
		int32_t yy=(int32_t)floor(y);
		if (draussenq>0) setP(xx,yy,awert);
		else {
			int gf=getP(x,y);
			if ( (gf != keepf) && (gf < BLOBOFFSET) ) {
				setP(xx,yy,awert);
				draussenq=1;
			}
		}

		y=vy+(1-t)*rvy;
		x=vx+(1-t)*rvx;
		xx=(int32_t)floor(x);
		yy=(int32_t)floor(y);
		if (draussenz>0) setP(xx,yy,awert);
		else {
			int32_t gf=getP(x,y);
			if ( (gf != keepf) && (gf < BLOBOFFSET) ) {
				setP(xx,yy,awert);
				draussenz=1;
			}
		}
	}
}

int32_t main(int32_t argc,char** argv) {
	// if there's no = sign in any of the argv
	// the old command line parameter structure is used
	// in the other case, the new one
	
	printf("  FUNC=string\n  CMD=calc\n  or CMD=basin\n  or CMD=period\n  LEN=n / c=re,im or re,re,im,im / A=re,im / RVCG=n\n");
	
	int32_t cmd=CMD_CALC;
	
	// standard
	getBoundingBoxfA=getBoundingBoxfA_z2c;
	_FUNC=FUNC_Z2C;
	MAXCOLOREXPECTED=1;
	bitpower=2;
	RANGE0=-2;
	RANGE1=2;
	SCREENWIDTH=1048;
	REVCGBITS=6;
	seedC0re=seedC1re=-1.0;
	seedC0im=seedC0im=0.0;
	FAKTORAre=FAKTORAim=0.0;
	
	// command line structure

	REVCGBITS=6;
	SCREENWIDTH=1024;
	
	for(int i=1;i<argc;i++) {
		upper(argv[i]);
		if (strstr(argv[i],"FUNC=")==argv[i]) {
			if (!strcmp(&argv[i][5],"Z2C")) {
				getBoundingBoxfA=getBoundingBoxfA_z2c;
				_FUNC=FUNC_Z2C;
				MAXCOLOREXPECTED=1;
				bitpower=2;
			} else
			if (!strcmp(&argv[i][5],"Z2AZC")) {
				getBoundingBoxfA=getBoundingBoxfA_z2azc;
				_FUNC=FUNC_Z2AZC;
				MAXCOLOREXPECTED=1;
				bitpower=2;
			} else
			if (!strcmp(&argv[i][5],"Z3AZC")) {
				getBoundingBoxfA=getBoundingBoxfA_z3azc;
				_FUNC=FUNC_Z3AZC;
				MAXCOLOREXPECTED=2;
				bitpower=3;
			} else
			if (!strcmp(&argv[i][5],"Z4AZC")) {
				getBoundingBoxfA=getBoundingBoxfA_z4azc;
				_FUNC=FUNC_Z4AZC;
				MAXCOLOREXPECTED=3;
				bitpower=4;
			} else
			if (!strcmp(&argv[i][5],"Z5AZC")) {
				getBoundingBoxfA=getBoundingBoxfA_z5azc;
				_FUNC=FUNC_Z5AZC;
				MAXCOLOREXPECTED=4;
				bitpower=5;
			} else
			if (!strcmp(&argv[i][5],"Z6AZC")) {
				getBoundingBoxfA=getBoundingBoxfA_z6azc;
				_FUNC=FUNC_Z6AZC;
				MAXCOLOREXPECTED=5;
				bitpower=6;
			}
		} else
		if (strstr(argv[i],"CMD=")==argv[i]) {
			if (!strcmp(&argv[i][4],"BASIN")) {
				cmd=CMD_BASIN;
			} else
			if (!strcmp(&argv[i][4],"PERIOD")) {
				cmd=CMD_PERIOD;
			} 
		} else if (strstr(argv[i],"C=")==argv[i]) {
			double r0,r1,i0,i1; // not NTYP
			// command line parameters are always considered double no matter the datatype used
			if (sscanf(&argv[i][2],"%lf,&lf,&lf,&lf",&r0,&r1,&i0,&i1) == 4) {
				seedC0re=floor(r0*DENOM225)/DENOM225;
				seedC1re=floor(r1*DENOM225)/DENOM225;
				seedC0im=floor(i0*DENOM225)/DENOM225;
				seedC1im=floor(i1*DENOM225)/DENOM225;
			} else
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				seedC0re=seedC1re=floor(r0*DENOM225)/DENOM225;
				seedC0im=seedC1im=floor(i0*DENOM225)/DENOM225;
			}
		} 
		else if (strstr(argv[i],"A=")==argv[i]) {
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
				RANGE0=-a;
				RANGE1= a;
			}
		}
	} // i
	
	COMPLETE0=RANGE0;
	COMPLETE1=RANGE1;
	
	if (REVCGBITS < 4) REVCGBITS=4;

	REVCGBLOCKWIDTH=(1 << REVCGBITS);
	if (SCREENWIDTH >= REVCGBLOCKWIDTH) {
		REVCGmaxnumber=SCREENWIDTH >> REVCGBITS;
	} else REVCGmaxnumber=1;
	REVCGmaxnumberQ=REVCGmaxnumber*REVCGmaxnumber;
	
	if (SCREENWIDTH < 256) SCREENWIDTH=256;

	MEMWIDTH=(SCREENWIDTH >> 4); // 16 Pixel per 32 bit integer
	SPLITAFTER=(int)floor( (double)(1 << 30) / MEMWIDTH);
	
	test_for_sufficient_precision_bits();
	
	REFINEMENTLEVEL=(int)ceil(log(SCREENWIDTH)/log(2.0));
	
	char fn[1024],tmp2[1024],tmp3[1024];
	fn[0]=0;
	
	if (_FUNC==FUNC_Z2C) {
		sprintf(fn,"_L%02i_%sz2c_%s.bmp",
		REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2));
	} else
	if (_FUNC==FUNC_Z2AZC) {
		sprintf(fn,"_L%02i_%sz2azc_%s_%s.bmp",
		REFINEMENTLEVEL,
		NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
	} else
	if (_FUNC==FUNC_Z3AZC) {
		sprintf(fn,"_L%02i_%sz3azc_%s_%s.bmp",
		REFINEMENTLEVEL,
		NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
	} else
	if (_FUNC==FUNC_Z4AZC) {
		sprintf(fn,"_L%02i_%sz4azc_%s_%s.bmp",
		REFINEMENTLEVEL,
		NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
	} else
	if (_FUNC==FUNC_Z5AZC) {
		sprintf(fn,"_L%02i_%sz5azc_%s_%s.bmp",
		REFINEMENTLEVEL,
		NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
	} else
	if (_FUNC==FUNC_Z6AZC) {
		sprintf(fn,"_L%02i_%sz6azc_%s_%s.bmp",
		REFINEMENTLEVEL,
		NNTYPSTR,seedCstr(tmp2),FAKTORAstr(tmp3));
	}
	
	printf("file principal part %s\n",fn);
	
	data5=new Data5;

	scaleRangePerPixel = (RANGE1-RANGE0);
	scaleRangePerPixel /= SCREENWIDTH;
	scalePixelPerRange = SCREENWIDTH;
	scalePixelPerRange /= (RANGE1-RANGE0);;
	
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

	// main routine
	// has to be done for CALC, BASIN, PERIOD
	compute(); 

	// storing image(s)
	printf("saving image ...\n");
	data5->saveBitmap4(fn);
	// storing a trustworthily downscaled image
	if (SCREENWIDTH >= 16384) {
		printf("downscaling 16-fold in a trustworthy manner ...\n");
		data5->saveBitmap4_trustworthily_downscaled_16fold(fn);
	}
	// storing raw data
	printf("saving raw data ...\n");
	data5->saveRaw(fn);
	
	// data is now computed or loaded
	if (cmd==CMD_BASIN) {
		basin(fn);
	} else
	if (cmd==CMD_PERIOD) {
		periodicity(fn);
	} 
	
	
	delete data5;
	
    return 0;
}

