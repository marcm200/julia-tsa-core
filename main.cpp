/*
	My (Marc Meidlinger, July 2019)	implementation of the ingenious 
	trustworthy Julia set algorithm from the article:
	
	"Images of Julia sets that you can trust" by Luiz Henrique de Figueiredo, Diego Nehab,
	Jorge Stolfi, Joao Batista Oliveira from 2013
*/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"


// const definitions

const int32_t MAXPARENT=1024;

// pixel/square colors in 2bit format
const uint32_t SQUARE_GRAY=0b00;
const uint32_t SQUARE_WHITE=0b01;
const uint32_t SQUARE_BLACK=0b10;
const uint32_t SQUARE_GRAY_POTENTIALLY_WHITE=0b11;

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


// structs

// palette-entry for bitmaps
struct RGB4 {
	uint8_t R,G,B,alpha;
};

// a rectangle in the complex plane - used for a square and its bounding box
struct PlaneRect {
	double x0,x1,y0,y1;
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
	// the parents as a fixed-sized-array/ (future plans: memory manager)
	Parent parent[MAXPARENT];
		
	RevCGBlock();
	void addParent(const int32_t,const int32_t);
};

// pixel-coordinate rectangle
struct ScreenRect {
	int32_t x0,x1,y0,y1;
};

struct Gray_in_row {
	int32_t g0,g1;
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


// globals

double seedCre,seedCim; 
Data5 *data5;
double scaleRangePerPixel,scalePixelPerRange;
int64_t countsquares_white,countsquares_gray;
int64_t countsquares_black,countsquares_graypotw;
int32_t SPLITAFTER;
// region in the complex plane where all the currently still gray squares reside
double planegrayx0,planegrayx1;
double planegrayy0,planegrayy1;
// region in screen coordinates
int32_t encgrayx0,encgrayx1;
int32_t encgrayy0,encgrayy1;
// width of screen in pixel, must be a power of 2
int32_t SCREENWIDTH;
// number of 32bit-integers a screen row holds, equal to (SCREENWIDTH >> 4)
int32_t MEMBREITE;
int32_t RANGE0,RANGE1;
// complte holds the same values as range
double COMPLETE0,COMPLETE1;
// for the low-resolution reverse cell graph working on (usually) 64x64 pixel squares or bigger
int32_t REVCGBITS,REVCGBLOCKWIDTH;
int32_t REVCGmaxnumber,REVCGmaxnumberQ;


// forward declarations

// main function to compute the Julia set
void compute(void);
// defines the function used for iteration
inline void getBoundingBoxfA(PlaneRect&,PlaneRect&);
// constructing the low-level reverse cell-graph
void construct_static_reverse_cellgraph(void);
// squares that directly hit the special exterior outside COMPLETE
void find_special_exterior_hitting_squares(void);
// white and potentially-white cells wiul
void propagate_white(void);
void coloring_interior_black(void);
inline int32_t scrcoord_as_lowerleft(const double);
void copy_pixel_to_2x2grid(const uint32_t,uint32_t*);

void write2(FILE*,const uint8_t,const uint8_t);
void write4(FILE*,const uint8_t,const uint8_t,const uint8_t,const uint8_t);
inline double minimumD(const double,const double);
inline double maximumD(const double,const double);
inline double minimumD(const double,const double,const double,const double);
inline double maximumD(const double,const double,const double,const double);


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

void RevCGBlock::addParent(const int32_t ax,const int32_t ay) {
	if (howmany >= MAXPARENT) {
		fprintf(stderr,"Error. addParent: too many parents\n");
		exit(99);
	}
	
	// if box lies outside the screen, it will never be marked as to be visitied, so the
	// block can be ignored
	if (
		(ax < 0) ||
		(ay < 0) ||
		(ax >= REVCGmaxnumber) ||
		(ay >= REVCGmaxnumber) 
	) return; 
	
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
	uint32_t slen=SPLITAFTER*MEMBREITE;
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
		
		uint32_t *p=new uint32_t[allok*MEMBREITE];
		if (!p) {
			fprintf(stderr,"Memory error. data5 consstructor\n");
			exit(99);
		}
		
		// setting the points to the array
		// portion that belongs to the row
		for(int32_t i=0;i<SPLITAFTER;i++) {
			pixelrow[y]=&p[i*MEMBREITE];
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

inline int32_t scrcoord_as_lowerleft(const double a) {
	// calculating the screen coordinte of the pixel that contains the coordinate
	// if the coordinate lies on an edge/corner (and belongs to more than one pixel)
	// the pixel where it lies on the left,bottom edge/corner is returned
	return (int)floor( (a - COMPLETE0) * scalePixelPerRange );
}

// maximum of 4 doubles. often used by the boundingbox function
inline double maximumD(const double a,const double b,const double c,const double d) {
	double m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}

inline double minimumD(const double a,const double b,const double c,const double d) {
	double m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	return m;
}

inline double minimumD(const double a,const double b) {
	if (a < b) return a;
	return b;
}

inline double maximumD(const double a,const double b) {
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
inline void getBoundingBoxfA(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+seedCre;
	fA.x1=maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+seedCre;
	fA.y0=2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedCim;
	fA.y1=2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedCim;
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
	}
	
	PlaneRect A,bbxfA;
	
	A.y1=COMPLETE0;
	const double DD=REVCGBLOCKWIDTH*scaleRangePerPixel;
	int32_t parentx,parenty;
	
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
			
			for(int32_t by=scr.y0;by<=scr.y1;++by) {
				int32_t yoffset=by*REVCGmaxnumber;
				for(int32_t bx=scr.x0;bx<=scr.x1;++bx) {
					data5->revcgYX[yoffset+bx].addParent(parentx,parenty);
				}
			}
		} // x
	} // y
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
	
	const double DD16SCALE=16.0*scaleRangePerPixel;
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
			loops_till_saving_raw_data=10;
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

int32_t main(int32_t argc,char** argv) {
	// argv: scrwith seedreal seedimag rangeright tilebits
	if ( (argc<=1) || (sscanf(argv[1],"%i",&SCREENWIDTH) != 1) ) SCREENWIDTH=8192;
	if (SCREENWIDTH < 256) SCREENWIDTH=256;
	if ( (argc<=2) || (sscanf(argv[2],"%lf",&seedCre) != 1) ) seedCre=(double)12967936.0/33554432.0;
	if ( (argc<=3) || (sscanf(argv[3],"%lf",&seedCim) != 1) ) seedCim=(double)4710400.0/33554432.0;
	if ( (argc<=4) || (sscanf(argv[4],"%i",&RANGE1) != 1) ) RANGE1=2;
	if (RANGE1 < 2) RANGE1=2;
	if ( (argc<=5) || (sscanf(argv[5],"%i",&REVCGBITS) != 1) ) {
		if (SCREENWIDTH <= 65536) REVCGBITS=6;
		if (SCREENWIDTH == (2*65536)) REVCGBITS=7;
		else REVCGBITS=8;
	}
	if (REVCGBITS < 6) REVCGBITS=6;

	RANGE0=-RANGE1;
	COMPLETE0=RANGE0;
	COMPLETE1=RANGE1;
	REVCGBLOCKWIDTH=(1 << REVCGBITS);
	if (SCREENWIDTH >= REVCGBLOCKWIDTH) {
		REVCGmaxnumber=SCREENWIDTH >> REVCGBITS;
	} else REVCGmaxnumber=1;
	REVCGmaxnumberQ=REVCGmaxnumber*REVCGmaxnumber;

	MEMBREITE=(SCREENWIDTH >> 4); // 16 Pixel per 32 bit integer
	SPLITAFTER=(int)floor( (double)(1 << 30) / MEMBREITE);

	data5=new Data5;

	scaleRangePerPixel = (RANGE1-RANGE0);
	scaleRangePerPixel /= SCREENWIDTH;
	scalePixelPerRange = SCREENWIDTH;
	scalePixelPerRange /= (RANGE1-RANGE0);;

	// main routine
	compute(); 

	// storing image(s)
	printf("saving image ...\n");
	data5->saveBitmap4("_tsa_juliaset");
	// storing a trustworthily downscaled image
	if (SCREENWIDTH >= 16384) {
		printf("downscaling 16-fold in a trustworthy manner ...\n");
		data5->saveBitmap4_trustworthily_downscaled_16fold("_tsa_juliaset");
	}
	// storing raw data
	printf("saving raw data ...\n");
	data5->saveRaw("_tsa_juliaset");
	
	delete data5;
	
    return 0;
}

