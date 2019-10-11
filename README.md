# julia-tsa-core
Julia sets with a mathematical guarantee

Implementation of the excellent articles:
```
"Images of Julia sets that you can trust"
by Luiz-Henrique de Figueiredo, Diego Nehab, Jofge Stolfi, Joao Batista Oliveira-
from 2013

"Rigorous bounds for polynomial Julia sets"
by Luiz-Henrique de Figueiredo, Diego Nehab, Jofge Stolfi, Joao Batista Oliveira-
```

## Purpose of this C++-code:

Command line tool to compute polynomial Julia sets of the form z^n+A*z+c for degree 2 to 6 with a mathematical guarentee, their basins of attraction and the periodic cycles.

### Disclaminer

Although I tried my best to ensure correctness of the implementation and therefore the results,
the code comes with no warranty.

## README is organized as:

0. Quick start 
1. Background information on the algorithm and the current implementation
2. Computing some example sets explaining the features of the code
3. Command-line parameters
4. Contact and links


## (0) Quick starting example:

![example Julia sets](./_juliasets.gif)

Compile main.cpp with any suitable C compiler supporting 64-bit integers (best with mathematical optimizations
for speed enabled). 

Start executable (assuming it is called juliatsacore_d.exe throughout README - d stands for the C++ double datatype,
other types see (1)) by entering in a command prompt:

`juliatsacore_d.exe func=z3azc len=10 cmd=calc a=1.25,0 c=0,0.025 range=4` 

The file `_L10_z3azc_c_ia_0_0_x_0.025_0.025_A_1.25_0.bmp_Y00X00.bmp` contains a 4-bit bitmap with the computed set.

Note the black interior, the white exterior and the gray cells containing the boundary as described in the articles.

A very interesting example - __float128 precision needed to be mathematically accurate - for basins of attraction is:

`juliatsacore_qd.exe func=z6azc cmd=basin len=12 a=1.0625,0 c=0,0.0234375 revcg=4 range=4`

Or start the batch files `_basin.bat` or `_periodicity.bat` - the resulting bitmaps have suffix `period` or `basin`.



## (1) Background

Contrary to the article, I did not implement a quad-tree representation. The code works
on a contiguous array (or several parts of 2 GB length) to get the fastest access times.
Each screen pixel represents a square in the complex number plane and needs 2 bits in memory for
color storage, grouped as 16 pixels in a 32-bit integer. 

For the final image, the following claims hold:

A black pixel means ALL complex numbers within the underlying square (including the edges
and corners) are starting points for bounded orbits and belong to the interior of the filled-in Julia set.

A white pixel means ALL complex numbers within the underlying square (including the edges
and corners) are starting points for escaping orbits and belong to the exterior.

Only gray pixel can (but do not have to) contain complex numbers that belong to the Julia set boundary itself 
and might or might not contain numbers belonging to the exterior and/or interior. In a larger image those gray squares 
will be split up and might then be colored to a definite interior or exterior.

During computation there are 4 colors: *black* and *white* as above, *gray* meaning pixel is unclear, and *gray-but-potentially-white*.
The latter meaning that at least one path from that square leads to the exterior. 

The software works in a main loop
in the routine `propagate_white` that goes over the image again and again until no color change occurs anymore, then all gray pixels are bounded and can be colored black. Furthermore if computation of an image
is continued in a larger screen width, the software also recognizes gray cells as being black if
all paths lead directly in one iteration to other black cells - by backpropagation all
iteration lengths are recognized.

The software takes command line arguments (in no particular order) and constructs the desired Julia set.
Ever so often temporary raw data is written to the hard drive, so the computation can be stopped
by simply closing the command prompt window.

Memory overhead to speeden up the computation includes that every row has a specific min and max
x coordinate value where the gray for that row resides. Those values will be adjusted if pixels are
judged as interior or exterior. At the beginning a static reverse cell graph in a low resolution
(usually 16x16 to 256x256 pixels were put together in a tile) and the preimages of every tile are
computed. When a gray pixel changes its color the preimages' tiles are set to *have to be visited*
and checked in one of the next rounds.

The desired and necessary C++ data type for the computation of the bounding box can be commented in or out
at the start of the source code:

`#define _DOUBLE`
`#define _LONGDOUBLE`
`#define _QUADMATH`

and then recompiling the software.

The main routine is `compute()`, the main struct is `data5` and the most important variables are `SCREENWIDTH`
and `seedC0re`, `seedC0im`, `seedC1re`, `seedC1im` and `FAKTORAre`, `FAKTORAim`.

The software first tries to find a file named `_in.raw_header`. If that exists
the raw data (in the files named `_in.raw_0001`, `_in_raw_0002` and so on) is loaded and 
computation will be resumed, i.e. starting with searching for new white cells by the function `propagate_white()`.

A special case is when the read-in data was computed with half the current command line screen
width. In that case the read-in image is doubled, every pixel is copied to a 2x2 grid keeping
its color and computation resumes in the higher screen width - simulating the refinement process
as described in the article. This might save computation time when screen width is around 64k or higher by using
previously computed information (for 16k or less I just compute the image from scratch). Unformation of potentially-white cells is ignoried during that blow-up process.

The software does not perform many error checks and was designed mainly for speed and complete memory usage.
If an error occurs, it prompts a message and exits immediately "dirty", leaving garbage collection
to the command-line window.

The output of the code are several files:

`_L10_z3azc_c_ia_0_0_x_0.025_0.025_A_1.25_0.bmp_Y00X00.bmp`
	which is a 4-bit bitmap representing the image (at most 2 GB in size, if the resulting file would be bigger,
	the image is tiled into appropriate parts (hence the YX-part of the name), starting
	from bottom to top and left to right)

`..._twd16.bmp`
	if the screen width is 2^14 or more, the current data is 16-fold downsized in a trustworthy manner 
	- meaning, if a 16x16 grid is all in one color, the resulting 1 pixel in the downscaled image is set 
	to that color, otherwise it will be gray.
	That image is just saved without file size check for now as an 8-bit bitmap.

files named `...raw_NNNN` with NNNN being a number from 1 to as much as needed (again,
	2 GB as maximum file size) and a file `...raw_header`. Those store the raw data so that
	the computation can be resumed or a bigger screen width can be computed using already generated information - if files are renamed appropriately to _in.raw_NNNN etc.

the `_temp.raw_*` files as mentioned above, in case the computation is interrupted by user, software or
	hardware crash. Those can be renamed to _in.raw_NNNN / header also and used to restart computation.

A bitmap ending with `...basins.bmp` (if command specified) contains the basin of attractions differently colored, similarily
files ending in `...periodicity.bmp`show the periodic cycle, its immediate basins and the direction in
which the cycle is traversed - from thicker to thinner blue lines.

Number representability is a prerequisite for the polynomial coefficents given in the command line, 
any value p/q where q is a power
of 2 can be accurately represented. Fractional numbers provided for c and A are internally made
representable by performing: floor(number provided * 2^25) / 2^25.

The software performs a preliminary check if the currently
used C++ data type (double, long double, __float128) is offering enough bits for the current resolution,
function and complex number range.

In the quadratic case, a bailout of 2 (range=2) is mathematically sufficient. For higher 
order polynomials, the value must be adapted to accomodate for larger shapes. The complex
plane represented on the screen goes from -bailout to +bailout in both axis. Integer bailout is
used so that a pixel has a plane width that can always be accurately represented (2*bailout / SCREENWIDTH).

The first article above, paragraph 7, **Extension to higher-degree polynomials** depicts
an estimate for general polynomials. For the functions provided in the software a value 2+|c|+|A| is sufficient, usually for the values used, 4 is taken. Too high a bailout is working, too low will compromise the mathematical guarantee.


## (2) Examples and features

#### (a) Computing small images from scratch

Enter one of the following lines in a command prompt.

`juliatsacore_d.exe cmd=calc func=z2c c=-1,0 len=10`

`juliatsacore_qd.exe func=z6azc cmd=basin len=11 a=1.0625,0 c=0,0.0234375 revcg=4 range=4`


#### (b) Resuming computation

Start a computation with the following command:

`juliatsacore_d.exe func=z2c cmd=calc c=-1,0 len=14`

At some point the message *saving raw data as temporary* appears. If the computation
continues after that, stop the command-line window before the image is calculated
completely.

The software saved the files `_temp.raw_header` and `_temp.raw_0001`.
Rename those to `_in.raw_header` and `_in.raw_0001`.

And start the computation again with the same command. The software then continues judging
pixels using all the information that has been calculated thus far.

Delete the `_in_raw*` files if no longer needed as the software always first tries to use those files if present.

#### (c) Increasing an image to a higher screen width

Compute a small version of the basilica z^2-1 (deleting `_in.raw*` files beforehand).

`juliatsacore_d.exe len=11 cmd=calc c=-1,0 func=z2c`

Rename the final files `_L11_...raw_0001` to `_in.raw_0001` and `_L11_....raw_header` to `_in.raw_header`.

Then start the software again with double the screen width, i.e. len=12:

`juliatsacore_d.exe c=-1,0 func=z2c cmd=calc len=12`

The software uses the already computed image to build a bigger version. Afterwards delete `_in.raw*` files.


## (3) Command-line parameters;

`FUNC=string` (if not provided, standard value is Z2C)
The desired polynomial to use. Implemented are z2azc, z3azc ... until z6azc for the polynomials z^n+A*z+c, where
the degree n is from 2 to 6 and the factor A is a complex number. For performance reasons a special function
FUNC=z2c for z^2+c is also implemented which uses a more optimized version of computing the bounding box (time benefit not measured).

`CMD=string` (standard value CALC)
1. CMD=CALC: The software computes the set, saves the image and the final data.
2. CMD=BASIN: Determining the basin(s) of attraction. Like with any start it first
searches for the files _in.raw... (described above), if not present, the Julia set is calculcated before
the basins of attraction are being searched for and colored differently. The final image result is in a bitmap with ending `_basins.bmp? 
3. CMD=PERIOD: The software tries to find the immediate basins of attracting periodic points if present
and colors them. If no ...in.raw files are present, the set is calculated from scratch as for BASIN.

`LEN=integer` (standard value 10)
The screen width is set to 2^integer pixel.
Images must be at least 2^8 pixels and can go up to
2^31 in principle. The largest I computed thus far is, however, 2^19 pixels in width.

'C=double1,double2`
Sets the seed value: double1 as real part, double2 as imaginary part.

`C=double1,double2,double3,double4`
Sets the seed value to a complex interval: The real part being [double1..double2], the imaginary part
being [double3..double4]

`A=double1,double2`
Sets the degree-1 coefficiant accordingly.

<b>Note</b> that real numbers given in the command line as parameters are always treated as C++ double, no matter 
what underlying data type is used in the binary.

`RANGE=integer` (standard value 2)
The complex plane where the whole Julia set is definitely contained is set to -integer .. +integer in
both axis.

`REVCG=integer` (standard value 6)
The reverse cell graph to speeden computation up uses groups of 2^integer x 2^integer pixels.
The larger the integer, the less space the reverse cell graph needs and more memory is usable for
the image itself. Whether the reverse graph speeds the current computation up or not is dependent
on the function and the coefficients. 
Numbers below 4 are not possible. 

<b>Note</b> if available memory is not sufficient for image size and size of reverse cell graph,
the program will terminate with a memory bad_alloc error message. Increasing the REVCG parameter
will reduce memory usage.


## (4) Contact

Please direct any comments to:

marcm200@freenet.de

forum: https://fractalforums.org/fractal-mathematics-and-new-theories/28/julia-sets-true-shape-and-escape-time/2725

Marc Meidlinger, July-October 2019

