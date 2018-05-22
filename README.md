Contents
========

1.  [Description](#description)
2.  [License](#license)
3.  [Requirements](#requirements)
4.  [Installation](#installation)
5.  [Instructions](#instructions)
6.  [Contact](#contact)
7.  [Contributing](#contributing)

Description
===========

| ![image](https://bytebucket.org/kuraiev/halftones/raw/a4b3975053b1a37a8e8c085e040e50708e9ab890/test/lena1.jpg =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/a4b3975053b1a37a8e8c085e040e50708e9ab890/test/halftone_shiaufan.png =512x512) |
|:-:|:-:|
| Original | Halftone via Error diffusion (Shiau-Fan dithering) |

Halftones is a Python library used to compute several halftones types and them inverses. This code was developed in part based on
the following papers:

1. Freitas, Pedro Garcia, et al. "[Error Concealment Using a Halftone
Watermarking Technique](http://dx.doi.org/10.1109/SIBGRAPI.2012.50)." Graphics, Patterns and Images (SIBGRAPI), 2012
25th SIBGRAPI Conference on. IEEE, 2012.

2. Freitas, Pedro Garcia, Mylene CQ Farias, and Aletéia PF de Araújo. "[Fast Inverse Halftoning Algorithm for Ordered Dithered Images](http://dx.doi.org/10.1109/SIBGRAPI.2011.14)." Graphics, Patterns and Images (Sibgrapi), 2011 24th SIBGRAPI Conference on. IEEE, 2011.


License
=======

Halftones is released under [GNU GPL version
2.](http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt)


Requirements
============

Halftones requires Python 2.6 (or above), Scipy 0.10.0 (or above) compiled with Fortran support (so, a Fortran compiler must be installed and configured in your system), Scikit-learn and Scikit-image. This project was also tested on Python 3 and it works fine.

Installation
===========

1. Clone the last version
> `# hg clone https://bitbucket.org/kuraiev/halftones`
2. Go to the directory where 'setup.py' os located
> `cd halftones`
3. Install using distributils
> `# sudo python setup.py install`

Instructions
============

Here is some examples of halftones computed using this library.

Ordered Dithering
-----------------

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_ordered_dither_balanced_centered_points.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_ordered_dither_central_white_points.png =512x512) |
|:-:|:-:|
| Balanced centered point | Central white point |

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_ordered_dither_clustered_dots.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_ordered_dither_diagonal_matrix.png =512x512) |
|:-:|:-:|
| Clustered dots | Diagonal ordered matrix with balanced centered points |


Error Diffusion
---------------

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_atkinson.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_burkes.png =512x512) |
|:-:|:-:|
| Atkinson's | Burkes |

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_floyd_steinberg.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_jarvis.png =512x512) |
|:-:|:-:|
| Floyd-Steinberg | Jarvis |

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_stucki.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_sierra3.png =512x512) |
|:-:|:-:|
| Stucki | Sierra (3x3) |

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_sierra2.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_sierra_simple.png =512x512) |
|:-:|:-:|
| Sierra (2x2) | Sierra (simple) |

Inverse Halftoning
------------------

| ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/halftone_ordered_dither_dispersed_dots.png =512x512) | ![image](https://bytebucket.org/kuraiev/halftones/raw/039ce90f32056a056f116a50e40c342ceae46989/test/inverse_jarvis.png =512x512) |
|:-:|:-:|
| Dispersed dots halftoning | Inverse halftoning (grayscale) |


The repository 'test' contains some examples of how to use the functions.

Contact
=======

Please send all comments, questions, reports and suggestions (especially if 
you would like to contribute) to **sawp@sawp.com.br**

Contributing
============

If you would like to contribute with new algorithms, increment the code
performance, documentation or another kind of modifications, please
contact me. The only requirements are: keep the code compatible with
PEP8 standardization and licensed by GPLv2.

References
==========

[1]  Foley, J.D. and A. van Dam, Fundamentals of Interactive Computer
     Graphics, Addison-Wesley, Reading, MA, 1982.

          This is a standard reference for many graphic techniques which has
          not declined with age.  Highly recommended.  This edition is out
          of print but can be found in many university and engineering
          libraries.  NOTE: This book has been updated and rewritten, and
          this new version is currently in print as:

     Foley, J.D., A. van Dam, S.K. Feiner, and J.F. Hughes;  Computer
     Graphics: Principles and Practice. Addison-Wesley, Reading, MA, 1990.

          This rewrite omits some of the more technical data of the 1982
          edition, but has been updated to include information on error-
          diffusion and the Floyd-Steinberg filter.  Currently on computer
          bookstore shelves and rather expensive (around $75 list price).

[2]  Bayer, B.E., "An Optimum Method for Two-Level Rendition of Continuous
     Tone Pictures," IEEE International Conference on Communications,
     Conference Records, 1973, pp. 26-11 to 26-15.  

          A short article proving the optimality of Bayer's pattern in the
          dispersed-dot ordered dither.  

[3]  Ulichney, R., Digital Halftoning, The MIT Press, Cambridge, MA, 1987.

          This is the best book I know of for describing the various black
          and white dithering methods.  It has clear explanations (a little
          higher math may come in handy) and wonderful illustrations.  It
          does not contain any code, but don't let that keep you from
          getting this book.  Computer Literacy normally carries it but the
          title is often sold out.  

          [MFM note:  I can't describe how much information I got from this
          book!  Several different writers have praised this reference to
          the skies, and I can only concur.  Some of it went right over my
          head -- it's heavenly for someone who is thrilled by Fourier
          analysis -- but the rest of it is a clear and excellent treatment
          of the subject.  I had to request it on an interlibrary loan, but
          it was worth the two weeks' wait and the 25 cents it cost me for
          the search.  University or engineering libraries would be your
          best bet, as would technical bookstores.]

[4]  Floyd, R.W. and L. Steinberg, "An Adaptive Algorithm for Spatial Gray
     Scale."  SID 1975, International Symposium Digest of Technical Papers,
     vol 1975m, pp. 36-37.

          Short article in which Floyd and Steinberg introduce their filter. 
          
[5]  Daniel Burkes is unpublished, but can be reached at this address:

          Daniel Burkes
          TerraVision, Inc.
          2351 College Station Road, Suite 563
          Athens, GA  30305

     or via CIS at UID# 72077,356.  The Burkes error filter was submitted to
     the public domain on September 15, 1988 in an unpublished document,
     "Presentation of the Burkes error filter for use in preparing
     continuous-tone images for presentation on bi-level devices."  The file
     BURKES.ARC, in LIB 15 (Publications) of the CIS Graphics Support Forum,
     contains this document as well as sample images.

[6]  Jarvis, J.F., C.N. Judice, and W.H. Ninke, "A Survey of Techniques for
     the Display of Continuous Tone Pictures on Bi-Level Displays," Computer
     Graphics and Image Processing, vol. 5, pp. 13-40, 1976.

[7]  Stucki, P., "MECCA - a multiple-error correcting computation algorithm
     for bilevel image hardcopy reproduction."  Research Report RZ1060, IBM
     Research Laboratory, Zurich, Switzerland, 1981.

[8]  Heckbert, P. "Color Image Quantization for Frame Buffer Display." 
     Computer Graphics (SIGGRAPH 82), vol. 16, pp. 297-307, 1982.

[9]  Frankie Sierra is unpublished, but can be reached via CIS at UID#
     76356,2254.  Pictorial presentations of his filters can be found in LIB
     17 (Developer's Den) of the CIS Graphics Support Forum as the files
     DITER1.GIF, DITER2.GIF, DITER6.GIF, DITER7.GIF, DITER8.GIF, and
     DITER9.GIF.

[10] J.F.R. "Frank" Slinkman is unpublished, but can be reached via CIS at
     UID# 72411,650.  The file NUDTHR.ARC in LIB 17 (Developer's Den) of the
     CIS Graphics Support Forum contains his document "New Dithering Method
     for Non-Square Pixels" as well as sample images and encoding program.

[11] Lawrence Gozum is unpublished, but can be reached via CIS at UID#
     73437,2372.  His document "Notes of IDTVGA Dithering Method" can be
     found in LIB 17 (Developer's Den) of the CIS Graphics Support Forum as
     the file IDTVGA.TXT.

[12] Robert M. Crawford is unpublished, but can be reached via CIS at UID#
     76356,741.  The file DGIF.ZIP in LIB 17 (Developer's Den) of the CIS
     Graphics Support Forum contains documentation, sample images, and demo
     program.

[13] Knuth, D.E., "Digital Halftones by Dot Diffusion." ACM Transactions on
Graphics, Vol. 6, No. 4, October 1987, pp 245-273.

[14] Rogers, D.F., Procedural Elements for Computer Graphics, McGraw-Hill, New
York, 1985.

[15] Kuto, S., "Continuous Color Presentation Using a Low-Cost Ink Jet Printer,"
Proc. Computer Graphics Tokyo 84, 24-27 April, 1984, Tokyo, Japan.

[16] Mitchell, W.J., R.S. Liggett, and T. Kvan, The Art of Computer Graphics 
Programming, Van Nostrand Reinhold Co., New York, 1987.

[17] Pavlidis, T., Algorithms for Graphics and Image Processing, Computer Science
Press, Rockville, MD, 1982.

[18] Shiau, Jeng-nan, and Zhigang Fan. ``Set of easily implementable 
     coefficients in error diffusion with reduced worm artifacts.''
     Electronic Imaging: Science & Technology. International Society
     for Optics and Photonics, 1996.
