# -*- coding: utf-8 -*-
#
# halft_np.ones.py - Modules for halftone
#
#
# Autor: Pedro Garcia Freitas [sawp@sawp.com.br]
# License: GPL V.2
#
# see: http://www.sawp.com.br
#
# Sep 2011
#
import halftones.compiled as _compiled
import numpy as _np
import scipy.misc as _spm


def __get_ordered_threshold(type="clustered dot"):
    """
    Compute the threshold matrix

    Parameters
    ----------
    type: string with threshold type name

    Return
    -------
    matrix: threshold matrix
    """
    if type.lower() == "diagonal ordered matrix":
        S1 = _np.array([[13, 9, 5, 12],
                        [6, 1, 0, 8],
                        [10, 2, 3, 4],
                        [14, 7, 11, 15]])
        S2 = _np.array([[18, 22, 26, 19],
                        [25, 30, 31, 23],
                        [21, 29, 28, 27],
                        [17, 24, 20, 16]])
        S21 = (_np.concatenate((S1, S2), 1), _np.concatenate((S2, S1), 1))
        M = _np.concatenate(S21)
    elif type.lower() == "clustered dots":
        M = _np.array([[34, 29, 17, 21, 30, 35],
                       [28, 14, 9, 16, 20, 31],
                       [13, 8, 4, 5, 15, 19],
                       [12, 3, 0, 1, 10, 18],
                       [27, 7, 2, 6, 23, 24],
                       [33, 26, 11, 22, 25, 32]])
    elif type.lower() == "central white point":
        M = _np.array([[34, 25, 21, 17, 29, 3],
                       [30, 13, 9, 5, 12, 24],
                       [18, 6, 1, 0, 8, 20],
                       [22, 10, 2, 3, 4, 16],
                       [26, 14, 7, 11, 15, 28],
                       [35, 31, 19, 23, 27, 32]])
    elif type.lower() == "balanced centered point":
        M = _np.array([[30, 22, 16, 21, 33, 35],
                       [24, 11, 7, 9, 26, 28],
                       [13, 5, 0, 2, 14, 19],
                       [15, 3, 1, 4, 12, 18],
                       [27, 8, 6, 10, 25, 29],
                       [32, 20, 17, 23, 31, 34]])
    elif type.lower() == "dispersed dots":
        M = _np.array([[32, 16, 20, 34, 18, 22],
                       [12, 0, 4, 14, 2, 6],
                       [28, 8, 24, 30, 10, 26],
                       [35, 19, 23, 33, 17, 21],
                       [15, 3, 7, 13, 1, 5],
                       [31, 11, 27, 29, 9, 25]])
    return M


def __ordered_dithering(I, type="clustered dot"):
    """
    Encode the image using a ordered dithering type

    Parameters
    ----------
    im: input grayscale image
    type: string with threshold type (diagonal ordered matrix, clustered dots
          central white point, balanced centered point, dispersed dots)

    Return
    -------
    binary: output image
    """
    def __round(i, j):
        pos = _np.ceil(float(i) / float(j))
        return pos.astype(_np.int64)
    S = __get_ordered_threshold(type)
    si = I.shape
    ss = S.shape
    ts = (__round(i,j) for i, j in zip(si, ss))
    Simg = _np.tile(S, ts)
    Simg = Simg[:si[0], :si[1]]
    Simg = Simg + 0.5
    N = S.max() - S.min() + 2
    D = 255.0 / (N - 1)
    Q = I.astype(_np.uint8) / D
    O = Q > Simg
    return 255 * O


def ordered_dithering_diagonal_ordered_matrix(grayscale):
    """ Diagonal ordered matrix with balanced centered points

    dithered = ordered_dithering_diagonal_ordered(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    return __ordered_dithering(grayscale, "diagonal ordered matrix")


def ordered_dithering_clustered_dots(im):
    """ Clustered dots

    dithered = ordered_dithering_clustered_dots(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    return __ordered_dithering(im, "clustered dots")


def ordered_dithering_central_white_point(im):
    """ Central white point ordered dithering

    dithered = ordered_dithering_central_white_point(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    return __ordered_dithering(im, "central white point")


def ordered_dithering_balanced_centered_point(im):
    """ Balanced centered point ordered dithering

    dithered = ordered_dithering_balanced_centered_point(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    return __ordered_dithering(im, "balanced centered point")


def ordered_dithering_dispersed_dots(im):
    """ Dispersed dots (Bayer with 3x3 threshold Matrix)

    dithered = ordered_dithering_dispersed_dots(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    return __ordered_dithering(im, "dispersed dots")


def ordered_combinatorial3(im):
    """
    Encode the image using combinatorial dithering for halftone with 2^3 levels

    Parameters
    ----------
    im: input image

    Return
    -------
    binary: output image
    """
    masked = (8.0 / 255.0) * _spm.imfilter(im, 'edge_enhance')
    masked = _np.floor(masked).astype('uint8')
    masked[masked > 7] = 7
    binary = _compiled.halftone.ordered_comb3_iterator(masked)
    return binary


def ordered_combinatorial2(im):
    """
    Encode the image using combinatorial dithering for halftone with 2^2 levels

    Parameters
    ----------
    im: input image

    Return
    -------
    binary: output image
    """
    masked = (4.0 / 255.0) * _spm.imfilter(im, 'edge_enhance')
    masked = _np.floor(masked).astype('uint8')
    masked[masked > 3] = 3
    binary = _compiled.halftone.ordered_comb2_iterator(masked)
    return binary


def ordered_combinatorial4(im):
    """
    Encode the image using combinatorial dithering for halftone with 2^4 levels

    Parameters
    ----------
    im: input image

    Return
    -------
    binary: output image
    """
    masked = (16.0 / 255.0) * _spm.imfilter(im, 'edge_enhance')
    masked = _np.floor(masked).astype('uint8')
    masked[masked > 15] = 15
    binary = _compiled.halftone.ordered_comb4_iterator(masked)
    return binary


def ordered_dithering_generalized_bayer(grayscale, pattern_size=3):
    """ Modified Ordered Dithering Algorithm.

    dithered = ordered_dithering_generalized_bayer(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale
    pattern_size: dot-pattern size (default = 3) (must be even)

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    P.Garcia, M.Farias, and A.Araujo, ``Fast Inverse Halftone for
    Ordered Dithering'' in SIBIGRAPI, 2011.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2011 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    bayesian_mask = __generate_bayesian_matrix(pattern_size)
    masks = __generate_ordered_dithering(bayesian_mask)
    binary = __generalized_bayer(grayscale, masks)
    return binary


def __generate_bayesian_matrix(n):
    """
    Return dot-pattern matrix of indicies for ordered dithering of 2^n order

    Parameters
    ----------
    * n: integer for order of mask (must be even)

    Return
    ----------
    * mask: bayesian matrix with dot patter order to generate ordered dithering

    Written by Pedro Garcia Freitas <sawp@sawp.com.br>
    Copyright 2011 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    mask = _np.array([[0, 2], [3, 1]])
    mask_order = 2
    while mask_order < n:
        u = _np.ones((mask_order, mask_order))
        m11 = 4 * mask + 2 * u
        m12 = 4 * mask
        m21 = 4 * mask + u
        m22 = mask + u
        top = _np.concatenate((m11, m12), axis=1)
        down = _np.concatenate((m21, m22), axis=1)
        mask = _np.concatenate((top, down), axis=0)
        mask_order = 2 * mask_order
    return mask


def __generate_ordered_dithering(bayesian_mask):
    """
    Return a set of ordered dithering dot patterns

    Parameters
    ----------
    * bayesian_mask: numpy matrix with mask order

    Return
    ----------
    * dot_pattern: numpy matrix all dot patterns

    Written by Pedro Garcia Freitas <sawp@sawp.com.br>
    Copyright 2011 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (m, n) = bayesian_mask.shape
    total = m * n - 1
    dot_pattern = [_np.zeros((m, m))]
    for i in range(total):
        last = dot_pattern[i]
        new = (bayesian_mask == i)
        dot_pattern.append(last + new)
    dot_pattern = _np.array(dot_pattern)
    return dot_pattern


def __generalized_bayer(image, masks):
    """
    Convert a grayscale image in binary halftone image with n levels.

    Parameters
    ----------
    * im: numpy matrix, original grayscale image
    * masks: numpy array, set with all dot patterns used in halftone

    Return
    ----------
    * binary: numpy matrix, dithered binary image

    Written by Pedro Garcia Freitas <sawp@sawp.com.br>
    Copyright 2011 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    # create and rescale new image to fit levels
    (levels, m, n) = masks.shape
    (h, l) = image.shape
    (step_x, step_y) = masks[0].shape
    masked = (levels / 255.0) * _spm.imfilter(image, 'edge_enhance')
    masked = _np.floor(masked).astype('uint8')
    masked[masked >= levels] = levels - 1
    binary = _np.zeros((m * h, n * l))

    # generate the halftoned image_path
    k = 0
    r = 0
    for i in range(h):
        for j in range(l):
            mask = int(masked[i, j])
            selected = masks[mask]
            xs = i + k
            xf = i + k + step_x
            ys = j + r
            yf = j + r + step_y
            binary[xs:xf, ys:yf] = selected[:, :]
            r = r + step_x - 1
        r = 0
        k = k + step_y - 1
    return binary


def error_diffusion_jarvis(grayscale):
    """ Original Jarvis' Dithering Algorithm.

    dithered = error_diffusion_jarvis(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    J.Jarvis, C.Judice, and W.Ninke, ``A survey of techniques for the
    display of continuous tone pictures on bilevel displays,''
    em Comp. Graph and Image  Proc., vol.~5, pp.~13--40, 1976.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2010 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.jarvis_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_stucki(grayscale):
    """ Original Stucki's Dithering Algorithm.

    dithered = error_diffusion_stucki(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    Stucki, P., ``MECCA - a multiple-error correcting computation algorithm
    for bilevel image hardcopy reproduction.''  Research Report RZ1060, IBM
    Research Laboratory, Zurich, Switzerland, 1981.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.stucki_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_shiaufan(grayscale):
    """ Original Shiau-Fan's Dithering Algorithm.

    dithered = error_diffusion_shiaufan(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    Shiau, Jeng-nan, and Zhigang Fan. ``Set of easily implementable 
    coefficients in error diffusion with reduced worm artifacts.''
    Electronic Imaging: Science & Technology. International Society
    for Optics and Photonics, 1996.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.shiaufan_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_atkinson(grayscale):
    """ Original Atkinson's Dithering Algorithm.

    dithered = error_diffusion_atkinson(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.atkinson_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_sierra_simple(grayscale):
    """ Original Frankie Sierra's Dithering Algorithm.

    dithered = error_diffusion_sierra_simple(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    Frankie Sierra is unpublished, but can be reached via CIS at UID#
    76356,2254.  Pictorial presentations of his filters can be found in LIB
    17 (Developer's Den) of the CIS Graphics Support Forum as the files
    DITER1.GIF, DITER2.GIF, DITER6.GIF, DITER7.GIF, DITER8.GIF, and
    DITER9.GIF.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 1))
    infer = threshold * _np.ones((1, N + 2))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.sierra24a_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_sierra3(grayscale):
    """ Original Sierra's Algorithm.

    dithered = error_diffusion_sierra3(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    Frankie Sierra is unpublished, but can be reached via CIS at UID#
    76356,2254.  Pictorial presentations of his filters can be found in LIB
    17 (Developer's Den) of the CIS Graphics Support Forum as the files
    DITER1.GIF, DITER2.GIF, DITER6.GIF, DITER7.GIF, DITER8.GIF, and
    DITER9.GIF.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.sierra3_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_sierra2(grayscale):
    """ Original Sierra's Dithering Algorithm.

    dithered = error_diffusion_sierra2(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    Frankie Sierra is unpublished, but can be reached via CIS at UID#
    76356,2254.  Pictorial presentations of his filters can be found in LIB
    17 (Developer's Den) of the CIS Graphics Support Forum as the files
    DITER1.GIF, DITER2.GIF, DITER6.GIF, DITER7.GIF, DITER8.GIF, and
    DITER9.GIF.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.sierra2_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_burkes(grayscale):
    """ Original Burke's Dithering Algorithm.

    dithered = error_diffusion_burkes(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    The Burkes error filter was submitted to the public domain on
    September 15, 1988 in an unpublished document,
    ``Presentation of the Burkes error filter for use in preparing
    continuous-tone images for presentation on bi-level devices.''
    The file BURKES.ARC, in LIB 15 (Publications) of the CIS Graphics
    Support Forum, contains this document as well as sample images.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 2))
    infer = threshold * _np.ones((2, N + 4))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.burkes_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered


def error_diffusion_floyd_steinberg(grayscale):
    """ Original Floyd-Steinberg Dithering Algorithm.

    dithered = error_diffusion_floyd_steinberg(grayscale)

    Parameters
    ----------
    grayscale: numpy ndarray in grayscale

    Return
    ----------
    dithered: numpy ndarray, dithered image

    References
    ----------
    Floyd, R.W. and L. Steinberg, ``An Adaptive Algorithm for Spatial Gray
    Scale.''  SID 1975, International Symposium Digest of Technical Papers,
    vol 1975m, pp. 36-37.

    Written by Pedro Garcia Freitas [sawp@sawp.com.br]
    Copyright 2014 by Pedro Garcia Freitas

    see: http://www.sawp.com.br
    """
    (M, N) = grayscale.shape
    threshold = 127.5
    supers = threshold * _np.ones((M, 1))
    infer = threshold * _np.ones((1, N + 2))
    tmp = _np.concatenate((supers, grayscale, supers), axis=1)
    tmp = _np.concatenate((tmp, infer))
    dithered = _compiled.halftone.floyd_steinberg_iterator(tmp)
    dithered = dithered.astype('uint8')
    return dithered
