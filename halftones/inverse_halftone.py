# -*- coding: utf-8 -*-
#
# inverse_halftone.py - Modules for inverse halftone.
#
#
# Autor: Pedro Garcia Freitas [sawp@sawp.com.br]
# License: GPL V.2
#
# see: http://www.sawp.com.br
#
# Sep 2014
#
import numpy as _np
import scipy.misc as _spm
from skimage.restoration import denoise_tv_chambolle as _tv
import halftones.compiled as _compiled


def inverse_halftone_ordered_combinatorial3(binary):
    '''
    Decode the halftone combinatorial dithering to a graylevel image.

    Parameters
    ----------
    binary: input binary image

    Return
    -------
    gray: output grayscale image
    '''
    binary = binary.astype('bool')
    invh = _compiled.inverse_halftone.inverse_ordered_comb3_iterator
    gray = invh(binary)
    grayscale = 255.0 * (gray - gray.min()) / (gray.max() - gray.min())
    return grayscale.astype('uint8')


def inverse_halftone_ordered_combinatorial2(binary):
    '''
    Decode the halftone combinatorial dithering to a graylevel image.

    Parameters
    ----------
    binary: input binary image

    Return
    -------
    gray: output grayscale image
    '''
    binary = binary.astype('bool')
    invh = _compiled.inverse_halftone.inverse_ordered_comb2_iterator
    gray = invh(binary)
    grayscale = 255.0 * (gray - gray.min()) / (gray.max() - gray.min())
    return grayscale.astype('uint8')


def inverse_halftone_ordered_combinatorial4(binary):
    '''
    Decode the halftone combinatorial dithering to a graylevel image.

    Parameters
    ----------
    binary: input binary image

    Return
    -------
    gray: output grayscale image
    '''
    def what_level(lvl):
        level = __generate_combinated_levels4()
        l = 0
        lev = lvl.tolist()
        for i in level:
            if lev == i:
                break
            l += 1
        return l

    def get_color_from_level(lvl):
        intervalsize = 255.0 / 16.0
        base = lvl * intervalsize
        top = base + intervalsize
        color = base + (top - base) * _np.random.rand()
        return color

    (h, l) = binary.shape
    gray = _np.zeros((h, l / 4))
    for line in range(h):
        for col in range(0, l, 4):
            bits = binary[line, col:(col + 4)]
            lvl = what_level(bits)
            color = get_color_from_level(lvl)
            gray[line, col / 4] = color
    return gray


def __generate_combinated_levels4():
    level = [[0, 0, 0, 0]]
    level += [[0, 0, 0, 1]]
    level += [[0, 0, 1, 0]]
    level += [[0, 0, 1, 1]]
    level += [[0, 1, 0, 0]]
    level += [[0, 1, 0, 1]]
    level += [[0, 1, 1, 0]]
    level += [[0, 1, 1, 1]]
    level += [[1, 0, 0, 0]]
    level += [[1, 0, 0, 1]]
    level += [[1, 0, 1, 0]]
    level += [[1, 0, 1, 1]]
    level += [[1, 1, 0, 0]]
    level += [[1, 1, 0, 1]]
    level += [[1, 1, 1, 0]]
    level += [[1, 1, 1, 1]]
    return level


def __normalize(im):
    oldmax, oldmin = im.max(), im.min()
    oldrange = oldmax - oldmin
    newmin, newmax = 0.0, 255.0
    newrange = newmax - newmin
    scaled = (im - oldmin) / oldrange
    normalized = newrange * scaled + newmin
    return normalized


def __denoise(im):
    ff = _np.fft.fft2(im)
    keep_fraction = 0.22
    r, c = ff.shape
    ff[r*keep_fraction:r*(1-keep_fraction)] = 0
    ff[:, c*keep_fraction:c*(1-keep_fraction)] = 0
    restored = _np.fft.ifft2(ff).real
    return __normalize(restored)


def inverse_ordered_dithering_generalized(binary, masksize=3):
    binary = binary.astype('bool')
    invh = _compiled.inverse_halftone.inverse_ordered_dithering_iterator
    gray = invh(binary, masksize)
    grayscale = 255.0 * (gray - gray.min()) / (gray.max() - gray.min())
    denoised = __denoise(grayscale)
    smoothed = _tv(denoised, weight=10)
    return smoothed.astype('uint8')


def inverse_fbih(binary):
    binary = binary.astype('bool')
    invh = _compiled.fbih.inverse_halftone
    gray = invh(binary)
    grayscale = 255.0 * (gray - gray.min()) / (gray.max() - gray.min())
    return grayscale.astype('uint8')
