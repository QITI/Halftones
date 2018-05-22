import halftones

from scipy.misc import *

gray = imread('lena1.jpg', True)
# halftones
jarvis = halftones.halftone.error_diffusion_jarvis(gray)
floyd_steinberg = halftones.halftone.error_diffusion_floyd_steinberg(gray)
stucki = halftones.halftone.error_diffusion_stucki(gray)
burkes = halftones.halftone.error_diffusion_burkes(gray)
sierra3 = halftones.halftone.error_diffusion_sierra3(gray)
sierra2 = halftones.halftone.error_diffusion_sierra2(gray)
sierra_simple = halftones.halftone.error_diffusion_sierra_simple(gray)
atkinson = halftones.halftone.error_diffusion_atkinson(gray)
shiaufan = halftones.halftone.error_diffusion_shiaufan(gray)
combinatorial3x3 = halftones.halftone.ordered_combinatorial3(gray)
combinatorial2x2 = halftones.halftone.ordered_combinatorial2(gray)
combinatorial4x4 = halftones.halftone.ordered_combinatorial4(gray)
OD = halftones.halftone.ordered_dithering_generalized_bayer(gray, 3)
diagonal_matrix = halftones.halftone.ordered_dithering_diagonal_ordered_matrix(gray)
clustered_dots = halftones.halftone.ordered_dithering_clustered_dots(gray)
central_white_points = halftones.halftone.ordered_dithering_central_white_point(gray)
balanced_centered_points = halftones.halftone.ordered_dithering_balanced_centered_point(gray)
dispersed_dots = halftones.halftone.ordered_dithering_dispersed_dots(gray)

# inverse halftones
inverseJarvis = halftones.inverse_halftone.inverse_fbih(jarvis)
inverseComb2x2 = halftones.inverse_halftone.inverse_halftone_ordered_combinatorial2(combinatorial2x2)
inverseComb3x3 = halftones.inverse_halftone.inverse_halftone_ordered_combinatorial3(combinatorial3x3)
inverseComb4x4 = halftones.inverse_halftone.inverse_halftone_ordered_combinatorial4(combinatorial4x4)
inverseOD = halftones.inverse_halftone.inverse_ordered_dithering_generalized(OD, 3)

# save some figures
imsave('halftone_ordered_dither_bayer.png', OD)
imsave('halftone_jarvis.png', jarvis)
imsave('halftone_stucki.png', stucki)
imsave('halftone_burkes.png', burkes)
imsave('halftone_sierra3.png', sierra3)
imsave('halftone_sierra2.png', sierra2)
imsave('halftone_sierra_simple.png', sierra_simple)
imsave('halftone_atkinson.png', atkinson)
imsave('halftone_shiaufan.png', shiaufan)
imsave('halftone_floyd_steinberg.png', floyd_steinberg)
imsave('halftone_ordered_dither_diagonal_matrix.png', diagonal_matrix)
imsave('halftone_ordered_dither_clustered_dots.png', clustered_dots)
imsave('halftone_ordered_dither_central_white_points.png', central_white_points)
imsave('halftone_ordered_dither_balanced_centered_points.png', balanced_centered_points)
imsave('halftone_ordered_dither_dispersed_dots.png', dispersed_dots)

imsave('inverse_ordered_dither.png', inverseOD)
imsave('inverse_jarvis.png', inverseJarvis)
