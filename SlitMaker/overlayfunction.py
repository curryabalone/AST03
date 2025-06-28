import matplotlib.pyplot as plt 
import numpy as np 

# Here is where I get all the input variables from:
# name = string of whatever the object name is. Strictly speaking isn't needed for the actual overlay
# slit_length = the length in keck pixels (scale is 0.1185 arcsecs/pixel for DEIMOS) the slit lengthwise (i.e. parallel to the longest sides). Found in 1d header item SLITX1-SLITX0 
# slit_PA = rotation of the slit off the vertical axis (i.e. straight up is 0 degrees/radians). Found in 1d header item SLITPA
# target_x = ra in degrees of the target. Get it either from the hubble catalog or the 1d header item RA_OBJ, although it will need to be converted into degrees
# target_y = dec in degrees of the target. Get it either from the hubble catalog or the 1d header item DEC_OBJ, although it will need to be converted into degrees
# slit_width = Width of the slit in arcseconds. As far as I know for DEIMOS it is always 0.8, but maybe double check it. Unfortunately not found in the 1d header from what I can tell.
# length_pos = the distance in keck pixels (scale is 0.1185 arcsecs/pixel for DEIMOS) of the target from the bottom of the slit lengthwise. Found in 1d header by formula (hdul1d[b].header['OBJPOS'] + hdul1d[r].header['OBJPOS'])/2, where b and r are the indicies for the blue and red halves (in two separate headers in the file, use the Horne halves if possible)
# wid_offset = Offset in the widthwise direction (i.e. parallel to the short ends of the slit). Supposed to be zero. Until we figure out how to determine it empirically from the A-band you can safely make the assumption that this remains true (tl;dr wid_off is always 0 for now)

def overlay(slit_length, slit_PA, target_x, target_y, slit_width, length_pos, wid_offset):
    radius = 5
    keckscale = 0.1185 #arcsec per pixel in keck. Unfortunately this info doesn't seem to be in the header, but apparently is the scale the DEIMOS system universally operates on.
    
    
    # overplotting the radius circle and the target position -- will add slit later
    (xcen, ycen) = (target_x, target_y)

    def getslitcoords(slit_length, slit_PA, slit_width, len_pos, wid_off):
        rect_dx = 0.5 * slit_width * np.array([1., 1., -1., -1., 1.]) -wid_off #width is in arcsecs already
        rect_dy = slit_length * keckscale * np.array([0, 1, 1, 0, 0]) -(len_pos*keckscale) #into arcsecs
        rect = np.stack([rect_dx, rect_dy], axis=0)
        rect /= 3600 # convert arcsec to degree #this provides a vertical rectangle shifted according to the length- and width-wise offset of the target
        print(rect)
        cospa = np.cos(np.deg2rad(slit_PA))
        sinpa = np.sin(np.deg2rad(slit_PA))
        rot = np.array([[cospa, sinpa],
                        [-sinpa, cospa]])  # rotates sky into mask-aligned
        rect_dxi, rect_deta = np.dot(rot.T, rect)

        return rect_dxi, rect_deta

    #slit boundaries positions

    slitx,slity = getslitcoords(slit_length, slit_PA, slit_width,length_pos,wid_offset) #both in degrees

    dx_slit = slitx + xcen 
    dy_slit = slity + ycen #this should shift it according to the location of the target in the array

    return dx_slit, dy_slit