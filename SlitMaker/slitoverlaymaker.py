import matplotlib.pyplot as plt 
import numpy as np 
from astropy.io import fits 
import utils

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
    (xcen, ycen) = (deltax/2, deltax/2) #in your case, this should be the center of the section of the sky you are looking at is in delta_ra/dec(to simplify, you could just set the middle to be (0,0), but I don't know what the (0,0) of your examples was supposed to be)

    x,y = () #in this slot should be the delta_ra and delta_dec (respectively) of whatever the target point is.

    xadjust = x - xcen
    yadjust = y - ycen #this tells the offset from the center of the array of the target

    xcen += xadjust
    ycen += yadjust #shifts the relative center of the image to the target (at least in terms of coordinates)


    apcircle = plt.Circle((xcen, ycen), radius, color='m', lw=3, fill=False)
    ax1.add_artist(apcircle) #this just draws a circle around the target of the radius defined at the top of this function. Modify as needed to suit your preferences.
    
    #plot the slit

    def getslitcoords(slit_length, slit_PA, slit_width, len_pos, wid_off):
        rect_dx = 0.5 * slit_width * np.array([1., 1., -1., -1., 1.]) -wid_off #width is in arcsecs already
        rect_dy = slit_length * keckscale * np.array([0, 1, 1, 0, 0]) -(len_pos*keckscale) #into arcsecs
        rect = np.stack([rect_dx, rect_dy], axis=0)
        rect /= 3600 # convert arcsec to degree #this provides a vertical rectangle shifted according to the length- and width-wise offset of the target
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

def slit_center(slit_PA, target_x, target_y, slit_length, R1, R2):
	#find theta in terms of the PA; deg
	if (slit_PA >= 0) & (slit_PA <= 90):
		theta = slit_PA
		xsign = 1
		ysign = 1
	elif (slit_PA > 90) & (slit_PA <=180):
		theta = 180 - slit_PA
		xsign = -1
		ysign = 1
	elif (slit_PA > 180) & (slit_PA < 270):
		theta = slit_PA - 180
		xsign = -1
		ysign = -1
	elif (slit_PA == 270) | (slit_PA == -90):
		theta = 90
		xsign = -1
		ysign = -1
	elif (slit_PA > 270) & (slit_PA < 360):
		theta = 360 - slit_PA
		xsign = -1
		ysign = 1
	elif slit_PA == -270:
		theta = 90
		xsign = 1
		ysign = 1
	elif slit_PA == -180:
		theta = 0
		xsign = -1
		ysign = -1
	elif (slit_PA < -270) & (slit_PA > -360):
		theta = 360 + slit_PA
		xsign = 1
		ysign = 1
	elif (slit_PA < -180) & (slit_PA > -270):
		theta = abs(180 + slit_PA)
		xsign = 1
		ysign = -1
	elif (slit_PA < -90) & (slit_PA >-180):
		theta = 180 + slit_PA
		xsign = -1
		ysign = -1
	elif (slit_PA < 0) & (slit_PA > -90):
		theta = abs(slit_PA)
		xsign = -1
		ysign = 1

	#calculate the center of the slit
	xc = target_x + 0.5 * (slit_length - R1 - R2) * (0.1185 / 60**2) * xsign * abs(np.sin(np.deg2rad(theta))) 
	yc = target_y + 0.5 * (slit_length - R1 - R2) * (0.1185 / 60**2) * ysign * abs(np.cos(np.deg2rad(theta))) 

	return xc, yc

def make_plot(APID, slit_length, slit_PA, target_x, target_y, slit_width, R1, R2, mask, slit, radiusc, radiusir):
    
    #plot the HST images 
    file814 = r'C://Users//mkoga//Kogan Research//hststamps//{}_F814W.fits'.format(APID)
    hdul1 = fits.open(file814)
    imdata1 = hdul1[0].data
    file475 = r'C://Users//mkoga//Kogan Research//hststamps//{}_F475W.fits'.format(APID)
    hdul2 = fits.open(file475)
    imdata2 = hdul2[0].data
    img1 = imdata1+imdata2
    file110 = r'C://Users//mkoga//Kogan Research//hststamps//{}_F110W.fits'.format(APID)
    hdul3 = fits.open(file110)
    imdata3 = hdul3[0].data
    file160 = r'C://Users//mkoga//Kogan Research//hststamps//{}_F160W.fits'.format(APID)
    hdul4 = fits.open(file160)
    imdata4 = hdul4[0].data
    #set up the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    deltax = hdul1[0].header['NAXIS1']
    for ax in (ax1, ax2):
        ax.set_xticks([])
        ax.set_yticks([])
        ax2.set_xlim(0, deltax) 
        ax2.set_ylim(deltax, 0)
    fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    fig.suptitle(str(APID), y=1.06, fontsize=45)
    ax1.text(5, -5, 'F475W + F814W', fontsize=30)
    ax2.text(240, -5, 'F110W + F160W', fontsize=30)
    
    #finish plotting image
    img1 = imdata1+imdata2
    img2 = imdata3+imdata4
    hstscale = hdul1[0].header['D001SCAL'] #0.05 arcsec per pixel in hst
    keckscale = 0.1185 #arcsec per pixel in keck
    ax1.imshow(img1, extent=[0, deltax, deltax, 0])
    ax2.imshow(img2, extent=[0,deltax,deltax,0],cmap='gray')
	
    # overplotting the radius circle and the target position -- will add slit later
    (xcen, ycen) = (deltax/2, deltax/2)
    dx = xcen
    dy = ycen
    apcirclec = plt.Circle((ycen, ycen), radiusc, color='m', lw=3, fill=False)
    ax1.add_artist(apcirclec)
    apcircleir = plt.Circle((ycen, ycen), radiusir, color='m', lw=3, fill=False)
    ax2.add_artist(apcircleir)
    
    #plot the slit
    slit_x, slit_y = slit_center(slit_PA, target_x, target_y, slit_length, R1, R2) #in degrees
    
    def getslitcoords(slit_length, slit_x, slit_y, slit_PA, slit_width):
        rect_dx = 0.5 * slit_width * np.array([1., 1., -1., -1., 1.]) #width is already in arcsec
        rect_dy = 0.5 * slit_length * keckscale * np.array([1., -1., -1., 1., 1.]) #into arcsecs
        rect = np.stack([rect_dx, rect_dy], axis=0)
        rect /= 3600 # convert arcsec to degree
        cospa = np.cos(np.deg2rad(slit_PA))
        sinpa = np.sin(np.deg2rad(slit_PA))
        rot = np.array([[cospa, -sinpa],
                        [sinpa, cospa]])  # rotates sky into mask-aligned
        rect_dxi, rect_deta = np.dot(rot.T, rect)
        rect_ra, rect_dec = utils.xieta2radec(rect_dxi, rect_deta, slit_x, slit_y) #make sure degs is ok
        return rect_ra, rect_dec

    #slit boundaries positions
    slitx = getslitcoords(slit_length, slit_x, slit_y, slit_PA, slit_width)[0]
    slity = getslitcoords(slit_length, slit_x, slit_y, slit_PA, slit_width)[1] #both in degrees
    # print(f' the slits are at {slitx} and {slity}')
    # print(f' target is at ({target_x},{target_y})')
    dx_slit = (target_x - slitx)  * 3600.0 / hstscale * 0.75 + xcen #3600 is to convert to arcsec, scale is to convert to pixels; 0.75 makes it not crooked
    dy_slit = (target_y - slity)  * 3600.0 / hstscale  + ycen
    # print(f' the new coords are at {dx_slit} and {dy_slit}')
    #plotting
    ax1.plot(dx_slit, dy_slit, lw=3, c='green')
    ax2.plot(dx_slit, dy_slit, lw=3, c='green')
    # plt.savefig('/Volumes/Titan/clusters/plots/overlays/{}_{}_{}.png'.format(mask, slit, ID), bbox_inches='tight', dpi=200)
    plt.close()
    
def make_plot_single(imageName, slit_length, slit_PA, target_ra, target_dec, slit_width, R1, R2, radiusc):
    
    # plot hst images
    hdul1 = fits.open(imageName)
    imdata1 = hdul1[0].data
    
    #set up the plot
    fig, ax1 = plt.subplots(1)
    deltax = hdul1[0].header['NAXIS1']
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlim(0, deltax)
    ax1.set_ylim(deltax, 0)
    fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    plt.title(str(imageName), y=1.06, fontsize=45)
    
    #finish plotting image
    img1 = imdata1
    hstscale = 0.05 #0.05 arcsec per pixel in hst
    keckscale = 0.1185 #arcsec per pixel in keck
    ax1.imshow(img1, extent=[0, deltax, deltax, 0])
    
    # overplotting the radius circle and the target position -- will add slit later
    (xcen, ycen) = (deltax/2, deltax/2)
    dx = xcen
    dy = ycen
    apcircle = plt.Circle((ycen, ycen), radiusc, color='m', lw=3, fill=False)
    ax1.add_artist(apcircle)
    
    #plot the slit
    slit_x, slit_y = slit_center(slit_PA, target_x, target_y, slit_length, R1, R2) #in degrees
    
    def getslitcoords(slit_length, slit_x, slit_y, slit_PA, slit_width):
        rect_dx = 0.5 * slit_width * np.array([1., 1., -1., -1., 1.]) #width is in arcsecs already
        rect_dy = 0.5 * slit_length * keckscale * np.array([1., -1., -1., 1., 1.]) #into arcsecs
        rect = np.stack([rect_dx, rect_dy], axis=0)
        rect /= 3600 # convert arcsec to degree
        # print(f'the unrotated rectangle, in pixels, is at {rect*3600/hstscale}')
        cospa = np.cos(np.deg2rad(slit_PA))
        sinpa = np.sin(np.deg2rad(slit_PA))
        rot = np.array([[cospa, -sinpa],
                        [sinpa, cospa]])  # rotates sy into mask-aligned
        rect_dxi, rect_deta = np.dot(rot.T, rect)
        rect_ra, rect_dec = utils.xieta2radec(rect_dxi, rect_deta, slit_x, slit_y) #make sure degs is ok
        return rect_ra, rect_dec

    #slit boundaries positions
    slitx = getslitcoords(slit_length, slit_x, slit_y, slit_PA, slit_width)[0]
    slity = getslitcoords(slit_length, slit_x, slit_y, slit_PA, slit_width)[1] #both in degrees
    # print(f' the slits are at {slitx} and {slity}')
    # print(f' target is at ({target_x},{target_y})')
    dx_slit = (target_x - slitx)  * 3600.0 / hstscale * 0.75 + xcen #3600 is to convert to arcsec, scale is to convert to pixels; 0.75 makes it not crooked
    dy_slit = (target_y - slity)  * 3600.0 / hstscale  + ycen
    # print(f' the new coords are at {dx_slit} and {dy_slit}')
    #plotting
    ax1.plot(dx_slit, dy_slit, lw=3, c='green')
    plt.show()
    # plt.savefig('/Volumes/Titan/clusters/plots/overlays/{}_{}_{}.png'.format(mask, slit, ID), bbox_inches='tight', dpi=200)
    plt.close()