import os
import astropy.units as u
import numpy as np
import numpy.ma as ma
import math
import sys
import matplotlib.pyplot as plt
from photutils import centroid_sources
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from regions import CirclePixelRegion, PixCoord
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.coordinates import ICRS,  FK5, SkyCoord
from astroquery.gaia import Gaia
from matplotlib.colors import LogNorm
CRVAL1 =[]
CRVAL2 = []
CD1_1 =[]
CD1_2 =[]
CD2_1 =[]
CD2_2 =[]
VALS = []
Big_Fudge = 1.5 #Arbitrary constant to scale large, bright sources. 
small_fudge = .2 #Arbitrary constant to sclae smaller, dim sources.

def get_wcs(images): #extracts necessary WCS info from science images.
    L = int(len(images))
    weight = np.ones(L)
    weight[0] = 10
    for W in images:
        CRVAL1.append(W.header['CRVAL1'])
        CRVAL2.append(W.header['CRVAL2'])
        CD1_1.append(W.header['CD1_1'])
        CD1_2.append(W.header['CD1_2'])
        CD2_1.append(W.header['CD2_1'])
        CD2_2.append(W.header['CD2_2'])
    NAXIS = images[0].header['NAXIS']
    NAXIS1 = images[0].header['NAXIS1']
    NAXIS2 = images[0].header['NAXIS2']
    CRPIX1 = images[0].header['CRPIX1']
    CRPIX2 = images[0].header['CRPIX2']
    CTYPE1 = images[0].header['CTYPE1']
    CTYPE2 = images[0].header['CTYPE2']
    CUNIT1 = images[0].header['CUNIT1']
    CUNIT2 = images[0].header['CUNIT2'] 
    VAL1 = np.average(CRVAL1, weights=weight) #need to add weighting back into this but letting it go until can get things working again.
    VAL2 = np.average(CRVAL2, weights=weight)
    VALS.append(VAL1)
    VALS.append(VAL2)
    CD11 = np.average(CD1_1, weights=weight)
    CD12 = np.average(CD1_2, weights=weight)
    CD21 = np.average(CD2_1, weights=weight)
    CD22 = np.average(CD2_2, weights=weight)
    total_wcs = WCS(naxis = int(NAXIS))
    total_wcs.wcs.crpix = [CRPIX1, CRPIX2]
    total_wcs.wcs.crval = [VAL1, VAL2]
    total_wcs.wcs.cunit = [CUNIT1, CUNIT2]
    total_wcs.wcs.ctype = [CTYPE1, CTYPE2]
    total_wcs.wcs.cd = [(CD11, CD12), (CD21, CD22)]
    total_wcs.array_shape = [NAXIS2, NAXIS1]
    return total_wcs
            
        
        
def reference(im): #runs GAIA cone search
   ra = VALS[0]
   dec = VALS[1]
   coord = SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'icrs')
   radius = u.Quantity(.3, u.deg) #.3 deg is good enough for normal sized images can be adjusted.
   Gaia.ROW_LIMIT = -1 # -1 returns all objects in search.
   j = Gaia.cone_search_async(coord, radius)
   r = j.get_results()
   print('done searching')
   return r
   


def mask(image, scis): # WILL REQUIRE TWO ARGS: image to be masked (i.e. template image) and list of science images passed from main.
    data = image
    w = get_wcs(scis)
    empty_mask = ma.zeros([data.shape[0]+100, data.shape[1]+100]) #adds padding to the image so that all sources in the frame can be masked
    #print(w)
    f = reference(image)
    for source in f:
        if  source['ra_error'] > .3 and  source['dec_error'] > .3: #filters out ra  &  dec error greater than .3 degrees
            continue
        print('STEP1')
        if source['phot_g_mean_mag'] >18 and source['phot_bp_mean_mag']  > 18 and source['phot_rp_mean_mag'] > 18:
            continue #filters  out sources with  RGB mag  greater than 18 -> sources dimmer than mag 18 (should be or?)
        print('STEP2')
        point = SkyCoord(ra = source['ra'],  dec= source['dec'], unit = (u.deg, u.deg))
        center = skycoord_to_pixel(point, w, origin =1, mode = 'all') #This and previous line convert GAIA output source coordinates to pixel value
        print(center)
        if (center[0]  <= 0) or (center[0] >= 3054):
            continue #makes sure sources are in frame
        print('STEP3')
        if (center[1] <= 0) or (center[1]  >= 2042): 
            continue
        print('STEP4')
        x_cent = float(center[0])
        y_cent = float(center[1])
        x, y = centroid_sources(data, x_cent, y_cent, box_size = 30)
        if math.isnan(x[0]) == True: #sometimes centroid sources returns unusable NAN values, this will reset them to GAIA centers if NAN.
             x[0] = center[0]
        if math.isnan(y[0]) == True:
             y[0] = center[1]
        bk_list = []
        i = 0
        while i != 100:
            x_back = np.random.random(100) * 60 #use random smaples instead
            y_back = np.random.random(100) * 60
            back = np.stack((x_back, y_back), axis = -1)
            avlist = []
            for pnt in back:
                if (x[0]+pnt[0]-30 <= 0) or (x[0]+pnt[0]-30 >= data.shape[1]): #making  sure to  exclude points outside the frame. Come up  with better way for   background sampling? Photutils????
                    continue
                if (y[0]+pnt[1]-30 <= 0) or (y[0]+pnt[1]-30 >= data.shape[0]):
                    continue
                avlist.append(data[int((y[0]-30)+pnt[1]), int((x[0]-30)+pnt[0])])
            bk_list.append(np.sum(avlist)/(int(len(avlist))))
            i +=1
        AV = np.sum(bk_list)/(int(len(bk_list)))
        C = 1
        if data[int(y[0]), int(x[0])] > 10000:
            C = Big_Fudge
        if data[int(y[0]), int(x[0])] - AV  <  20:
            C = small_fudge
        print('Step5')
        Rs = 4    
        R1 = 0
        while data[int(y[0]), int(x[0]+R1)] > AV:
            R1 += .1
            if (x[0]+R1) > data.shape[1]:
                R1 =0
                Rs -= 1
                break
            if R1 >= 17.5:
                break
        R2 = 0
        while data[int(y[0]+R2), int(x[0])] > AV:
            R2 += .1
            if (y[0]+R2) > data.shape[0]:
                R2 =0
                Rs -= 1
                break
            if R2 >= 17.5:
                break
        R3 = 0
        while data[int(y[0]), int(x[0]-R3)] > AV:
            R3 += .1
            if (x[0]-R3) < 0:
                R3 =0
                Rs -= 1
                break
            if R3 >= 17.5:
                break
        R4 = 0
        while data[int(y[0]-R4), int(x[0])] > AV:
            R4 += .1
            if (y[0]-R4) > 0:
                R4 =0
                Rs -= 1
                break
            if R4 >= 17.5:
                break
        R = float(((R1+R2+R3+R4)/Rs)*C)
        #print(R)
        center_2 = PixCoord(x[0]+50, y[0]+50)
        circle = CirclePixelRegion(center_2, R)
        mask = circle.to_mask()
        #data[mask.bbox.slices] *= 1 - mask.data
        empty_mask[mask.bbox.slices] += mask
    print('Step6')
    dlt = np.arange(0, 50, 1)
    np.delete(empty_mask, [i for i in dlt], axis= 0)
    np.delete(empty_mask, [i+data.shape[0] for i in dlt], axis= 0)
    np.delete(empty_mask, [i for i in dlt], axis= 1)
    np.delete(empty_mask, [i+data.shape[1] for i in dlt], axis= 1)
    #print(empty_mask.data)
    print('Step7')
    plt.imshow(empty_mask.data, norm=LogNorm(vmin=10, vmax=100))
    lg = np.log(empty_mask.data+.001)
    plt.imsave(fname='new_mask.png', arr=lg, cmap='viridis', dpi = 300) 
    #plt.show()
    #temp_im = fits.PrimaryHDU(data).writeto('/home/luca/Desktop/mask_test')
    #return data
                                
#mask(img)                         
