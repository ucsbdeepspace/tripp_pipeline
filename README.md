# sdi_pipeline
SDI Pipeline 2.0.0
The purpose of the pipeline is to detect transients. The aim is to make this user-friendly, widely applicable, and efficient. To that end we have a few steps:

1. Create a CAT file if it doesn't exist already: This is done using Source Extractor Python (SEP). This is a lightweight, python-native version of Sextractor. The parameters we get out of it are calculated in the same way SEXtractor calculates them. We are limited in the information we get out of SEP though. For example, we only get xy information and not RA DEC information, even though the information to do the conversion is available.

2. Align the images using astroalign: The images taken over an observation period will be positioned differently, resulting in different image frames. To be able to ensure that the stars that our pipeline finds are consistently the same stars, we align the images. This is done using astroalign, a lightweight python package that DOES NOT use WCS information to transform images. The image with the lowest SNR is selected as the reference image to which all other images are aligned. This is also used sometimes to stack images sometimes, though I do not know that we do that in our pipeline.
Makes a list of the brightest images
Create a 2d-tree and performs a nearest neighbor search
Then finds the 4 closest neighbors to the star (parameters used are scale and rotation and translationally invariant)
Find the best triangles from those methods
Test the linear transformation for the rest of the data and if it fits more than 80% of them then they use the model. 

3. Combine/stack the images and perform difference image analysis using Bramich or Adaptive Bramich: Combine works by taking the median value of the pixels in the aligned images. This combined image serves as the template image. We then subtract each aligned image from the template image to produce residuals. There are currently three options as to how to do this: np.subtract, a matrix element subtraction pixel by pixel, Bramich, and Adaptive Bramich. The latter two algorithms fall under the ois python package umbrella. The algorithm finds kernels using a least squares approach. Then it performs kernel subtraction rather than pixel to pixel, preserving features and accounting for background variations. This method lets us find the variable stars (difference image analysis)

4. Extract variable stars: Once variable stars are found, we use SEP on the Residual images to extract astrometry about the sources remaining. SEP can only extract information about the position of the star on the image. No other information useful to us is reasonably preserved in the residuals. We use these positions to identify the transient images.


What the work has been for the past year (pipeline restructure and why):

The way that the earlier pipeline was constructed, the data structures outputted were broken, making the identification of transients by cross-correlation very difficult. Along with this, certain scripts purposes were misunderstood, so some algorithms were being misused. This time last year, the pipeline began its renovation, and has only completed testing as of a month ago. The way the new pipeline structures and transfers its data is by saving the data to the fits files. This means that the output after the four steps of the pipeline consists of the following HDUs:

1. ALGN - The aligned pixel data, and we are working to ensure that the headers of these correspond to the reference image headers, and that this is a reasonable thing to do
2. SUB - The residuals pixel data, and the header is also the reference image header.
3. CAT - The catalog of all of the sources found by sep in the unprocessed images.
4. XRT - The catalog of the variable sources
5. REF - The catalog of stars that are both in our image and in the Gaia database.
This process of renovation, addition of additional features such as display, collate, ref etc. has taken us the whole year.
