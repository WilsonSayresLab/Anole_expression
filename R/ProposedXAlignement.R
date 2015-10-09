####################################################################################################
# Make a line chart illustrating the alignment bewtween proposed X-chromosome and GG15
####################################################################################################

#----------
#1. X starting positions by scaffold (in Mb): 
#----------

#Divide the starting loci for the anole and chicken by 1,000,000 to obtain loci in Mb. Add the following to the anole locus start for each scaffold to obtain linear starting points. Each is a rounded value of the scaffold length added to the previous value (obtained from UCSC):

Scaffold	Start	Color 		Length

GL343282.1	0	green4		1.8
GL343364.1	1.8	blue3		1.1
LGb		2.9	firebrick4	3.3
GL343550.1	6.2	blueviolet	.6
GL343423.1	6.8	tomato3		.9			
GL343913.1	7.7	darkgreen	.2
GL343947.1	7.9	navy		.2
GL343338.1	8.1	red3		1.3
GL343417.1	9.4	purple4		.9

# Set the following chicken chromosomes to fixed values ot separate them from GG 15:

GG 1	13.5
GG 2 	14
GG 19	14.5

#----------
#2. Copy columns in the following order in Excel (characters not in <> need to be copied to every row in the column):
#----------

segments(<propXstart>,<y1>,<gg_start_Mb>,<y2>,col=<color>)

#Copy columns into a text editor and replace tabs with spaces. The resulting text should resemble the following R script.

#----------
#3. Plot in R: 
#----------

# Add chromosome image, title ("Alignment between Proposed X-Chromosome and Chicken Chromosomes"), x-axis label ("Order on Anole Scaffolds"), and scaffold labels in GIMP.

{plot(16.1,2, type="n",xaxt='n',xlim=c(0,15.6),ylim=c(0,2),yaxt="n") 
segments(0.077515,0,13.5,2,col="green4")
segments(0.884957,0,14.5,2,"green4")
segments(0.033185,0,13.5,2,"green4")
segments(0.68496,0,7.877611,2,"green4")
segments(0.739104,0,7.866185,2,"green4")
segments(0.469263,0,6.476769,2,"green4")
segments(0.449853,0,6.458296,2,"green4")
segments(0.182112,0,6.301498,2,"green4")
segments(0.112724,0,6.259423,2,"green4")
segments(0.107011,0,6.25485,2,"green4")
segments(0.037793,0,6.212456,2,"green4")
segments(0.500824,0,5.935811,2,"green4")
segments(1.639515,0,0.77501,2,"green4")
segments(1.536457,0,0.708057,2,"green4")
segments(1.495441,0,0.66412,2,"green4")
segments(1.463449,0,0.653952,2,"green4")
segments(1.413479,0,0.625943,2,"green4")
segments(1.339479,0,0.546959,2,"green4")
segments(1.294534,0,0.510691,2,"green4")
segments(1.263283,0,0.458204,2,"green4")
segments(1.243229,0,0.443711,2,"green4")
segments(1.201636,0,0.399359,2,"green4")
segments(1.200189,0,0.398134,2,"green4")
segments(1.189878,0,0.391656,2,"green4")
segments(1.188351,0,0.390543,2,"green4")
segments(1.145954,0,0.353466,2,"green4")
segments(1.138838,0,0.341573,2,"green4")
segments(0.950208,0,0.105047,2,"green4")
segments(0.932623,0,0.083607,2,"green4")
segments(0.864633,0,0.032278,2,"green4")
segments(8.945563,0,12.584397,2,"red3")
segments(8.997252,0,12.572529,2,"red3")
segments(9.002622,0,12.568401,2,"red3")
segments(9.113275,0,12.452471,2,"red3")
segments(9.224526,0,12.30982,2,"red3")
segments(8.853895,0,11.191801,2,"red3")
segments(8.829223,0,11.185341,2,"red3")
segments(8.332812,0,10.698914,2,"red3")
segments(8.893128,0,9.904857,2,"red3")
segments(8.816818,0,9.892017,2,"red3")
segments(8.762592,0,9.820949,2,"red3")
segments(8.731916,0,9.794821,2,"red3")
segments(8.629368,0,9.728774,2,"red3")
segments(8.356131,0,9.652897,2,"red3")
segments(8.367507,0,9.640424,2,"red3")
segments(8.378336,0,9.566084,2,"red3")
segments(8.463576,0,9.50431,2,"red3")
segments(8.494128,0,9.484624,2,"red3")
segments(8.556686,0,9.441167,2,"red3")
segments(8.551988,0,9.437097,2,"red3")
segments(8.530653,0,9.417065,2,"red3")
segments(8.513258,0,9.394644,2,"red3")
segments(8.101146,0,9.388357,2,"red3")
segments(8.126845,0,9.33479,2,"red3")
segments(8.171804,0,9.293941,2,"red3")
segments(8.190264,0,9.288062,2,"red3")
segments(8.19754,0,9.257597,2,"red3")
segments(8.254659,0,9.206727,2,"red3")
segments(8.266738,0,9.196692,2,"red3")
segments(8.281181,0,9.18589,2,"red3")
segments(8.290379,0,9.178197,2,"red3")
segments(8.303364,0,9.150903,2,"red3")
segments(1.800465,0,3.179902,2,"blue3")
segments(1.809969,0,2.778011,2,"blue3")
segments(1.884687,0,2.711273,2,"blue3")
segments(1.913177,0,2.622866,2,"blue3")
segments(1.952557,0,2.61499,2,"blue3")
segments(1.963123,0,2.568608,2,"blue3")
segments(2.031255,0,2.342102,2,"blue3")
segments(2.602105,0,1.502227,2,"blue3")
segments(2.658394,0,1.436912,2,"blue3")
segments(2.699884,0,1.41727,2,"blue3")
segments(2.827192,0,1.300825,2,"blue3")
segments(2.807484,0,1.293466,2,"blue3")
segments(2.85839,0,1.279112,2,"blue3")
segments(9.498659,0,11.149676,2,"purple4")
segments(9.483832,0,11.142536,2,"purple4")
segments(9.470279,0,11.128801,2,"purple4")
segments(9.453873,0,11.107873,2,"purple4")
segments(9.436616,0,11.060992,2,"purple4")
segments(9.517453,0,10.827349,2,"purple4")
segments(9.529011,0,10.814684,2,"purple4")
segments(9.536904,0,10.80452,2,"purple4")
segments(9.549208,0,10.803673,2,"purple4")
segments(9.566392,0,10.794681,2,"purple4")
segments(9.589952,0,10.779334,2,"purple4")
segments(9.597291,0,10.775401,2,"purple4")
segments(9.603202,0,10.698915,2,"purple4")
segments(10.148458,0,10.110679,2,"purple4")
segments(10.160187,0,10.095236,2,"purple4")
segments(10.185376,0,10.077413,2,"purple4")
segments(10.203944,0,10.066049,2,"purple4")
segments(9.832562,0,10.007213,2,"purple4")
segments(9.784204,0,9.992497,2,"purple4")
segments(9.778577,0,9.991442,2,"purple4")
segments(9.75216,0,9.972759,2,"purple4")
segments(9.408636,0,8.128303,2,"purple4")
segments(7.173378,0,10.066049,2,"tomato3")
segments(7.612599,0,6.176354,2,"tomato3")
segments(7.555331,0,6.11706,2,"tomato3")
segments(7.465191,0,6.093727,2,"tomato3")
segments(7.312853,0,6.044223,2,"tomato3")
segments(7.266272,0,6.028493,2,"tomato3")
segments(7.193055,0,5.994099,2,"tomato3")
segments(7.133744,0,5.949729,2,"tomato3")
segments(7.108248,0,5.935811,2,"tomato3")
segments(7.010271,0,5.890702,2,"tomato3")
segments(6.988503,0,5.878062,2,"tomato3")
segments(6.875273,0,5.81044,2,"tomato3")
segments(6.582213,0,5.705935,2,"blueviolet")
segments(6.458329,0,5.634231,2,"blueviolet")
segments(6.428074,0,5.612621,2,"blueviolet")
segments(6.325296,0,5.574193,2,"blueviolet")
segments(6.244269,0,5.469305,2,"blueviolet")
segments(6.202036,0,5.446352,2,"blueviolet")
segments(6.238769,0,15.5,2,"blueviolet")
segments(7.704282,0,6.530372,2,"darkgreen")
segments(7.770283,0,6.520471,2,"darkgreen")
segments(7.790259,0,6.500745,2,"darkgreen")
segments(7.825073,0,6.496565,2,"darkgreen")
segments(7.83725,0,5.432308,2,"darkgreen")
segments(7.900001,0,8.849198,2,"navy")
segments(7.927833,0,8.58082,2,"navy")
segments(7.960391,0,8.481375,2,"navy")
segments(7.977406,0,8.434483,2,"navy")
segments(7.986484,0,8.378113,2,"navy")
segments(8.005099,0,8.326506,2,"navy")
segments(6.161512,0,12.207814,2,"firebrick4")
segments(5.732962,0,11.571879,2,"firebrick4")
segments(5.194739,0,11.407922,2,"firebrick4")
segments(5.239025,0,11.34425,2,"firebrick4")
segments(5.296371,0,11.323464,2,"firebrick4")
segments(5.313952,0,11.287023,2,"firebrick4")
segments(5.603148,0,11.205594,2,"firebrick4")
segments(5.618498,0,11.201187,2,"firebrick4")
segments(5.651633,0,9.892017,2,"firebrick4")
segments(5.661403,0,7.969262,2,"firebrick4")
segments(5.17913,0,5.425109,2,"firebrick4")
segments(5.166525,0,5.407225,2,"firebrick4")
segments(5.12729,0,5.368868,2,"firebrick4")
segments(5.055741,0,5.281004,2,"firebrick4")
segments(5.048701,0,5.275356,2,"firebrick4")
segments(5.036997,0,5.264075,2,"firebrick4")
segments(4.957694,0,5.226402,2,"firebrick4")
segments(4.937808,0,5.208507,2,"firebrick4")
segments(4.931204,0,5.201856,2,"firebrick4")
segments(4.886348,0,5.080304,2,"firebrick4")
segments(4.777079,0,5.05346,2,"firebrick4")
segments(4.769884,0,5.03837,2,"firebrick4")
segments(4.735711,0,5.002524,2,"firebrick4")
segments(4.699187,0,4.977536,2,"firebrick4")
segments(4.682913,0,4.956025,2,"firebrick4")
segments(4.671294,0,4.948665,2,"firebrick4")
segments(4.659453,0,4.941748,2,"firebrick4")
segments(4.650187,0,4.938113,2,"firebrick4")
segments(4.558933,0,4.855484,2,"firebrick4")
segments(4.547252,0,4.833491,2,"firebrick4")
segments(4.495171,0,4.734252,2,"firebrick4")
segments(4.318065,0,4.565314,2,"firebrick4")
segments(4.213349,0,4.490046,2,"firebrick4")
segments(4.198637,0,4.485385,2,"firebrick4")
segments(3.883503,0,4.260707,2,"firebrick4")
segments(3.470551,0,3.796426,2,"firebrick4")
segments(3.426496,0,3.750315,2,"firebrick4")
segments(3.189226,0,3.532964,2,"firebrick4")
segments(3.073812,0,3.411518,2,"firebrick4")
segments(5.069271,0,3.274476,2,"firebrick4")
segments(2.984719,0,3.219108,2,"firebrick4")
}
