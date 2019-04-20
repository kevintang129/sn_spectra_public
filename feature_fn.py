# current
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#helper functions
import helper as hp
import os
# from linear_fit import linear_fit
from scipy.signal import argrelextrema, savgol_filter
from sklearn import linear_model


class linear_fit:
# class that takes in x,y values and is able to find the wavelength of the feature 
# through fitting and subtracting the continuum from the spectrum. Then a polynomial
# fitting is done of the remaining feature in order to find the minimum.

# x, y - input data to be plotted
# mini - index of minimum value found
# left_edge - left bound of continuum linear fit
# right edge - right bound of continuum linear fit
	def __init__(self, x, y, mini, left_edge, right_edge):
		self.x = x
		self.y = y
		self.min_idx = mini
		self.left_edge = left_edge
		self.right_edge = right_edge
		self.linear_edge = int(args.linear_region // (x[1]-x[0]))

	def fit_continuum(self):
	#Fits a line across the two bounds of the feature in order to find the continuum. Returns
	#the line's y values, slope, and y-intercept.

		x_fit = np.array(x[lin_high-self.linear_edge:lin_high+self.linear_edge] + 
			x[lin_low-self.linear_edge:lin_low+self.linear_edge])
		y_fit = np.array(y[lin_high-self.linear_edge:lin_high+self.linear_edge] + 
			y[lin_low-self.linear_edge:lin_low+self.linear_edge])

		reg = linear_model.LinearRegression()
		reg.fit(x_fit.reshape(-1,1),y_fit.reshape(-1,1))
		self.y_pred = reg.predict(np.array(x).reshape(-1,1))
		self.m = reg.coef_[0]
		self.b = reg.intercept_
		return (self.y_pred, self.m, self.b)


	def find_feature(self):
	#Subtracts the continuum, then takes the remaining feature and applies a polynomial regression 
	#to find the minimum wavelength, or the wavelength of the feature. Returns the subtracted spectrum,
	#the range of wavelengths the feature is fit over, and the fitted values.
		y_sub = [y - (self.m*x+self.b) for x,y in zip(x, y)]

		if lin_low > lin_high: #sees which bound is the higher bound actually
			low = lin_high-1
			high = lin_low+1
		else:
			low = lin_low+1
			high = lin_high-1
		# print(x[high])
		self.poly_edge = int(args.feature_region * min(self.min_idx-low, high-self.min_idx))

		# symmetrizes the polynomial fitting if points are skewed
		right_wv_diff = x[high]-x[min_idx]
		left_wv_diff =x[min_idx]-x[low]
		if right_wv_diff - left_wv_diff >= 50.0:
			high = 2*min_idx - low
		if left_wv_diff - right_wv_diff >= 50.0:
			low = 2*min_idx - high

		z = np.polyfit(np.squeeze(x[low+self.poly_edge:high-self.poly_edge]), 
			np.squeeze(y_sub[low+self.poly_edge:high-self.poly_edge]),args.order)
		fit = np.poly1d(z)
		data = np.linspace(x[low],x[high],high-low-2*self.poly_edge)
		return (y_sub, data, fit, low, high)

# Code for accepting command line arguments
# filename - name of file 
# bin_factor - amount of binning on the data
# smooth - order of smoothing polynomial
# order - order of final polynomial fitting on feature
# linear_region - +- number of angstroms from linear fit bounds to form linear fitting region
# feature_region - subtracting how many angstroms from upper bounds for final polynomial fit region
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-f", '--filename', type=str, required=True)
parser.add_argument("-r", '--redshift',type=float, default=0)
parser.add_argument("-b", '--bin_low', type=int, default=1)
parser.add_argument("-d", '--bin_high', type = int, default=5)
parser.add_argument("-s", '--smooth', type=int, default=3)
parser.add_argument("-o", '--order', type=int, default=5)
parser.add_argument("-l", '--linear_region', type=int, default = 20)
parser.add_argument("-p", '--feature_region', type=float, default = 0.3)
args = parser.parse_args()

sn = args.filename.split('-')
data = "data/"+sn[0]+'/'+args.filename
db = pd.read_csv(data, delim_whitespace = True, names = ['wavelength', 'values'])
average_v = 0
for b in range(args.bin_low, args.bin_high+1):
	m = db.values
	# print(m)
	x = np.delete(m,[1,2],axis=1) 
	y = np.delete(m,[0,2],axis=1) 
	x = x.ravel()
	y = y.ravel()
	x = [wv/(1+args.redshift) for wv in x]
	x = hp.rebin(x,b) #rebin to reduce noise
	y = hp.rebin(y,b)
	y_smooth = savgol_filter(y, 51, args.smooth) #smooth the curve
	
	#Finding the bounds for the linear fit. Starts at min wv in region, then scales up until 
	#finds a value that is greater. Compare with surounding values to check if actually a max or just noise.
	min_wv, min_val, min_idx = hp.find_min(x,y_smooth)
	# print(min_wv)
	step = x[1]-x[0]
	edge = int(10*b // step)
	j,k = min_idx, min_idx
	while j < len(y_smooth)-1:
		if y_smooth[j] > y_smooth[j+1]:
			if((y_smooth[j] > y_smooth[j-edge]) & (y_smooth[j] > y_smooth[j+edge])):
				break
		j+=1
	lin_high = j
	
	while k > 0:
		if y_smooth[k] > y_smooth[k+1]:
			if((y_smooth[k] > y_smooth[k-edge]) & (y_smooth[k] > y_smooth[k+edge])):
				break
		k-=1
	lin_low = k
	
	#Create a linear_fit class that can subtract the continuum and find the feature.
	line = linear_fit(x, y, min_idx, lin_low, lin_high)
	y_pred = line.fit_continuum()[0]
	sub_data, fit_wv_region, poly_fitted, low_bound_fit, high_bound_fit = line.find_feature()
	fitted_values = poly_fitted(fit_wv_region)

	for i,(a,b) in enumerate(zip(fit_wv_region, fitted_values)):
		if b < fitted_values[i-1] and b < fitted_values[i+1]:
			min_wv = a
		if i == len(fitted_values)-2:
			break
	
	min_diff = 999999
	k=0
	# for w in x[low_bound_fit:high_bound_fit]:

	for w in x[low_bound_fit:high_bound_fit]:
		if abs(min_wv-w) < min_diff:
			min_diff = abs(min_wv-w)
			min_w = w
			min_wi = k
		k+=1

	err_point = 0
	err_1 = 0
	err_2 = 0
	for c1, wv1 in enumerate(x[min_wi+low_bound_fit:high_bound_fit]):
		# print(wv1)
		# print(y_comp[min_wi+c1], poly_fitted(wv1))
		if sub_data[low_bound_fit+min_wi+c1] > poly_fitted(wv1):
			if err_point == 3:
				# print(x[min_wi:high_bound_fit][min_wi+c1])
				err_1 = wv1-min_wv
				break
			else:
				# print(x[min_wi:high_bound_fit][min_wi+c1])
				err_point += 1 
	
	err_point = 0
	lower_end = x[low_bound_fit:min_wi+low_bound_fit]
	lower_end.reverse()
	for c2, wv2 in enumerate(lower_end):
		# print(wv2, y_comp[min_wi-c2], poly_fitted(wv2))
		if sub_data[low_bound_fit+min_wi-c2] > poly_fitted(wv2):
			# print(wv2)
			# print(lower_end[c2])
			if err_point == 3:
				# print(x[min_wi:high_bound_fit][min_wi+c2])
				err_2 = min_wv-wv2
				break
			else:
				# print(x[min_wi:high_bound_fit][min_wi+c2])
				err_point += 1 
	error = max(err_1, err_2)
	vel = hp.find_vel(min_wv,6355)
	print('%f +- %f' % (min_wv,error))
	print('%f' % vel + " km/s")
	average_v += vel
	
	#plotting the results
	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()
	ax1.plot(x,y)
	ax1.plot(x, y_smooth)
	ax1.plot(x, y_pred)

	ax2.plot(fit_wv_region, fitted_values)
	ax2.scatter(x[low_bound_fit+line.poly_edge:high_bound_fit-line.poly_edge], 
		sub_data[low_bound_fit+line.poly_edge:high_bound_fit-line.poly_edge])
	
	ax2.scatter(x[low_bound_fit:low_bound_fit+line.poly_edge], 
		sub_data[low_bound_fit:low_bound_fit+line.poly_edge], c = "r")
	
	ax2.scatter(x[high_bound_fit-line.poly_edge:high_bound_fit], 
		sub_data[high_bound_fit-line.poly_edge:high_bound_fit], c = 'g')
	
	# ax2.scatter(x[min_idx], line.find_feature()[0][min_idx], c = "r")
	ax2.scatter(x[low_bound_fit:high_bound_fit][min_wi],sub_data[low_bound_fit + min_wi], c = 'r')
	# ax1.set_title(args.filepath)
	dir_path = os.path.dirname(os.path.realpath(__file__))
	# if not os.path.exists(dir_path+'/plots/'):
	#     os.makedirs(directory+'/plots/')
	print(dir_path)
	# fig1.savefig(dir_path+'/plots/'+args.filename[:-4]+'_bin_%.0f.png' %b)
	# fig2.savefig(dir_path+'/data/'+sn[0]+'/'+spectra.name[:-4]+'_bin_'+str(b)+'_feature.png' %b)
	# fig1.savefig('spectra/plots/'+args.filename[:-3]+'_bin_%f.png' %b)
	# fig2.savefig('spectra/plots/'+args.filename[:-3]+'_bin_%f_feature.png' %b)
	plt.show()	
	plt.close('all')

	# date = spectra.name.split('-')[1].split('.')[0]
	# data_sheet = dir_path+'/data/'+args.filepath+'/'+spectra.name[:-4]+'.csv'
	# header = ["bin_factor", "velocity"]
	# if os.path.exists(data_sheet):
	# 	with open(data_sheet,'rb') as file1:
 # 	  		existingLines = [line for line in csv.reader(file1, delimiter=',')]
 # 	  else:
	# 	with open (data_sheet,'a') as filedata:                            
 # 	    		writer = csv.DictWriter(filedata, delimiter=',', fieldnames=header)
 # 	    		writer.writerow(data) 
average_v = average_v / 5