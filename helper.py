def rebin(vect, a):
# IDL rebin translated to python
# ; this does basically the same thing as IDL rebin, but you don't have
# ; to worry about the length of the old vector being an integer
# ; multiple of the length of the new vector: the end elements of the
# ; old array are just ignored
  
# ; the calling sequence is somewhat different in that you pass my rebin
# ; the old vector plus the number of pixels to average over

#a is the number you want to divide the length of the array by 
#e.g. myrebin (vect_length_8, 2) will return a vector of length 4
	npix = len(vect)
	i = 0
	newvect = []
	while i < npix:
		temp = vect[i]
		for j in range(i, i+a-1):
			if j >= npix:
				return newvect
			temp = temp+vect[j]
		newvect.append(temp/(1.0*a))
		i = i+a
	return newvect

def find_min(x,y):
#Finds the minimum wavelength in a specific region of the spectra. Finds starting point
	mini, maxi, min_wv, i = 9999999, 0, 0, 0
	print(x[-1],y[-1])
	for wv, val in zip(x,y): #finds the minimum wavelength within this range
		if ((wv >= 5800) & (wv <= 6600)):
			if val < mini:
				min_idx = i
				mini = val
				min_wv = wv
			elif wv > 6600:
				break
		i+=1
	return (min_wv, mini, min_idx)

#returns velocity in km/s using relativistic redshift
def find_vel(obs, lab):
	return 2.99792458e5 * (obs**2-lab**2)/(obs**2+lab**2)
