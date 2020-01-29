#!/usr/bin/env python3

import sys
import os
import glob

catlookup = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_BY_CATALOG.txt"

#convlist = glob.glob("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/*/*FULLGENOME*csv")
convlist = glob.glob("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/*/*S.txt")
#convlist = glob.glob("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/FULL_ENV/DISTANCES/*.txt")


with open(catlookup,"r") as catfile:
	catlines = catfile.readlines()

dicti = {}

for line in catlines:
	cat,lat,long = line.split("\t")
	latlong = lat+" "+long.strip()
	new_cats = dicti.get(latlong, "None")
	if new_cats == "None":
		new_cats = [cat]
	else:
		new_cats.append(cat)
	dicti[latlong] = new_cats

for key,val in dicti.items():
	#print(key)
	val2 = ",".join(val)
	dicti[key] = val2

for filetoconvert in convlist:
	
	print(filetoconvert)
	
	with open(filetoconvert,"r") as convfile:
		conv = convfile.read()
	
	#lines = conv.split("\n")
	
	for key,val in dicti.items():
		temp = conv.replace(key,val)
		conv = temp
	
	with open(filetoconvert+".converted","w") as outfile:
		outfile.write(conv)

