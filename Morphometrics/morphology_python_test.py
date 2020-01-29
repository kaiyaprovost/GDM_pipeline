#import csv 
#import math 

file_string = "/Users/kprovost/Documents/Classes/Stephanie/Pipilo_fuscus_measurements_3october2017.csv"

#print("reading in file")

## reading in the file 
with open(file_string,"r") as filename:
	#wholefile = filename.read()
	linesfile = filename.readlines() ## usually better 
	
#print("read in file")
	
#filename = open(....)
#filename.close(....)

bill_length = [] ## designate empty list
## empty string = "" 
## empty tuple = ()
## empty dictionary = {}

## for every line in the file we read, do something
#for line in linesfile[1:6]: ## this does first five lines AFTER the header
for line in linesfile[1:]: ## this does all lines AFTER the header
	#print(line) ## prints whole line
	subset = line[0:20] ## gets first 20 characters
	#print(subset) ## prints the first 20 characters
	splitline = line.split(",") ## separate text by commas
	#print(splitline) 
	## bill length is column 8
	bill_text = splitline[8]
	print(bill_text)
	
	try:
		bill = float(bill_text)
		bill_length.append(bill)
	except:
		## do nothing 
		pass
	
	#bill = float(bill_text)
	#print(bill)
	#bill_length.append(bill)
	
print(bill_length)

def mean(list):
	## take the mean of a list 
	sumlist = sum(list)
	length = len(list)
	mean = sumlist/length 
	return(mean)
	#return(sum(list)/len(list))

mean_bill = mean(bill_length)
print(mean_bill)

#print("end split lines")