a={}
f2 = open("all.txt", 'r')## takes the obo parsed file . format is doid "\t" disease name or umls ids
for line in f2:
	y=line.strip().split("\t")
##creating dictionary with first column as key and next column as values
	if y[0] in a:
		a[y[0]].append(y[1])
	else:
		a[y[0]] = [y[1]]
#print (a) ## check whether dict formation is correct
for key, values in a.items(): 
	print (key+"\t"+"\t".join(values)) ## print tab separated values. here we have multiple values for one key.
