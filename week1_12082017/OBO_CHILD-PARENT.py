from __future__ import print_function, division
import json
from collections import defaultdict
fname = "doid-merged.obo"
term_head = "[Term]"
#Keep the desired object data here
all_objects = {}
def add_object(d):
	if "is_obsolete" in d:
		return

	#Gather desired data into a single list,
	# and store it in the main all_objects dict
	key = d["id"][0]
	#is_a = d["is_a"]
	#alt_id = d["alt_id"]
	namespace = d["name"]
	xref = d["xref"]
	#namespace = d["namespace"]
	#relationship = d["relationship"]
	#print (relationship)
	#is_obsolete = d["is_obsolete"]    
	#consider = d["consider"]
	#replaced_by = d ["replaced_by"]

	#Remove the next line if you want to keep the is_a description info
	#is_a = [s.partition(' ! ')[0] for s in is_a]
	#relationship = [p.partition (' ! ')[0] for p in relationship]
	#relationship1 = [q.partition (' ')[2] for q in relationship]
	#print (relationship1)
	#print (type(is_a))
	#all_objects[key] = d["name"] + d["namespace"] + is_a + relationship + is_obsolete + consider + replaced_by
	#all_objects[key] = is_obsolete + consider
	#all_objects[key] =  is_a + relationship
	all_objects[key] =  namespace + xref
#A temporary dict to hold object data
current = defaultdict(list)
with open(fname) as f:
	#Skip header data
	for line in f:
		if line.rstrip() == term_head:
			break
	for line in f:
		line = line.rstrip()
		if not line:
			#ignore blank lines
			continue
		if line == term_head:
			#end of term
			add_object(current)           
			current = defaultdict(list)
		else:
			#accumulate object data into current
			key, _, val = line.partition(": ")
			current[key].append(val)

if current:
	add_object(current)    
#print("\nall_objects =")
#print(all_objects)

#print(json.dumps(all_objects, indent = 4, sort_keys=True))
#fid=open("obo_final.txt",'w')
for i in all_objects.keys():
 #  print(i) 
	met=all_objects.get(i)
#   print(met)
#   print (i+"\t"+"\t".join(met))
	for k in range(len(met)):
		print (i + "\t" + (met[k]))

