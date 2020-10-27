import sys
import os
import subprocess

#input file list

datalistname = str(sys.argv[3]);
idatalist = open(datalistname, "r")


idatafiles = [f.rstrip() for f in idatalist if ".root" in f ]
print(idatafiles)

numFiles = int(sys.argv[1]) 
ifile = 1

nthreads = sys.argv[2]

execpath = sys.argv[4]


if (numFiles == 0 ):
	numFiles = len(idatafiles)
print( "from runmacro numFiles == ", numFiles," nthreads == ",nthreads)
cmd = execpath+"/compiledThreads "+str(nthreads)+" "

for f in idatafiles:
	if (ifile <= numFiles):
		pathsplit = f.split("/")
		namesplit = pathsplit[-1].split(".")
		tag = namesplit[0]
	
		cmd = cmd + "\""+f+"\" "
	ifile = ifile + 1

print cmd
os.system(cmd)	
