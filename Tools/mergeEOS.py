import sys
from sys import argv
import os
sys.path.append('cfgs/')
import subprocess
import threading, Queue, time
verbose = False
users = {
	#"jan":["srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/jschulte/SingleMuon/"]
	#"jan":["srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/jschulte/"]
	"jan":["/eos/cms/store/group/phys_muon/jschulte/13TeV/MC"]
	#"jan":["gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/SingleMuon/"]

}

def __getCommandOutput2(command):
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        if int(err) == 256:
            print 'Path in %s does not exist!' % (command)
            return ' '
        else:
            raise RuntimeError, '%s failed w/ exit code %d' % (command, err)
    return data


def getPathList(path):

		
	command = 'eos ls ' + path

    	output = __getCommandOutput2(command).splitlines()
	result = []

	for name in output:
		#print name
		#if '_' in name:
		result.append(path+'/'+name)

	return result

def getFileList(path,result):

		
	command = 'eos ls ' + path

    	output = __getCommandOutput2(command).splitlines()
	print output
	if len(output) > 0:
		if output[0] == path:
			result = output
			return result
		
		for name in output:
			if not 'root' in name and not 'failed' in name and not "171217" in name:
				result = getFileList(path+'/'+name,result)
			elif not "failed" in name and not "DQM" in name and not "171217" in name:
				result.append(path+'/'+name)

#		if not "root" in output:
#			for newPath in output:
#				if not "failed" in newPath:
#					result = getFileList(path+"/"+newPath,result)			
#		else:
#			for index, file in enumerate(output):
#				output[index] = path+"/"+file
#			result = result + output
	return result



def main():
        import argparse
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    	parser.add_argument("--merge",dest="merge", action="store_true", default=False, help='merge expected limits')
    	parser.add_argument("--tag",dest="tag", default="", help='tag')
        args = parser.parse_args()
	tag = args.tag	
	global outDir
	

	
	path = users['jan'][0]+"/"+args.tag
	print "Searching for results in %s"%path

	paths = getPathList(path)
	print paths
	for path in paths:
		#path = paths[0]
		sample = path.split('/')[-1]
		files =  getFileList(path,[])
	
		if len(files) > 0:	
			print "merging %d files for %s"%(len(files),sample)
			command = ["hadd","-f","muonNtuple_%s_%s.root"%(sample,tag)]			
			command += files
			subprocess.call(command,stdout=open(os.devnull, 'wb'))
		else:
			print "no output for sample ",  tag+"/"+sample
	
main()
