import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import json


#Read in snif file

def snif_to_mosfit(sfile):
	snif_file = open(sfile)
	snif_str = snif_file.read()
	snif_data = json.loads(snif_str)
	sn_name = list(snif_data.keys())[0]
	snif_model = snif_data[sn_name]['models']
	print('python3 -m mosfit -e '+sn_name, end =" ")
	print('-i 1 -N 100 -F ', end =" ")

	for thing in snif_model['parameters']:
		val = snif_model['parameters'][thing]['value']
		print(thing+' '+str(val)+ ' ', end=" ")
	print()

	return 0

def main(argv):
	inputfile = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hs:",["snif="])
	except getopt.GetoptError:
		print('test.py -s <snif_file> -m <mosfit_file> -o <output_file>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('convert.py -i <sniffile> -m <mosfitfile> -o <outputfile>')
			sys.exit()
		elif opt in ("-s", "--snif"):
			sniffile = arg
	snif_to_mosfit(sniffile)
if __name__ == "__main__":
	main(sys.argv[1:])