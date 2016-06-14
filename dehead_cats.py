import os
from os.path import join
import argparse

def main(cat_dir):

	for cat in os.listdir(cat_dir):
	    if cat.startswith("step"):
	        os.system("tail -n +2 " + join(cat_dir, cat) + " > " + join(cat_dir, "a" + cat))
        	os.system("mv " + join(cat_dir, "a" + cat) + " " + join(cat_dir, cat))
        #os.system("sleep 1")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'cat-dir',
		help='steps of random catalog to be de-column-headed')
	args = parser.parse_args()

	main(args.cat_dir)
