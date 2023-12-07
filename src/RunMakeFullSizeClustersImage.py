import MakeFullSizeClustersImage_pkg as MakeCluster
import argparse

## --------------------------------------------------------------------------------
## Argument parsing
## --------------------------------------------------------------------------------

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="annData *.h5ad file")
parser.add_argument("-o", "--out", help="Output filename", nargs='?')

# Read arguments from command line
args = parser.parse_args()


print("\n\n")
print("Input image file: {0}".format(args.file))
print("\n")

if args.out == None:
    outfile = args.file + '.FullSizeImage_LeidenClusters.tif'
else:
    outfile = args.out

print("\n\n")
print("Output image file: {0}".format(outfile))
print("\n")
cluster = MakeCluster.Cluster(file=args.file, out_file=outfile)
cluster.process()
OUTIMG = cluster.vectorized_method()
write_image = cluster.write_image(outimg=OUTIMG)