import argparse

import xml.etree.ElementTree as ET
import geojson

import numpy as np
import anndata as ad
import scanpy as sc

## --------------------------------------------------------------------------------
##
## --------------------------------------------------------------------------------

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="AnnData File")
parser.add_argument("-c", "--coords", help="Coordinates File (GeoJSON file format; *.json)")
parser.add_argument("-o", "--out", help="Output filename", nargs='?', const='')

# Read arguments from command line
args = parser.parse_args()

print("\n\n")
print("Input image file: {0}".format(args.file))
print("\n")
print("Annotations coordinates file: {0}".format(args.coords))
print("\n")

with open(args.coords) as f:
    gj = geojson.load(f)
FeatureCollection = gj['features']

CropRegions = list(filter(lambda region: region["properties"]["classification"]["name"] == "CropRegion", FeatureCollection))

print("Found {0} annotations with 'CropRegion' label. They are: [left, top, width, heigth]".format(len(CropRegions)))

crop_region_specs = []
for Region in CropRegions:
    left = min([i[0] for i in Region["geometry"]["coordinates"][0]])
    top = min([i[1] for i in Region["geometry"]["coordinates"][0]])
    width = max([i[0] for i in Region["geometry"]["coordinates"][0]]) - left
    height = max([i[1] for i in Region["geometry"]["coordinates"][0]]) - top

    crop_region_specs.append([left,top,width,height])

print(crop_region_specs)


## --------------------------------------------------------------------------------
## Read in AnnData and assign annotations to observations based on
## whether the observation spatial coordinates are contained with polygons defined
## by annotation coordinates
## --------------------------------------------------------------------------------

adata = ad.read_h5ad(args.file)

cr = 1    
for CropR in crop_region_specs:

    left = CropR[0]
    top = CropR[1]
    width = CropR[2]
    height = CropR[3]


    
print("\n\n\nDone!")

