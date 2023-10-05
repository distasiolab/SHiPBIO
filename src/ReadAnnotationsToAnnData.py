import argparse

import xml.etree.ElementTree as ET
import geojson

from numba import jit


import numpy as np
import pandas as pd
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




@jit(nopython=True)
def PointsInPolygon(points, poly):
    x,y = points[:,0], points[:,1]
    n = len(poly)
    inside = np.zeros(len(x),np.bool_)
    p2x = 0.0
    p2y = 0.0
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        idx = np.nonzero((y > min(p1y,p2y)) & (y <= max(p1y,p2y)) & (x <= max(p1x,p2x)))[0]
        if len(idx):    # <-- Fixed here. If idx is null skip comparisons below.
            if p1y != p2y:
                xints = (y[idx]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            if p1x == p2x:
                inside[idx] = ~inside[idx]
            else:
                idxx = idx[x[idx] <= xints]
                inside[idxx] = ~inside[idxx]    

        p1x,p1y = p2x,p2y
    return inside 


with open(args.coords) as f:
    gj = geojson.load(f)
FeatureCollection = gj['features']


Annotations = []

for F in FeatureCollection:
    try:
        # Make sure the right keys are there
        _ = F['properties']['classification']['name']
        Annotations.append(F)
    except KeyError:
        pass
    
AnnotationNames = list(set([A["properties"]["classification"]["name"] for A in Annotations]))

print("Found {0} annotations.".format(len(Annotations)))
print("Found {0} UNIQUE annotations:".format(len(AnnotationNames)))
print(AnnotationNames)


## --------------------------------------------------------------------------------
## Read in AnnData and assign annotations to observations based on
## whether the observation spatial coordinates are contained with polygons defined
## by annotation coordinates
## --------------------------------------------------------------------------------

adata = ad.read_h5ad(args.file)

X_origin = [min(adata.obsm['X_spatial'][:,0]),min(adata.obsm['X_spatial'][:,1])]

for Annotation in Annotations:
    print(Annotation["properties"]["classification"]["name"])
    InPolygon = PointsInPolygon(np.array(adata.obsm['X_spatial'] - X_origin), np.array(Annotation['geometry']['coordinates'][0]))
    print(np.sum(InPolygon))
    
print("\n\n\nDone!")

