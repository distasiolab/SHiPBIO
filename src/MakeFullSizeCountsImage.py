import argparse

import numpy as np
import anndata as ad
import scanpy as sc

from PIL import Image, ImageFilter, ImageEnhance
from matplotlib import cm

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
import argparse
import os
import xml.etree.ElementTree as ET
import geojson

from numba import jit


import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

## --------------------------------------------------------------------------------
## Input arguments
## --------------------------------------------------------------------------------

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="AnnData File")
parser.add_argument("-c", "--coords", help="Coordinates File (GeoJSON file format; *.json)")
parser.add_argument("-o", "--out", help="Output filename")

# Read arguments from command line
args = parser.parse_args()

if args.out is None:
    outfilename = os.path.splitext(args.file)[0] + "_annotated.h5ad"
else:
    outfilename = args.out

print("\n")
print("Input image file: {0}".format(args.file))
print("Annotations coordinates file: {0}".format(args.coords))
print("Output AnnData file: {0}".format(outfilename))

## --------------------------------------------------------------------------------
## Method Definitions
## --------------------------------------------------------------------------------

def dim(a):
    if not type(a) == list:
        return []
    return [len(a)] + dim(a[0])


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


## --------------------------------------------------------------------------------
## Main Method
## --------------------------------------------------------------------------------

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

print("\n")
print("Found {0} annotations.".format(len(Annotations)))
print("Found {0} UNIQUE annotations: ".format(len(AnnotationNames))+" ".join(AnnotationNames))

## --------------------------------------------------------------------------------
## Read in AnnData and assign annotations to observations based on
## whether the observation spatial coordinates are contained with polygons defined
## by annotation coordinates
## --------------------------------------------------------------------------------

adata = ad.read_h5ad(args.file)

for AN in AnnotationNames:
    adata.obs[AN] = False

X_origin = [min(adata.obsm['X_spatial'][:,0]),min(adata.obsm['X_spatial'][:,1])]

print("\n")
for Annotation in Annotations:

    print(len(dim(Annotation['geometry']['coordinates'][0])))
    if len(dim(Annotation['geometry']['coordinates'][0])) > 2:
        coords = np.array(Annotation['geometry']['coordinates'][0][0])
    else:
        coords = np.array(Annotation['geometry']['coordinates'][0])

    print(coords)
        
    InPolygon = PointsInPolygon(np.array(adata.obsm['X_spatial'] - X_origin), coords)
    adata.obs[Annotation["properties"]["classification"]["name"]][np.where(InPolygon)[0]] = True
    print("{0} contains {1} points.".format(Annotation["properties"]["classification"]["name"],np.sum(InPolygon)))


adata.write(outfilename)
print("\n")
print("Wrote AnnData file: {0}".format(outfilename))
    
print("\n\n\nDone!")

print("\n")

if args.out == None:
    outfile =  args.file + '.FullSizeImage_TotalCounts.tif'
else:
    outfile = args.out
print("\n\n")
print("Output image file: {0}".format(outfile))
print("\n")

## --------------------------------------------------------------------------------
## Run
## --------------------------------------------------------------------------------


print('Loading...')

ad_filename = args.file
adata = ad.read_h5ad(ad_filename)

# QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


## Filtering
print('Filtering...')
#sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


print('Building Output Image...')

#fig = sc.pl.embedding(adata, basis='spatial', color='leiden', return_fig=True)


minX = np.min(adata.obsm['X_spatial'][:,0])
maxX = np.max(adata.obsm['X_spatial'][:,0])
minY = np.min(adata.obsm['X_spatial'][:,1])
maxY = np.max(adata.obsm['X_spatial'][:,1])


downsample = 1

X,Y = np.meshgrid(np.arange(minX, maxX, downsample), np.arange(minY, maxY, downsample))

@np.vectorize
def BuildCountsImage(x,y):
    # Find all spots within distance 0.5 of the grid location
    spots = np.where( np.logical_and( np.abs(adata.obsm['X_spatial'][:,0] - x) < downsample, np.abs(adata.obsm['X_spatial'][:,1] - y) < downsample ))
    spots = spots[0]
    if spots.size > 0:
        # Get their cluster labels
        return np.mean(adata.obs['total_counts'])
    else:
        return np.nan
    

OUTIMG = BuildCountsImage(X,Y)
OUTIMG = OUTIMG.astype(float)

OUTIMG_norm = (OUTIMG - np.nanmin(OUTIMG))/(np.nanmax(OUTIMG) - np.nanmin(OUTIMG))
im = Image.fromarray(np.uint8(cm.gist_rainbow(OUTIMG_norm)*255)).convert('RGB')
imc = ImageEnhance.Contrast(im.filter(ImageFilter.BoxBlur(radius=1)).filter(ImageFilter.GaussianBlur(radius=2))).enhance(50)

# Save output image file
imc.save(outfile)

print('Done!')

