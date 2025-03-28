import argparse
import os
import geojson
import numpy as np
import anndata as ad

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

def PointsInPolygon(points, poly):
    x, y = points[:,0], points[:,1]
    n = len(poly)
    inside = np.zeros(len(x), np.bool_)
    p2x, p2y = 0.0, 0.0
    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        idx = np.nonzero((y > min(p1y, p2y)) & (y <= max(p1y, p2y)) & (x <= max(p1x, p2x)))[0]
        if len(idx):  # If idx is non-empty, proceed with the calculation
            if p1y != p2y:
                xints = (y[idx] - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1x == p2x:
                inside[idx] = ~inside[idx]
            else:
                idxx = idx[x[idx] <= xints]
                inside[idxx] = ~inside[idxx]    
        p1x, p1y = p2x, p2y
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
print("Found {0} UNIQUE annotations: ".format(len(AnnotationNames)) + " ".join(AnnotationNames))

## --------------------------------------------------------------------------------
## Read in AnnData and assign annotations to observations based on
## whether the observation spatial coordinates are contained with polygons defined
## by annotation coordinates
## --------------------------------------------------------------------------------

adata = ad.read_h5ad(args.file)

for AN in AnnotationNames:
    adata.obs[AN] = False

X_origin = [min(adata.obsm['X_spatial'][:, 0]), min(adata.obsm['X_spatial'][:, 1])]

print("\n")
for Annotation in Annotations:
    
    if Annotation['geometry']['type'] == 'Polygon':
        # Single polygon, just use the coordinates directly
        coords = np.array(Annotation['geometry']['coordinates'][0]) 
        InPolygon = PointsInPolygon(np.array(adata.obsm['X_spatial'] - X_origin), coords)
        adata.obs[Annotation["properties"]["classification"]["name"]][np.where(InPolygon)[0]] = True
    elif Annotation['geometry']['type'] == 'MultiPolygon':
        # Multiple polygons MultiPolygon, iterate over each
        for poly in Annotation['geometry']['coordinates']:
            coords = np.array(poly[0])  
            InPolygon = PointsInPolygon(np.array(adata.obsm['X_spatial'] - X_origin), coords)
            adata.obs[Annotation["properties"]["classification"]["name"]][np.where(InPolygon)[0]] = True

    print("{0} contains {1} points.".format(Annotation["properties"]["classification"]["name"], adata.obs[Annotation["properties"]["classification"]["name"]].sum()))

adata.write(outfilename)
print("\n")
print("Wrote AnnData file: {0}".format(outfilename))

print("\n\n\nDone!")

