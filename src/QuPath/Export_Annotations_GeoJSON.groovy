



def annotations = getAnnotationObjects()
def entry = getProjectEntry() 

String imageLocation = getCurrentImageData().getServer().getURIs()[0].getPath()
def outfile = imageLocation.concat('_Annotations.json')

// 'FEATURE_COLLECTION' is standard GeoJSON format for multiple objects
exportObjectsToGeoJson(annotations, outfile, "FEATURE_COLLECTION")

print('Exported '.concat(annotations.size().toString()).concat(' Annotation(s) to: ').concat(outfile))
print('Done')
