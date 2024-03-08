// ----------------------------------------------------------------------------------------
// GetDetectionParents.groovy
// Marcello DiStasio
// March, 2024
//
// QuPath script that writes to a CSV file the list of parents for each detection object in 
// an image
// ----------------------------------------------------------------------------------------



// ----------------------------------------------------------------------------------------
// Method defintions
// ----------------------------------------------------------------------------------------
def get_parentlist(it) {
    ArrayList<String> parents = new ArrayList<String>();
    
    x = it;
    while (x.getLevel() != 1) {
        parents.add(x.getParent().getName())
        x = x.getParent()    
    }
    
    for (int i = 0; i < parents.size() / 2; i++) {
            String temp = parents.get(i);
            parents.set(i, parents.get(parents.size() - i - 1));
            parents.set(parents.size() - i - 1, temp);
        }
 
    return(parents)
}


// ----------------------------------------------------------------------------------------
// Main method
// ----------------------------------------------------------------------------------------
resolveHierarchy()
detections = getDetectionObjects()

// Set up file
String imageLocation = getCurrentImageData().getServer().getURIs()[0].getPath()
def outputFileName = imageLocation.concat('_Detection_Annotation_Hierarchy.csv')


FileWriter writer = new FileWriter(outputFileName);

// Construct output array
ArrayList<String> datalines = new ArrayList<String>();
for (int i = 0; i < detections.size(); i++) {
    ArrayList<String> line = new ArrayList<String>();
    line.add(detections[i].getID().toString());
    line.addAll(get_parentlist(detections[i]));
    String csv_line = String.join(",", line);
    writer.append(csv_line);
    writer.append('\n')
}
writer.close();


print('Wrote '.concat(outputFileName))
print('Done!')
