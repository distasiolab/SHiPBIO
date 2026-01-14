import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.PathDetectionObject
import qupath.lib.objects.PathClass

// -------------------------
// SETTINGS
// -------------------------

// Convert only detections inside selected annotations?
boolean ONLY_WITHIN_SELECTED_ANNOTATIONS = false

// Optional: require detections to have this class to be converted (null = no filter)
String REQUIRED_DETECTION_CLASS = null  // e.g. "Positive" or null

// ---- Classification behavior for NEW annotations ----
// Choose ONE of these modes by setting CLASS_MODE:
//
// "COPY"  -> annotation gets the same class as the detection
// "FIXED" -> annotation gets FIXED_ANNOTATION_CLASS
// "MAP"   -> annotation class determined by CLASS_MAP (fallback: MAP_DEFAULT_CLASS)
String CLASS_MODE = "FIXED"

// For FIXED mode:
String FIXED_ANNOTATION_CLASS = "Tumor"

// For MAP mode (keys are detection class names, values are annotation class names):
def CLASS_MAP = [
    "Positive": "Tumor",
    "Negative": "Stroma",
    "Unclassified": "Other"
]
String MAP_DEFAULT_CLASS = "Other"   // used if detection class not in map or null

// Delete detections after converting?
boolean DELETE_DETECTIONS_AFTER = false


// -------------------------
// HELPERS
// -------------------------

PathClass getOrCreateClass(String name) {
    if (name == null) return null
    def pc = getPathClass(name)
    if (pc != null) return pc
    // In QuPath, classes are usually created implicitly when referenced
    // getPathClass(name) often returns a PathClass even if newly referenced,
    // but if your version doesn't, we still handle it by returning getPathClass(name) again.
    return getPathClass(name)
}

PathClass resolveAnnotationClass(PathDetectionObject det) {
    def detClass = det.getPathClass()
    def detName = detClass == null ? null : detClass.toString()

    if (CLASS_MODE == "COPY") {
        return detClass
    }
    if (CLASS_MODE == "FIXED") {
        return getOrCreateClass(FIXED_ANNOTATION_CLASS)
    }
    if (CLASS_MODE == "MAP") {
        def targetName = (detName != null && CLASS_MAP.containsKey(detName)) ? CLASS_MAP[detName] : MAP_DEFAULT_CLASS
        return getOrCreateClass(targetName)
    }

    // If unknown mode, don't set a class
    return null
}


// -------------------------
// MAIN
// -------------------------

def hierarchy = getCurrentHierarchy()
if (hierarchy == null) {
    print "No hierarchy found (is an image open?)"
    return
}

// Candidate detections
Collection candidates
def selectedAnnotations = getSelectedObjects().findAll { it.isAnnotation() }

if (ONLY_WITHIN_SELECTED_ANNOTATIONS) {
    if (selectedAnnotations.isEmpty()) {
        print "ONLY_WITHIN_SELECTED_ANNOTATIONS is true, but no annotations are selected."
        return
    }
    candidates = []
    selectedAnnotations.each { ann ->
        candidates.addAll(hierarchy.getObjectsForROI(PathDetectionObject, ann.getROI()))
    }
} else {
    candidates = hierarchy.getDetectionObjects()
}

if (candidates == null || candidates.isEmpty()) {
    print "No detections found."
    return
}

// Optional filter: required detection class
PathClass requiredDetClass = null
if (REQUIRED_DETECTION_CLASS != null) {
    requiredDetClass = getPathClass(REQUIRED_DETECTION_CLASS)
    if (requiredDetClass == null) {
        print "Required detection class '${REQUIRED_DETECTION_CLASS}' not found."
        return
    }
}

// Convert
def detections = candidates.findAll { it instanceof PathDetectionObject }
if (requiredDetClass != null) {
    detections = detections.findAll { it.getPathClass() == requiredDetClass }
}

print "Detections selected for conversion: ${detections.size()}"

def newAnnotations = []
detections.each { PathDetectionObject det ->
    def roi = det.getROI()
    if (roi == null) return

    def ann = new PathAnnotationObject(roi)

    def annClass = resolveAnnotationClass(det)
    if (annClass != null) {
        ann.setPathClass(annClass)
    }

    newAnnotations << ann
}

hierarchy.addPathObjects(newAnnotations)

if (DELETE_DETECTIONS_AFTER) {
    hierarchy.removeObjects(detections, true)
}

fireHierarchyUpdate()
print "Created ${newAnnotations.size()} annotations (CLASS_MODE=${CLASS_MODE})."
