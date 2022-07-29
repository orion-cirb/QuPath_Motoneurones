// Imports
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.RoiTools
import qupath.lib.roi.ROIs
import static qupath.lib.gui.scripting.QPEx.*
import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.objects.*


// Init project
setImageType('Fluorescence')
def project = getProject()
def imageDir = new File(project.getImageList()[0].getUris()[0]).getParent()

// Create results file and write headers
def resultsDir = buildFilePath(imageDir, '/Results')
if (!fileExists(resultsDir)) mkdirs(resultsDir)
def resultsFile = new File(buildFilePath(resultsDir, 'Results.csv'))
resultsFile.createNewFile()
def resHeaders = 'Image name\tNeuron ID\tStart slice\tStop slice\tVolume (um^3)\tBackground mean AF546 intensity\tNeuron mean AF546 intensity\tNucleus mean AF546 intensity\tCytoplasm mean AF546 intensity\n'
resultsFile.write(resHeaders)

// Build CellPose models
def cellposeNeurons = Cellpose2D.builder('cyto2')
        .pixelSize(0.25)          // Resolution for detection
        .channels('AF488-T2')    // Select detection channel(s)
//        .tileSize(2048)
//        .preprocess(ImageOps.Filters.median(1))           // List of preprocessing ImageOps to run on the images before exporting them
//        .maskThreshold(-0.2)            // Threshold for the mask detection, defaults to 0.0
//        .flowThreshold(0.5)             // Threshold for the flows, defaults to 0.4
        .diameter(100)           // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
        .excludeEdges()                  // Clears objects touching the edge of the image (Not of the QuPath ROI)
//        .clusterDBSCAN()               // Use DBSCAN clustering to avoid over-segmenting long object
        .measureShape()                  // Add shape measurements
        .measureIntensity()              // Add cell measurements (in all compartments)
        .classify(PathClassFactory.getPathClass('Neuron', makeRGB(0,255,0)))
        .useGPU()                        // Optional: Use the GPU if configured, defaults to CPU only
        .build()

def cellposeNuclei = Cellpose2D.builder('cyto2')
        .pixelSize( 0.25)        // Resolution for detection
        .channels( 'DAPI-T3')   // Select detection channel(s)
//        .tileSize(2048)
//        .preprocess( ImageOps.Filters.median(1))           // List of preprocessing ImageOps to run on the images before exporting them
//        .maskThreshold(-0.2)           // Threshold for the mask detection, defaults to 0.0
//        .flowThreshold(0.5)            // Threshold for the flows, defaults to 0.4
        .diameter(25)            // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
        .excludeEdges()                  // Clears objects touching the edge of the image (Not of the QuPath ROI)
//        .clusterDBSCAN()               // Use DBSCAN clustering to avoid over-segmenting long object
        .measureShape()                  // Add shape measurements
        .measureIntensity()              // Add cell measurements (in all compartments)
        .classify(nucleusPathClass = PathClassFactory.getPathClass('Nucleus', makeRGB(0,0,255)))
        .useGPU()                       // Optional: Use the GPU if configured, defaults to CPU only
        .build()

def sortDetectionsBySlice(detections, nbSlices) {
    def sortedDetections = []
    for (def i=0; i < nbSlices; i++) {
        sortedDetections << detections.findAll{it.getParent().getName() == 'Slice ' + i}
    }
    return sortedDetections
}

def stitch3D(neurons) {
    def neuronID = 0;
    for (n in neurons[0]) {
        n.setName("Neuron " + neuronID)
        n.getMeasurementList().putMeasurement("Neuron ID", neuronID)
        n.getMeasurementList().putMeasurement("Z slice", 0)
        neuronID++
    }

    for (def i=0; i < neurons.size()-1; i++) {
        for (n2 in neurons[i+1]) {
            def n2ROI = n2.getROI()
            for (n1 in neurons[i]) {
                def n1ROI = n1.getROI()
                if (RoiTools.areaContains(n1ROI, n2ROI.getCentroidX(), n2ROI.getCentroidY())) {
                    def id = n1.getMeasurementList().getMeasurementValue("Neuron ID")
                    n2.setName("Neuron " + (int) id)
                    n2.getMeasurementList().putMeasurement("Neuron ID", id)
                    break
                }
            }
            if (! n2.getMeasurementList().containsNamedMeasurement("Neuron ID")) {
                n2.setName("Neuron " + neuronID)
                n2.getMeasurementList().putMeasurement("Neuron ID", neuronID)
                neuronID++
            }
            n2.getMeasurementList().putMeasurement("Z slice", i+1)
        }
    }
    fireHierarchyUpdate()
    return neuronID
}

def getNucleiInNeurons(nuclei, neurons) {
    def filteredNuclei = []
    for (def i=0; i < neurons.size(); i++) {
        for (nucleus in nuclei[i]) {
            def roiNucleus = nucleus.getROI()
            for (neuron in neurons[i]) {
                def roiNeuron = neuron.getROI()
                def intersection = RoiTools.intersection([roiNucleus, roiNeuron])
                if (intersection.getArea() == roiNucleus.getArea()) {
                    def id = neuron.getMeasurementList().getMeasurementValue("Neuron ID")
                    nucleus.setName("Nucleus " + (int) id)
                    nucleus.getMeasurementList().putMeasurement("Neuron ID", id)
                    nucleus.getMeasurementList().putMeasurement("Z slice", neuron.getMeasurementList().getMeasurementValue("Z slice"))
                    neuron.addPathObject(nucleus)
                    filteredNuclei << PathObjectTools.transformObject(nucleus, null, true)
                    break
                }
            }
        }
    }
    resolveHierarchy()
    return(filteredNuclei)
}

def getNeuronsWithNuclei(neurons, neuronIDMax) {
    def filteredNeurons = []
    for(def id = 0; id < neuronIDMax; id++) {
        def neuron = neurons.flatten().findAll{it.getName() == ("Neuron " + id)}
        def hasNucleus = false
        for (n in neuron) {
            if (n.hasChildren()) {
                hasNucleus = true
                break
            }
        }
        if ((hasNucleus == true) && (neuron.size() > 1)) {
            for (n in neuron) {
                filteredNeurons << n
            }
        }
    }


    return filteredNeurons
}

// Get mean intensity of a population of objects
def saveNeuronsParameters(imgName, resultsFile, neurons, neuronIDMax, backgrounds, pixelCal) {
    def toMeasure = 'AF546-T1: Mean'

    for(def id = 0; id < neuronIDMax; id++) {
        def neuron = neurons.flatten().findAll{ it.getName() == ("Neuron " + id)}
        def nbNeuronRois = neuron.size()
        if (nbNeuronRois > 0) {
            def slices = []
            def bgMean = 0
            def neuronMean = 0
            def neuronVolume = 0
            def nucleusMean = 0
            def nbNucleusRois = 0
            def cytoplasmMean = 0
            for (n in neuron) {
                slices << n.getMeasurementList().getMeasurementValue('Z slice')
                def bg = backgrounds[n.getMeasurementList().getMeasurementValue('Z slice')].getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: ' + toMeasure)
                bgMean += bg
                neuronMean += n.getMeasurementList().getMeasurementValue(toMeasure) - bg

                def area = n.getMeasurementList().getMeasurementValue('Area µm^2')
                neuronVolume += area * pixelCal.ZSpacingMicrons
                def cytoplasmSum = n.getMeasurementList().getMeasurementValue(toMeasure) * area
                def areaCytoplasm = area

                if (n.hasChildren()) {
                    def nucleus = n.getChildObjects()[0]
                    nucleusMean += nucleus.getMeasurementList().getMeasurementValue(toMeasure) - bg
                    nbNucleusRois++
                    cytoplasmSum -= nucleus.getMeasurementList().getMeasurementValue(toMeasure) * nucleus.getMeasurementList().getMeasurementValue('Area µm^2')
                    areaCytoplasm -= nucleus.getMeasurementList().getMeasurementValue('Area µm^2')
                }
                cytoplasmMean += cytoplasmSum / areaCytoplasm - bg
            }
            bgMean /= nbNeuronRois
            neuronMean /= nbNeuronRois
            nucleusMean /= nbNucleusRois
            cytoplasmMean /= nbNeuronRois

            //Neuron - nucleus mean AF546 intensity
            def results = imgName + '\t' + id + '\t' + (int)slices.min() + '\t' + (int)slices.max() + '\t' + neuronVolume + '\t' + bgMean + '\t' + neuronMean + '\t' + nucleusMean + '\t' + cytoplasmMean + '\n'
            resultsFile << results
        }
    }
}

// Save annotations
def saveAnnotations(imgName) {
    def path = buildFilePath(imgName + '.annot')
    def annotations = getAnnotationObjects()
    new File(path).withObjectOutputStream {
        it.writeObject(annotations)
    }
    println('Annotations saved')
}

// Loop over images in project
for (entry in project.getImageList()) {
    def imageData = entry.readImageData()
    def server = imageData.getServer()
    def cal = server.getPixelCalibration()
    def imgName = entry.getImageName()
    println ''
    println ''
    println '--------- ANALYZING IMAGE ' + imgName + ' ---------'

    setBatchProjectAndImage(project, imageData)
    def annotations = getAnnotationObjects()
    def background = annotations.findAll{it.getName().contains('Background')}[0]
    if (background == null) {
        Dialogs.showErrorMessage("Problem", "Please provide background ROI in image " + imgName)
        continue
    }
    def backgroundX = background.getROI().getBoundsX()
    def backgroundY = background.getROI().getBoundsY()
    def backgroundHeight = background.getROI().getBoundsHeight()
    def backgroundWidth = background.getROI().getBoundsWidth()

    def backgrounds = []
    def rois = []
    for (def i=0; i < server.nZSlices(); i++) {
        def bg = new PathAnnotationObject(ROIs.createRectangleROI(backgroundX, backgroundY, backgroundWidth, backgroundHeight, ImagePlane.getPlane(i, 0)))
        addObject(bg)
        backgrounds << bg

        def roi = new PathAnnotationObject(ROIs.createRectangleROI(0, 0, server.getWidth(), server.getHeight(), ImagePlane.getPlane(i, 0)))
        roi.setName("Slice " + i)
        addObject(roi)
        rois << roi
    }
    deselectAll()
    selectObjects(backgrounds)
    runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": true,  ' +
            '"channel2": true,  "channel3": true, "doMean": true,  "doStdDev": false,  "doMinMax": false,  "doMedian": false,  "doHaralick": false}')

    println '--------- Detecting neurons ---------'
    deselectAll()
    selectObjects(rois)
    cellposeNeurons.detectObjects(imageData, getSelectedObjects())
    def neurons = sortDetectionsBySlice(getDetectionObjects(), server.nZSlices())
    def neuronIDMax = stitch3D(neurons)

    println '--------- Detecting nuclei ---------'
    deselectAll()
    selectObjects(rois)
    cellposeNuclei.detectObjects(imageData, getSelectedObjects())
    def nuclei = sortDetectionsBySlice(getDetectionObjects(), server.nZSlices())
    def filteredNuclei = getNucleiInNeurons(nuclei, neurons)
    def filteredNeurons = getNeuronsWithNuclei(neurons, neuronIDMax)

    saveNeuronsParameters(imgName, resultsFile, filteredNeurons, neuronIDMax, backgrounds, cal)
    selectObjects(backgrounds)
    clearSelectedObjects()
    clearDetections()
    addObjects(filteredNeurons)
    resolveHierarchy()
    saveAnnotations(buildFilePath(resultsDir, imgName))
}
println '--------- Analysis done! ---------'


