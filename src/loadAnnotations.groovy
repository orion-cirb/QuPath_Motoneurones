import org.apache.commons.io.FilenameUtils
import static qupath.lib.scripting.QP.*

// Load annotations from file
def imgDir = new File(project.getImageList()[0].getUris()[0]).getParent()
def imgName = getCurrentImageData().getServer().getMetadata().getName()

// Delete all annotations
clearAllObjects()

// Find annotations files for current image
def p = ~/${imgName}.*\.annot/
def resultsDir = new File(buildFilePath(imgDir+'/Results'))
resultsDir.eachFileMatch(p) {file ->
    new File(file.path).withObjectInputStream {
        def annotations = it.readObject()
        print('Adding annotation ' + annotations.toString())
        addObjects(annotations)
    }
}
resolveHierarchy()
