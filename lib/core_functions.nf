//
// This file holds several functions specific to the workflow/mhc_hammer.nf in the MHChammer pipeline
//

def filterByMetadataField(fieldName, fieldValue) {
    { channelItem -> channelItem[0].get(fieldName) == fieldValue }
}

def check_hlahd_results(file_tuple) {
    def fileContent = file_tuple[1].readLines()
    def counter = 0
    for (line in fileContent) {
        if (line.contains('not typed')) {
            counter++
        }
    }
    return counter == fileContent.size()
}

def check_tumour_results(input_file) {
    def fileContent = input_file.readLines()
    def counter = 0
    for (line in fileContent) {
        if (line.contains('hla_alleles.csv')) {
            counter++
        }
    }
    return counter == fileContent.size()
}

def mutations_detected_check(input_channel) {
    def mutationFiles = input_channel[4]
    def mutationCount = mutationFiles.size()
    if (mutationCount == 0) {
        return false
    } else {
        return true
    }
}