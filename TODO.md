# TODO
 
## Issues
* Flanking sequences are not called correctly
    * Causes spacers to be called as flankers
* Speed degrades rapidly with the number of patterns identified
    * Probable cause is the singleton finder
* In some circumstances closely related groups are not getting split apart
* The final output coverage of the spacer does not always correspond to the number of sources
* Spacers go missing in the graph from Crass, but they are present if you perform a regular overlap assembly
    * This is not a sequencing error issue as even when searching with mismatches, there is not sign of the spacers
* Spacer strings only get created during the initial pass of the graph building.  However, spacer objects can also be made from the p-graph
    * This means that spacers may exist in the graph, but not in the output sequences.  Maybe the cause of the above issue
* If a slave DR is a perfect palindrome of the master DR the offset will be undefined
    

## Wanted Features
* Multi-threading the search algorithms
* Using paired read information in the search algorithms
* Ability to run Crass on genomes or assembled contigs
    * Output a gff3 formatted file
* Screen for known eukaryotic microsatellites for datasets that may be host contaminated

## Improvements
* Use read information better when building the graph
* Allow for multiple DR types to be in a single sequence. Very important if Crass is to work on genomes 
