# Lineage-Mapper

![status](https://img.shields.io/badge/status-passed-blue.svg)

@Author: Ali Hashmi

The Wolfram Language implementation of the overlap-based cell tracker is based on an article by `Chalfoun et al, Scientific Reports, 2016` (https://www.nature.com/articles/srep36984). The algorithm can robustly track cells, corrects for any incorrect mergers/fusions of cells in an automated fashion and detects possible mitotic events using stringent criterions for mother/daughter cells. The algorithm used for tracking is not specific to any segmentation scheme, which indeed is a great thing and it makes lineage mapper more robust than many tracking schemes out there ! (Read the article for more detail). In short, the tracking scheme uses a minimization scheme over a cost matrix (Hungarian Algorithm). In this implementation, I have used a graph theoretic approach to generate the same results. 


***Note to reader: This particular implemenation enhances Lineage Mapper by introducing additional functionalities. The implementation may seem quite cryptic unless you have a good understanding of and familiarity with functional, pattern matching and rule-based programming as well as the standard and non-standard evaluation sequence in the Wolfram Language. So if you wish to use Lineage Mapper or want some functionality incorporated, please do not hesistate to ask. For users interested in the authors' version, download FIJI or ImageJ plugin from:***
https://pages.nist.gov/Lineage-Mapper/


##### `To do:`
###### * 1. make edge tracking faster



Shown below are a few capabilities of "Lineage Mapper" and some of the additional functionality

#### 1. `Segment and Track a group of cells`

* ##### example 1

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/alain.gif)

* ##### example 2

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitsmask.gif) | ![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitsmasksegtracked.gif)

#### 2. `View the displacements of cells`

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/centroiddispmap.png)

#### 3. `Get Information about cell neighbours for a particular time`

##### **here the first column depicts the cell of interest and the second column the labels of its neighbours**

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/cellneighbours.png)


#### 4. `View individual cell or a set of cells in isolation`

##### we can isolate cell(s) from the mask for better viewing 

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/seesingleormultiplecells.png)


#### 5. `Tracking an Edge`

##### ****shown in red****

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitedgetrack.gif)


#### 6. `Incorrect Fusions are dealt with`

##### Sometimes during segmentation neighbouring cells can fuse together to form a single cell. This behaviour is in most cases not desired. Lineage Mapper can break such clusters of cells.

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/trackedandmaskcorrected.png)

#### 7. `Generate Lineage Table/Tree and BirthDeath Table`

##### **we can get the parent-daughter associations in form of either a table or a graph**

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/lineageTree%26table.png)


##### **Note: the death toll is comprised of only those cells that are unassigned from the previous frame. This can change if one uses the centroid information of the tracked cells. Notably, cells that do not move much across multiple frames can be considered dead. Information about apoptosis can be obtained separately (using LineageMapperTools.m) **

#### 8. `Fusions of cells into clusters`

##### cells can be tracked while they merge or fuse together to form colonies.

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/fusions1.png)

##### A fusion graph can be generated (the graph/tree delineates cells that fuse together to form new colonies/entities).

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/fusionstree.png)

##### Lineage tree/table for mitotic events and the birth-death table.

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/fusions2.png)


#### 9. `Cell as Regions for Computations`

##### cells as image can be converted to regions which can be further used for computations (solving PDEs etc..)

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/solveequations%20over%20a%20cell.png)

