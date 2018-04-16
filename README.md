# Lineage-Mapper
@Author: Ali Hashmi

The Wolfram Language implementation of the overlap-based cell tracker is based on an article by `Chalfoun et al, Scientific Reports, 2016` (https://www.nature.com/articles/srep36984). The algorithm can robustly track cells, correct for any incorrect mergers/fusions of cells in an automated fashion and detect possible mitotic events using stringent criterions for mother/daughter cells. The algorithm used for tracking is not specific to the segmentation scheme, which indeed is a great thing and makes lineage mapper more robust than many tracking schemes out there ! (Read the article for more detail). In short, the tracking scheme uses a minimization scheme over a cost matrix (Hungarian Algorithm). In this implementation, I have used a graph theoretic approach to generate the same results. 


***Note to reader: This particular implemenation enhances Lineage Mapper by introducing additional functionalities. The implementation may seem quite cryptic unless you have a good understanding of and familiarity with functional, pattern matching and rule-based programming as well as the standard and non-standard evaluation sequence in the Wolfram Language. So if you wish to use Lineage Mapper or want some functionality incorporated, please do not hesistate to ask. Moreover, as of now, I have not implemented the case when fusions or mergers of cells is desired by the users, since fusion of cells is usually not the case of interest. For users interested in the latter, download FIJI or ImageJ plugin from:***
https://pages.nist.gov/Lineage-Mapper/


##### `To do:`
###### * 1. make edge tracking faster



Here are a few capabilities of "Lineage Mapper" and some the functionality on top

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

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/tablesandtree.png)


##### **Note: the death toll for now is comprised of only those cells that are unassigned from the previous frame. This can change if one uses the centroid informations of the tracked cells. Notably, cells that do not move much across multiple frames can be considered dead.**


#### 8. `Cell as Regions for Computations`

##### cells as image can be converted to regions which can be further used for computations (solving PDEs etc..)

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/solveequations%20over%20a%20cell.png)

