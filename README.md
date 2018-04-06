# Lineage-Mapper

The Wolfram Language implementation of the overlap-based cell tracker is based on an article by `Chalfoun et al, Scientific Reports, 2016`. The algorithm can robustly track cells, correct for any incorrect mergers/fusions of cells in an automated fashion and detect possible mitotic events using stringent criterions for mother/daughter cells. The algorithm used for tracking is not specific to the segmentation scheme, which indeed is a great thing and makes lineage mapper more robust than many tracking schemes out there ! (Read the article for more detail)

****Note: The implementation is not easy to read and understand unless you have a good deal of understanding of and familiarity with functional, pattern matching and rule-based programming and the standard and non-standard evaluation sequence in the Wolfram Language. So if you wish to use it or add some functionality, please do not hesistate to ask. Moreover, as of now, I have not implemented the case when fusions or mergers of cells is desired by the users, since fusion of cells is usually not the case of interest****

Here are a few capabilities that i have introduced on top of "Lineage Mapper"

## 1. `Segment and Track a group of cells`

![alt text]


![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitsmask.gif) | ![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitsmasksegtracked.gif)

## 2. `View the displacements of cells`

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/centroiddispmap.png)

## 3. `Get Information about cell neighbours for a particular time`

*****here the first column depicts the cell of interest and the second column the labels of its neighbours*****

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/cellneighbours.png)


## 4. `View individual cell or a set of cells in isolation`

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/seesingleormultiplecells.png)


## 5. `Tracking an Edge`

****shown in red****

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitedgetrack.gif)


