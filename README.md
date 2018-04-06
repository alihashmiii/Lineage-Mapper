# Lineage-Mapper
The Wolfram Language implementation of the overlap-based cell tracker is based on an article by Chalfoun et al, Scientific Reports, 2016. The algorithm can robustly track cells, correct for any incorrect mergers/fusions of cells in an automated fashion and detect possible mitotic events using stringent criterions for mother/daughter cells. Read the article for more detail. The algorithm used for tracking is not specific to the segmentation scheme, which indeed is a great thing !

**Note: I have not implemented the case when fusions or mergers of cells is desired, since it is usually not the case of interest**

Here are a few capabilities that i have introduced on top of "Lineage Mapper"

1. Segment and Track a group of cells

binarized mask (fed as input)

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitsmask.gif)

after segmentation and tracking

![alt text](https://github.com/alihashmiii/Lineage-Mapper/blob/master/uploadReadMe/benoitsmasksegtracked.gif)
