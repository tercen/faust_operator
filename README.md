# FAUST operator

##### Description

FAUST: Full annotation using shape-constrained trees.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the clusters (e.g. cells) 
`y-axis`| is the value of measurement signal of the channel/marker

Output relations|.
---|---
`umapX`| numeric, coordinate on first UMAP dimension
`umapY`| numeric, coordinate on second UMAP dimension
`faustLabels`| factor, FAUST annotation
`Computed Table 3`| Diagnostic histograms.

##### Details

The operator is a wrapper for the `faust` function of the `faust` [R package](https://github.com/RGLab/FAUST/).

##### See Also

[flowsom_operator](https://github.com/tercen/flowsom_operator)
