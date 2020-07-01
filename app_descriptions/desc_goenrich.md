#### <a name='goenrich'></a>GO Enrichment:

The **enrichGO** method (from the `clusterProfiler` package) performs an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs.
The list of genes affected in the selected condition is used as a gene set, with all genes as a background. <br>
Columns in the resulting table show the following information: <br>

- *ID*: The ID for the GO term (with a button to open the corresponding entry to the AmiGO database)
- *Description*: The description of the GO term
- *GeneRatio*: The ratio between the number of genes associated with the term VS the ones not associated, in the subset of the affected genes for the condition
- *BgRatio*: The ratio between the number of genes associated with the term VS the ones not associated, in the set of background genes
- *pvalue*: The p-value for over-representation of the GO term in the set of affected genes versus the background
- *p.adjust*: The adjusted p-value, after correction via the Benjamini-Hochberg method (to control the FDR, False Discovery Rate)
- *qvalue*: The q-value for the gene set, as a means to control the positive False Discovery Rate (pFDR)
- *geneID*: The list of affected genes associated with the term, separated by the "/" character symbol
- *Count*: The number of affected genes associated with the GO term
