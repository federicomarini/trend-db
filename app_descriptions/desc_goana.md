#### <a name='goana'></a>GO Enrichment:

The **GOana** method (from the `limma` package) performs an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs.
The list of genes affected in the selected condition is used as a gene set with all genes as a background. <br>
Columns in the resulting table show the following information: <br>

- *GO*: the ID for the GO term (with a button to open the corresponding entry to the AmiGO database)
- *Term*: The description of the GO term
- *Ont*: The ontology that the GO term belongs to ('BP', 'CC' or 'MF')
- *N*: the number of genes in the GO term 
- *DE*: the number of genes in the provided gene set (the ones affected in the TREND-DB condition)
- *P.DE*: the p-value for over-representation of the GO term in the set

Click the **GOana** button to create the table.
