#### <a name='goana'></a>GO Enrichment:

The **Goana** method (from the `limma` package) performs an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs.
The list of genes affected in the selected condition is used as a gene set with all genes as a background. <br>
Columns in the resulting table show the following: <br>
- *Term* GO Term 
- *Ont* ontology that the GO term belongs to ('BP', 'CC' or 'MF')
- *N* number of genes in the GO term 
- *DE* number of genes in the provided gene set
- *P.DE* p-value for over-representation of the GO term in the set

Click the **Goana** button to create the table.
