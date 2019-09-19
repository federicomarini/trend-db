This page contains descriptions for the different features of the TREND-DB. 

#### Data Preview

The **Data Preview** tab shows an overview of all shortening index values of each gene across all conditions. 
With the '**Filter by shortening index**' slider, the table can be subsetted to only show genes with a shortening index value (absolute) above the threshold in at least one condition. 
Upon selecting a row (= gene) in the table, a summary table of the selected gene is printed, showing the conditions in which the gene is affected with the respective shortening indices and p-values.

#### Main View

The **Main View** tab is divided into **Condition View** and **Gene View**.

TODO: describe the trendnetwork

In the **Condition View**, users can select a condition to view a table containing all genes that are affected in this condition.
By clicking the **Goana** button, an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs is performed and the results are shown in a table. 
The list of genes affected in the selected condition is used as a gene set with all genes as a background.

In the **Gene View**, users can select a gene to print a short summary and a table containing all conditions where the selected gene is affected.

#### Gene Plot

In the **Gene Plot** tab, users can plot the selected gene (using the `Gviz` package). 
Aside from the reference track of the selected gene and its chromosome, the control track as well as the selected condition are shown.
The default condition is the condition with the maximum SI for this gene, it can however be changed using a drop-down menu.

Adjusting the **Number of similarly affected** slider will print genes that are affected in a similar way as the selected gene (with the amount of them determined by the slider).
Similarly affected genes are obtained by computing the distance matrix (using the shortening indices of all genes across conditions) and looking for genes with the smallest distances.

#### Genome Browser

The **Genome Browser** tab includes an instance of the UCSC Genome Browser, where users can take a more in-depth look at the data. 

#### Contact

- Federico Marini <a href='mailto:marinif@uni-mainz.de' target='_top'>(marinif@uni-mainz.de)</a>
- Denise Scherzinger <a href='mailto:denscher@uni-mainz.de' target='_top'>(denscher@uni-mainz.de)</a>
