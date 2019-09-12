<div><br>This page contains descriptions for the different features of the TREND-DB. </div><br>
      <h4>Data Preview</h4>
      <div>The <strong>Data Preview</strong> tab shows an overview of all shortening index values of each gene across all conditions. With the
      '<strong>Filter by shortening index</strong>' slider, the table can be subsetted to only show genes with a shortening index value (absolute) above the threshold
      in at least one condition. Upon selecting a row (= gene) in the table, a summary table of the selected gene is printed, showing the conditions in which the gene is affected
      with the respective shortening indices and p-values.</div> <br> <br>
      <h4>Main View</h4>
      <div>The <strong>Main View</strong> tab is divided into <strong>Condition View</strong> and <strong>Gene View</strong>.
      In the <strong>Condition View</strong>, users can select a condition to view a table containing all genes that are affected in this condition.
      By clicking the <strong>Goana</strong> button, an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs
      is performed and the results are shown in a table.  The list of genes affected in the selected condition is used as a gene set with all genes as a background.
      In the <strong>Gene View</strong>, users can select a gene to print a short summary and a table containing all conditions where the selected gene is affected.</div><br><br>
      <h4>Gene Plot</h4>
      <div>In the <strong>Gene Plot</strong> tab, users can plot the selected gene (using the <strong>Gviz</strong> package). Aside from the reference track of the selected gene and its chromosome, the control track as well as the selected condition are shown.
      The default condition is the condition with the maximum SI for this gene, it can however be changed using a drop-down menu.
      Adjusting the <strong> Number of similarily affected</strong> slider will print genes that are affected in a similar way as the selected gene (with the amount of them determined by the slider).
      Similarily affected genes are obtained by computing the distance matrix (using the shortening indices of all genes across conditions) and looking for genes with the smallest distances.</div><br><br>
      <h4>Genome Browser</h4>
      <div> The <strong>Genome Browser</strong> tab includes an instance of the UCSC Genome Browser, where users can take a more in-depth look at the data. </div><br><br>
      <h4>Contact</h4>
      <div>Denise Scherzinger <a href='mailto:denscher@uni-mainz.de' target='_top'>(denscher@uni-mainz.de)</a></div>