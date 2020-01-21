---


---

<blockquote>
<h1 id="normalisation">Normalisation</h1>
</blockquote>
<p><a href="https://osca.bioconductor.org/normalization.html">https://osca.bioconductor.org/normalization.html</a></p>
<p>Systematic differences in sequencing coverage between libraries is a problem in single-cell RNA sequencing.</p>
<p>Causes:</p>
<ul>
<li>technical differences in cDNA capture</li>
<li>PCR amplification efficiency across cells.</li>
</ul>
<p>Normalization removes these differences so they dont interfere with comparisons between cells &gt; observed heterogeneity or differential expression within the cell population is driven by biology and not technical biases.</p>
<p>Difference between <strong>normalization</strong> and <strong>batch correction</strong>:</p>
<ul>
<li>Normalization occurs regardless of the batch structure and only <strong>considers technical biases</strong></li>
<li>Batch correction only occurs <strong>across batches</strong> and must consider both technical biases and biological differences.</li>
</ul>
<p>Technical biases affect genes in a similar manner (e.g., length, GC content).<br>
Biological differences between batches is highly unpredictable. Important to avoid conflating “normalized” and “batch-corrected” data, as these refer to different things.</p>
<p><strong>Scaling normalization</strong> is the simplest class of normalization strategies:</p>
<ul>
<li>Divide all counts for each cell by a cell-specific scaling factor = “<strong>size factor</strong>”  (Anders and Huber  <a href="https://osca.bioconductor.org/normalization.html#ref-anders2010differential">2010</a>).</li>
<li>Assumption = any cell-specific bias (e.g., capture or amplification efficiency) <strong>affects all genes equally</strong> via scaling of the expected mean count for that cell.</li>
<li>Size factor for each cell represents the estimate of the relative bias in that cell. Dividing a cells counts by its size factor removes that bias.</li>
<li>Resulting “<strong>normalized expression values</strong>” are used in clustering and dimensionality	 reduction.</li>
</ul>
<h2 id="library-size-normalization">1. Library size normalization</h2>
<p>Library size normalization is the simplest strategy for performing scaling normalization.<br>
Library size = the total sum of counts across all genes for each cell. Library size scales with any cell-specific biases.<br>
“Library size factor” for each cell is directly proportional to its library size.<br>
Proportionality constant ensures the mean size factor across all cells is equal to 1. The normalized expression values are then on the same scale as the original counts, which is useful for interpretation - especially when dealing with transformed data.</p>
<p>Library size factors assumes there is no “imbalance” in the differentially expressed genes between any pair of cells. Any upregulation for a subset of genes is cancelled out by the same magnitude of downregulation in a different subset of genes. This ensures that the library size is an unbiased estimate of the relative cell-specific bias by avoiding composition effects  (Robinson and Oshlack  <a href="https://osca.bioconductor.org/normalization.html#ref-robinson2010scaling">2010</a>). However, <strong>balanced DE is not generally present in scRNA-seq</strong>, so library size normalization may not yield accurate normalized expression values.</p>
<p>In practice, normalization accuracy is not a major consideration for exploratory scRNA-seq. Composition biases do not usually affect the separation of clusters, only the magnitude - and to a lesser extent, direction - of the log-fold changes between clusters or cell types.<br>
<strong>Library size normalization is usually sufficient</strong> where the aim is to identify clusters and the top markers that define each cluster.</p>
<h2 id="normalization-by-deconvolution-composition-bias">2. Normalization by deconvolution (Composition Bias)</h2>
<p>Composition biases are present when any unbalanced differential expression exists between samples.<br>
E.g two cells where a single gene  XX  is upregulated in one cell  AA  compared to the other cell  BB. This upregulation means that either (i) more sequencing resources are devoted to  XX  in  AA, thus decreasing coverage of all other non-DE genes when the total library size of each cell is experimentally fixed (e.g., due to library quantification); or (ii) the library size of  AA  increases when  XX  is assigned more reads or UMIs, increasing the library size factor and yielding smaller normalized expression values for all non-DE genes. In both cases, the <strong>net effect is that non-DE genes in  AA  will incorrectly appear to be downregulated compared to  BB.</strong></p>
<p>Removal of composition biases is well-studied in bulk RNAseq. Normalization is performed with <code>estimateSizeFactorsFromMatrix()</code>  in  <em><a href="https://bioconductor.org/packages/3.11/DESeq2">DESeq2</a></em>   or with  <code>calcNormFactors()</code>  in  <em><a href="https://bioconductor.org/packages/3.11/edgeR">edgeR</a></em> . These <strong>assume that most genes are not DE between cells.</strong> Systematic differences in counts across non-DE genes between two cells represents bias that is used to compute a size factor for its removal.</p>
<p>Single-cell data can be problematic for these bulk normalization methods due to the <strong>dominance of low and zero counts</strong>. To overcome this, we <strong>pool counts from many cells</strong> to increase the size of the counts for accurate size factor estimation  (Lun, Bach, and Marioni  <a href="https://osca.bioconductor.org/normalization.html#ref-lun2016pooling">2016</a>).  <strong>Pool-based size factors</strong> are then “<strong>deconvolved</strong>” into cell-based factors for normalization of each cell’s expression profile. This is performed using  <code>calculateSumFactors()</code>  in  <em><a href="https://bioconductor.org/packages/3.11/scran">scran</a></em>.</p>
<p>Pre-clustering step  <code>quickCluster()</code> : <strong>cells in each cluster are normalized separately</strong> and the size factors are rescaled to be comparable across clusters. This avoids the assumption that most genes are non-DE across the entire population - only a non-DE majority is required between pairs of clusters, which is a weaker assumption for highly heterogeneous populations. <code>quickCluster()</code>  uses an approximate algorithm for PCA relying on stochastic initialization so we need to set the random seed (via  <code>set.seed()</code>) for reproducibility.</p>
<p><strong>Deconvolution size factors</strong> exhibit cell type-specific deviations from the <strong>library size factors</strong>. This is due to composition biases from strong differential expression between cell types. Use of the deconvolution size factors adjusts for these biases to improve normalization accuracy for downstream applications.</p>
<p>Accurate normalization is most important for <strong>estimation and interpretation of per-gene statistics</strong>. Composition biases compromise DE analyses by systematically shifting the log-fold changes. However, <strong>library size normalization is best for cell-based analyses such as clustering</strong>. The presence of composition biases implies strong differences in expression profiles, so changing the normalization strategy is unlikely to affect the outcome of clustering.</p>
<h2 id="normalization-by-spike-ins">3. Normalization by spike-ins</h2>
<p>Spike-in normalization assumes that the same amount of spike-in RNA was added to each cell. Systematic differences in the coverage of the spike-in transcripts can only be due to cell-specific biases, e.g., in capture efficiency or sequencing depth. To remove these biases, we <strong>equalize spike-in coverage across cells by scaling with “spike-in size factors”.</strong> Spike-in normalization requires no assumption about the biology (i.e., the absence of many DE genes).</p>
<p>It assumes that the spike-in transcripts were (i) added at a constant level to each cell, and (ii) respond to biases in the same relative manner as endogenous genes.</p>
<p>Practically, spike-in normalization should be <strong>used if differences in the total RNA content of individual cells are of interest</strong>. An increase in a cells overall amount of endogenous RNA will not increase its spike-in size factor. This ensures that the effects of total RNA content on expression across the population are not removed upon scaling. Other normalization methods simply interpret any change in total RNA content as part of the bias and remove it.</p>
<p>Whether or not total RNA content is relevant – and the choice of normalization strategy – depends on the biological hypothesis. In most cases, changes in total RNA content are not interesting and can be normalized out by applying the library size or deconvolution factors. However, this is not the case with cell cycle activity or T cell activation. Spike-in normalization will preserve RNA content differences such that any changes in expression between biological groups have the correct sign.</p>
<p><strong>However!</strong>  Regardless of whether we care about total RNA content, it is <strong>critical that spike-in transcripts are normalized using spike-in size factors</strong>. Size factors computed from counts for endogenous genes should not be applied to the spike-in transcripts, because size factors captures differences in total RNA content that are not experienced by spike-ins. <strong>Normalizing spike-in counts with gene-based size factors leads to over-normalization and incorrect quantification.</strong></p>
<h2 id="applying-size-factors">Applying size factors</h2>
<h3 id="scaling-and-log-transforming">Scaling and log-transforming</h3>
<p>Once we have computed the size factors, we use the  <code>logNormCounts()</code>  function from  <em><a href="https://bioconductor.org/packages/3.11/scater">scater</a></em>  to compute normalized expression values for each cell. This is done by dividing the count for each gene/spike-in transcript with the appropriate size factor for that cell. The function also log-transforms the normalized values, creating a new assay called  <code>"logcounts"</code>. (Technically, these are “log-transformed normalized expression values”, but that’s too much of a mouthful to fit into the assay name.) These log-values will be the basis of our downstream analyses in the following chapters.</p>
<pre><code>set.seed(100)
clust.zeisel &lt;- quickCluster(sce.zeisel) 
sce.zeisel &lt;- computeSumFactors(sce.zeisel, cluster=clust.zeisel, min.mean=0.1)
sce.zeisel &lt;- logNormCounts(sce.zeisel)
assayNames(sce.zeisel)
</code></pre>
<pre><code>## [1] "counts"    "logcounts"
</code></pre>
<p>The log-transformation is useful as differences in the log-values represent log-fold changes in expression. This is important in downstream procedures based on Euclidean distances, which includes many forms of clustering and dimensionality reduction. By operating on log-transformed data, we ensure that these procedures are measuring distances between cells based on log-fold changes in expression. Or in other words, which is more interesting - a gene that is expressed at an average count of 50 in cell type  AA  and 10 in cell type  BB, or a gene that is expressed at an average count of 1100 in  AA  and 1000 in  BB? Log-transformation focuses on the former by promoting contributions from genes with strong relative differences.</p>
<p>When log-transforming, we typically add a pseudo-count to avoid undefined values at zero. Larger pseudo-counts will effectively shrink the log-fold changes between cells towards zero for low-abundance genes, meaning that downstream high-dimensional analyses will be driven more by differences in expression for high-abundance genes. Conversely, smaller pseudo-counts will increase the relative contribution of low-abundance genes. Common practice is to use a pseudo-count of 1, for the simple pragmatic reason that it preserves sparsity in the original matrix (i.e., zeroes in the input remain zeroes after transformation). This works well in all but the most pathological scenarios  (A. Lun  <a href="https://osca.bioconductor.org/normalization.html#ref-lun2018overcoming">2018</a>).</p>
<p>Incidentally, the addition of the pseudo-count is the motivation for the centering of the size factors at unity. This ensures that both the pseudo-count and the normalized expression values are on the same scale; a pseudo-count of 1 can be interpreted as an extra read or UMI for each gene. In practical terms, centering means that the shrinkage effect of the pseudo-count diminishes as sequencing depth improves. This correctly ensures that estimates of the log-fold change in expression (e.g., from differences in the log-values between groups of cells) become increasingly accurate with deeper coverage. In contrast, if we applied a constant pseudo-count to some count-per-million-like measure, accuracy of the subsequent log-fold changes would never improve regardless of how much additional sequencing we performed.</p>
<h3 id="downsampling-and-log-transforming">7.5.2  Downsampling and log-transforming</h3>
<p>In rare cases, direct scaling of the counts is not appropriate due to the effect described by  A. Lun (<a href="https://osca.bioconductor.org/normalization.html#ref-lun2018overcoming">2018</a>). Briefly, this is caused by the fact that the mean of the log-normalized counts is not the same as the log-transformed mean of the normalized counts. The difference between them depends on the mean and variance of the original counts, such that there is a systematic trend in the mean of the log-counts with respect to the count size. This typically manifests as trajectories correlated strongly with library size even after library size normalization, as shown in Figure  <a href="https://osca.bioconductor.org/normalization.html#fig:cellbench-lognorm-fail">7.5</a>  for synthetic scRNA-seq data generated with a pool-and-split approach  (Tian et al.  <a href="https://osca.bioconductor.org/normalization.html#ref-tian2019benchmarking">2019</a>).</p>
<pre><code># TODO: move to scRNAseq.
library(BiocFileCache)
bfc &lt;- BiocFileCache(ask=FALSE)
qcdata &lt;- bfcrpath(bfc, "https://github.com/LuyiTian/CellBench_data/blob/master/data/mRNAmix_qc.RData?raw=true")

env &lt;- new.env()
load(qcdata, envir=env)
sce.8qc &lt;- env$sce8_qc

# Library size normalization and log-transformation.
sce.8qc &lt;- logNormCounts(sce.8qc)
sce.8qc &lt;- runPCA(sce.8qc)
gridExtra::grid.arrange(
    plotPCA(sce.8qc, colour_by=I(factor(sce.8qc$mix))),
    plotPCA(sce.8qc, colour_by=I(librarySizeFactors(sce.8qc))),
    ncol=2
)
</code></pre>
<p><img src="https://osca.bioconductor.org/P2_W03.normalization_files/figure-html/cellbench-lognorm-fail-1.png" alt="PCA plot of all pool-and-split libraries in the SORT-seq CellBench data, computed from the log-normalized expression values with library size-derived size factors. Each point represents a library and is colored by the mixing ratio used to construct it (left) or by the size factor (right)."></p>
<p>Figure 7.5: PCA plot of all pool-and-split libraries in the SORT-seq CellBench data, computed from the log-normalized expression values with library size-derived size factors. Each point represents a library and is colored by the mixing ratio used to construct it (left) or by the size factor (right).</p>
<p>As the problem arises from differences in the sizes of the counts, the most straightforward solution is to downsample the counts of the high-coverage cells to match those of low-coverage cells. This uses the size factors to determine the amount of downsampling for each cell required to reach the 1st percentile of size factors. (The small minority of cells with smaller size factors are simply scaled up. We do not attempt to downsample to the smallest size factor, as this would result in excessive loss of information for one aberrant cell with very low size factors.) We can see that this eliminates the library size factor-associated trajectories from the first two PCs, improving resolution of the known differences based on mixing ratios (Figure  <a href="https://osca.bioconductor.org/normalization.html#fig:cellbench-lognorm-downsample">7.6</a>). The log-transformation is still necessary but no longer introduces a shift in the means when the sizes of the counts are similar across cells.</p>
<pre><code>sce.8qc2 &lt;- logNormCounts(sce.8qc, downsample=TRUE)
sce.8qc2 &lt;- runPCA(sce.8qc2)
gridExtra::grid.arrange(
    plotPCA(sce.8qc2, colour_by=I(factor(sce.8qc2$mix))),
    plotPCA(sce.8qc2, colour_by=I(librarySizeFactors(sce.8qc2))),
    ncol=2
)
</code></pre>
<p><img src="https://osca.bioconductor.org/P2_W03.normalization_files/figure-html/cellbench-lognorm-downsample-1.png" alt="PCA plot of pool-and-split libraries in the SORT-seq CellBench data, computed from the log-transformed counts after downsampling in proportion to the library size factors. Each point represents a library and is colored by the mixing ratio used to construct it (left) or by the size factor (right)."></p>
<p>Figure 7.6: PCA plot of pool-and-split libraries in the SORT-seq CellBench data, computed from the log-transformed counts after downsampling in proportion to the library size factors. Each point represents a library and is colored by the mixing ratio used to construct it (left) or by the size factor (right).</p>
<p>While downsampling is an expedient solution, it is statistically inefficient as it needs to increase the noise of high-coverage cells in order to avoid differences with low-coverage cells. It is also slower than simple scaling. Thus, we would only recommend using this approach after an initial analysis with scaled counts reveals suspicious trajectories that are strongly correlated with the size factors. In such cases, it is a simple matter to re-normalize by downsampling to determine whether the trajectory is an artifact of the log-transformation.</p>

