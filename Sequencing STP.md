---


---

<blockquote>
<h1 id="sequencing-stp">Sequencing STP</h1>
</blockquote>
<h1 id="step-1-cell-counter">Step 1: Cell Counter</h1>
<p>Sequencing library is provided to Milly at the Sequencing STP in an Eppendorf Tube. That sample is treated as a single batch, irrespective of multiple sequencing runs.</p>
<p>The sample is run through an Automated Cell Counter which provides:</p>
<ol>
<li>Viability (trichome stain) - A minimum viability of 70% is required to proceed to sequencing.</li>
<li>Cell Count - ideally 10,000 cells per sample. Expect ~40-50% cells to be recovered in sequencing (i.e. ~ 4,500 cells)</li>
</ol>
<h1 id="step-2-10x-preparation">Step 2: 10x preparation</h1>
<p>Sample is loaded into gel beads with a lysis reagent &amp; emulsion.<br>
UMI library barcodes are added<br>
PCR amplification<br>
QC step provides:</p>
<ol>
<li>Library size</li>
<li>Molarity (molecule number)<br>
Pooling step creates a library concentration of ideally 4 nM</li>
</ol>
<h1 id="step-3-sequencing">Step 3: Sequencing</h1>
<p>Cost is ~£2,500 per sample</p>
<p>Initial trial sequencing run at LOW sequencing depth:</p>
<ul>
<li>gives proxy for total cell number</li>
<li>it is a poor proxy if there are lots of dead cells as the ambient RNA is falsely counted as cells (cell number is overestimated)</li>
</ul>
<p>A second sequencing run is performed to achieve ideally 50,000 reads per cell. If this isnt achieved then subsequent sequencing runs are performed.</p>
<h1 id="step-4-cell-ranger-qc">Step 4: Cell Ranger QC</h1>
<p>The sequencing STP runs Cell Ranger on each run separately. This is a QC measure for them to give cell number &amp; reads per cell.<br>
This should not be used for downstream analysis as Cell Ranger <code>count</code> needs to be run on all fastq files from the same sample (i.e. merge all sequencing runs since they are from the same batch).</p>
<h1 id="step-5-fastq-output">Step 5: Fastq output</h1>
<p>FASTQ files are put into a folder on CAMP e.g<br>
/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190806_K00102_0374_AH5KN3BBXY/fastq/TAH421A1</p>
<p>SC19137 = project ID<br>
190806_K00102_0374_AH5KN3BBXY = sequencing run ID<br>
TAH421A1 = library ID</p>
<p>The Fastq output has the following structure:<br>
TAH421A1_S4_L001_I1_001.fastq.gz</p>
<p>TAH421A1 = Library ID<br>
S4 = Position for library demultiplexing (used by STP but irrelevent for downstream analysis)<br>
L001 = Lane ID<br>
R1 = Read Pair (R1, R2). Also I1 - this is very short nucleotide sequencing needed for cell ranger demultiplexing.<br>
001.fastq.gz = ending on all fastq files</p>
<p>Cell Ranger <code>count</code> should then be run on all fastq files with the same Library ID in a single step.</p>

