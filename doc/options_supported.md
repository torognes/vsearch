# Supported command line options in VSEARCH

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/usearch_algo.html)

Options related to chimera filtering not included. Other options may also be missing.

Both single (-) and double dash (--) in front of options are allowed.

Defaults are indicated in parentheses.

Command line example:

`./vsearch-0.0.6-macosx-x86_64 --usearch_global q1000.fas --db big_file.fas --strand plus --id 0.9 --alnout alnout.aln`


## Basic options

- [x] `--usearch_global <filename>`

- [x] `--db <filename>`

- [x] `--strand <plus|both>`

- [x] `--threads <integer>` (Option recognized but has no effect)

- [ ] `--cluster_smallmem <filename>`

- [ ] `--cluster_fast <filename>`

- [x] `--sizein <int>`


## Output options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/output_files.html)

- [x] `--alnout <filename>`

- [x] `--rowlen <integer(64)>`

- [x] `--userout <filename>`

- [ ] `--fastapairs <filename>`

- [x] `--blast6out <filename>`

- [x] `--matched <filename>`

- [x] `--notmatched <filename>`

- [x] `--dbmatched <filename>`

- [x] `--dbnotmatched <filename>`

- [ ] `--maxhits <integer (0)>`

- [ ] `--top_hits_only`

- [ ] `--output_no_hits`

- [x] `--sizeout`

- [ ] `--centroids <filename>`

- [ ] `--consout <filename>`

- [ ] `--clusters <filename>`

- [ ] `--msaout <filename>`

- [x] `--uc <filename>`

- [x] `--output <filename>`


## Accept options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/accept_options.html)

- [x] `--id <real>`

- [ ] `--evalue <real>`

- [ ] `--query_cov <real>`

- [ ] `--target_cov <real>`

- [ ] `--idprefix <string>`

- [ ] `--idsuffix <string>`

- [ ] `--minqt <real>`

- [ ] `--maxqt <real>`

- [ ] `--minsl <real>`

- [ ] `--maxsl <real>`

- [ ] `--leftjust`

- [ ] `--rightjust`

- [x] `--self`

- [ ] `--selfid`

- [ ] `--maxid <real>`

- [ ] `--minsizeratio <real>`

- [ ] `--maxsizeratio <real>`

- [ ] `--maxdiffs <integer>`

- [ ] `--maxsubs <integer>`

- [ ] `--maxgaps <integer>`

- [ ] `--mincols <integer>`

- [ ] `--maxqsize <integer>`

- [ ] `--mintsize <integer>`

- [ ] `--mid <real>`


## Termination options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/termination_options.html)

- [x] `--maxaccepts <integer (1)>`

- [x] `--maxrejects <integer (32)>`


## Weak hits options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/weak_hits.html)

- [x] `--weak_id <real>`

- [ ] `--weak_evalue <real>`


## Indexing options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/indexing_options.html)

- [x] `--wordlength <integer (8)>`

- [ ] `--pattern <bitpattern (11111111)>`

- [ ] `--alpha <nt|aa>`

- [ ] `--dbstep <integer (1)>`

- [ ] `--dbaccelpct <real (1.0)>`

- [ ] `--dbmask <method (fastnucleo)>`

- [ ] `--slots <integer (65536)>`


## Masking options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/masking_options.html)

- [ ] `--qmask <fastamino|fastnucleo|seg|dust|none|soft>`

- [ ] `--dbmask <fastamino|fastnucleo|seg|dust|none|soft>`

- [ ] `--hardmask`


## Alignment parameters

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/aln_params.html)

Only integer scores and penalties are allowed. The defaults are multiplied by two to compensate for this.

- [x] `--match <integer (2)>`

- [x] `--mismatch <integer (-4)>`

- [x] `--gapopen <integer[LRIEQT]/... (20I/2E)>`

- [x] `--gapext <integer[LRIEQT]/... (2I/1E)>`

- [ ] `--matrix <filename (blosum62 (aa), +1/-2 (nt))>`

- [ ] `--lopen <real (10.0)`

- [ ] `--lext <real (1.0)`


## Alignment heuristics options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/aln_heuristics.html)

- [x] `--fulldp`

- [ ] `--hspw <integer (5 (nt), 3 (aa))>`

- [ ] `--xdrop_u <real (16.0)>`

- [ ] `--xdrop_g <real (32.0)>`

- [ ] `--xdrop_nw <real (16.0)>`

- [ ] `--minhsp <integer (16)>`

- [ ] `--band <integer (16)>`


## Karlin Altschul statistics parameters

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/karlin_altschul.html)

- [ ] `--ka_gapped_lambda <real>`

- [ ] `--ka_ungapped_lambda <real>`

- [ ] `--ka_gapped_k <real>`

- [ ] `--ka_ungapped_k <real>`

- [ ] `--ka_dbsize <integer>`


## Special options

[Link to official USEARCH documentation](http://www.drive5.com/usearch/manual/opt_maxseqlength.html)

- [x] `--maxseqlength <integer>`


## Clustering options

- [ ] `--cons_truncate`


## Other options

- [x] `--derep_fulllength`

- [x] `--sortbysize`

- [x] `--sortbylength`

- [x] `--userfields <string>`

- [x] `--uc_allhits`

- [x] `--rowlen <integer>`

- [x] `--notrunclabels`

- [x] `--minuniquesize <integer>`

- [x] `--minsize <integer>`

- [x] `--maxsize <integer>`

- [x] `--minseqlength <integer>`

- [x] `--maxseqlength <integer>`

- [x] `--relabel <string>`

- [x] `--sizeout`

- [x] `--topn <integer>`
