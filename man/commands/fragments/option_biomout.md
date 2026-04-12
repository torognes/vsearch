`--biomout` *filename*
: Write an OTU table to *filename* in the biom version 1.0 JSON file
  format. The OTUs are represented by the cluster centroids. Sample
  identifiers are extracted from sequence headers (`;sample=abc123;`
  or `;barcodelabel=abc123;` patterns, or the initial part of the
  header). OTU identifiers are extracted from centroid headers
  (`;otu=def789;` pattern, or the initial part of the header, or via
  relabelling options). Taxonomy information is extracted from centroid
  headers (`;tax=...;` pattern) if available.
