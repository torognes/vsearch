`--tabbedout` *filename*
: Write chimera detection results to *filename* as an eighteen-column
  tab-separated file, with one row per chimera. Columns are:

    1.  score: dummy value, always set to 99.9999
    2.  query label
    3.  parent A label
    4.  parent B label
    5.  parent C label ("*" if there are only two parents)
    6.  QModel: maximum global similarity percentage (always 100.0%)
    7.  QA: global similarity percentage with parent A
    8.  QB: global similarity percentage with parent B
    9.  QC: global similarity percentage with parent C (0.00 if only two parents)
    10. QT: highest similarity percentage with any parent
    11. left yes: ignored, always set to zero
    12. left no: ignored, always set to zero
    13. left abstain: ignored, always set to zero
    14. right yes: ignored, always set to zero
    15. right no: ignored, always set to zero
    16. right abstain: ignored, always set to zero
    17. dummy value, always set to 0.00
    18. chimeric status, always set to Y (only chimeras are reported)
