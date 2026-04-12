`--uchimeout` *filename*
: Write chimera detection results to *filename* using an 18-field,
  tab-separated uchime-like format. Output row order may vary when
  using multiple threads. Use `--uchimeout5` for a format compatible
  with usearch version 5 and earlier. The 18 fields are:

    1.  score: higher score means a more likely chimeric alignment.
    2.  Q: query sequence label.
    3.  A: parent A sequence label.
    4.  B: parent B sequence label.
    5.  T: top parent sequence label (parent most similar to the query).
    6.  idQM: percentage of similarity between query (Q) and the model
        (M) constructed as a part of parent A and a part of parent B.
    7.  idQA: percentage of similarity between query (Q) and parent A.
    8.  idQB: percentage of similarity between query (Q) and parent B.
    9.  idAB: percentage of similarity between parent A and parent B.
    10. idQT: percentage of similarity between query (Q) and top parent (T).
    11. LY: yes votes in the left part of the model.
    12. LN: no votes in the left part of the model.
    13. LA: abstain votes in the left part of the model.
    14. RY: yes votes in the right part of the model.
    15. RN: no votes in the right part of the model.
    16. RA: abstain votes in the right part of the model.
    17. div: divergence, defined as (idQM - idQT).
    18. YN: query is chimeric (Y), not chimeric (N), or borderline (?).
