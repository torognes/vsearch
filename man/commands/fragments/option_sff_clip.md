`--sff_clip`
: Clip sequences at both ends as indicated by the clipping coordinates
  stored in the SFF file (`clip_qual_left`, `clip_qual_right`,
  `clip_adapter_left`, `clip_adapter_right`). Without this option, no
  clipping is performed and bases that would have been clipped are
  written in lower case, while the remaining bases are in upper case.
