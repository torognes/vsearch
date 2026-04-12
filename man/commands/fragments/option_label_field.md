`--label_field` *string*
: Restrict word matching (with `--label_word` or `--label_words`) to a
  named field in the header. A field is a semicolon-delimited
  key=value pair in the header. For example, with
  `--label_field "abc"` and `--label_word "123"`, the header
  `seq1;abc=123` matches because the field `abc` equals `123`.
