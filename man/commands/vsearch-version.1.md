% vsearch-version(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-version --- write version information


# SYNOPSIS

| **vsearch** **\-v|\-\-version**


# DESCRIPTION

The vsearch command `--version` writes version information to the
*standard error* `stderr(3)`. A citation for the vsearch publication,
and the support status for gzip- and bzip2-compressed input files are
written to the *standard output* `stdout(3)`.


# OPTIONS

## mandatory options

None


## core options

None


## secondary options

#(./fragments/option_log.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Use the `--version` command to get vsearch's version number:
:

```sh
vsearch \
    --version 2>&1 | \
    grep -oE "v[0-9]+[.][0-9]+[.][0-9]+"
```


# SEE ALSO

[`vsearch(1)`](./vsearch.1.md)

#(./fragments/footer.md)
