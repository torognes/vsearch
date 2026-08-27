# Expand the "#(path)" include directives found in the markdown manual
# sources: such a line is replaced by the contents of the file it names
# (paths are relative to the including file, hence the callers change
# directory first).  Anything following the closing parenthesis is
# dropped, and, as with the perl one-liner this replaces, the included
# text is followed by an empty line -- the newline of the directive
# line itself.  A path that cannot be read expands to nothing.
#
# awk is already required to build vsearch (automake uses it), perl is
# not, and perl is absent from the base system of some of the platforms
# we support.  Every markdown source ends with a newline, which is what
# lets the concatenation below reproduce the file byte for byte.

function expand(path,   buffer, line) {
    buffer = ""
    while ((getline line < path) > 0)
        buffer = buffer line "\n"
    close(path)
    return buffer
}

{
    # the directive must open the line, and the path must not be empty
    if (substr($0, 1, 2) == "#(") {
        closing = 0
        for (i = length($0); i > 3; i--)
            if (substr($0, i, 1) == ")") {
                closing = i
                break
            }
        if (closing > 3) {
            printf "%s\n", expand(substr($0, 3, closing - 3))
            next
        }
    }
    print
}
