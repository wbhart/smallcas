# Mathematical algorithm files are checked here.
# Dispatch files and memory.c are exempt by project rule.
# Blank lines and lines containing only a closing brace do not count.
# The limit is 59 substantive lines, set by direct subresultant xgcd propagation.

function limited(line, file)
{
    if (file ~ /\/zz\.c$/ && line ~ /^sc_value \*sc_zz_/)
        return 1
    if (file ~ /\/zz_poly\.c$/ &&
        (line ~ /^sc_value \*sc_zz_poly_/ || line ~ /^int sc_zz_poly_/))
        return 1
    if (file ~ /\/poly\.c$/ && line ~ /^sc_value \*sc_polynomial_ring\(/)
        return 1
    return 0
}

function function_name(line, name)
{
    name = line
    if (name ~ /^sc_value \*/ )
        sub(/^sc_value \*/, "", name)
    else if (name ~ /^int /)
        sub(/^int /, "", name)
    else
        sub(/^.*[[:space:]]+/, "", name)
    sub(/^\*/, "", name)
    sub(/\(.*/, "", name)
    return name
}

function counted(line)
{
    if (line ~ /^[[:space:]]*$/)
        return 0
    if (line ~ /^[[:space:]]*}[[:space:]]*$/)
        return 0
    return 1
}

limited($0, FILENAME) {
    active = 1
    start = FNR
    name = function_name($0)
    depth = 0
    seen = 0
    lines = 0
}

active {
    if (counted($0))
        lines++
    text = $0
    opens = gsub(/{/, "{", text)
    closes = gsub(/}/, "}", text)
    depth += opens - closes
    if (opens > 0)
        seen = 1
    if (seen && depth == 0) {
        if (lines > 59) {
            print FILENAME ":" start ": " name " has " lines " counted lines"
            bad = 1
        }
        active = 0
    }
}

END {
    exit bad
}
