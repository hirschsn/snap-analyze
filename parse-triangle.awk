#!/usr/bin/awk -f
# For license details see LICENSE.
#
# Parses bond information from sana-bond and calculates a histogram
# of triangles that have at least one broken angular bond.

# "(x)" -> "x"
function unparen(s,   spl) { # spl is local
    split(s, spl, "[\(\)]")
    return spl[2]
}

BEGIN {
    fmtstr[0] = "%i %lf%% %i %i %lf %lf\n"
    fmtstr[1] = "%3i %5.1lf%% (%3i/%3i) Average bond length: %lf Average deriv. from phi_0: %lf\n"
    # Use -v formatted_output=1 to view output formatted.
    if (!formatted_output)
        formatted_output = 0

    # Use -v "broken" to display information about broken bonds
    if (!wanted)
        wanted = "valid"
}

{
    state = $15
    a[state]++

    bond_id = unparen($4)

    if (state == wanted) {
        v[bond_id]++
        len[bond_id] += unparen($16)
        # phi - phi_0
        abw[bond_id] += $9
    }
    t[bond_id]++
}

END {
    for (i in a)
        printf("%10s %i\n", i, a[i]);

    printf("\nHistogram of still angular bonds for which the opposite pair bond is %s:\n", wanted)
    for (i in v) {
        # This if is to easily allow "for (i in t)" to be printed
        if (v[i] > 0) {
            avglen = len[i] / v[i]
            avgabw = abw[i] / v[i]
        } else {
            avglen = 0.0
            avgabw = 0.0
        }
        printf(fmtstr[formatted_output], i, v[i]/t[i]*100, v[i], t[i], avglen, avgabw);
    }
}

