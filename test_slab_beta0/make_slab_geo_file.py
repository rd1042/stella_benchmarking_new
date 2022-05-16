""" """


import re
import numpy as np

if __name__ == "__main__":
    print("Hello world!")
    slab_geo_template = "stella_sims/input_cbc.geometry"
    slab_geo_outfile = "stella_sims/unsheared_slab.geometry"

    template_file = open(slab_geo_template, "r")
    template_lines = template_file.readlines()
    template_file.close()

    slab_geo_text=""
    nchars = 12 # each value string must be exactly 12 characters long.
    # The first few lines are the same
    for template_line in template_lines[0:4]:
        slab_geo_text += template_line

    # Go through all the other lines, except the last line which is empty
    for template_line in template_lines[4:-1]:
        values = re.split("\s+", template_line.strip())
        # values = [alpha, zed, zeta, bmag, gradpar, gds2, gds21, gds22, gds23, gds24, gbdrift, cvdrift, gbdrift0, bmag_psi0, btor]
        # Because of how fortran reads the file, each value string must be exactly 12 characters long.
        slab_geo_line=values[0].rjust(nchars) + values[1].rjust(nchars) + values[2].rjust(nchars)
        # Replace bmag, which is values[3]
        slab_geo_line+= "0.1000E+01".rjust(nchars)
        # gradpar unchanged
        slab_geo_line+= values[4].rjust(nchars)
        ## Replace gd21, gds22, gds23, which are values[5, 6, 7]
        ## See E. Highcock thesis for these quantities
        # gds2 = 1
        slab_geo_line+= "0.1000E+01".rjust(nchars)
        # gds22 = gds23 = 0
        slab_geo_line+= "0.0000E+00".rjust(nchars) + "0.0000E+00".rjust(nchars)
        # Not sure about gds23, gds24 (values[8, 9]), but only used if neoclassical terms included
        slab_geo_line+= values[8].rjust(nchars) + values[9].rjust(nchars)
        # Replace gbdrift, cvdrift, gbdrift0, which are values[10, 11, 12]
        slab_geo_line+= "0.0000E+00".rjust(nchars) + "0.0000E+00".rjust(nchars) + "0.0000E+00".rjust(nchars)
        # stella doesn't actually try to read bmag_psi0, btor, so don't bother specifying them
        slab_geo_line+="\n"
        slab_geo_text+=slab_geo_line

    slab_geo_file = open(slab_geo_outfile, "w")
    slab_geo_file.write(slab_geo_text)
    slab_geo_file.close()
