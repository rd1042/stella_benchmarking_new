""" """

import re
import numpy as np
import sys

def compare_golder():
    """ """
    def extract_golder(txtfile):
        """ """
        myfile=open(txtfile, "r")
        filetext = myfile.read()
        myfile.close()
        entries = re.split("\s+", filetext.strip())
        golder_list = []
        for entry in entries:
            try:
                # Remove the brackets
                if entry[0] == "(" and entry[-1] ==")":
                    entry=entry[1:-1]
                else:
                    print("expected brackets not present. entry = ", entry)
                    sys.exit()
                reim = re.split(",", entry)
                golder_list.append(float(reim[0]) + 1j*float(reim[1]))
            except ValueError:
                print("ERROR!")
                print("entry = ", entry)
                sys.exit()
        return golder_list

    nisl_txtfile = "sims/nisl/nisl_tmp_3.txt"
    leapfrog_txtfile = "sims/leapfrog/leapfrog_tmp_3.txt"

    golder_nisl = extract_golder(nisl_txtfile)
    golder_leapfrog = extract_golder(leapfrog_txtfile)

    golder_diff = np.array(golder_nisl) - np.array(golder_leapfrog)
    print("golder_nisl = ", golder_nisl[:10])
    print("golder_leapfrog = ", golder_leapfrog[:10])
    print("golder_diff = ", golder_diff)
    print("max(abs(golder_diff)) = ", np.max(abs(golder_diff)))


    return

if __name__ == "__main__":
    compare_golder()
