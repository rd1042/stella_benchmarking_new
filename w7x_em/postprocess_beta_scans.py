""" """

from save_pickles_from_stella_scans import postprocess_beta_scan
import make_beta_scans as make_scans

if __name__ == "__main__":
    print("Hello world")
    # postprocess_beta_scan(make_scans.folder_name_dkh_explicit_betalr)
    postprocess_beta_scan(make_scans.folder_name_dkh_explicit_betalr_fbpar0)
