""" """

from save_pickles_from_stella_scans import postprocess_fprim_tprim_ky_scan
import make_fprim_tprim_ky_scans as make_scans

if __name__ == "__main__":
    print("Hello world")
    postprocess_fprim_tprim_ky_scan(make_scans.folder_name_1)
    postprocess_fprim_tprim_ky_scan(make_scans.folder_name_2)
