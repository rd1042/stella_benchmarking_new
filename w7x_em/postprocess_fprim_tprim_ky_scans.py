""" """

from save_pickles_from_stella_scans import postprocess_fprim_tprim_ky_scan
import make_fprim_tprim_ky_scans as make_scans

if __name__ == "__main__":
    print("Hello world")
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_1)
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_2)
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_expl_higher_ky_em)
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_expl_higher_ky_es)
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_impl_higher_ky_em)
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_impl_higher_ky_es)

    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_impl_em_good_resolution)
    # postprocess_fprim_tprim_ky_scan(make_scans.folder_name_impl_em_good_resolution_beta001)
    postprocess_fprim_tprim_ky_scan(make_scans.folder_name_impl_em_good_resolution_beta003)
