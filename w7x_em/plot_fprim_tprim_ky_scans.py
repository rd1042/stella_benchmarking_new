""" """

import numpy as np
import pickle
import matplotlib.pyplot as plt
from plot_2d_utils import uniquearrays2meshgrids
import math
import sys

import make_fprim_tprim_ky_scans as make_scans

def for_thesis_make_fprim_tprim_ky_scan_kjm3():
    """For fancy plotting """

    def make_plots(freq_ax, gamma_ax, freq_cbax, gamma_cbax, fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid,
                   gammaom_fprim_tprim_meshgrid, freq_levels, gamma_levels, freq_ticks, gamma_ticks
                   ):
        """ """
        # gamma_levels=20
        # freq_levels=20
        gammaom_contours = gamma_ax.contourf(fprim_mesh, tprim_mesh, gammaom_fprim_tprim_meshgrid, gamma_levels, cmap="inferno", extend="max")
        # freq_lim = math.ceil(max(abs(freq_fprim_tprim_meshgrid.min()), freq_fprim_tprim_meshgrid.max())*10)/10
        # freq_spacing = 0.1 # 0.02
        # freq_levels = np.arange(-freq_lim, freq_lim, freq_spacing)
        freq_contours = freq_ax.contourf(fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid, freq_levels, cmap="PiYG", extend="both")

        fig.colorbar(gammaom_contours, cax=gamma_cbax, ticks=gamma_ticks )
        fig.colorbar(freq_contours, cax=freq_cbax, ticks=freq_ticks)


        return

    def init_plots():

        fig = plt.figure(figsize=[13, 8])
        frac = (13.0/8.0)
        # Axis placement parameters.
        vspace_min = 0.1
        ax_hspace_min = 0.08
        cb_hspace = 0.01
        plot_dim_y_try = 0.38
        plot_dim_x_try = plot_dim_y_try/frac
        right = 0.93
        left = 0.05
        top = 0.99
        bottom = 0.1
        cb_width = 0.015
        plot_width_with_min_spacing = (right - left - 2*ax_hspace_min - 3*cb_hspace - 3*cb_width)/3
        plot_dim_x = plot_width_with_min_spacing
        plot_dim_y = plot_dim_x * frac

        col1_lhs = left
        col1_cb_lhs = left + plot_dim_x + cb_hspace
        col2_lhs = col1_cb_lhs + cb_width + ax_hspace_min
        col2_cb_lhs = col2_lhs + plot_dim_x + cb_hspace
        col3_lhs = col2_cb_lhs + cb_width + ax_hspace_min
        col3_cb_lhs = col3_lhs + plot_dim_x + cb_hspace
        # col2_cb_lhs = 0.94
        #row3_bottom = 0.08
        row2_bottom = bottom # row3_bottom + plot_dim_y + vspace
        row1_bottom = row2_bottom + plot_dim_y + vspace_min


        ########################################################################
        # Create axes.
        # [ax1, ax2] = [gammaom, freq]
        # [ax3, ax4] = [gammap2, phase velocity]
        ########################################################################
        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax2 = fig.add_axes([col1_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax3 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax4 = fig.add_axes([col2_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax5 = fig.add_axes([col3_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax6 = fig.add_axes([col3_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax2 = fig.add_axes([col1_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax3 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax4 = fig.add_axes([col2_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax5 = fig.add_axes([col3_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax6 = fig.add_axes([col3_cb_lhs, row2_bottom, cb_width, plot_dim_y])

        return (fig, ax1, ax2, ax3, ax4, ax5, ax6, cbax1, cbax2, cbax3, cbax4, cbax5, cbax6)

    def finish_plots(freq_ticklabels_list, gamma_ticklabels_list):
        # Formatting, axis labels etc.

        x_ticklabelfontsize = 14
        x_labelfontsize = 20
        my_xticklength = 3
        my_xtickwidth = 1
        y_labelfontsize = x_labelfontsize
        y_ticklabelfontsize = x_ticklabelfontsize
        cb_ticklabelsize = 14
        cbax_title_fontsize = 20
        ax_title_fontsize = 20

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
            #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
            ax.tick_params("x", which="major", top=True, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.set_xlim([np.min(unique_fprim), np.max(unique_fprim)])
            ax.set_ylim([np.min(unique_tprim), np.max(unique_tprim)])
            ax.set_xlabel(r"$a/L_n$", fontsize=x_labelfontsize)
            ax.set_ylabel(r"$a/L_T$", fontsize=y_labelfontsize)
            ax.set_xticks([0, 2, 4, 6, 8])
            ax.set_xticklabels(["0", "2", "4", "6", "8"], fontsize=x_ticklabelfontsize)
            ax.set_yticks([0, 2, 4, 6, 8])
            ax.set_yticklabels(["0", "2", "4", "6", "8"], fontsize=y_ticklabelfontsize)

        ax1.set_title(r"$\tilde{k}_y=1.5$", fontsize=ax_title_fontsize)
        ax3.set_title(r"$\tilde{k}_y=3$", fontsize=ax_title_fontsize)
        ax5.set_title(r"$\tilde{k}_y=4.5$", fontsize=ax_title_fontsize)

        for cbax in [cbax1, cbax3, cbax5]:
            cbax.set_title(r"$\omega$", fontsize=cbax_title_fontsize)
        for cbax in [cbax2, cbax4, cbax6]:
            cbax.set_title(r"$\gamma$", fontsize=cbax_title_fontsize)

        #plt.show()

        cbax1.set_yticklabels(freq_ticklabels_list[0], fontsize=cb_ticklabelsize)
        cbax3.set_yticklabels(freq_ticklabels_list[1], fontsize=cb_ticklabelsize)
        cbax5.set_yticklabels(freq_ticklabels_list[2], fontsize=cb_ticklabelsize)
        cbax2.set_yticklabels(gamma_ticklabels_list[0], fontsize=cb_ticklabelsize)
        cbax4.set_yticklabels(gamma_ticklabels_list[1], fontsize=cb_ticklabelsize)
        cbax6.set_yticklabels(gamma_ticklabels_list[2], fontsize=cb_ticklabelsize)
        # Uncomment to save the figure.
        plt.savefig(str("images/kjm3_fprim_tprim_beta0.01") + ".png")
        plt.close()

        return


    folder_1 = make_scans.folder_name_impl_em_good_resolution_beta001
    pickle_longname = folder_1 + "/omega_fprim_tprim_ky_array.pickle"
    file = open(pickle_longname, "rb")
    [pickle_string, unique_fprim, unique_tprim, unique_ky,
        gammaom_fprim_tprim_ky_array, gammasafe_fprim_tprim_ky_array,
        freq_fprim_tprim_ky_array] = pickle.load(file)
    file.close()
    print("unique_ky = ", unique_ky)
    freq_fprim_tprim_ky_array = -freq_fprim_tprim_ky_array # Swap sign because of equilibrium
    (fig, ax1, ax2, ax3, ax4, ax5, ax6, cbax1, cbax2, cbax3, cbax4, cbax5, cbax6) = init_plots()
    ky_15_idx = unique_ky.index(1.5)
    ky_3_idx = unique_ky.index(3)
    ky_45_idx = unique_ky.index(4.5)

    freq_axes = [ax1, ax3, ax5]
    freq_cbaxes = [cbax1, cbax3, cbax5]
    gamma_axes = [ax2, ax4, ax6]
    gamma_cbaxes = [cbax2, cbax4, cbax6]

    nfreq = 50
    ngamma = 50
    freq_levels_15 = np.linspace(-2.5, 2.5, nfreq)
    freq_ticks_15 = [-2, -1, 0, 1, 2]
    freq_ticklabels_15 = ["-2", "-1", "0", "1", "2"]
    freq_levels_3 = np.linspace(-10, 10, nfreq)
    freq_ticks_3 = [-10, -5, 0, 5, 10]
    freq_ticklabels_3 = ["-10", "-5", "0", "5", "10"]
    freq_levels_45 = np.linspace(-10, 10, nfreq)
    freq_ticks_45 = [-10, -5, 0, 5, 10]
    freq_ticklabels_45 = ["-10", "-5", "0", "5", "10"]
    gamma_levels_15 = np.linspace(0, 0.6, ngamma)
    gamma_ticks_15 = [0, 0.25, 0.5]
    gamma_ticklabels_15 = ["0", "0.25", "0.5"]
    gamma_levels_3 = np.linspace(0, 0.9, ngamma)
    gamma_ticks_3 = [0, 0.3, 0.6, 0.9]
    gamma_ticklabels_3 = ["0", "0.3", "0.6", "0.9"]
    gamma_levels_45 = np.linspace(0, 2.7, ngamma)
    gamma_ticks_45 = [0, 1, 2]
    gamma_ticklabels_45 = ["0", "1", "2"]
    freq_levels_list = [freq_levels_15, freq_levels_3, freq_levels_45]
    gamma_levels_list = [gamma_levels_15, gamma_levels_3, gamma_levels_45]
    freq_ticks_list = [freq_ticks_15, freq_ticks_3, freq_ticks_45]
    gamma_ticks_list = [gamma_ticks_15, gamma_ticks_3, gamma_ticks_45]
    freq_ticklabels_list = [freq_ticklabels_15, freq_ticklabels_3, freq_ticklabels_45]
    gamma_ticklabels_list = [gamma_ticklabels_15, gamma_ticklabels_3, gamma_ticklabels_45]

    for counter, ky_idx in enumerate([ky_15_idx,
                                      ky_3_idx,
                                      ky_45_idx]):
        fprim_mesh, tprim_mesh, [freq_meshgrid,
            gammaom_meshgrid] = uniquearrays2meshgrids(unique_fprim,
                        unique_tprim, [freq_fprim_tprim_ky_array[:,:,ky_idx],
                            gammasafe_fprim_tprim_ky_array[:,:,ky_idx],
                        ])
        make_plots(freq_axes[counter], gamma_axes[counter],
                    freq_cbaxes[counter],
                    gamma_cbaxes[counter],
                    fprim_mesh, tprim_mesh, freq_meshgrid,
                    gammaom_meshgrid,
                    freq_levels_list[counter],
                    gamma_levels_list[counter],
                    freq_ticks_list[counter], gamma_ticks_list[counter]
                    )

    finish_plots(freq_ticklabels_list, gamma_ticklabels_list)

    return

def for_thesis_make_fprim_tprim_ky_scan_w003():
    """For fancy plotting """

    def make_plots(freq_ax, gamma_ax, freq_cbax, gamma_cbax, fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid,
                   gammaom_fprim_tprim_meshgrid, freq_levels, gamma_levels, freq_ticks, gamma_ticks
                   ):
        """ """
        # gamma_levels=20
        # freq_levels=20
        gammaom_contours = gamma_ax.contourf(fprim_mesh, tprim_mesh, gammaom_fprim_tprim_meshgrid, gamma_levels, cmap="inferno", extend="max")
        # freq_lim = math.ceil(max(abs(freq_fprim_tprim_meshgrid.min()), freq_fprim_tprim_meshgrid.max())*10)/10
        # freq_spacing = 0.1 # 0.02
        # freq_levels = np.arange(-freq_lim, freq_lim, freq_spacing)
        freq_contours = freq_ax.contourf(fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid, freq_levels, cmap="PiYG", extend="both")

        fig.colorbar(gammaom_contours, cax=gamma_cbax, ticks=gamma_ticks )
        fig.colorbar(freq_contours, cax=freq_cbax, ticks=freq_ticks)


        return

    def init_plots():

        fig = plt.figure(figsize=[13, 8])
        frac = (13.0/8.0)
        # Axis placement parameters.
        vspace_min = 0.1
        ax_hspace_min = 0.08
        cb_hspace = 0.01
        plot_dim_y_try = 0.38
        plot_dim_x_try = plot_dim_y_try/frac
        right = 0.93
        left = 0.05
        top = 0.99
        bottom = 0.1
        cb_width = 0.015
        plot_width_with_min_spacing = (right - left - 2*ax_hspace_min - 3*cb_hspace - 3*cb_width)/3
        plot_dim_x = plot_width_with_min_spacing
        plot_dim_y = plot_dim_x * frac

        col1_lhs = left
        col1_cb_lhs = left + plot_dim_x + cb_hspace
        col2_lhs = col1_cb_lhs + cb_width + ax_hspace_min
        col2_cb_lhs = col2_lhs + plot_dim_x + cb_hspace
        col3_lhs = col2_cb_lhs + cb_width + ax_hspace_min
        col3_cb_lhs = col3_lhs + plot_dim_x + cb_hspace
        # col2_cb_lhs = 0.94
        #row3_bottom = 0.08
        row2_bottom = bottom # row3_bottom + plot_dim_y + vspace
        row1_bottom = row2_bottom + plot_dim_y + vspace_min


        ########################################################################
        # Create axes.
        # [ax1, ax2] = [gammaom, freq]
        # [ax3, ax4] = [gammap2, phase velocity]
        ########################################################################
        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax2 = fig.add_axes([col1_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax3 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax4 = fig.add_axes([col2_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax5 = fig.add_axes([col3_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax6 = fig.add_axes([col3_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax2 = fig.add_axes([col1_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax3 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax4 = fig.add_axes([col2_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax5 = fig.add_axes([col3_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax6 = fig.add_axes([col3_cb_lhs, row2_bottom, cb_width, plot_dim_y])

        return (fig, ax1, ax2, ax3, ax4, ax5, ax6, cbax1, cbax2, cbax3, cbax4, cbax5, cbax6)

    def finish_plots(freq_ticklabels_list, gamma_ticklabels_list):
        # Formatting, axis labels etc.

        x_ticklabelfontsize = 14
        x_labelfontsize = 20
        my_xticklength = 3
        my_xtickwidth = 1
        y_labelfontsize = x_labelfontsize
        y_ticklabelfontsize = x_ticklabelfontsize
        cb_ticklabelsize = 14
        cbax_title_fontsize = 20
        ax_title_fontsize = 20

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
            #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
            ax.tick_params("x", which="major", top=True, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.set_xlim([np.min(unique_fprim), np.max(unique_fprim)])
            ax.set_ylim([np.min(unique_tprim), np.max(unique_tprim)])
            ax.set_xlabel(r"$a/L_n$", fontsize=x_labelfontsize)
            ax.set_ylabel(r"$a/L_T$", fontsize=y_labelfontsize)
            ax.set_xticks([0, 2, 4, 6, 8])
            ax.set_xticklabels(["0", "2", "4", "6", "8"], fontsize=x_ticklabelfontsize)
            ax.set_yticks([0, 2, 4, 6, 8])
            ax.set_yticklabels(["0", "2", "4", "6", "8"], fontsize=y_ticklabelfontsize)

        ax1.set_title(r"$\tilde{k}_y=1.5$", fontsize=ax_title_fontsize)
        ax3.set_title(r"$\tilde{k}_y=3$", fontsize=ax_title_fontsize)
        ax5.set_title(r"$\tilde{k}_y=4.5$", fontsize=ax_title_fontsize)

        for cbax in [cbax1, cbax3, cbax5]:
            cbax.set_title(r"$\omega$", fontsize=cbax_title_fontsize)
        for cbax in [cbax2, cbax4, cbax6]:
            cbax.set_title(r"$\gamma$", fontsize=cbax_title_fontsize)

        #plt.show()

        cbax1.set_yticklabels(freq_ticklabels_list[0], fontsize=cb_ticklabelsize)
        cbax3.set_yticklabels(freq_ticklabels_list[1], fontsize=cb_ticklabelsize)
        cbax5.set_yticklabels(freq_ticklabels_list[2], fontsize=cb_ticklabelsize)
        cbax2.set_yticklabels(gamma_ticklabels_list[0], fontsize=cb_ticklabelsize)
        cbax4.set_yticklabels(gamma_ticklabels_list[1], fontsize=cb_ticklabelsize)
        cbax6.set_yticklabels(gamma_ticklabels_list[2], fontsize=cb_ticklabelsize)
        # Uncomment to save the figure.
        plt.savefig(str("images/w003_fprim_tprim_beta003") + ".png")
        plt.close()

        return

    folder_1 = make_scans.w003_impl_em_lower_resolution_beta003

    pickle_longname = folder_1 + "/omega_fprim_tprim_ky_array.pickle"
    file = open(pickle_longname, "rb")
    [pickle_string, unique_fprim, unique_tprim, unique_ky,
        gammaom_fprim_tprim_ky_array, gammasafe_fprim_tprim_ky_array,
        freq_fprim_tprim_ky_array] = pickle.load(file)
    file.close()
    print("unique_ky = ", unique_ky)
    (fig, ax1, ax2, ax3, ax4, ax5, ax6, cbax1, cbax2, cbax3, cbax4, cbax5, cbax6) = init_plots()
    ky_15_idx = unique_ky.index(1.5)
    ky_3_idx = unique_ky.index(3)
    ky_45_idx = unique_ky.index(4.5)

    freq_axes = [ax1, ax3, ax5]
    freq_cbaxes = [cbax1, cbax3, cbax5]
    gamma_axes = [ax2, ax4, ax6]
    gamma_cbaxes = [cbax2, cbax4, cbax6]

    nfreq = 50
    ngamma = 50
    freq_levels_15 = np.linspace(-2.5, 2.5, nfreq)
    freq_ticks_15 = [-2, -1, 0, 1, 2]
    freq_ticklabels_15 = ["-2", "-1", "0", "1", "2"]
    freq_levels_3 = np.linspace(-10, 10, nfreq)
    freq_ticks_3 = [-10, -5, 0, 5, 10]
    freq_ticklabels_3 = ["-10", "-5", "0", "5", "10"]
    freq_levels_45 = np.linspace(-10, 10, nfreq)
    freq_ticks_45 = [-10, -5, 0, 5, 10]
    freq_ticklabels_45 = ["-10", "-5", "0", "5", "10"]
    gamma_levels_15 = np.linspace(0, 0.6, ngamma)
    gamma_ticks_15 = [0, 0.25, 0.5]
    gamma_ticklabels_15 = ["0", "0.25", "0.5"]
    gamma_levels_3 = np.linspace(0, 0.9, ngamma)
    gamma_ticks_3 = [0, 0.3, 0.6, 0.9]
    gamma_ticklabels_3 = ["0", "0.3", "0.6", "0.9"]
    gamma_levels_45 = np.linspace(0, 2.7, ngamma)
    gamma_ticks_45 = [0, 1, 2]
    gamma_ticklabels_45 = ["0", "1", "2"]
    freq_levels_list = [freq_levels_15, freq_levels_3, freq_levels_45]
    gamma_levels_list = [gamma_levels_15, gamma_levels_3, gamma_levels_45]
    freq_ticks_list = [freq_ticks_15, freq_ticks_3, freq_ticks_45]
    gamma_ticks_list = [gamma_ticks_15, gamma_ticks_3, gamma_ticks_45]
    freq_ticklabels_list = [freq_ticklabels_15, freq_ticklabels_3, freq_ticklabels_45]
    gamma_ticklabels_list = [gamma_ticklabels_15, gamma_ticklabels_3, gamma_ticklabels_45]

    for counter, ky_idx in enumerate([ky_15_idx,
                                      ky_3_idx,
                                      ky_45_idx]):
        fprim_mesh, tprim_mesh, [freq_meshgrid,
            gammaom_meshgrid] = uniquearrays2meshgrids(unique_fprim,
                        unique_tprim, [freq_fprim_tprim_ky_array[:,:,ky_idx],
                            gammasafe_fprim_tprim_ky_array[:,:,ky_idx],
                        ])
        make_plots(freq_axes[counter], gamma_axes[counter],
                    freq_cbaxes[counter],
                    gamma_cbaxes[counter],
                    fprim_mesh, tprim_mesh, freq_meshgrid,
                    gammaom_meshgrid,
                    freq_levels_list[counter],
                    gamma_levels_list[counter],
                    freq_ticks_list[counter], gamma_ticks_list[counter]
                    )

    finish_plots(freq_ticklabels_list, gamma_ticklabels_list)

    return

def for_thesis_make_fprim_tprim_ky_scan_w015_027():
    """For fancy plotting """

    def make_plots(freq_ax, gamma_ax, freq_cbax, gamma_cbax, fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid,
                   gammaom_fprim_tprim_meshgrid, freq_levels, gamma_levels, freq_ticks, gamma_ticks
                   ):
        """ """
        # gamma_levels=20
        # freq_levels=20
        gammaom_contours = gamma_ax.contourf(fprim_mesh, tprim_mesh, gammaom_fprim_tprim_meshgrid, gamma_levels, cmap="inferno", extend="max")
        # freq_lim = math.ceil(max(abs(freq_fprim_tprim_meshgrid.min()), freq_fprim_tprim_meshgrid.max())*10)/10
        # freq_spacing = 0.1 # 0.02
        # freq_levels = np.arange(-freq_lim, freq_lim, freq_spacing)
        freq_contours = freq_ax.contourf(fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid, freq_levels, cmap="PiYG", extend="both")

        fig.colorbar(gammaom_contours, cax=gamma_cbax, ticks=gamma_ticks )
        fig.colorbar(freq_contours, cax=freq_cbax, ticks=freq_ticks)


        return

    def init_plots():

        fig = plt.figure(figsize=[10, 10])
        frac = (10.0/10.0)
        # Axis placement parameters.
        vspace_min = 0.1
        ax_hspace_min = 0.08
        cb_hspace = 0.01
        plot_dim_y_try = 0.38
        plot_dim_x_try = plot_dim_y_try/frac
        right = 0.93
        left = 0.08
        top = 0.99
        bottom = 0.1
        cb_width = 0.015
        plot_width_with_min_spacing = (right - left - ax_hspace_min - 2*cb_hspace - 2*cb_width)/2
        plot_dim_x = plot_width_with_min_spacing
        plot_dim_y = plot_dim_x * frac

        col1_lhs = left
        col1_cb_lhs = left + plot_dim_x + cb_hspace
        col2_lhs = col1_cb_lhs + cb_width + ax_hspace_min
        col2_cb_lhs = col2_lhs + plot_dim_x + cb_hspace
        col3_lhs = col2_cb_lhs + cb_width + ax_hspace_min
        col3_cb_lhs = col3_lhs + plot_dim_x + cb_hspace
        # col2_cb_lhs = 0.94
        #row3_bottom = 0.08
        row2_bottom = bottom # row3_bottom + plot_dim_y + vspace
        row1_bottom = row2_bottom + plot_dim_y + vspace_min


        ########################################################################
        # Create axes.
        # [ax1, ax2] = [gammaom, freq]
        # [ax3, ax4] = [gammap2, phase velocity]
        ########################################################################
        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax2 = fig.add_axes([col1_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax3 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax4 = fig.add_axes([col2_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax2 = fig.add_axes([col1_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax3 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax4 = fig.add_axes([col2_cb_lhs, row2_bottom, cb_width, plot_dim_y])

        return (fig, ax1, ax2, ax3, ax4, cbax1, cbax2, cbax3, cbax4)

    def finish_plots(freq_ticklabels_list, gamma_ticklabels_list):
        # Formatting, axis labels etc.

        x_ticklabelfontsize = 14
        x_labelfontsize = 20
        my_xticklength = 3
        my_xtickwidth = 1
        y_labelfontsize = x_labelfontsize
        y_ticklabelfontsize = x_ticklabelfontsize
        cb_ticklabelsize = 14
        cbax_title_fontsize = 20
        ax_title_fontsize = 20

        for ax in [ax1, ax2, ax3, ax4]:
            #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
            #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
            ax.tick_params("x", which="major", top=True, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.set_xlim([np.min(unique_fprim), np.max(unique_fprim)])
            ax.set_ylim([np.min(unique_tprim), np.max(unique_tprim)])
            ax.set_xlabel(r"$a/L_n$", fontsize=x_labelfontsize)
            ax.set_ylabel(r"$a/L_T$", fontsize=y_labelfontsize)
            ax.set_xticks([0, 2, 4, 6, 8])
            ax.set_xticklabels(["0", "2", "4", "6", "8"], fontsize=x_ticklabelfontsize)
            ax.set_yticks([0, 2, 4, 6, 8])
            ax.set_yticklabels(["0", "2", "4", "6", "8"], fontsize=y_ticklabelfontsize)

        ax1.set_title(r"w015", fontsize=ax_title_fontsize)
        ax3.set_title(r"w027", fontsize=ax_title_fontsize)

        for cbax in [cbax1, cbax3]:
            cbax.set_title(r"$\omega$", fontsize=cbax_title_fontsize)
        for cbax in [cbax2, cbax4]:
            cbax.set_title(r"$\gamma$", fontsize=cbax_title_fontsize)

        #plt.show()

        cbax1.set_yticklabels(freq_ticklabels_list[0], fontsize=cb_ticklabelsize)
        cbax3.set_yticklabels(freq_ticklabels_list[1], fontsize=cb_ticklabelsize)
        cbax2.set_yticklabels(gamma_ticklabels_list[0], fontsize=cb_ticklabelsize)
        cbax4.set_yticklabels(gamma_ticklabels_list[1], fontsize=cb_ticklabelsize)
        plt.savefig("draft_fprim_tprim_w015_w027.png")
        plt.close()

        return

    pickle_longname_015 = folder_2 + "/omega_fprim_tprim_ky_array.pickle"
    pickle_longname_027 = folder_3 + "/omega_fprim_tprim_ky_array.pickle"

    nfreq = 50
    ngamma = 50
    freq_levels_w015 = np.linspace(-10, 10, nfreq) # np.linspace(-2.5, 2.5, nfreq)
    freq_ticks_w015 = [-10, -5, 0, 5, 10] # [-2, -1, 0, 1, 2]
    freq_ticklabels_w015 =  ["-10", "-5", "0", "5", "10"] # ["-2", "-1", "0", "1", "2"]
    freq_levels_w027 = np.linspace(-10, 10, nfreq) # np.linspace(-2.5, 2.5, nfreq)
    freq_ticks_w027 = [-10, -5, 0, 5, 10] # [-2, -1, 0, 1, 2]
    freq_ticklabels_w027 = ["-10", "-5", "0", "5", "10"] # ["-2", "-1", "0", "1", "2"]
    gamma_levels_w015 = np.linspace(0, 0.9, ngamma) # np.linspace(0, 0.6, ngamma)
    gamma_ticks_w015 = [0, 0.3, 0.6, 0.9] # [0, 0.25, 0.5]
    gamma_ticklabels_w015 = ["0", "0.3", "0.6", "0.9"] # ["0", "0.25", "0.5"]
    gamma_levels_w027 = np.linspace(0, 0.9, ngamma) # np.linspace(0, 0.6, ngamma)
    gamma_ticks_w027 = [0, 0.3, 0.6, 0.9] # [0, 0.25, 0.5]
    gamma_ticklabels_w027 = ["0", "0.3", "0.6", "0.9"] # ["0", "0.25", "0.5"]

    #freq_levels_3 = np.linspace(-10, 10, nfreq)
    # freq_ticks_3 = [-10, -5, 0, 5, 10]
    # freq_ticklabels_3 = ["-10", "-5", "0", "5", "10"]
    # gamma_levels_3 = np.linspace(0, 0.9, ngamma)
    # gamma_ticks_3 = [0, 0.3, 0.6, 0.9]
    # gamma_ticklabels_3 = ["0", "0.3", "0.6", "0.9"]

    freq_levels_list = [freq_levels_w015, freq_levels_w027]
    gamma_levels_list = [gamma_levels_w015, gamma_levels_w027]
    freq_ticks_list = [freq_ticks_w015, freq_ticks_w027]
    gamma_ticks_list = [gamma_ticks_w015, gamma_ticks_w027]
    freq_ticklabels_list = [freq_ticklabels_w015, freq_ticklabels_w027]
    gamma_ticklabels_list = [gamma_ticklabels_w015, gamma_ticklabels_w027]

    (fig, ax1, ax2, ax3, ax4, cbax1, cbax2, cbax3, cbax4) = init_plots()

    freq_axes = [ax1, ax3]
    freq_cbaxes = [cbax1, cbax3]
    gamma_axes = [ax2, ax4]
    gamma_cbaxes = [cbax2, cbax4]

    for eq_idx, pickle_longname in enumerate([pickle_longname_015, pickle_longname_027]):
        file = open(pickle_longname, "rb")
        [pickle_string, unique_fprim, unique_tprim, unique_ky,
            gammaom_fprim_tprim_ky_array, gammasafe_fprim_tprim_ky_array,
            freq_fprim_tprim_ky_array] = pickle.load(file)
        file.close()

        # ky_idx = unique_ky.index(1.5)
        ky_idx = unique_ky.index(3)
        # ky_45_idx = unique_ky.index(4.5)

        fprim_mesh, tprim_mesh, [freq_meshgrid,
            gammaom_meshgrid] = uniquearrays2meshgrids(unique_fprim,
                        unique_tprim, [freq_fprim_tprim_ky_array[:,:,ky_idx],
                            gammasafe_fprim_tprim_ky_array[:,:,ky_idx],
                        ])
        ## Reverse the sign of the frequency - I think these equilibria are
        ## "the other way around" compared to kjm3
        freq_meshgrid = - freq_meshgrid
        make_plots(freq_axes[eq_idx], gamma_axes[eq_idx],
                    freq_cbaxes[eq_idx],
                    gamma_cbaxes[eq_idx],
                    fprim_mesh, tprim_mesh, freq_meshgrid,
                    gammaom_meshgrid,
                    freq_levels_list[eq_idx],
                    gamma_levels_list[eq_idx],
                    freq_ticks_list[eq_idx], gamma_ticks_list[eq_idx]
                    )

    finish_plots(freq_ticklabels_list, gamma_ticklabels_list)

    return

def make_fprim_tprim_ky_scan(folder_longname, save_name_prefix="draft"):
    """For looking at the data, not fancy plotting """

    def make_plots(ky_val, fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid,
                   gammaom_fprim_tprim_meshgrid,
                   gammasafe_fprim_tprim_meshgrid,
                   save_name):
        """Create nice-looking plots for the 3 quantities."""

        fig = plt.figure(figsize=[12, 12])

        # Axis placement parameters.
        plot_dim_x = 0.38
        plot_dim_y = plot_dim_x
        col1_lhs = 0.05
        col2_lhs = 0.54
        col1_cb_lhs = 0.45
        col2_cb_lhs = 0.94
        row1_bottom = 0.54
        row2_bottom = 0.08
        cb_width = 0.015

        ########################################################################
        # Create axes.
        # [ax1, ax2] = [gammaom, freq]
        # [ax3, ax4] = [gammap2, phase velocity]
        ########################################################################
        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax2 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax3 = fig.add_axes([col1_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax2 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax3 = fig.add_axes([col1_cb_lhs, row2_bottom, cb_width, plot_dim_y])

        ########################################################################
        # Put the data on the axes. We plot gammap2 first, as it is less likely
        # to have anaolous values than gammaom (prevent skewing the contours with
        # anomolous values).
        ########################################################################
        gammap2_contours = ax3.contourf(fprim_mesh, tprim_mesh, gammasafe_fprim_tprim_meshgrid, 20, cmap="inferno")
        # if (np.isfinite(gammasafe_fprim_tprim_meshgrid)).all():
        #     gammaom_contours = ax1.contourf(fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid, gammap2_contours.levels, cmap="inferno")
        # else:
        gammaom_contours = ax1.contourf(fprim_mesh, tprim_mesh, gammaom_fprim_tprim_meshgrid, 20, cmap="inferno")

        freq_lim = math.ceil(max(abs(freq_fprim_tprim_meshgrid.min()), freq_fprim_tprim_meshgrid.max())*10)/10

        freq_spacing = 0.02
        freq_contours = ax2.contourf(fprim_mesh, tprim_mesh, freq_fprim_tprim_meshgrid, np.arange(-freq_lim, freq_lim, freq_spacing), cmap="PiYG")

        fig.colorbar(gammaom_contours, cax=cbax1)
        fig.colorbar(freq_contours, cax=cbax2)
        fig.colorbar(gammap2_contours, cax=cbax3)

        # Label the plot with kx, ky
        fig.text(0.7, 0.3, (r"$k_y = $" + str(ky_val)), fontsize=18)

        # Formatting, axis labels etc.
        for ax in [ax1, ax2, ax3]:
            #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
            #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
            ax.set_xlim([np.min(unique_fprim), np.max(unique_fprim)])
            ax.set_ylim([np.min(unique_tprim), np.max(unique_tprim)])
            ax.set_xlabel("fprim")
            ax.set_ylabel("tprim")

        # Titles
        ax1.set_title(r"$\gamma$ (.omega file)");# ax4.set_title(r"$\gamma$")
        ax2.set_title(r"$\omega$ (.omega file)"); #ax5.set_title(r"$\omega$")
        ax3.set_title(r"$\gamma$ ($\phi^2$ fit)"); #ax6.set_title(r"$k_y$")

        plt.savefig(save_name)
        plt.close()

        # Uncomment to save the figure.
        # plt.savefig(str(target_ky) + ".png")
        # plt.close()
        return

    pickle_longname = folder_longname + "/omega_fprim_tprim_ky_array.pickle"
    file = open(pickle_longname, "rb")
    [pickle_string, unique_fprim, unique_tprim, unique_ky,
        gammaom_fprim_tprim_ky_array, gammasafe_fprim_tprim_ky_array,
        freq_fprim_tprim_ky_array] = pickle.load(file)
    file.close()

    for ky_idx, ky_val in enumerate(unique_ky):
        fprim_mesh, tprim_mesh, [freq_meshgrid, gammaom_meshgrid,
            gammasafe_meshgrid] = uniquearrays2meshgrids(unique_fprim,
                        unique_tprim, [freq_fprim_tprim_ky_array[:,:,ky_idx],
                            gammaom_fprim_tprim_ky_array[:,:,ky_idx], gammasafe_fprim_tprim_ky_array[:,:,ky_idx],
                        ])
        ky_str = "ky{:.3f}".format(ky_val)
        save_name = save_name_prefix + "_" + ky_str + ".eps"
        make_plots(ky_val, fprim_mesh, tprim_mesh, freq_meshgrid,
                    gammaom_meshgrid,
                    gammasafe_meshgrid,
                    save_name
                    )
    return

if __name__ == "__main__":
    print("Hello world")

    # make_fprim_tprim_ky_scan(make_scans.folder_name_expl_ky05, save_name_prefix="w7x_expl_em.eps")
    # make_fprim_tprim_ky_scan(make_scans.folder_name_impl_ky05, save_name_prefix="w7x_expl_em.eps")
    # make_fprim_tprim_ky_scan(make_scans.folder_name_expl_higher_ky_em, save_name_prefix="w7x_expl_em.eps")
    # make_fprim_tprim_ky_scan(make_scans.folder_name_expl_higher_ky_es, save_name_prefix="w7x_expl_es.eps")
    # make_fprim_tprim_ky_scan(make_scans.folder_name_impl_higher_ky_em, save_name_prefix="w7x_impl_em.eps")
    # make_fprim_tprim_ky_scan(make_scans.folder_name_impl_higher_ky_es, save_name_prefix="w7x_impl_es.eps")
    # for_thesis_make_fprim_tprim_ky_scan_kjm3()
    for_thesis_make_fprim_tprim_ky_scan_w003()
