""" """


import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import math
from scipy.interpolate import interp2d


default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']


nonlinear_file = "sims/leapfrog/input_leapfrog_restarted.nonlinear_quantities"
naky = 5
nakx = 7
ny = 14
nx = 10

dx = 1.9230769230769229
dy = 1.3463968515384828
dt = 0.134
xmax = (nx) * dx
ymax = (ny) * dy

ky = np.array([0.0000000000000000, 0.33333333333333331, 0.66666666666666663, 1.0000000000000000, 1.3333333333333333])
kx = np.array([0.0000000000000000, 0.32672563597333848, 0.65345127194667696,
                0.98017690792001544, -0.98017690792001544, -0.65345127194667696, -0.32672563597333848])

x_grid = np.arange(0, xmax, dx )
y_grid = np.arange(0, ymax, dy )
x_grid_upsampled = np.arange(0, xmax, dx/2 )
y_grid_upsampled = np.arange(0, ymax, dy/2 )
x_idxs = np.arange(0, nx, 1, dtype="int")
y_idxs = np.arange(0, ny, 1, dtype="int")
print("x_grid, xmax = ", x_grid, xmax)
## Get the idxs in 2D
x_idxs_2d = np.zeros((ny,nx), dtype="int")
y_idxs_2d = np.zeros((ny,nx), dtype="int")
x_grid_2d = np.zeros((ny,nx))
y_grid_2d = np.zeros((ny,nx))
x_grid_2d_upsampled = np.zeros((2*ny,2*nx))
y_grid_2d_upsampled = np.zeros((2*ny,2*nx))

very_small_dt = 1e-8

for yidx in range(0, ny):
    x_idxs_2d[yidx, :] = x_idxs
    x_grid_2d[yidx, :] = x_grid

for xidx in range(0, nx):
    y_idxs_2d[:,xidx] = y_idxs
    y_grid_2d[:,xidx] = y_grid

for yidx in range(0, 2*ny):
    x_grid_2d_upsampled[yidx, :] = x_grid_upsampled

for xidx in range(0, 2*nx):
    y_grid_2d_upsampled[:,xidx] = y_grid_upsampled
#
# print("x_grid = ", x_grid)
# print("y_grid = ", y_grid)
# print("x_grid_2d = ", x_grid_2d)
# print("y_grid_2d = ", y_grid_2d)

def get_array_for_real_plaintext(plaintext_block):
    """ """
    #print("plaintext_block = ", plaintext_block)
    # Get rid of "quantity = " term
    cropped_plaintext = re.split(":", plaintext_block.strip())[1]
    lines_str_list = re.split("\n", cropped_plaintext.strip())

    data_array = np.zeros((ny, nx))
    for iy in range(0, ny):
        entry_strs = re.split("\s+", lines_str_list[iy].strip())
        for ix in range(0, nx):
            data_array[iy,ix] = float(entry_strs[ix])

    #print("data_array = ", data_array)

    return data_array

def get_array_for_complex_plaintext(plaintext_block):
    """ """
    cropped_plaintext = re.split(":", plaintext_block.strip())[1]
    # Each line is a single iy i.e. coresponding to a single y value
    lines_str_list = re.split("\n", cropped_plaintext.strip())

    data_array = np.zeros((naky, nakx), dtype="complex")
    for iky in range(0, naky):
        entry_strs = re.split("\s+", lines_str_list[iky].strip())
        for ikx in range(0, nakx):
            [real_str, imag_str] = re.split(",", entry_strs[ikx])
            real_val = float(real_str[1:])
            imag_val = float(imag_str[:-1])
            data_array[iky,ikx] = real_val + 1j*imag_val


    return data_array

def get_arrays_from_nonlinear_data():
    """ """
    myfile = open(nonlinear_file, "r")
    data = myfile.read()
    myfile.close()
    data_blocks = re.split("XXXXXXXXXXXXXXXXXXXXXX", data.strip())
    golder_block = data_blocks[0]
    golderyx_block = data_blocks[1]
    dgold_dy_block = data_blocks[2]
    dgold_dx_block = data_blocks[3]
    vchiold_y_block = data_blocks[4]
    vchiold_x_block = data_blocks[5]

    golderyx_array = get_array_for_real_plaintext(golderyx_block)
    dgold_dy_array = get_array_for_real_plaintext(dgold_dy_block)
    dgold_dx_array = get_array_for_real_plaintext(dgold_dx_block)
    vchiold_y_array = get_array_for_real_plaintext(vchiold_y_block)
    vchiold_x_array = get_array_for_real_plaintext(vchiold_x_block)
    golder_array = get_array_for_complex_plaintext(golder_block)

    return [golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
        vchiold_y_array, vchiold_x_array]

def create_upsampled_grid(data_array):
    """Given a data array, return an upsampled version of the array.
    Original array samples on an x, y grid:
    *   *   *   *   *

    *   *   *   *   *

    *   *   *   *   *

    We want to sample at midpoints:
    & & & & & & & & & &
    * & * & * & * & * &
    & & & & & & & & & &
    * & * & * & * & * &
    & & & & & & & & & &
    * & * & * & * & * &

    where & denotes a new sampling point. The final & in x, y should be an
    interpolation of the final * and the first * (periodic).
    """
    (ny, nx) = data_array.shape
    data_array_upsampled = np.zeros((2*ny, 2*nx))

    ## Probably can do this with scipy's interpolate library, but
    ## want to ensure BCs are correct.
    total_upsampled_gridpoints = 4 * nx*ny
    upsampled_x_idx = 0
    upsampled_y_idx = 0

    for upsampled_xidx in range(0, 2*nx):
        for upsampled_yidx in range(0, 2*ny):
            ## 4 cases to consider:
            # (1) x_idx even, y_idx even. Use the original data point
            # (2) x_idx even, y_idx odd. Upsampled value is 1/2 y_up, 1/2 y_down
            # (3) x_idx odd, y_idx even. Upsampled value is 1/2 x_left, 1/2 x_right
            # (4) x_idx odd, y_idx odd. Upsampled value is 1/4 (x_left,y_up),
            #                            1/4 (x_left,y_down), 1/4 (x_right,y_up), 1/4 (x_right,y_down)
            if (upsampled_xidx%2 == 0) and (upsampled_yidx%2 == 0):
                data_array_upsampled[upsampled_yidx, upsampled_xidx] = data_array[int(upsampled_yidx/2), int(upsampled_xidx/2)]
            elif (upsampled_xidx%2 == 0) and (upsampled_yidx%2 != 0):
                yidx_down = math.floor(upsampled_yidx/2)
                # %(ny) means that the final upsampled_yidx  uses yidx_up = 0
                yidx_up = math.ceil(upsampled_yidx/2)%(ny)
                data_array_upsampled[upsampled_yidx, upsampled_xidx] = (0.5*data_array[yidx_down, int(upsampled_xidx/2)]
                                                        +   0.5*data_array[yidx_up, int(upsampled_xidx/2)])
            elif (upsampled_xidx%2 != 0) and (upsampled_yidx%2 == 0):
                xidx_left = math.floor(upsampled_xidx/2)
                # %(nx) means that the final upsampled_xidx  uses xidx_right = 0
                xidx_right = math.ceil(upsampled_xidx/2)%(nx)
                data_array_upsampled[upsampled_yidx, upsampled_xidx] = (0.5*data_array[int(upsampled_yidx/2), xidx_left]
                                                        +   0.5*data_array[int(upsampled_yidx/2), xidx_right])

            elif (upsampled_xidx%2 != 0) and (upsampled_yidx%2 != 0):
                xidx_left = math.floor(upsampled_xidx/2)
                # %(nx) means that the final upsampled_xidx  uses xidx_right = 0
                xidx_right = math.ceil(upsampled_xidx/2)%(nx)
                yidx_down = math.floor(upsampled_yidx/2)
                # %(ny) means that the final upsampled_yidx  uses yidx_up = 0
                yidx_up = math.ceil(upsampled_yidx/2)%(ny)
                data_array_upsampled[upsampled_yidx, upsampled_xidx] = (0.25*data_array[yidx_down, xidx_left]
                                                        +   0.25*data_array[yidx_up, xidx_left]
                                                        +   0.25*data_array[yidx_down, xidx_right]
                                                        +   0.25*data_array[yidx_up, xidx_right])

    return data_array_upsampled

def nisl_step_ritchie(golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
                vchiold_y_array, vchiold_x_array, with_diagnostics=False):
    """The NISL step according to Ritchie; here we don't calculate the (approximate)
    departure point, but instead (attempt to) find the trajectory from t^n which arrives
    closest to our gridpoint at t^(n+1)"""

    def update_p_and_q(p_array, q_array, vchiold_x_array_upsampled, vchiold_y_array_upsampled,
                        yidx_for_upsampled_array, xidx_for_upsampled_array):
        """Update guess for p, q, and the idxs corresponding to the upsampled
        array."""
        p_array = np.rint(2 * dt/dx * vchiold_x_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array]).astype("int")
        q_array = np.rint(2 * dt/dy * vchiold_y_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array]).astype("int")
        xidx_for_upsampled_array = (2*x_idxs_2d - p_array)%(2*nx)
        yidx_for_upsampled_array = (2*y_idxs_2d - q_array)%(2*ny)
        return p_array, q_array, yidx_for_upsampled_array, xidx_for_upsampled_array

    #######################################################################
    # gnew(x_i, y_j) = golder[x_i - p_ij*dx, y_j - * q_ij*dy]
    #                       + rhs_ij
    #  with
    #
    #  rhs_ij = - (vresidual_x_ij*dgold_dx[x_i - p_ij*dx/2, y_j - * q_ij*dy/2]
    #            + vresidual_y_ij*dgold_dy[x_i - p_ij*dx/2, y_j - * q_ij*dy/2])
    #
    # vresidual_x_ij = vchiold_x[x_i - p_ij*dx/2, y_j - * q_ij*dy/2]
    #                    - p_ij * dx/(2*dt)
    # vresidual_y_ij = vchiold_y[x_i - p_ij*dx/2, y_j - * q_ij*dy/2]
    #                    - q_ij * dy/(2*dt)
    #
    #######################################################################
    ### OLD IDEA
    # Estimate the departure points based on vchiold
    # #####################################################################
    # Idea: Perform
    # (a) p_ij = q_ij = 0
    # (b) p_ij = NINT(2 dt/dx vchiold_x[x_i - p*dx/2, y_j - q*dy/2])
    #     q_ij = NINT(2 dt/dy vchiold_y[x_i - p*dx/2, y_j - q*dy/2])
    # (c) Repeat (b) a few times
    #
    # To get vchiold_x,y[(x_i - p*dx/2), (y_j - q*dy/2)], need to upsample
    # vchiold_x,y
    # #####################################################################

    # p=0 means shift non-upsampled xidx by 0, and shift upsampled xidx by 0
    # p=1 means shift non-upsampled xidx by 1, and shift upsampled xidx by 1
    #     (because upsampled x is being shifted by p*dx/2, but is sampled
    #      every dx/2)

    p_array = np.zeros((ny, nx), dtype="int")
    q_array = np.zeros((ny, nx), dtype="int")

    # Upsample - double the number of points in the vchi, gold, golder grids.
    # golderyx_array_upsampled = create_upsampled_grid(golderyx_array)

    dgold_dy_array_upsampled = create_upsampled_grid(dgold_dy_array)
    dgold_dx_array_upsampled = create_upsampled_grid(dgold_dx_array)
    vchiold_x_array_upsampled = create_upsampled_grid(vchiold_x_array)
    vchiold_y_array_upsampled = create_upsampled_grid(vchiold_y_array)
    # To get the arrival points:
    # xnew = (xold + (2*dt*vchiold_x))%xmax
    # ynew = (yold + (2*dt*vchiold_y))%ymax

    # As a diagnostic, plot the regular grid and the arrival locations.
    if with_diagnostics:
        marker_size = 20.

        xnew = (x_grid_2d_upsampled - (dt*vchiold_x_array_upsampled))%xmax
        ynew = (y_grid_2d_upsampled - (dt*vchiold_y_array_upsampled))%ymax

        fig = plt.figure(figsize=[12, 8])
        ax1 = fig.add_subplot(111)
        ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=marker_size, label="upsampled grid")
        ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
        ax1.scatter(xnew.flatten(), ynew.flatten(), s=marker_size, label="arrival points")
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        ax1.grid(True)
        ax1.set_xlim([-1, 23])
        ax1.legend(loc="upper right")
        plt.show()
        # sys.exit()

    upsampled_xidxs = np.arange(0, 2*nx, 1, dtype="int")
    upsampled_yidxs = np.arange(0, 2*ny, 1, dtype="int")
    # usampled_xidxs_2d = np.zeros((2*ny, 2*nx))
    # usampled_yidxs_2d = np.zeros((2*ny, 2*nx))

    # Update p, q;
    # (b) p_ij = NINT(2 dt/dx vchiold_x[x_i - p*dx/2, y_j - q*dy/2])
    #     q_ij = NINT(2 dt/dy vchiold_y[x_i - p*dx/2, y_j - q*dy/2])

    # The idxs are different between the "normal" (i.e. non-upsampled) arrays
    # and the upsampled arrays. Get the correct upsampled idx for
    # (x_i - p*dx/2), (y_j - q*dy/2) here.
    #
    # Each (x,y) gridpoint has an integer value of phalf and qhalf; from this
    # we want to get the idxs fro the upsampled quantities; the t^n quantities
    # have twice as many gridpoints, so for a particular (xidx, yidx, phalf, qhalf)
    # the corresponding idxs of the upsampled data points are
    # xidx_for_upsampled_data = 2*(xidx - pidx/2) = 2*xidx - pidx
    # yidx_for_upsampled_data = 2*(yidx - qidx/2) = 2*yidx - qidx

    xidx_for_upsampled_array = (2*x_idxs_2d - p_array)%(2*nx)
    yidx_for_upsampled_array = (2*y_idxs_2d - q_array)%(2*ny)
    # We've got the idxs - update p, q
    for counter in range(0, 1000, 1):
        #print("counter = ", counter)
        p_array, q_array, yidx_for_upsampled_array, xidx_for_upsampled_array = update_p_and_q(
                    p_array, q_array, vchiold_x_array_upsampled, vchiold_y_array_upsampled,
                    yidx_for_upsampled_array, xidx_for_upsampled_array)
        # print("xidx_for_upsampled_array[yidx, xidx], yidx_for_upsampled_array[yidx, xidx] = ",
        #         xidx_for_upsampled_array[yidx, xidx], yidx_for_upsampled_array[yidx, xidx])
        # print("p_array[yidx, xidx], q_array[yidx, xidx] = ", p_array[yidx, xidx], q_array[yidx, xidx])

    xidx_for_norm_array = (x_idxs_2d - p_array)%nx
    yidx_for_norm_array = (y_idxs_2d - q_array)%ny
    if with_diagnostics:
        print("p_array = ", p_array)  # For the stella-given quantities, vchi small so p=q=0 everywhere.
        print("q_array = ", q_array)

    # Calculate the residual velocities
    vchiresidual_x = vchiold_x_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array] - p_array*dx/(2*dt)
    vchiresidual_y = vchiold_y_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array] - q_array*dy/(2*dt)

    # print("vchiold_x, p, q, p_array*dx/(2*dt), vchiresidual_x = ",
    #         vchiold_x_array_upsampled[yidx_for_upsampled_array[yidx,xidx], xidx_for_upsampled_array[yidx,xidx]],
    #         p_array[yidx, xidx], q_array[yidx, xidx], p_array[yidx, xidx]*dx/(2*dt),
    #         vchiresidual_x[yidx, xidx])

    # for yidx in y_idxs:
    #     for xidx in x_idxs:
    #
    #         print("vchiold_x, p, q, p_array*dx/(2*dt), vchiresidual_x = ",
    #                 vchiold_x_array_upsampled[yidx_for_upsampled_array[yidx,xidx], xidx_for_upsampled_array[yidx,xidx]],
    #                 p_array[yidx, xidx], q_array[yidx, xidx], p_array[yidx, xidx]*dx/(2*dt),
    #                 vchiresidual_x[yidx, xidx])
    Courant_num_array = (vchiold_x_array*dt/dx + vchiold_y_array*dt/dy)
    Courant_residual_array = (vchiresidual_x*dt/dx + vchiresidual_y*dt/dy)
    if with_diagnostics:
        print("max Courant no = ", np.max(abs(Courant_num_array)))
        print("max residual Courant no = ", np.max(abs(Courant_residual_array)))
    #print("dx/dt, max(vchiresidual_x), dy/dt, max(vchiresidual_y) = ", dx/dt, np.max(vchiresidual_x), dy/dt, np.max(vchiresidual_y))
    # Calculate rhs_ij
    #  rhs_ij = - (vresidual_x_ij*dgold_dx[x_i - p_ij*dx/2, y_j - * q_ij*dy/2]
    #            + vresidual_y_ij*dgold_dy[x_i - p_ij*dx/2, y_j - * q_ij*dy/2])
    rhs_array = - (vchiresidual_x * dgold_dy_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array]
                    + vchiresidual_y * dgold_dx_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array] )

    # Calculate gnew
    gnewyx_array = golderyx_array[yidx_for_norm_array, xidx_for_norm_array] + rhs_array

    if with_diagnostics:
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.imshow(golderyx_array)
        ax2.imshow(gnewyx_array)
        ax1.set_ylabel("x idx")
        ax2.set_ylabel("x idx")
        ax1.set_xlabel("y idx")
        ax2.set_xlabel("y idx")
        plt.show()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(y_grid, golderyx_array[:,5])
        ax1.plot(y_grid, gnewyx_array[:,5])
        print("x_grid[5] = ", x_grid[5] )
        ax1.set_ylabel("g")
        ax1.set_xlabel("y")
        plt.show()

    return gnewyx_array

def benchmark_rk2_vs_nisl_different_velocity_fields_spatially_constant_temporally_varying():
    """ """

    [golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
        vchiold_y_array, vchiold_x_array] = get_arrays_from_nonlinear_data()

    ## Let's try artificially increasing vchi to exceed the CFL condition.
    scaling_fac = 1
    vchiold_x_array = vchiold_x_array * scaling_fac
    vchiold_y_array = vchiold_y_array * scaling_fac

    # ALternative scheme - spatially constant advecting velocity.
    v_x = dx/dt * 0.123456
    constant_vchi_x = np.ones((ny, nx)) * v_x
    v_y = dy/dt * 0.434123
    constant_vchi_y = np.ones((ny, nx)) * v_y

    # Another alternative - sinusoidally varying velocity
    vx_sine_amp = 1.0
    vy_sine_amp = 1.0
    vx_sinusoidal = np.zeros((ny, nx))
    vy_sinusoidal = np.zeros((ny, nx))
    for yidx in range(0, ny):
        vx_sinusoidal[yidx,:] = vx_sine_amp*np.cos(x_grid*2*np.pi/xmax)
    for xidx in range(0, nx):
        vy_sinusoidal[:, xidx] = vy_sine_amp*np.cos(y_grid*2*np.pi/ymax)

    # Another alternative - sinusoidally varying velocity, but the "other way" arounf
    # Appears to work better (g conserves better)
    vx_sine_amp_other_way = 1.0
    vy_sine_amp_other_way = 1.0
    vx_sinusoidal_other_way = np.zeros((ny, nx))
    vy_sinusoidal_other_way = np.zeros((ny, nx))
    for xidx in range(0, nx):
        vx_sinusoidal_other_way[:,xidx] = vx_sine_amp_other_way*np.cos(y_grid*2*np.pi/ymax)
    for yidx in range(0, ny):
        vy_sinusoidal_other_way[yidx,:] = vy_sine_amp_other_way*np.cos(x_grid*2*np.pi/xmax)

    # Another alternative - sinusoidally varying velocity, but in both dimensions
    vx_sine_amp_both = 1.0
    vy_sine_amp_both = 1.0
    vx_sinusoidal_both = np.zeros((ny, nx))
    vy_sinusoidal_both = np.zeros((ny, nx))
    for xidx in range(0, nx):
        for yidx in range(0, ny):
            vx_sinusoidal_both[yidx,xidx] = vx_sine_amp_both*np.cos(y_grid[yidx]*2*np.pi/ymax)*np.cos(x_grid[xidx]*2*np.pi/xmax)
            vy_sinusoidal_both[yidx,xidx] = vy_sine_amp_both*np.cos(y_grid[yidx]*2*np.pi/ymax)*np.cos(x_grid[xidx]*2*np.pi/xmax)
    take_many_steps(golderyx_array, constant_vchi_y, constant_vchi_x, nstep=1000, schemes=["RK2", "NISL"],
                    time_varying_velocity=True, velocity_omega=1)
    take_many_steps(golderyx_array, vchiold_y_array, vchiold_x_array, nstep=500, schemes=["RK2", "NISL"],
                    time_varying_velocity=True, velocity_omega=1)
    take_many_steps(golderyx_array, vy_sinusoidal, vx_sinusoidal, nstep=1000, schemes=["RK2", "NISL"],
                    time_varying_velocity=True, velocity_omega=1)
    take_many_steps(golderyx_array, vy_sinusoidal_both, vx_sinusoidal_both, nstep=1000, schemes=["RK2", "NISL"],
                    time_varying_velocity=True, velocity_omega=1)
    take_many_steps(golderyx_array, vy_sinusoidal_other_way, vx_sinusoidal_other_way, nstep=1000, schemes=["RK2", "NISL"],
                    time_varying_velocity=True, velocity_omega=1)
    take_many_steps(golderyx_array, np.zeros((ny, nx)), vchiold_x_array, nstep=1000, schemes=["RK2", "NISL"],
                    time_varying_velocity=True, velocity_omega=1)
    print("Finished")

    return

def get_approx_departure_point(yidx, xidx, vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                n_timesteps=2):
    """Find the approximate departure point which takes us from an (in general,
    off-grid) point at t^{n-1} to the gridpoint labelled by (xidx, yidx) at
    t^{n+1}.

    A complication is that we want not only the values of (x,y) within our
    periodic grid, but we'll also need the advecting velocity. To this end, we want to
    get the depatrure point on the non-periodic grid, rather than the periodic
    grid.  """
    upsampled_xidx = 2*xidx
    upsampled_yidx = 2*yidx

    # location of the gridpoint (=arrival point)
    x_nonperiodic = x_grid[xidx]; y_nonperiodic=y_grid[yidx]
    # Normalised values; these are integer at cell boundaries.
    time_remaining = n_timesteps*dt
    xhistory = []
    yhistory = []
    xhistory.append(x_nonperiodic); yhistory.append(y_nonperiodic)

    # Velocities in the initial cell - put here rather than in the loop,
    # because need to update velocities carefully (the "divergent velocities")
    # problem.
    # NB these are the velocities of the cells, but we'll actually be moving
    # back in time.
    u_x = vchiold_x_array_upsampled[upsampled_yidx, upsampled_xidx]
    u_y = vchiold_y_array_upsampled[upsampled_yidx, upsampled_xidx]

    max_iterations = max((100*np.mean(abs(vchiold_x_array_upsampled))*dt/dx) , (100*np.mean(abs(vchiold_y_array_upsampled))*dt/dy))
    #print("max_iterations = ", max_iterations)
    counter=0
    while time_remaining > 0 and counter < max_iterations:
        counter +=1
        ### TODO: Find something sensible to do if the velocities are
        ### divergent (pushing between boxes)
        # Velocity with which we move (although we'll actually be moving
        # back in time.)

        # Take a step in (x,y) based on the velocity at [x, y, t^n].
        xnorm = (2*x_nonperiodic/dx - 0.5); ynorm = (2*y_nonperiodic/dy - 0.5)
        # Our boundaries are nonperiodic.
        if u_x > 0 :
            xboundary = ((np.floor(xnorm) + 0.5) * dx/2)
        else:
            xboundary = ((np.ceil(xnorm)  + 0.5) * dx/2)
        if u_y > 0 :
            yboundary = ((np.floor(ynorm) + 0.5) * dy/2)
        else:
            yboundary = ((np.ceil(ynorm)  + 0.5) * dy/2)

        # print("upsampled_xidx, upsampled_yidx = ", upsampled_xidx, upsampled_yidx )
        # print("xnorm, u_x, xboundary, ynorm, u_y, yboundary = ", xnorm, u_x, xboundary, ynorm, u_y, yboundary)
        # Calculate the time required to reach the nearest boundaries in y and x
        # Careful though - if our velocities are too small we'll get <>dist_dt=inf
        min_vx_magnitude = (dx/dt)*1e-10   # We know it's very unlikely that the trajectory with this velocity will make it into the next cell
        min_vy_magnitude = (dy/dt)*1e-10   # We know it's very unlikely that the trajectory with this velocity will make it into the next cell

        if abs(u_y) < min_vy_magnitude:
            ydist_dt = 10*dt    # This ensures ydist_dt > time_remaining
        else:
            ydist_dt = -(yboundary - y_nonperiodic)/u_y
        if abs(u_x) < min_vx_magnitude:
            xdist_dt = 10*dt # This ensures xdist_dt > time_remaining
        else:
            xdist_dt = -(xboundary - x_nonperiodic)/u_x

        #print("ydist_dt, xdist_dt = ", ydist_dt, xdist_dt )

        if (ydist_dt > time_remaining) and (xdist_dt > time_remaining):
            #print("STOPPING before next boundary")
            # Update location
            y_nonperiodic = (y_nonperiodic - u_y * time_remaining)
            x_nonperiodic = (x_nonperiodic - u_x * time_remaining)
            time_remaining = 0
        else:
            if ydist_dt < xdist_dt:
                #print("Hit y!")
                # We've hit the boundary in y - update the non-periodic
                # position of the particle and the periodic position of the
                # particle.
                y_nonperiodic = (yboundary - u_y*very_small_dt )   # slightly overstep, so we're definitely in the next cell
                x_nonperiodic = (x_nonperiodic - u_x * (ydist_dt + very_small_dt))

                # Update the values of u_x, u_y
                if u_y > 0:
                    # +ve u_y, so going back in time we're going in the -ve y direction; our new cell is
                    # below the old cell
                    upsampled_yidx = (upsampled_yidx - 1)%(2*ny)
                else:
                    # -ve u_y, so going back in time we're going in the +ve y direction; our new cell is
                    # below the old cell
                    upsampled_yidx = (upsampled_yidx + 1)%(2*ny)

                # Update the velocities
                u_x = vchiold_x_array_upsampled[upsampled_yidx, upsampled_xidx]
                # Update u_y, but if the sign is different, we're going to
                # bounce back and forth (this indicates the velocity falling
                # to zero somewhere between the 2 cell centres). To avoid the "bouncing",
                # set velocity to zero.
                u_y_new = vchiold_y_array_upsampled[upsampled_yidx, upsampled_xidx]
                if (u_y * u_y_new) < 0:
                    # Opposite signs, so set u_y=0
                    u_y=0
                else:
                    u_y = u_y_new

                # Update time_remaining. Include the "very small dt" contribution
                time_remaining = time_remaining - (ydist_dt + very_small_dt)
            else:
                #print("Hit x!")
                # Hit the boundary in x
                x_nonperiodic = (xboundary - u_x*very_small_dt)    # slightly overstep, so we're definitely in the next cell
                y_nonperiodic = (y_nonperiodic - u_y * (xdist_dt + very_small_dt))

                # Update the values of u_x, u_y
                if u_x > 0:
                    # +ve u_x, so going back in time we're going in the -ve x direction; our new cell is
                    # to the left of the old cell
                    upsampled_xidx = (upsampled_xidx - 1)%(2*nx)
                else:
                    # -ve u_y, so going back in time we're going in the +ve y direction; our new cell is
                    # below the old cell
                    upsampled_xidx = (upsampled_xidx + 1)%(2*nx)

                # Update velocities
                u_y = vchiold_y_array_upsampled[upsampled_yidx, upsampled_xidx]
                # Update u_x, but if the sign is different, we're going to
                # bounce back and forth (this indicates the velocity falling
                # to zero somewhere between the 2 cell centres). To avoid the "bouncing",
                # set velocity to zero.
                u_x_new = vchiold_x_array_upsampled[upsampled_yidx, upsampled_xidx]
                if (u_x * u_x_new) < 0:
                    # Opposite signs, so set u_y=0
                    u_x=0
                else:
                    u_x = u_x_new

                # Update time_remaining. Include the "very small dt" contribution
                time_remaining = time_remaining - (xdist_dt + very_small_dt)

        xhistory.append(x_nonperiodic); yhistory.append(y_nonperiodic)

    ########################################################################
    ##### DIAGNOSTIC PLOTS #################################################
    ########################################################################
    # def basic_diagnostic_plot_trajectories():
    #     """The "vanilla" plot - show paths and gridpoints. """
    #     marker_size = 20.
    #     fig = plt.figure(figsize=[12, 8])
    #     ax1 = fig.add_subplot(111)
    #     ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=marker_size, label="upsampled grid")
    #     ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
    #     ax1.scatter(xhistory, yhistory, s=marker_size, label="trajectory")
    #     for hist_idx in range(0, len(xhistory)-1):
    #         ax1.plot([xhistory[hist_idx], xhistory[hist_idx+1]], [yhistory[hist_idx], yhistory[hist_idx+1]])
    #     ax1.set_xlabel(r"$x$")
    #     ax1.set_ylabel(r"$y$")
    #     ax1.grid(True)
    #     #ax1.set_xlim([-1, 23])
    #     ax1.legend(loc="upper right")
    #     plt.show()

    ########################################################################

    return y_nonperiodic, x_nonperiodic, yhistory, xhistory

def plot_grids_and_trajectories_for_thesis():
    """ """

    [golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
        vchiold_y_array, vchiold_x_array] = get_arrays_from_nonlinear_data()
    scaling_factor = 30 # 30 for thesis
    vchiold_y_array = vchiold_y_array * scaling_factor
    vchiold_x_array = vchiold_x_array * scaling_factor
    vchiold_x_array_upsampled = create_upsampled_grid(vchiold_x_array)
    vchiold_y_array_upsampled = create_upsampled_grid(vchiold_y_array)

    approx_departure_points_x = np.zeros((ny, nx))
    approx_departure_points_y = np.zeros((ny, nx))
    xhist_list = []
    yhist_list = []
    # Find the approximate departure point
    for yidx in range(0, ny):
        for xidx in range(0, nx):
            y, x, yhist, xhist = get_approx_departure_point(yidx, xidx,
                                        vchiold_y_array_upsampled, vchiold_x_array_upsampled)
            approx_departure_points_x[yidx, xidx] = x
            approx_departure_points_y[yidx, xidx] = y
            xhist_list.append(xhist)
            yhist_list.append(yhist)
    Courant_num_array = (vchiold_x_array*dt/dx + vchiold_y_array*dt/dy)
    print("max(abs(courant_number)) =  ", np.max(abs(Courant_num_array)))
    fig = diagnostic_plot_trajectories(vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                    approx_departure_points_y, approx_departure_points_x,
                                    yhist_list, xhist_list, which_plot="0")
    plt.savefig("images/trajectory_plot_for_thesis_0.eps")
    plt.close()
    fig = diagnostic_plot_trajectories(vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                    approx_departure_points_y, approx_departure_points_x,
                                    yhist_list, xhist_list, which_plot="1")
    plt.savefig("images/trajectory_plot_for_thesis_1.eps")
    plt.close()
    fig = diagnostic_plot_trajectories(vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                    approx_departure_points_y, approx_departure_points_x,
                                    yhist_list, xhist_list, which_plot="2")
    plt.savefig("images/trajectory_plot_for_thesis_2.eps")
    plt.close()
    fig = diagnostic_plot_trajectories(vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                    approx_departure_points_y, approx_departure_points_x,
                                    yhist_list, xhist_list, which_plot="3")
    plt.savefig("images/trajectory_plot_for_thesis_c3p47_3.eps")
    plt.close()

    fig = diagnostic_plot_ritchie(vchiold_y_array_upsampled, vchiold_x_array_upsampled, which_plot="1")
    plt.savefig("images/trajectory_plot_for_thesis_c3p47_rithcie_1.eps")
    fig = diagnostic_plot_ritchie(vchiold_y_array_upsampled, vchiold_x_array_upsampled, which_plot="2")
    plt.savefig("images/trajectory_plot_for_thesis_c3p47_rithcie_2.eps")
    plt.close()
    return

def diagnostic_plot_trajectories(#x_grid_upsampled, y_grid_upsampled,
                                vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                approx_departure_points_y, approx_departure_points_x,
                                yhist_list, xhist_list,
                                #nx, ny, dx, dy,
                                which_plot="1"
                                ):
    """ For a single "step" (i.e. a single velocity field), make one of 3 possible plots:
    (0) A grid of cells showing gridpoints + upsampled gridpoints
    (1) A grid of cells showing cell boundaries + velocities
    (2) A grid of cells showing gridpoints + cell boundaries + trajectories.
    (3) A grid of cells showing gridpoints + depature points"""

    x_grid_upsampled_boundaries = x_grid_upsampled + dx/4
    y_grid_upsampled_boundaries = y_grid_upsampled + dy/4

    ################ Initialise plot #########################
    arrow_head_width = 0.1
    left = 0.13
    right = 0.97
    top = 0.97
    bottom = 0.1
    width = right - left
    height = top - bottom

    xlabel_fontsize = 40
    ylabel_fontsize = xlabel_fontsize
    xticklabel_fontsize = 30
    yticklabel_fontsize = xticklabel_fontsize
    legend_fontsize = 20

    grid_markersize=80
    foot_markersize=80
    upsampled_grid_markersize = 60

    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes((left, bottom, width, height))

    def make_gridpoints():
        """ """
        ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=grid_markersize, label="(x,y) grid", zorder=10)
        return

    def make_upsampled_gridpoints():
        """ """
        ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=upsampled_grid_markersize, label="(x,y) upsampled grid", zorder=00)
        return

    def make_grids():
        """ """
        ################ Make grids ##############################
        # To make the horizontal lines: make 1 horizontal line per
        # y_grid_upsampled_boundaries, starting at x=0 and ending at max(x_grid_upsampled_boundaries)
        horizontal_lines = []
        for diag_yval in y_grid_upsampled_boundaries:
            horizontal_line_xvals = [0, max(x_grid_upsampled_boundaries)]
            horizontal_line_yvals = [diag_yval, diag_yval]
            horizontal_lines.append([horizontal_line_xvals, horizontal_line_yvals])

        # To make the vertical lines: make 1 vertical line per
        # x_grid_upsampled_boundaries, starting at y=0 and ending at max(y_grid_upsampled_boundaries)
        vertical_lines = []
        for diag_xval in x_grid_upsampled_boundaries:
            vertical_line_xvals = [diag_xval, diag_xval]
            vertical_line_yvals = [0, max(y_grid_upsampled_boundaries)]
            vertical_lines.append([vertical_line_xvals, vertical_line_yvals])

        ax1.plot(horizontal_lines[0][0], horizontal_lines[0][1], ls="--", c="gray", label="cell boundary")
        for horizontal_line in horizontal_lines:
            ax1.plot(horizontal_line[0], horizontal_line[1], ls="--", c="gray")
        for vertical_line in vertical_lines:
            ax1.plot(vertical_line[0], vertical_line[1], ls="--", c="gray")

        ######################################################################

        return

    def make_velocities():
        """ """
        ################ Make velocity arrows ##############################
        # Normalise velocities such that the largest veloccity occupies a cell length/height.
        # That is, want max(unorm_x) = dx/2 and max(unorm_y) = dy/2
        x_scaling_fac = dx/2 / np.max(abs(vchiold_x_array_upsampled))
        unorm_x = vchiold_x_array_upsampled * x_scaling_fac
        y_scaling_fac = dy/2 / np.max(abs(vchiold_y_array_upsampled))
        unorm_y = vchiold_y_array_upsampled * y_scaling_fac
        # Want to represent these velocities with an arrow, which is centered on
        # the gridpoint. So the starting point of the arrow should be [x - unorm_x/2, y - unorm_y/2]
        arrows = [] # Each item in arrows is a list describing a single arrow; [x, y, delta_x, delta_y]
        for diag_upsampled_xidx in range(0, 2*nx):
            for diag_upsampled_yidx in range(0, 2*ny):
                arrow_x = x_grid_2d_upsampled[diag_upsampled_yidx, diag_upsampled_xidx] - unorm_x[diag_upsampled_yidx, diag_upsampled_xidx]/2
                arrow_y = y_grid_2d_upsampled[diag_upsampled_yidx, diag_upsampled_xidx] - unorm_y[diag_upsampled_yidx, diag_upsampled_xidx]/2
                arrow_dx = unorm_x[diag_upsampled_yidx, diag_upsampled_xidx]
                arrow_dy = unorm_y[diag_upsampled_yidx, diag_upsampled_xidx]
                arrows.append([arrow_x, arrow_y, arrow_dx, arrow_dy])

        ax1.plot([-102, -101], [-101, -102], color="blue", marker=">", markersize=14,
                  label="velocity")
        for arrow in arrows:
            ax1.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)

        ####################################################################

    def make_trajectories():
        """ """
        ################## Make trajectories ###############################
        for trajectory_idx in range(len(xhist_list)):
            xhistory = xhist_list[trajectory_idx]
            yhistory = yhist_list[trajectory_idx]

            #ax1.scatter(xhistory, yhistory, s=marker_size, label="trajectory")
            for hist_idx in range(0, len(xhistory)-1):
                ax1.plot([xhistory[hist_idx], xhistory[hist_idx+1]], [yhistory[hist_idx], yhistory[hist_idx+1]],
                         c="red")
        ax1.plot([-100, -99], [-100, -99], label="trajectory",
                 c="red")

        ####################################################################

    def make_departure_points():
        """ """
        ################## Departure points ################################
        ax1.scatter(approx_departure_points_x, approx_departure_points_y, marker="X", c="green", label="trajectory foot", s=foot_markersize)
        ####################################################################


    if which_plot == "0":
        make_gridpoints()
        make_upsampled_gridpoints()
        ax1.legend(loc="upper center", fontsize=legend_fontsize, ncol=2).set_zorder(30)
    elif which_plot == "1":
        make_grids()
        make_velocities()
        ax1.legend(loc="upper center", fontsize=legend_fontsize, ncol=2).set_zorder(30)
    elif which_plot == "2":
        make_gridpoints()
        make_grids()
        make_trajectories()
        ax1.legend(loc="upper center", fontsize=legend_fontsize, ncol=3).set_zorder(30)
    elif which_plot == "3":
        make_gridpoints()
        make_departure_points()
        ax1.legend(loc="upper center", fontsize=legend_fontsize, ncol=2).set_zorder(30)
    else:
        print("which_plot not recognised!")

    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    ax1.set_xlim([-2.5, 20.5])
    ax1.set_ylim([-1, 20.5])
    ax1.set_xlabel(r"$\tilde{x}$", fontsize=xlabel_fontsize)
    ax1.set_ylabel(r"$\tilde{y}$", fontsize=ylabel_fontsize)
    ax1.set_xticks([0, 5, 10, 15, 20])
    ax1.set_xticklabels([r"$0$", r"$5$", r"$10$", r"$15$", r"$20$"], fontsize=xticklabel_fontsize)
    ax1.set_yticks([0, 5, 10, 15])
    ax1.set_yticklabels([r"$0$", r"$5$", r"$10$", r"$15$"], fontsize=xticklabel_fontsize)


    return fig

def diagnostic_plot_ritchie(#x_grid_upsampled, y_grid_upsampled,
                                vchiold_y_array_upsampled, vchiold_x_array_upsampled,
                                #nx, ny, dx, dy,
                                which_plot="1"
                                ):
    """ For a single "step" (i.e. a single velocity field), make one of 3 possible plots:
    (1) upsampled gridpoints + trajectories according to Ritchie
    (2) gridpoints + departure points + arrival points
    """
    x_grid_upsampled_boundaries = x_grid_upsampled + dx/4
    y_grid_upsampled_boundaries = y_grid_upsampled + dy/4

    ################ Initialise plot #########################
    arrow_head_width = 0.15
    left = 0.13
    right = 0.97
    top = 0.97
    bottom = 0.1
    width = right - left
    height = top - bottom

    xlabel_fontsize = 40
    ylabel_fontsize = xlabel_fontsize
    xticklabel_fontsize = 30
    yticklabel_fontsize = xticklabel_fontsize
    legend_fontsize = 20

    grid_markersize=80
    foot_markersize=80
    upsampled_grid_markersize = 60

    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes((left, bottom, width, height))

    def make_gridpoints():
        """ """
        ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=grid_markersize, label="(x,y) grid", zorder=10)
        return

    def make_upsampled_gridpoints():
        """ """
        ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=upsampled_grid_markersize, label="(x,y) upsampled grid", zorder=00)
        return

    def make_trajectories():
        """ """
        ################## Make trajectories ###############################
        arrows = []

        for xidx in range(0, 2*nx):
            for yidx in range(0, 2*ny):
                arrow_x = x_grid_2d_upsampled[yidx, xidx] - vchiold_x_array_upsampled[yidx, xidx]*dt
                arrow_y = y_grid_2d_upsampled[yidx, xidx] - vchiold_y_array_upsampled[yidx, xidx]*dt
                arrow_dx = vchiold_x_array_upsampled[yidx, xidx]*2*dt
                arrow_dy = vchiold_y_array_upsampled[yidx, xidx]*2*dt
                arrows.append([arrow_x, arrow_y, arrow_dx, arrow_dy])

        ax1.plot([-102, -101], [-101, -102], color="blue", marker=">", markersize=14,
                  label=r"$v_{E\times B}.2\Delta\tilde{t}$")
        for arrow in arrows:
            ax1.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)


        ####################################################################

    def make_departure_and_arrival_points():
        """ """
        departure_x = []
        departure_y = []
        arrival_x = []
        arrival_y = []
        for xidx in range(0, 2*nx):
            for yidx in range(0, 2*ny):
                departure_x.append(x_grid_2d_upsampled[yidx, xidx] - vchiold_x_array_upsampled[yidx, xidx]*dt)
                departure_y.append(y_grid_2d_upsampled[yidx, xidx] - vchiold_y_array_upsampled[yidx, xidx]*dt)
                arrival_x.append(x_grid_2d_upsampled[yidx, xidx] + vchiold_x_array_upsampled[yidx, xidx]*dt)
                arrival_y.append(y_grid_2d_upsampled[yidx, xidx] + vchiold_y_array_upsampled[yidx, xidx]*dt)
        ################## Departure points ################################
        ax1.scatter(departure_x, departure_y, marker="X", c="red", label="trajectory foot", s=foot_markersize)
        ax1.scatter(arrival_x, arrival_y, marker="X", c="green", label="trajectory arrival point", s=foot_markersize)
        ####################################################################

    if which_plot == "1":
        #make_upsampled_gridpoints()
        make_trajectories()
        ax1.legend(loc="upper center", fontsize=legend_fontsize, ncol=2).set_zorder(30)
    elif which_plot == "2":
        make_gridpoints()
        make_departure_and_arrival_points()
        ax1.legend(loc="upper center", fontsize=legend_fontsize, ncol=3).set_zorder(30)
    else:
        print("which_plot not recognised!")

    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    ax1.set_xlim([-2.5, 20.5])
    ax1.set_ylim([-1, 20.5])
    ax1.set_xlabel(r"$\tilde{x}$", fontsize=xlabel_fontsize)
    ax1.set_ylabel(r"$\tilde{y}$", fontsize=ylabel_fontsize)
    ax1.set_xticks([0, 5, 10, 15, 20])
    ax1.set_xticklabels([r"$0$", r"$5$", r"$10$", r"$15$", r"$20$"], fontsize=xticklabel_fontsize)
    ax1.set_yticks([0, 5, 10, 15])
    ax1.set_yticklabels([r"$0$", r"$5$", r"$10$", r"$15$"], fontsize=xticklabel_fontsize)


    return fig

if __name__ == "__main__":
    plot_grids_and_trajectories_for_thesis()
