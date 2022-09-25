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

def get_dgdx_and_dgdy(g):
    """Given g(y,x), calculate dg/dx and dg/dy using a space-centered second order
    scheme."""

    gleft = np.roll(g, 1, axis=1); gright = np.roll(g, -1, axis=1)
    dgdx = (gright - gleft)/(2*dx)
    gdown = np.roll(g, 1, axis=0); gup = np.roll(g, -1, axis=0)
    dgdy = (gup - gdown)/(2*dy)

    return dgdy, dgdx

def nisl_step_finding_departure_point(golderyx_array, dgold_dy_array, dgold_dx_array,
                vchiold_y_array, vchiold_x_array, n_timesteps=2):
    """Take a NISL step, but finding the "actual" (approximate) departure point.
    Guaranteed to have good stability properties """

    very_small_dt = 1e-8

    p_array = np.zeros((ny, nx), dtype="int")
    q_array = np.zeros((ny, nx), dtype="int")

    # Upsample - double the number of points in the vchi, gold, golder grids.
    # golderyx_array_upsampled = create_upsampled_grid(golderyx_array)

    dgold_dy_array_upsampled = create_upsampled_grid(dgold_dy_array)
    dgold_dx_array_upsampled = create_upsampled_grid(dgold_dx_array)
    vchiold_x_array_upsampled = create_upsampled_grid(vchiold_x_array)
    vchiold_y_array_upsampled = create_upsampled_grid(vchiold_y_array)
#################
# [Plots in here I think
####################
    def get_approx_departure_point(yidx, xidx):
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
        def basic_diagnostic_plot_trajectories():
            """The "vanilla" plot - show paths and gridpoints. """
            marker_size = 20.
            fig = plt.figure(figsize=[12, 8])
            ax1 = fig.add_subplot(111)
            ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=marker_size, label="upsampled grid")
            ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
            ax1.scatter(xhistory, yhistory, s=marker_size, label="trajectory")
            for hist_idx in range(0, len(xhistory)-1):
                ax1.plot([xhistory[hist_idx], xhistory[hist_idx+1]], [yhistory[hist_idx], yhistory[hist_idx+1]])
            ax1.set_xlabel(r"$x$")
            ax1.set_ylabel(r"$y$")
            ax1.grid(True)
            #ax1.set_xlim([-1, 23])
            ax1.legend(loc="upper right")
            plt.show()

        def diagnostic_plot_trajectories():
            """ A more complicated plot - show paths and gridpoints, and boundaries and
            cell velocities."""
            x_grid_upsampled_boundaries = x_grid_upsampled + dx/4
            y_grid_upsampled_boundaries = y_grid_upsampled + dy/4

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


            marker_size = 20.
            arrow_head_width = 0.1
            fig = plt.figure(figsize=[12, 8])
            ax1 = fig.add_subplot(111)
            ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=marker_size, label="upsampled grid")
            ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
            ax1.scatter(xhistory, yhistory, s=marker_size, label="trajectory")
            for horizontal_line in horizontal_lines:
                ax1.plot(horizontal_line[0], horizontal_line[1], ls="--", c="gray")
            for vertical_line in vertical_lines:
                ax1.plot(vertical_line[0], vertical_line[1], ls="--", c="gray")
            for arrow in arrows:
                ax1.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)
            for hist_idx in range(0, len(xhistory)-1):
                ax1.plot([xhistory[hist_idx], xhistory[hist_idx+1]], [yhistory[hist_idx], yhistory[hist_idx+1]])
            ax1.set_xlabel(r"$x$")
            ax1.set_ylabel(r"$y$")
            #ax1.grid(True)
            #ax1.set_xlim([-1, 23])
            ax1.legend(loc="upper right")
            plt.show()

            return

        diagnostic_plot_trajectories()
        ########################################################################

        return y_nonperiodic, x_nonperiodic

    approx_departure_points_x = np.zeros((ny, nx))
    approx_departure_points_y = np.zeros((ny, nx))

    # Find the approximate departure point
    for yidx in range(0, ny):
        for xidx in range(0, nx):
            y, x = get_approx_departure_point(yidx, xidx)
            approx_departure_points_x[yidx, xidx] = x
            approx_departure_points_y[yidx, xidx] = y

    ############################################################################
    ### DIAGNSOTIC PLOT
    ############################################################################

    def make_diagnostic_plot_departure_points():
        """ """
        marker_size = 20.
        arrow_head_width = 0.1
        fig = plt.figure(figsize=[12, 8])
        ax1 = fig.add_subplot(111)
        ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=marker_size, label="upsampled grid")
        ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
        ax1.scatter(approx_departure_points_x.flatten(), approx_departure_points_y.flatten(), s=marker_size,label="departure point")
        # ax1.scatter(xhistory, yhistory, s=marker_size, label="trajectory")
        # for horizontal_line in horizontal_lines:
        #     ax1.plot(horizontal_line[0], horizontal_line[1], ls="--", c="gray")
        # for vertical_line in vertical_lines:
        #     ax1.plot(vertical_line[0], vertical_line[1], ls="--", c="gray")
        # for arrow in arrows:
        #     ax1.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)
        # for hist_idx in range(0, len(xhistory)-1):
        #     ax1.plot([xhistory[hist_idx], xhistory[hist_idx+1]], [yhistory[hist_idx], yhistory[hist_idx+1]])
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        ax1.grid(True)
        #ax1.set_xlim([-0.96, 25])
        #ax1.set_xlim([-1, 23])
        ax1.legend(loc="upper right")
        plt.show()

        return

    #make_diagnostic_plot_departure_points()
    ############################################################################

    # Now we've got approximate departure points, work out p and q.
    # We've got an approximate departure point; first calculate a velocity which
    # takes us from the approximate departure point to the gridpoint.
    velocity_array_x = (x_grid_2d - approx_departure_points_x)/(n_timesteps*dt)
    velocity_array_y = (y_grid_2d - approx_departure_points_y)/(n_timesteps*dt)

    p_array = np.rint((x_grid_2d - approx_departure_points_x)/dx).astype("int")
    q_array = np.rint((y_grid_2d - approx_departure_points_y)/dy).astype("int")
    #print("max (p, q) = ", max(np.max(abs(p_array)), np.max(abs(q_array))))
    # print("q_array = ", q_array)

    xidx_for_norm_array = (x_idxs_2d - p_array)%nx
    yidx_for_norm_array = (y_idxs_2d - q_array)%ny
    xidx_for_upsampled_array = (2*x_idxs_2d - p_array)%(2*nx)
    yidx_for_upsampled_array = (2*y_idxs_2d - q_array)%(2*ny)

    vchiresidual_x = velocity_array_x - p_array*dx/(n_timesteps*dt)
    vchiresidual_y = velocity_array_y - q_array*dy/(n_timesteps*dt)

    Courant_num_array = (vchiold_x_array*dt/dx + vchiold_y_array*dt/dy)
    Courant_residual_array = (vchiresidual_x*dt/dx + vchiresidual_y*dt/dy)
    # print("max Courant no = ", np.max(abs(Courant_num_array)))
    # print("max residual Courant no = ", np.max(abs(Courant_residual_array)))
    #print("dx/dt, max(vchiresidual_x), dy/dt, max(vchiresidual_y) = ", dx/dt, np.max(vchiresidual_x), dy/dt, np.max(vchiresidual_y))
    # Calculate rhs_ij
    #  rhs_ij = - (vresidual_x_ij*dgold_dx[x_i - p_ij*dx/2, y_j - * q_ij*dy/2]
    #            + vresidual_y_ij*dgold_dy[x_i - p_ij*dx/2, y_j - * q_ij*dy/2]) *dt (or 2dt)
    if n_timesteps == 1:
        ## The single-step version:
        rhs_array = - dt* (vchiresidual_x * dgold_dx_array[yidx_for_norm_array, xidx_for_norm_array]
                        + vchiresidual_y * dgold_dy_array[yidx_for_norm_array, xidx_for_norm_array] )
    elif n_timesteps == 2:
        # print("dgold_dx_array = ", dgold_dx_array)
        # print("dgold_dx_array[yidx_for_norm_array, xidx_for_norm_array] = ", dgold_dx_array[yidx_for_norm_array, xidx_for_norm_array])
        # sys.exit()
        ## The 3-step version
        rhs_array = - 2*dt*(vchiresidual_x * dgold_dx_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array]
                + vchiresidual_y * dgold_dy_array_upsampled[yidx_for_upsampled_array, xidx_for_upsampled_array] )
    else:
        print("In trouble! n_timesteps has unexpected value")
        print("n_timesteps = ", n_timesteps)
        sys.exit()

    # Calculate gnew
    gnewyx_array = golderyx_array[yidx_for_norm_array, xidx_for_norm_array] + rhs_array

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.plot(y_grid, golderyx_array[:,5])
    # ax1.plot(y_grid, gnewyx_array[:,5])
    # ax1.plot(y_grid, golderyx_array[:,9])
    # ax1.plot(y_grid, gnewyx_array[:,9])
    # print("x_grid[5] = ", x_grid[5] )
    # ax1.set_ylabel("g")
    # ax1.set_xlabel("y")
    # plt.show()

    return gnewyx_array

def take_many_nisl_steps(golderyx_array, vchiold_y_init, vchiold_x_init, nstep=100,
                        time_varying_velocity=False, velocity_omega=1):
    """Take a series of steps using NISL, collecting some diagnostic quantities
    for debugging/error-checking. We can't calculate vchi each step easily (would
    need to implement field equations) so keep velocity constant in time."""

    gzero = golderyx_array
    sum_g = []
    sum_gsquared = []
    sum_g.append(np.sum(gzero))
    sum_gsquared.append(np.sum(gzero*gzero))
    sum_over_x_gsquared_array = np.zeros((ny, nstep+1))
    sum_over_y_gsquared_array = np.zeros((nx, nstep+1))
    sum_over_x_gsquared_array[:,0] = np.sum(gzero*gzero, axis=1)
    sum_over_y_gsquared_array[:,0] = np.sum(gzero*gzero, axis=0)

    for istep in range(0, nstep):
        if time_varying_velocity:
            time = istep*dt
            vchiold_x_array = vchiold_x_init * np.cos(velocity_omega*time)
            vchiold_y_array = vchiold_y_init * np.cos(velocity_omega*time)
        else:
            vchiold_x_array = vchiold_x_init
            vchiold_y_array = vchiold_y_init
        ## Some diagnostic info
        if istep%50 == 0:
            print("istep = ", istep, "v_x[0,0] = ", vchiold_x_array[0,0])
        if istep == 0:
            print("Single-step iteration")
            ### Single-step scheme
            # Calcaulte dgolder/dx
            dgolder_dy,  dgolder_dx = get_dgdx_and_dgdy(golderyx_array)
            gnewyx_array = nisl_step_finding_departure_point(golderyx_array, dgolder_dy, dgolder_dx,
                                        vchiold_y_array, vchiold_x_array, n_timesteps=1)
            goldyx_array = gnewyx_array
        else:
            dgold_dy, dgold_dx = get_dgdx_and_dgdy(goldyx_array)
            gnewyx_array = nisl_step_finding_departure_point(golderyx_array, dgold_dy, dgold_dx,
                                                        vchiold_y_array, vchiold_x_array, n_timesteps=2)
            # Update arrays for g
            golderyx_array = goldyx_array
            goldyx_array = gnewyx_array
            # Also need to get dgold_dx
        sum_g.append(np.sum(gnewyx_array))
        sum_gsquared.append(np.sum(gnewyx_array*gnewyx_array))
        sum_over_x_gsquared_array[:,istep+1] = np.sum(gnewyx_array*gnewyx_array, axis=1)
        sum_over_y_gsquared_array[:,istep+1] = np.sum(gnewyx_array*gnewyx_array, axis=0)

    def diagnostic_plot_many_steps():
        """ """
        # fig = plt.figure()
        # ax1 = fig.add_subplot(111)
        # ax1.plot(y_grid, gzero[:,5])
        # ax1.plot(y_grid, gnewyx_array[:,5])
        # ax1.plot(y_grid, gzero[:,9])
        # ax1.plot(y_grid, gnewyx_array[:,9])
        # ax1.set_ylabel("g")
        # ax1.set_xlabel("y")
        # plt.show()

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(range(0, nstep+1), sum_g)
        ax2.plot(range(0, nstep+1), sum_gsquared)
        ax2.set_xlabel("time idx")
        ax1.set_ylabel(r"$\sum_{x,y} g(x,y)$")
        ax2.set_ylabel(r"$\sum_{x,y} g^2$")
        for ax in [ax1, ax2]:
            ax.grid(True)

        plt.show()

        return

    # diagnostic_plot_many_steps()

    return sum_g, sum_gsquared, gnewyx_array, sum_over_x_gsquared_array, sum_over_y_gsquared_array

def isl_step(golderyx_array,
                vchiold_y_array, vchiold_x_array, n_timesteps=2, interp_kind="linear"):
    """Take an ISL step, finding trajectory foot & h at the foot
    using linear interpolation."""

    # Upsample - double the number of points in the vchi, gold, golder grids.
    # golderyx_array_upsampled = create_upsampled_grid(golderyx_array)

    def get_interpolated_quantity_linear(x, y, f_array):
        """Given a quantity f evaluated on the x,y grid (f_array), find the value
        of the function at arbitrary location (x,y). Begin just by using
        linear interpolation (likely inaccurate but hopefully good stability
        properties) """

        ## x and y may be beyond-grid. We need to rectify this so that the
        ## idxs are within-grid, but the distances are unchanged (i.e. can
        ## be beyond-grid.)
        ## Find the 4 nearest gridpoints
        xidx_left = int(np.floor(x/dx)) ; xidx_right = int(np.ceil(x/dx))
        xidx_left_on_grid = xidx_left%nx ; xidx_right_on_grid = xidx_right%nx
        yidx_bottom = int(np.floor(y/dy)) ; yidx_top = int(np.ceil(y/dy))
        yidx_bottom_on_grid = yidx_bottom%ny ; yidx_top_on_grid = yidx_top%ny
        # print("xidx_left, xidx_right, yidx_bottom, yidx_top = ", xidx_left, xidx_right, yidx_bottom, yidx_top)
        # print("within-grid xidx_left, xidx_right, yidx_bottom, yidx_top = ",
        #         xidx_left_on_grid, xidx_right_on_grid, yidx_bottom_on_grid, yidx_top_on_grid)
        ### There's 2 situations which can occur, in both x and y:
        # (1) The 2 idxs (left/right or top/bottom) are different, because the value
        #     of x or y is non-integer
        # (2) The 2 idxs are the same, because the value of x or y is integer
        ## Find the value of the quantities at the 4 nearby corners  . . .
        ## Need to think about BCs . . . .
        f_lb = f_array[yidx_bottom_on_grid, xidx_left_on_grid]
        f_lt = f_array[yidx_top_on_grid, xidx_left_on_grid]
        f_rb = f_array[yidx_bottom_on_grid, xidx_right_on_grid]
        f_rt = f_array[yidx_top_on_grid, xidx_right_on_grid]

        ## And the weightings of these points . . .
        # if the 2 idxs are the same, set the weighting of one of them to zero
        if yidx_bottom == yidx_top:
            w_bottom = 0
        else:
            w_bottom = 1-abs(yidx_bottom - y/dy)
        w_top = 1-abs(yidx_top - y/dy)

        if xidx_left == xidx_right:
            w_left = 0
        else:
            w_left = 1-abs(xidx_left - x/dx)
        w_right = 1-abs(xidx_right - x/dx)

        #print("w_bottom, w_top, w_left, w_right = ", w_bottom, w_top, w_left, w_right)
        w_lb = w_left * w_bottom
        w_lt = w_left * w_top
        w_rb = w_right * w_bottom
        w_rt = w_right * w_top
        # print("f = ", f_lb, f_lt, f_rb, f_rt)
        # print("weights = ", w_lb, w_lt, w_rb, w_rt)
        f_interpolated = f_lb*w_lb + f_lt*w_lt + f_rb*w_rb + f_rt*w_rt
        #print("f_interpolated = ", f_interpolated)
        return f_interpolated

    def get_interpolated_quantity_linear_scipy(x, y, f_array):
        """Given a quantity f evaluated on the x,y grid (f_array), find the value
        of the function at arbitrary location (x,y) using scipy's interp2d, with
        linear interpolation """

        linear_fit = interp2d(x_grid, y_grid, f_array, bounds_error=True, kind="linear")
        x_within_grid = x%(xmax-dx) ; y_within_grid = y%(ymax-dy)
        f_interpolated = linear_fit(x_within_grid, y_within_grid)

        return f_interpolated

    def get_interpolated_quantity_cubic_scipy(x, y, f_array):
        """Given a quantity f evaluated on the x,y grid (f_array), find the value
        of the function at arbitrary location (x,y) using scipy's interp2d, with
        linear interpolation """

        linear_fit = interp2d(x_grid, y_grid, f_array, bounds_error=True, kind="cubic")
        x_within_grid = x%(xmax-dx) ; y_within_grid = y%(ymax-dy)
        f_interpolated = linear_fit(x_within_grid, y_within_grid)

        return f_interpolated

    def get_interpolated_quantity(x, y, f_array):
        """ """
        if interp_kind == "linear":
            return get_interpolated_quantity_linear(x, y, f_array)
        elif interp_kind == "linear_scipy":
            return get_interpolated_quantity_linear_scipy(x, y, f_array)
        elif interp_kind == "cubic_scipy":
            return get_interpolated_quantity_cubic_scipy(x, y, f_array)

    def get_approx_departure_point_and_g_departure(yidx, xidx):
        """Find the approximate departure point which takes us from an (in general,
        off-grid) point at t^{n-1} to the gridpoint labelled by (xidx, yidx) at
        t^{n+1}.
        """

        ### Get approx departure point, using the iterative guess approach

        # location of the gridpoint (=arrival point)
        x_nonperiodic = x_grid[xidx]; y_nonperiodic=y_grid[yidx]
        # Normalised values; these are integer at cell boundaries.
        x = x_grid[xidx] ; y = y_grid[yidx]
        # print("yidx, xidx, y, x = ", yidx, xidx, y, x )
        max_iterations = 2 ; counter = 0
        alpha_x = 0 ; alpha_y = 0
        while counter < max_iterations:
            counter +=1

            ## The algorithm to iteratively guess the displacement of the
            ## trajectory foot is:
            # alpha_x = dt*u_x[x-alpha_x, y-alpha_y, t^n]
            # alpha_y = dt*u_y[x-alpha_x, y-alpha_y, t^n]
            ## But also need to account for boundaries

            ## Find u_x, u_y at departure points
            u_x = get_interpolated_quantity(x-alpha_x, y-alpha_y, vchiold_x_array)
            u_y = get_interpolated_quantity(x-alpha_x, y-alpha_y, vchiold_y_array)

            ## Update alpha_x, alpha_y
            alpha_x = dt*u_x ; alpha_y = dt*u_y
            # print("u_x, u_y = ", u_x, u_y)
            # Our boundaries are nonperiodic.
            # if u_x > 0 :
            #     xboundary = ((np.floor(xnorm) + 0.5) * dx/2)
            # else:
            #     xboundary = ((np.ceil(xnorm)  + 0.5) * dx/2)
            # if u_y > 0 :
            #     yboundary = ((np.floor(ynorm) + 0.5) * dy/2)
            # else:
            #     yboundary = ((np.ceil(ynorm)  + 0.5) * dy/2)

        ## Once we've got alpha_x, alpha_y, need to find golder at (x-2*alpha_x, y-2*alpha_y)
        if n_timesteps == 1:
            x_departure = x - alpha_x ; y_departure = y - alpha_y
        elif n_timesteps ==2:
            x_departure = x - 2*alpha_x ; y_departure = y - 2*alpha_y
        else:
            print("Error! n_timesteps = ", n_timesteps)
            sys.exit()
        g_at_departure = get_interpolated_quantity(x_departure, y_departure, golderyx_array)
        # print("g_at_departure = ", g_at_departure)
        # sys.exit()

        return x_departure, y_departure, g_at_departure

    # print("vchiold_y_array = ", vchiold_y_array)
    # print("vchiold_x_array = ", vchiold_x_array)
    gnewyx_array = np.zeros((ny, nx))
    approx_departure_points_x = np.zeros((ny, nx))
    approx_departure_points_y = np.zeros((ny, nx))

    # Find the approximate departure point, and the value of the function g at
    # the departure points.
    for yidx in range(0, ny):
        for xidx in range(0, nx):
            y, x, g_at_departure = get_approx_departure_point_and_g_departure(yidx, xidx)
            approx_departure_points_x[yidx, xidx] = x
            approx_departure_points_y[yidx, xidx] = y
            gnewyx_array[yidx, xidx] = g_at_departure

    ############################################################################
    ### DIAGNSOTIC PLOT
    ############################################################################

    def make_diagnostic_plot_departure_points():
        """ """
        marker_size = 20.
        arrow_head_width = 0.1
        fig = plt.figure(figsize=[12, 8])
        ax1 = fig.add_subplot(111)
        ax1.scatter(x_grid_2d_upsampled.flatten(), y_grid_2d_upsampled.flatten(), marker="x", s=marker_size, label="upsampled grid")
        ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
        ax1.scatter(approx_departure_points_x.flatten(), approx_departure_points_y.flatten(), s=marker_size,label="departure point")
        # ax1.scatter(xhistory, yhistory, s=marker_size, label="trajectory")
        # for horizontal_line in horizontal_lines:
        #     ax1.plot(horizontal_line[0], horizontal_line[1], ls="--", c="gray")
        # for vertical_line in vertical_lines:
        #     ax1.plot(vertical_line[0], vertical_line[1], ls="--", c="gray")
        # for arrow in arrows:
        #     ax1.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)
        # for hist_idx in range(0, len(xhistory)-1):
        #     ax1.plot([xhistory[hist_idx], xhistory[hist_idx+1]], [yhistory[hist_idx], yhistory[hist_idx+1]])
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        ax1.grid(True)
        #ax1.set_xlim([-0.96, 25])
        #ax1.set_xlim([-1, 23])
        ax1.legend(loc="upper right")
        plt.show()

        return

    #make_diagnostic_plot_departure_points()
    ############################################################################

    # print("max(golderyx_array) = ", np.max(abs(golderyx_array)))
    # print("max(gnewyx_array) = ", np.max(abs(gnewyx_array)))
    return gnewyx_array

def take_many_isl_steps(golderyx_array, vchiold_y_init, vchiold_x_init, nstep=100,
                        time_varying_velocity=False, velocity_omega=1, interp_kind="linear"):
    """Take a series of steps using interpolating semi-lagrange. As a quick-and-dirty
    approach, we're going to use linear interpolation for both the departure
    point calculations and the calculation of h at the departure point."""

    gzero = golderyx_array
    sum_g = []
    sum_gsquared = []
    sum_g.append(np.sum(gzero))
    sum_gsquared.append(np.sum(gzero*gzero))
    sum_over_x_gsquared_array = np.zeros((ny, nstep+1))
    sum_over_y_gsquared_array = np.zeros((nx, nstep+1))
    sum_over_x_gsquared_array[:,0] = np.sum(gzero*gzero, axis=1)
    sum_over_y_gsquared_array[:,0] = np.sum(gzero*gzero, axis=0)

    for istep in range(0, nstep):
        if time_varying_velocity:
            time = istep*dt
            vchiold_x_array = vchiold_x_init * np.cos(velocity_omega*time)
            vchiold_y_array = vchiold_y_init * np.cos(velocity_omega*time)
        else:
            vchiold_x_array = vchiold_x_init
            vchiold_y_array = vchiold_y_init
        ## Some diagnostic info
        if istep%50 == 0:
            print("istep = ", istep, "v_x[0,0] = ", vchiold_x_array[0,0])
        if istep == 0:
            print("Single-step iteration")
            ### Single-step scheme
            gnewyx_array = isl_step(golderyx_array,
                                    vchiold_y_array, vchiold_x_array, n_timesteps=1, interp_kind=interp_kind)
            goldyx_array = gnewyx_array
        else:
            gnewyx_array = isl_step(golderyx_array,
                                    vchiold_y_array, vchiold_x_array, n_timesteps=2, interp_kind=interp_kind)
            # Update arrays for g
            golderyx_array = goldyx_array
            goldyx_array = gnewyx_array
            # Also need to get dgold_dx
        sum_g.append(np.sum(gnewyx_array))
        sum_gsquared.append(np.sum(gnewyx_array*gnewyx_array))
        sum_over_x_gsquared_array[:,istep+1] = np.sum(gnewyx_array*gnewyx_array, axis=1)
        sum_over_y_gsquared_array[:,istep+1] = np.sum(gnewyx_array*gnewyx_array, axis=0)

    return sum_g, sum_gsquared, gnewyx_array, sum_over_x_gsquared_array, sum_over_y_gsquared_array

def rk2_step(g, dg_dy, dg_dx, vy_array, vx_array):
    """RK2 method consists of:
    fa = fold - h/2 * v*df/dx |(told)  # This is an estimate of f at thalf=told+h/2
    fnew = fold - h * v*df/dx |(told+h/2)   # Taking a full step, but with df/dt evaluated at t=told+h/2
    """
    ghalf = g - (dt/2) * (dg_dy*vy_array + dg_dx*vx_array)
    dg_dy_half, dg_dx_half = get_dgdx_and_dgdy(ghalf)
    gnew = g - dt * (dg_dy_half*vy_array + dg_dx_half*vx_array)
    return gnew

def take_many_rk2_steps(goldyx_array, vchiold_y_init, vchiold_x_init, nstep=100,
                        time_varying_velocity=False, velocity_omega=1):
    """Take a series of steps using RK2, collecting some diagnostic quantities
    for debugging/error-checking. We can't calculate vchi each step easily (would
    need to implement field equations) so keep velocity constant in time."""

    gzero = goldyx_array
    sum_g = []
    sum_gsquared = []
    sum_g.append(np.sum(gzero))
    sum_gsquared.append(np.sum(gzero*gzero))
    sum_over_x_gsquared_array = np.zeros((ny, nstep+1))
    sum_over_y_gsquared_array = np.zeros((nx, nstep+1))
    sum_over_x_gsquared_array[:,0] = np.sum(gzero*gzero, axis=1)
    sum_over_y_gsquared_array[:,0] = np.sum(gzero*gzero, axis=0)
    Courant_num_array = (vchiold_x_init*dt/dx + vchiold_y_init*dt/dy)
    print("max Courant no = ", np.max(abs(Courant_num_array)))
    def diagnostic_plot_velocity_field():
        """ A plot showing cell velocities."""
        x_grid_boundaries = x_grid + dx/2
        y_grid_boundaries = y_grid + dy/2

        # To make the horizontal lines: make 1 horizontal line per
        # y_grid_upsampled_boundaries, starting at x=0 and ending at max(x_grid_upsampled_boundaries)
        horizontal_lines = []
        for diag_yval in y_grid_boundaries:
            horizontal_line_xvals = [0, max(x_grid_boundaries)]
            horizontal_line_yvals = [diag_yval, diag_yval]
            horizontal_lines.append([horizontal_line_xvals, horizontal_line_yvals])

        # To make the vertical lines: make 1 vertical line per
        # x_grid_upsampled_boundaries, starting at y=0 and ending at max(y_grid_upsampled_boundaries)
        vertical_lines = []
        for diag_xval in x_grid_boundaries:
            vertical_line_xvals = [diag_xval, diag_xval]
            vertical_line_yvals = [0, max(y_grid_boundaries)]
            vertical_lines.append([vertical_line_xvals, vertical_line_yvals])

        # Normalise velocities such that the largest velocity occupies a cell diagonal.
        # That is, want max(unorm_x) <= dx/2 and max(unorm_y) <= dy/2, but
        # want to scale both by the same amount so it's an accurate representation
        # of the velocity direction.
        x_scaling_fac = dx / np.max(abs(vchiold_x_array))
        y_scaling_fac = dy / np.max(abs(vchiold_y_array))
        scaling_fac = min(x_scaling_fac, y_scaling_fac)
        unorm_x = vchiold_x_init * scaling_fac
        unorm_y = vchiold_y_init * scaling_fac

        # Want to represent these velocities with an arrow, which is centered on
        # the gridpoint. So the starting point of the arrow should be [x - unorm_x/2, y - unorm_y/2]
        arrows = [] # Each item in arrows is a list describing a single arrow; [x, y, delta_x, delta_y]
        for diag_xidx in range(0, nx):
            for diag_yidx in range(0, ny):
                arrow_x = x_grid_2d[diag_yidx, diag_xidx] - unorm_x[diag_yidx, diag_xidx]/2
                arrow_y = y_grid_2d[diag_yidx, diag_xidx] - unorm_y[diag_yidx, diag_xidx]/2
                arrow_dx = unorm_x[diag_yidx, diag_xidx]
                arrow_dy = unorm_y[diag_yidx, diag_xidx]
                arrows.append([arrow_x, arrow_y, arrow_dx, arrow_dy])


        marker_size = 20.
        arrow_head_width = 0.1
        fig = plt.figure(figsize=[12, 8])
        ax1 = fig.add_subplot(111)
        ax1.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
        for horizontal_line in horizontal_lines:
            ax1.plot(horizontal_line[0], horizontal_line[1], ls="--", c="gray")
        for vertical_line in vertical_lines:
            ax1.plot(vertical_line[0], vertical_line[1], ls="--", c="gray")
        for arrow in arrows:
            ax1.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        #ax1.grid(True)
        #ax1.set_xlim([-1, 23])
        ax1.legend(loc="upper right")
        plt.show()

        return

    # diagnostic_plot_velocity_field()

    for istep in range(0, nstep):
        if time_varying_velocity:
            time = istep*dt
            vchiold_x_array = vchiold_x_init * np.cos(velocity_omega*time)
            vchiold_y_array = vchiold_y_init * np.cos(velocity_omega*time)
        else:
            vchiold_x_array = vchiold_x_init
            vchiold_y_array = vchiold_y_init

        ## Some diagnostic info
        if istep%50 == 0:
            print("istep = ", istep, "v_x[0,0] = ", vchiold_x_array[0,0])

        # Calcaulte dgolder/dx
        dgold_dy,  dgold_dx = get_dgdx_and_dgdy(goldyx_array)
        gnewyx_array = rk2_step(goldyx_array, dgold_dy, dgold_dx,
                                vchiold_y_array, vchiold_x_array)
        goldyx_array = gnewyx_array

        sum_g.append(np.sum(gnewyx_array))
        sum_gsquared.append(np.sum(gnewyx_array*gnewyx_array))
        sum_over_x_gsquared_array[:,istep+1] = np.sum(gnewyx_array*gnewyx_array, axis=1)
        sum_over_y_gsquared_array[:,istep+1] = np.sum(gnewyx_array*gnewyx_array, axis=0)
    def diagnostic_plot_many_steps():
        """ """
        # fig = plt.figure()
        # ax1 = fig.add_subplot(111)
        # ax1.plot(y_grid, gzero[:,5])
        # ax1.plot(y_grid, gnewyx_array[:,5])
        # ax1.plot(y_grid, gzero[:,9])
        # ax1.plot(y_grid, gnewyx_array[:,9])
        # ax1.set_ylabel("g")
        # ax1.set_xlabel("y")
        # plt.show()

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(range(0, nstep+1), sum_g)
        ax2.plot(range(0, nstep+1), sum_gsquared)
        ax2.set_xlabel("time idx")
        ax1.set_ylabel(r"$\sum_{x,y} g(x,y)$")
        ax2.set_ylabel(r"$\sum_{x,y} g^2$")
        for ax in [ax1, ax2]:
            ax.grid(True)

        plt.show()

        return

    # diagnostic_plot_many_steps()

    return sum_g, sum_gsquared, gnewyx_array, sum_over_x_gsquared_array, sum_over_y_gsquared_array

def take_many_steps(goldyx_array, vchiold_y_array, vchiold_x_array, nstep=100, schemes=["RK2"],
                        time_varying_velocity=False, velocity_omega=1, figtitle="Blank" ,
                        save_prefix="images/default", interp_kinds=["linear"]):
    """Take a series of steps using RK2, collecting some diagnostic quantities
    for debugging/error-checking. We can't calculate vchi each step easily (would
    need to implement field equations) so keep velocity constant in time."""

    time_array = np.arange(0, nstep+1) * dt

    Courant_num_array = (vchiold_x_array*dt/dx + vchiold_y_array*dt/dy)
    print("max Courant no = ", np.max(abs(Courant_num_array)))
    ## Calculate max. dv/dx,y * dt
    dvx_dx = (vchiold_x_array - np.roll(vchiold_x_array, -1, axis=1) )/ dx
    dvy_dx = (vchiold_y_array - np.roll(vchiold_y_array, -1, axis=1) )/ dx
    dvx_dy = (vchiold_x_array - np.roll(vchiold_x_array, -1, axis=0) )/ dy
    dvy_dy = (vchiold_y_array - np.roll(vchiold_y_array, -1, axis=0) )/ dy

    max_dvx_dx_dt = np.max(abs(dvx_dx)*dt)
    max_dvy_dx_dt = np.max(abs(dvy_dx)*dt)
    max_dvx_dy_dt = np.max(abs(dvx_dy)*dt)
    max_dvy_dy_dt = np.max(abs(dvy_dy)*dt)
    print("max dvx/dx * dt = ", max_dvx_dx_dt)
    print("max dvy/dx * dt = ", max_dvy_dx_dt)
    print("max dvx/dy * dt = ", max_dvx_dy_dt)
    print("max dvy/dy * dt = ", max_dvy_dy_dt)

    def diagnostic_plot_velocity_field(ax):
        """ A plot showing cell velocities."""
        x_grid_boundaries = x_grid + dx/2
        y_grid_boundaries = y_grid + dy/2

        # To make the horizontal lines: make 1 horizontal line per
        # y_grid_upsampled_boundaries, starting at x=0 and ending at max(x_grid_upsampled_boundaries)
        horizontal_lines = []
        for diag_yval in y_grid_boundaries:
            horizontal_line_xvals = [0, max(x_grid_boundaries)]
            horizontal_line_yvals = [diag_yval, diag_yval]
            horizontal_lines.append([horizontal_line_xvals, horizontal_line_yvals])

        # To make the vertical lines: make 1 vertical line per
        # x_grid_upsampled_boundaries, starting at y=0 and ending at max(y_grid_upsampled_boundaries)
        vertical_lines = []
        for diag_xval in x_grid_boundaries:
            vertical_line_xvals = [diag_xval, diag_xval]
            vertical_line_yvals = [0, max(y_grid_boundaries)]
            vertical_lines.append([vertical_line_xvals, vertical_line_yvals])

        # Normalise velocities such that the largest velocity occupies a cell diagonal.
        # That is, want max(unorm_x) <= dx/2 and max(unorm_y) <= dy/2, but
        # want to scale both by the same amount so it's an accurate representation
        # of the velocity direction.
        x_scaling_fac = dx / np.max(abs(vchiold_x_array))
        y_scaling_fac = dy / np.max(abs(vchiold_y_array))
        scaling_fac = min(x_scaling_fac, y_scaling_fac)
        unorm_x = vchiold_x_array * scaling_fac
        unorm_y = vchiold_y_array * scaling_fac

        # Want to represent these velocities with an arrow, which is centered on
        # the gridpoint. So the starting point of the arrow should be [x - unorm_x/2, y - unorm_y/2]
        arrows = [] # Each item in arrows is a list describing a single arrow; [x, y, delta_x, delta_y]
        for diag_xidx in range(0, nx):
            for diag_yidx in range(0, ny):
                arrow_x = x_grid_2d[diag_yidx, diag_xidx] - unorm_x[diag_yidx, diag_xidx]/2
                arrow_y = y_grid_2d[diag_yidx, diag_xidx] - unorm_y[diag_yidx, diag_xidx]/2
                arrow_dx = unorm_x[diag_yidx, diag_xidx]
                arrow_dy = unorm_y[diag_yidx, diag_xidx]
                arrows.append([arrow_x, arrow_y, arrow_dx, arrow_dy])


        marker_size = 20.
        arrow_head_width = 0.1
        ax.scatter(x_grid_2d.flatten(), y_grid_2d.flatten(), s=60, label="grid")
        for horizontal_line in horizontal_lines:
            ax.plot(horizontal_line[0], horizontal_line[1], ls="--", c="gray")
        for vertical_line in vertical_lines:
            ax.plot(vertical_line[0], vertical_line[1], ls="--", c="gray")
        for arrow in arrows:
            ax.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color="blue", length_includes_head = True, head_width=arrow_head_width)

        return


    sum_g_list = []
    sum_gsquared_list = []
    sum_gxsquared_list = []
    sum_gysquared_list = []
    g_array_final_list = []
    for scheme_idx, scheme in enumerate(schemes):
        if scheme == "RK2":
            sum_g, sum_gsquared, g_array_final, sum_over_x_gsquared_array, sum_over_y_gsquared_array = take_many_rk2_steps(goldyx_array, vchiold_y_array, vchiold_x_array, nstep=nstep,
                                        time_varying_velocity=time_varying_velocity, velocity_omega=velocity_omega,)
            sum_g_list.append(sum_g)
            sum_gsquared_list.append(sum_gsquared)
            sum_gxsquared_list.append(sum_over_x_gsquared_array)
            sum_gysquared_list.append(sum_over_y_gsquared_array)
            g_array_final_list.append(g_array_final)


        elif scheme == "NISL":
            sum_g, sum_gsquared, g_array_final, sum_over_x_gsquared_array, sum_over_y_gsquared_array = take_many_nisl_steps(goldyx_array, vchiold_y_array, vchiold_x_array, nstep=nstep,
                                        time_varying_velocity=time_varying_velocity, velocity_omega=velocity_omega,)
            sum_g_list.append(sum_g)
            sum_gsquared_list.append(sum_gsquared)
            g_array_final_list.append(g_array_final)
            sum_gxsquared_list.append(sum_over_x_gsquared_array)
            sum_gysquared_list.append(sum_over_y_gsquared_array)
        elif scheme == "ISL":
            sum_g, sum_gsquared, g_array_final, sum_over_x_gsquared_array, sum_over_y_gsquared_array = take_many_isl_steps(
                                        goldyx_array, vchiold_y_array, vchiold_x_array, nstep=nstep,
                                        time_varying_velocity=time_varying_velocity, velocity_omega=velocity_omega,
                                        interp_kind=interp_kinds[scheme_idx])
            sum_g_list.append(sum_g)
            sum_gsquared_list.append(sum_gsquared)
            g_array_final_list.append(g_array_final)
            sum_gxsquared_list.append(sum_over_x_gsquared_array)
            sum_gysquared_list.append(sum_over_y_gsquared_array)
        else:
            print("scheme not recognised, quitting!")
            print("scheme = ", scheme)
            sys.exit()

    fig = plt.figure(figsize=[18,6])
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(224)
    diagnostic_plot_velocity_field(ax1)

    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    ax1.legend(loc="upper right")
    plt.suptitle(figtitle)
    plt.savefig(save_prefix + "_a.png")
    plt.close()

    my_linewidth = 3
    linestyles = ["-", "--", ":"]
    cols = [default_cols[0], default_cols[1], default_cols[2]]

    fig = plt.figure(figsize=[18,10])
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224)

    for scheme_idx, scheme in enumerate(schemes):
        ax1.plot(time_array, sum_g_list[scheme_idx], label=scheme, ls=linestyles[scheme_idx], lw=my_linewidth)
        ax2.plot(time_array, sum_gsquared_list[scheme_idx], ls=linestyles[scheme_idx], lw=my_linewidth)

        sum_over_y_gsquared_array = sum_gysquared_list[scheme_idx]
        sum_over_x_gsquared_array = sum_gxsquared_list[scheme_idx]

        for xidx in range(0, nx):
            ax3.plot(time_array, sum_over_y_gsquared_array[xidx,:], c=cols[scheme_idx], ls=linestyles[scheme_idx], alpha=0.8)
        for xidx in range(0, ny):
            ax4.plot(time_array, sum_over_x_gsquared_array[xidx,:], c=cols[scheme_idx], ls=linestyles[scheme_idx], alpha=0.8)

    ax1.legend(loc="best")
    ax2.set_xlabel("time idx")
    ax4.set_xlabel("time idx")
    ax1.set_ylabel(r"$\sum_{x,y} g(x,y)$")
    ax2.set_ylabel(r"$\sum_{x,y} g^2$")
    ax3.set_ylabel(r"$\sum_{y} g^2$")
    ax4.set_ylabel(r"$\sum_{x} g^2$")

    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
    for ax in [ax3, ax4]:
        ax.set_yscale("log")
    plt.suptitle(figtitle)
    plt.savefig(save_prefix + "_b.png")
    plt.close()
    fig = plt.figure(figsize=[18,6])
    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)

    axes = [ax2, ax3, ax4]

    orig = ax1.imshow(goldyx_array)
    fig.colorbar(orig, ax=ax1)

    for scheme_idx, scheme in enumerate(schemes):
        imshow_object = axes[scheme_idx].imshow(g_array_final_list[scheme_idx])
        fig.colorbar(imshow_object, ax=axes[scheme_idx])
    ax1.set_ylabel("x idx")
    ax2.set_ylabel("x idx")
    ax1.set_xlabel("y idx")
    ax2.set_xlabel("y idx")
    plt.suptitle(figtitle)
    plt.savefig(save_prefix + "_c.png")
    plt.close()

    return

def benchmark_rk2_vs_nisl_different_velocity_fields_time_constant_spatilly_varying():
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
    take_many_steps(golderyx_array, constant_vchi_y*0.1, constant_vchi_x*0.1, nstep=1000, schemes=["RK2", "NISL"],
                    figtitle="v=constant, small", save_prefix="images/constant_v_small")
    # take_many_steps(golderyx_array, constant_vchi_y, constant_vchi_x, nstep=1000, schemes=["RK2", "NISL"],
    #                 figtitle="v=constant", save_prefix="images/constant_v")
    # take_many_steps(golderyx_array, vy_sinusoidal, vx_sinusoidal, nstep=1000, schemes=["RK2", "NISL"],
    #                 figtitle="v=sinusoidal (1)", save_prefix="images/v_sinusoidal_1")
    # take_many_steps(golderyx_array, vy_sinusoidal_both, vx_sinusoidal_both, nstep=1000, schemes=["RK2", "NISL"],
    #                 figtitle="v=sinusoidal (2)", save_prefix="images/v_sinusoidal_2")
    # take_many_steps(golderyx_array, vy_sinusoidal_other_way, vx_sinusoidal_other_way, nstep=1000, schemes=["RK2", "NISL"],
    #                 figtitle="v=sinusoidal (3)", save_prefix="images/v_sinusoidal_3")
    # take_many_steps(golderyx_array, np.zeros((ny, nx)), vchiold_x_array, nstep=1000, schemes=["RK2", "NISL"],
    #                 figtitle="vy=0, vx=v_stella", save_prefix="images/vx_chi_vy_0")
    # take_many_steps(golderyx_array, vchiold_y_array, vchiold_x_array, nstep=500, schemes=["RK2", "NISL"],
    #                 figtitle="v=v_stella", save_prefix="images/vstella")
    print("Finished")

    return

def benchmark_nisl_different_vmag():
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
    take_many_steps(golderyx_array, constant_vchi_y, constant_vchi_x, nstep=1000, schemes=["RK2", "NISL"])
    # take_many_steps(golderyx_array, constant_vchi_y*0.2, constant_vchi_x*0.2, nstep=1000, schemes=["RK2", "NISL"])
    # take_many_steps(golderyx_array, constant_vchi_y*0.05, constant_vchi_x*0.05, nstep=1000, schemes=["RK2", "NISL"])
    # take_many_steps(golderyx_array, constant_vchi_y*2, constant_vchi_x*2, nstep=1000, schemes=["NISL"])
    # take_many_steps(golderyx_array, constant_vchi_y*10, constant_vchi_x*10, nstep=1000, schemes=["NISL"])
    # take_many_steps(golderyx_array, constant_vchi_y*2, constant_vchi_x*2, nstep=1000, schemes=["NISL"], time_varying_velocity=True, velocity_omega=1)
    # take_many_steps(golderyx_array, constant_vchi_y*2, constant_vchi_x*2, nstep=1000, schemes=["NISL"], time_varying_velocity=True, velocity_omega=10)
    # take_many_steps(golderyx_array, vchiold_y_array, vchiold_x_array, nstep=500, schemes=["RK2", "NISL"], time_varying_velocity=True, velocity_omega=1)
    # take_many_steps(golderyx_array, vchiold_y_array*3, vchiold_x_array*3, nstep=500, schemes=["RK2", "NISL"], time_varying_velocity=True, velocity_omega=1)
    # take_many_steps(golderyx_array, vchiold_y_array*3, vchiold_x_array*3, nstep=500, schemes=["RK2", "NISL"], time_varying_velocity=True, velocity_omega=10)
    # take_many_steps(golderyx_array, vchiold_y_array, vchiold_x_array, nstep=500, schemes=["NISL"], time_varying_velocity=False)
    # take_many_steps(golderyx_array, vchiold_y_array*10, vchiold_x_array*10, nstep=500, schemes=["RK2", "NISL"], time_varying_velocity=False)
    # take_many_steps(golderyx_array, vchiold_y_array*40, vchiold_x_array*40, nstep=500, schemes=["NISL"], time_varying_velocity=False)
    # take_many_steps(golderyx_array, vchiold_y_array*40, vchiold_x_array*40, nstep=500, schemes=["NISL"], time_varying_velocity=True, velocity_omega=0.01)
    # take_many_steps(golderyx_array, vchiold_y_array*40, vchiold_x_array*40, nstep=500, schemes=["NISL"], time_varying_velocity=True, velocity_omega=1)
    # take_many_steps(golderyx_array, vchiold_y_array*40, vchiold_x_array*40, nstep=500, schemes=["NISL"], time_varying_velocity=True, velocity_omega=10)
    print("Finished")

    return

def benchmark_isl():
    """ """

    [golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
        vchiold_y_array, vchiold_x_array] = get_arrays_from_nonlinear_data()

    ## A knob which can increase vchi to exceed the CFL condition.
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
    take_many_steps(golderyx_array, constant_vchi_y, constant_vchi_x, nstep=100, schemes=["RK2", "NISL", "ISL"],
                    interp_kinds=["linear", "linear_scipy", "cubic_scipy"])
    print("Finished")

    return

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

def test_nisl_step_ritchie():
    """ """

    [golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
        vchiold_y_array, vchiold_x_array] = get_arrays_from_nonlinear_data()

    gnewyx_array_1 = nisl_step_ritchie(golder_array, golderyx_array, dgold_dy_array, dgold_dx_array,
                    vchiold_y_array, vchiold_x_array)
    gnewyx_array_2 = nisl_step_finding_departure_point(golderyx_array, dgold_dy_array, dgold_dx_array,
                    vchiold_y_array, vchiold_x_array, n_timesteps=1)
    print("gnewyx_array_2 - gnewyx_array_1 = ", gnewyx_array_2 - gnewyx_array_1)

if __name__ == "__main__":

    # benchmark_rk2_vs_nisl_different_velocity_fields_time_constant_spatilly_varying()
    #benchmark_isl()
    # benchmark_rk2_vs_nisl_different_velocity_fields_spatially_constant_temporally_varying()
    # benchmark_nisl_different_vmag()
    test_nisl_step_ritchie()
    # nisl_step_finding_departure_point(golderyx_array, dgold_dy_array, dgold_dx_array,
    #              vchiold_y_array, vchiold_x_array)
