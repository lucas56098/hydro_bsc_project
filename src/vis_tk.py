## Provides functions to easily plot the outputs of the simulation
# import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from scipy.optimize import curve_fit
from scipy.optimize import brentq


# processes mesh_file and returns seeds, polygons and the given quantities in an array
def process_file(file_name, sort_option = 'none', print_option = False):
    
    # read data into data frame
    df = pd.read_csv(file_name, delimiter='|', header=None, skiprows=1)
    Q = []
    seeds = []
    polygons = []
    data = []

    for index, row in df.iterrows():
        
        # get seedpoints
        seed_data = row[0]
        seed_points = np.array(seed_data.split(','), dtype=float)

        # get polygons for cells
        coords_data = row[1]
        poly = []
        for item in coords_data.split(';'):
            if item:
                points = item.split(",")
                points = np.array(points, dtype=float)
                poly.append(points)
        poly.append(poly[0])

        # get Q value or list of quantities
        q_data = row[2]
        if ',' in str(q_data):
            q_split = np.array(q_data.split(','), dtype=float)
        else:
            q_split = float(q_data)

        data.append((seed_points, poly, q_split))

    # sort data by x or y direction if needed
    if sort_option == 'x':
        data.sort(key=lambda x: x[0][0])
    elif sort_option == 'y':
        data.sort(key=lambda x: x[0][1])

    seeds, polygons, Q = zip(*data)

    if print_option:
        print('processed file')

    return np.array(seeds), list(polygons), np.array(Q)


# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2D(polygons, Q, cmap = 'viridis', vmin = 0, vmax = 1, edgecolor = 'face', cbar_label = 'Q_value', title = '', xlabel = "", ylabel = "", save = True, save_name = 'image2D', figsize = (12, 10), logscale = False, logmin = 1e-18, xlim = (0, 1), ylim = (0, 1)):

    # optionally prepare Q for logscale
    if logscale:
        Q = np.where(np.isinf(Q), -np.inf, np.log10(np.maximum(Q, np.zeros(len(Q)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(figsize = figsize)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    collection = PolyCollection(polygons, array=Q, cmap=cmap, norm=norm, edgecolor=edgecolor) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')
    
    # add collection to axis and set axis options
    ax.add_collection(collection)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # set colorbar
    cbar = fig.colorbar(collection, ax=ax)
    cbar.set_label(cbar_label)

    # optional set title
    if title != '':
        plt.title(title)

    # optional save plot
    if save:
        plt.savefig('figures/' + save_name + '.png')

    plt.show()

# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2Dx3(polygons0, polygons1, polygons2, Q0, Q1, Q2, cmap = 'viridis', vmin = 0, vmax = 1, edgecolor = 'face', cbar_label = 'Q_value', title = '', xlabel = "", ylabel = "", save = True, save_name = 'image2D', figsize = (12, 10), logscale = False, logmin = 1e-18, xlim = (0, 1), ylim = (0, 1)):

    # optionally prepare Q for logscale
    if logscale:
        Q = np.where(np.isinf(Q), -np.inf, np.log10(np.maximum(Q, np.zeros(len(Q)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(1, 3, figsize = figsize)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    collection = PolyCollection(polygons0, array=Q0, cmap=cmap, norm=norm, edgecolor=edgecolor) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')
    
    # add collection to axis and set axis options
    ax[0].add_collection(collection)
    ax[0].set_xlim(xlim[0], xlim[1])
    ax[0].set_ylim(ylim[0], ylim[1])
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[0].set_title("t = 0")

    collection2 = PolyCollection(polygons1, array=Q1, cmap=cmap, norm=norm, edgecolor=edgecolor) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    ax[1].add_collection(collection2)
    ax[1].set_xlim(xlim[0], xlim[1])
    ax[1].set_ylim(ylim[0], ylim[1])
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(ylabel)
    ax[1].set_title("t = 1")

    collection3 = PolyCollection(polygons2, array=Q2, cmap=cmap, norm=norm, edgecolor=edgecolor) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    ax[2].add_collection(collection3)
    ax[2].set_xlim(xlim[0], xlim[1])
    ax[2].set_ylim(ylim[0], ylim[1])
    ax[2].set_xlabel(xlabel)
    ax[2].set_ylabel(ylabel)
    ax[2].set_title("t = 1.6")

    for i in [0, 1, 2]:
        ax[i].set_xticks([])
        ax[i].set_yticks([])

    # set colorbar
    cbar = fig.colorbar(collection, ax=ax, orientation='horizontal', aspect = 40, pad = 0.02)
    cbar.set_label(cbar_label)
    
    # optional set title
    if title != '':
        plt.title(title)

    # optional save plot
    if save:
        #plt.savefig('figures/' + save_name + '.png')
        plt.savefig('figures/' + save_name + '.pdf')

    plt.show()


# function to do a animation of the mesh evoulution in 2D
def animation2D(file_name, frames, fps=30, animation_name='animation2D', cbar_label='Q-value', cmap='viridis', edgecolor='face', vmin=0, vmax=1, title='', xlabel = "", ylabel = "", lim=(0, 1), figsize=(12, 10), logscale=False, logmin=1e-18, quantity_index = 1, plot_seeds = False):
    
    # define plot, norm, and still empty collection
    fig, ax = plt.subplots(figsize=figsize)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    collection = PolyCollection([], array=[], cmap=cmap, norm=norm)
    collection.set_edgecolor(edgecolor)

    # add empty collection to axis and set axis options
    ax.add_collection(collection)
    ax.set_xlim(lim[0], lim[1])
    ax.set_ylim(lim[0], lim[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # optional set title
    if title != '':
        ax.set_title(title)

    # set colorbar
    cbar = fig.colorbar(collection, ax=ax)
    cbar.set_label(cbar_label)

    # add empty scatter plot to eventually later on show seeds
    scatter_plot = ax.scatter([], [], color = 'tab:red', s = 10)

    # animation update function
    def update(frame):

        # open file for this frame
        file_path = "files/" + file_name + str(frame) + ".csv"
        seeds, polygons, Q = process_file(file_path)
        
        # get Q and polygons for collection
        Q = Q[:, quantity_index]
        if logscale:
            Q = np.where(np.isinf(Q), -np.inf, np.log10(np.maximum(Q, np.zeros(len(Q)) + logmin)))
        collection.set_paths(polygons)
        collection.set_array(Q)

        # optional plot seeds
        if plot_seeds:
            scatter_plot.set_offsets(seeds)
        
        return collection,

    # update loop with progress bar
    with tqdm(total=len(frames), desc="Generating Animation") as pbar:
        def wrapped_update(frame):
            result = update(frame)
            pbar.update(1)
            return result

        ani = FuncAnimation(fig, wrapped_update, frames=frames, blit=True, repeat=False)
        ani.save('figures/' + animation_name + '.gif', fps=fps)
    
    plt.show()    


# function to do 1D-plot of vmesh and cartesian Q values
def plot_1D(filenames, labels, sort_option, title = '', x_label = '', y_label = 'Q-Value', xlim = (0,1), ylim = (0, 1), save_name = '', bin_size = 0, logscale = False, quantity_index = 1):
    
    # go through all filenames
    for i in range(len(filenames)):

        # open file
        seed, polygons, Q = process_file("files/" + filenames[i] + ".csv", sort_option)

        # plot Q[:, quantity_index] sorted by x or y
        if sort_option == 'y':
            plt.scatter(seed[:, 1], Q[:, quantity_index], marker = '.', label= labels[i])
        elif sort_option == 'x':
            plt.scatter(seed[:, 0], Q[:, quantity_index], marker = '.', label= labels[i])

        # optional binned plot in x direction
        if bin_size != 0:
            avg_seed = [np.mean(seed[i:i+bin_size, 0]) for i in range(0, int(len(seed[:, 1])), bin_size)]
            avg_Q = [np.mean(Q[i:i+bin_size, quantity_index]) for i in range(0, int(len(Q[:, quantity_index])), bin_size)]
            plt.plot(avg_seed, avg_Q, label = labels[i] + 'avg', color = 'black')

    # optional set title
    if title != '':
        plt.title(title)

    # set labels, limits, legend
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    plt.legend()

    # optional logscale
    if logscale:
        plt.yscale('log')

    # optional save to file
    if save_name != '':
        plt.savefig("figures/" + save_name + ".png")

    plt.show()


## analytic solution to 1D advection of step function
def analytic_Q(x, t, velocity, a, b):
    dx = t * velocity
    if (x >= dx + a and x<b+dx):
        return 1
    else:
        return 0


# function to make 1D animation of the mesh
def animation1D(filerange, filenames, labels, sort_option, quantity_index, fps = 30, title = '', x_label = '', y_label = 'Q-Value', 
                xlim = (0,1), ylim = (0, 1), save_name = 'animation1D', bin_size = 0, analytic_solution = "", velocity = 0.5, a = 0, b = 0.1, 
                hl = 2, hr= 1, g = 1, x0 = 0.5, gamma = 5.0/3.0, left_state = (1,1,0), right_state = (0.1, 0.125, 0)):                          # left/right state = (Pressure, Density, Velocity)
    
    # function to update the plot for given frame
    def update1D(frame):

        # clear plot
        plt.cla()

        # go through all filenames
        for i in range(len(filenames)):

            # load file
            seed, polygons, Q = process_file("files/" + filenames[i] + str(frame) + ".csv", sort_option[i])

            # sort and plot by x, y, diagonal or filterd diagonal
            if sort_option[i] == 'y':
                plt.scatter(seed[:, 1], Q[:, quantity_index[i]], label= labels[i], s=2)
            elif sort_option[i] == 'x':
                plt.scatter(seed[:, 0], Q[:, quantity_index[i]], label= labels[i], s=3)
            elif sort_option[i] == 'diagonal':
                # plot a projection onto diagonal here
                plt.scatter(((1/np.sqrt(2)) * (seed[:, 0] + seed[:, 1])), Q[:, quantity_index[i]], label= labels[i], s=0.8)
            elif sort_option[i] == 'diagonal_filter':
                # calculate distance of points to diagonal
                distance_to_diagonal = np.abs(seed[:, 0] - seed[:, 1]) / np.sqrt(2)
                
                # filter indices by distance limit
                filtered_indices = distance_to_diagonal < 0.05
                filtered_seed = seed[filtered_indices]
                filtered_Q = Q[filtered_indices]

                # plot filtered seeds projected onto diagonal
                proj_seed = (1/np.sqrt(2)) * (filtered_seed[:, 0] + filtered_seed[:, 1])
                plt.scatter(proj_seed, filtered_Q[:, quantity_index[0]], label=labels[i], s=0.8)

            # optional binned plot
            if bin_size != 0:
                avg_seed = [np.mean(seed[i:i+bin_size, 1]) for i in range(0, int(len(seed[:, 1])), bin_size)]
                avg_Q = [np.mean(Q[i:i+bin_size, quantity_index[i]]) for i in range(0, int(len(Q[:, quantity_index[i]])), bin_size)]
                plt.plot(avg_seed, avg_Q, label = labels[i] + '_avg', color = 'black')
        
            # optional add analytic solution for advecting step
            if analytic_solution == "adv_step":
                lsp = np.linspace(xlim[0], xlim[1], 1000)
                Q_analytic = [analytic_Q(x, Q[0, 0], velocity, a, b) for x in lsp]
                plt.plot(lsp, Q_analytic, label = labels[i] + '_analytic', color = 'red')

            # optional add analytic solution for shallow water dam break
            if analytic_solution == "swe_dam" and i == range(len(filenames))[-1]:

                # calculate wave speeds
                cl = np.sqrt(g*hl)
                cr = np.sqrt(g*hr)

                # root of this function gives cm
                def get_cm(cm):
                    return -8 * cr**2 * cm**2 * (cl - cm)**2 + (cm**2 - cr**2)**2 * (cm**2 + cr**2)
                
                # get cm and hm by finding the root
                cm = brentq(get_cm, min(cl, cr), max(cl, cr))
                hm = (cm**2)/g

                # given cm, hm now get analytical solution
                def get_h_at_t_and_x(x, t, x0, cl, cm, cr, hl, hm, hr, g):

                    # positions of the shocks and rarefications
                    xa = x0 - cl*t                                              # leftmost rarefication
                    xb = x0 + t*(2*cl - 3*cm)                                   # rightmos rarefication
                    xc = x0 + t * ((2 * cm**2 * (cl - cm))/(cm**2 - cr**2))     # right shock

                    # return h according to analytical solution
                    if x < xa:
                        return hl
                    elif xa < x and x < xb:
                        return (4)/(9*g) * (cl - ((x-x0)/(2*t)))**2
                    elif xb < x and x < xc:
                        return hm
                    elif xc < x:
                        return hr
                    return None

                # plot analytical solution given we now have a way to calculate it
                lsp = np.linspace(xlim[0], xlim[1], 1000)
                h_analytic = [get_h_at_t_and_x(x, Q[0, 0], x0, cl, cm, cr, hl, hm, hr, g) for x in lsp]
                plt.plot(lsp, h_analytic, label = 'analytic', color = 'red', linewidth = 0.5)

            if analytic_solution == "shock_tube" and i == range(len(filenames))[-1]:
                
                ### For the analytical solution to the shock tube we use a 
                ### package written by Jerko Škifić. Feel free to check it out: https://github.com/ibackus/sod-shocktube
                import sodshock  # Copyright (c) 2015 Jerko Škifić

                positions, regions, values = sodshock.solve(left_state, right_state, (xlim[0], xlim[1], x0), Q[0, 0], gamma, 1000)

                val_names = ['time - cannot plot analytical solution for time', 'rho', 'u', 'v  - cannot plot analytical solution for y-velocity', 'e  - cannot plot analytical solution for energy density', 'p'] # only rho, u, p will work though (only ones that we can directly compare)

                plt.plot(values['x'], values[val_names[quantity_index[0]]], color = 'red', label = 'analytic', linewidth = 1)

        # optional title
        if title != '':
            plt.title(title)

        # set labels, limits, legend
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])
        plt.legend(loc = 'upper right')


    # define plot
    fig, ax = plt.subplots()

    # update loop with progress bar
    with tqdm(total=len(filerange), desc='Generating Animation') as pbar:
        def update_and_progress(frame):
            update1D(frame)
            pbar.update(1)
    
        ani = FuncAnimation(fig, update_and_progress, frames=filerange, repeat=True)
        ani.save('figures/' + save_name + '.gif', fps=fps)
        plt.show()   


# function to do L1 error plots over time
def plot_L1_error_over_time(filenames = ["L1_error"], labels = ["L1_error"], title = 'L1 error over time', xlabel = "time in [a.u]", ylabel = 'L1 error', save_name = 'L1_over_time', axvlinepos=0):
    
    # go through all filenames
    for i in range(len(filenames)):

        # load time and L1 error from file using data frame
        df = pd.read_csv('files/' + filenames[i] + '.csv', decimal = ',', header=None)
        time = df[0].astype(float).values
        L1 = df[1].astype(float).values

        # plot L1 over time
        plt.plot(time, L1, label = labels[i])

    # title, labels, legend
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()

    # optional vline
    if axvlinepos != 0:
        plt.axvline(axvlinepos, color = "black")

    # save plot
    plt.savefig("figures/"+ save_name +".png")
    plt.show() 


# get relative to a "analytic" high res solution the L1 error
def get_L1_rel(analytic, test):

    sa, pa, qa = process_file("files/" + analytic + ".csv", "x")
    s, p, q = process_file("files/" + test + ".csv", "x")

    h_ana = qa[:, 1]
    h_test = q[:, 1]

    nr = int(np.log(len(h_ana)/len(h_test))/np.log(2))
    
    for i in range(nr):
        h_ana = np.mean(h_ana.reshape(-1, 2), axis = 1)

    return np.sum(np.abs(h_ana-h_test))/len(h_test)


# a/x^b fit function for L1 convergence over N
def fit_function(x, a, b):
    return a*x+b


# function to plot L1 error over N
def plot_L1_error_over_N(filenames, N_list, index = -1, dataname = "L1 error", relative_path = "", only_first_two = False):

    L1s = []

    if relative_path == "":
        # for all files store L1 at given index in L1s list
        for i in range(len(filenames)):
            df = pd.read_csv('files/' + filenames[i] + '.csv', decimal = ',', header=None)
            L1 = df[1].astype(float).values
            L1s.append(L1[index])
    else:
        for i in range(len(filenames)):
            L1s.append(get_L1_rel(relative_path, filenames[i]))

    N_list = np.array(N_list)
    L1s = np.array(L1s)

    # fit 1/x function onto L1 over N data
    x = N_list
    y = L1s
    lsp = np.linspace(min(np.log(x)), max(np.log(x)), 10)
    if only_first_two:
        lsp = np.linspace(min(np.log(x[0:4])), max(np.log(x[0:4])), 10)
        x = x[0:2]
        y = y[0:2]
    popt, pcov = curve_fit(fit_function, np.log(x), np.log(y))
    print(popt)

    # plot result
    plt.scatter(N_list, L1s, marker = 'o', label = dataname + f" order: {np.abs(popt[0]):.2f}")
    plt.plot(np.exp(lsp), np.exp(fit_function(lsp, *popt)), color = 'grey')
    plt.title("L1 error over resolution")
    plt.xlabel("N_x")
    plt.ylabel("L1 error")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig("figures/L1_error_over_N.png")
    #plt.show() 


# function to plot difference in total Q over time
def plot_Q_diff_over_time(filename = 'total_Q_diff', logscale = True, title = 'Change in total Q over time', xlabel = "time in [a.u]", ylabel = '|total_Q - total_Q_initial|', save_name = 'delta_Q_total', bar = 0):

    # load file into data frame and get time and Q_diff
    df = pd.read_csv('files/' + filename + '.csv', decimal=',', header=None)
    times = df[0].astype(float).values
    diff_Q = df[1].astype(float).values

    # plot change
    plt.plot(times, np.abs(diff_Q))

    # optional log scale
    if logscale:
        plt.yscale('log')

    # optional vlines
    if bar != 0:
        plt.axvline(bar, color = 'tab:orange')
    
    # title, label, save
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig("figures/"+ save_name +".png")
    plt.show()    
