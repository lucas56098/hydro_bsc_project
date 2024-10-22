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
import csv
import random


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

### 2D FV PLOTTING OPTIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2D(s, polygons, Q, cmap = 'viridis', vmin = 0, vmax = 1, edgecolor = 'face', cbar_label = 'Q_value', title = '', xlabel = "", ylabel = "", save = True, save_name = 'image2D', figsize = (12, 10), logscale = False, logmin = 1e-18, xlim = (0, 1), ylim = (0, 1), rasterize = True, plot_seeds = False, seed_size = 10, xticks = [0.2, 0.4, 0.6, 0.8], yticks = [0.2, 0.4, 0.6, 0.8], edgewith = 1, show_cbar = True, seedcolor = "tab:blue", dpi = 1000):

    # optionally prepare Q for logscale
    if logscale:
        Q = np.where(np.isinf(Q), -np.inf, np.log10(np.maximum(Q, np.zeros(len(Q)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(figsize = figsize)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    collection = PolyCollection(polygons, array=Q, cmap=cmap, norm=norm, edgecolor=edgecolor, linewidth = edgewith) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')
    
    collection.set_rasterized(rasterize)

    # add collection to axis and set axis options
    ax.add_collection(collection)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([''] * len(ax.get_xticks()))
    ax.set_yticklabels([''] * len(ax.get_yticks()))
    
    if plot_seeds:
        ax.scatter(s[:, 0], s[:, 1], color = seedcolor, s=seed_size)


    # set colorbar
    if show_cbar:
        cbar = fig.colorbar(collection, ax=ax)
        cbar.set_label(cbar_label)

    # optional set title
    if title != '':
        plt.title(title)

    # optional save plot
    if save:
        plt.savefig('figures/' + save_name + '.pdf', dpi = dpi, bbox_inches='tight')

    #plt.show()


# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2Dx3(polygons0, polygons1, polygons2, Q0, Q1, Q2, cmap = 'viridis', vmin = 0, vmax = 1, edgecolors = ['face', 'face', 'face'], cbar_label = 'Q_value', subtitles = ["", "", ""], xlabel = "", ylabel = "", save = True, save_name = 'image2D', figsize = (12, 10), logscale = False, logmin = 1e-18, xlim = (0, 1), ylim = (0, 1), rasterize = True, cbar_orientation = "horizontal", cbar_aspect = 40, cbar_padding = 0.02, edge_widths = [1, 1, 1], dpi = 1000):

    # optionally prepare Q for logscale
    if logscale:
        Q0 = np.where(np.isinf(Q0), -np.inf, np.log10(np.maximum(Q0, np.zeros(len(Q0)) + logmin)))
        Q1 = np.where(np.isinf(Q1), -np.inf, np.log10(np.maximum(Q1, np.zeros(len(Q1)) + logmin)))
        Q2 = np.where(np.isinf(Q2), -np.inf, np.log10(np.maximum(Q2, np.zeros(len(Q2)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(1, 3, figsize = figsize)
    plt.subplots_adjust(wspace=0.07)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    collection = PolyCollection(polygons0, array=Q0, cmap=cmap, norm=norm, edgecolor=edgecolors[0], linewidth = edge_widths[0]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')
    
    collection.set_rasterized(rasterize)

    # add collection to axis and set axis options
    ax[0].add_collection(collection)
    ax[0].set_xlim(xlim[0], xlim[1])
    ax[0].set_ylim(ylim[0], ylim[1])
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[0].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[0].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[0].set_xticklabels([''] * len(ax[0].get_xticks()))
    ax[0].set_yticklabels([''] * len(ax[0].get_yticks()))
    ax[0].set_title(subtitles[0])

    collection2 = PolyCollection(polygons1, array=Q1, cmap=cmap, norm=norm, edgecolor=edgecolors[1], linewidth = edge_widths[1]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    collection2.set_rasterized(rasterize)

    ax[1].add_collection(collection2)
    ax[1].set_xlim(xlim[0], xlim[1])
    ax[1].set_ylim(ylim[0], ylim[1])
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(ylabel)
    ax[1].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[1].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[1].set_xticklabels([''] * len(ax[1].get_xticks()))
    ax[1].set_yticklabels([''] * len(ax[1].get_yticks()))
    ax[1].set_title(subtitles[1])

    collection3 = PolyCollection(polygons2, array=Q2, cmap=cmap, norm=norm, edgecolor=edgecolors[2], linewidth = edge_widths[2]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    collection3.set_rasterized(rasterize)

    ax[2].add_collection(collection3)
    ax[2].set_xlim(xlim[0], xlim[1])
    ax[2].set_ylim(ylim[0], ylim[1])
    ax[2].set_xlabel(xlabel)
    ax[2].set_ylabel(ylabel)
    ax[2].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[2].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[2].set_xticklabels([''] * len(ax[2].get_xticks()))
    ax[2].set_yticklabels([''] * len(ax[2].get_yticks()))
    ax[2].set_title(subtitles[2])

    # set colorbar
    cbar = fig.colorbar(collection, ax=ax, orientation=cbar_orientation, aspect = cbar_aspect, pad = cbar_padding)
    cbar.set_label(cbar_label)
    
    # optional save plot
    if save:
        plt.savefig('figures/' + save_name + '.pdf', dpi = dpi, bbox_inches='tight')

    plt.show()


# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2Dx4(polygons0, polygons1, polygons2, polygons3, Q0, Q1, Q2, Q3, cmap = 'viridis', vmin = 0, vmax = 1, edgecolors = ['face', 'face', 'face', 'face'], cbar_label = 'Q_value', subtitles = ["", "", "", ""], xlabel = "", ylabel = "", save = True, save_name = 'image2D', figsize = (12, 10), logscale = False, logmin = 1e-18, xlim = (0, 1), ylim = (0, 1), rasterize = True, cbar_orientation = "horizontal", cbar_aspect = 40, cbar_padding = 0.02, edge_widths = [1, 1, 1, 1], xticks = [0.2, 0.4, 0.6, 0.8], yticks = [0.2, 0.4, 0.6, 0.8], dpi = 1000):

    # optionally prepare Q for logscale
    if logscale:
        Q = np.where(np.isinf(Q), -np.inf, np.log10(np.maximum(Q, np.zeros(len(Q)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(1, 4, figsize = figsize)
    plt.subplots_adjust(wspace=0.07)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    collection = PolyCollection(polygons0, array=Q0, cmap=cmap, norm=norm, edgecolor=edgecolors[0], linewidth = edge_widths[0]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')
    
    collection.set_rasterized(rasterize)

    #fig.suptitle("high resolution")

    # add collection to axis and set axis options
    ax[0].add_collection(collection)
    ax[0].set_xlim(xlim[0], xlim[1])
    ax[0].set_ylim(ylim[0], ylim[1])
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[0].set_xticks(xticks)
    ax[0].set_yticks(yticks)
    ax[0].set_xticklabels([''] * len(ax[0].get_xticks()))
    ax[0].set_yticklabels([''] * len(ax[0].get_yticks()))
    ax[0].set_title(subtitles[0])

    collection2 = PolyCollection(polygons1, array=Q1, cmap=cmap, norm=norm, edgecolor=edgecolors[1], linewidth = edge_widths[1]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    collection2.set_rasterized(rasterize)

    ax[1].add_collection(collection2)
    ax[1].set_xlim(xlim[0], xlim[1])
    ax[1].set_ylim(ylim[0], ylim[1])
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(ylabel)
    ax[1].set_xticks(xticks)
    ax[1].set_yticks(yticks)
    ax[1].set_xticklabels([''] * len(ax[1].get_xticks()))
    ax[1].set_yticklabels([''] * len(ax[1].get_yticks()))
    ax[1].set_title(subtitles[1])

    collection3 = PolyCollection(polygons2, array=Q2, cmap=cmap, norm=norm, edgecolor=edgecolors[2], linewidth = edge_widths[2]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    collection3.set_rasterized(rasterize)

    ax[2].add_collection(collection3)
    ax[2].set_xlim(xlim[0], xlim[1])
    ax[2].set_ylim(ylim[0], ylim[1])
    ax[2].set_xlabel(xlabel)
    ax[2].set_ylabel(ylabel)
    ax[2].set_xticks(xticks)
    ax[2].set_yticks(yticks)
    ax[2].set_xticklabels([''] * len(ax[2].get_xticks()))
    ax[2].set_yticklabels([''] * len(ax[2].get_yticks()))
    ax[2].set_title(subtitles[2])

    collection4 = PolyCollection(polygons3, array=Q3, cmap=cmap, norm=norm, edgecolor=edgecolors[3], linewidth = edge_widths[3]) # edgecolor = 'none' / 'face'
    print('finished PolyCollection')

    collection4.set_rasterized(rasterize)

    ax[3].add_collection(collection4)
    ax[3].set_xlim(xlim[0], xlim[1])
    ax[3].set_ylim(ylim[0], ylim[1])
    ax[3].set_xlabel(xlabel)
    ax[3].set_ylabel(ylabel)
    ax[3].set_xticks(xticks)
    ax[3].set_yticks(yticks)
    ax[3].set_xticklabels([''] * len(ax[3].get_xticks()))
    ax[3].set_yticklabels([''] * len(ax[3].get_yticks()))
    ax[3].set_title(subtitles[3])

    # set colorbar
    cbar = fig.colorbar(collection, ax=ax, orientation=cbar_orientation, aspect = cbar_aspect, pad = cbar_padding)
    cbar.set_label(cbar_label)
    
    # optional save plot
    if save:
        plt.savefig('figures/' + save_name + '.pdf', dpi = dpi, bbox_inches='tight')

    plt.show()


# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2Dx2(polygons0, polygons1, Q0, Q1, cmap='viridis', vmin=0, vmax=1, edgecolors = ['face', 'face'], cbar_label='Q_value', 
              subtitle1='', subtitle2='', xlabel="", ylabel="", save=True, save_name='image2D', figsize=(12, 10), logscale=False, 
              logmin=1e-18, xlim=(0, 1), ylim=(0, 1), rasterize = True, dpi = 1000, edgewidth = 0.5):

    # optionally prepare Q for logscale
    if logscale:
        Q0 = np.where(np.isinf(Q0), -np.inf, np.log10(np.maximum(Q0, np.zeros(len(Q0)) + logmin)))
        Q1 = np.where(np.isinf(Q1), -np.inf, np.log10(np.maximum(Q1, np.zeros(len(Q1)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(1, 2, figsize=figsize)
    plt.subplots_adjust(wspace=0.05)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    collection0 = PolyCollection(polygons0, array=Q0, cmap=cmap, norm=norm, edgecolor=edgecolors[0], linewidth = edgewidth)
    print('finished PolyCollection')

    collection0.set_rasterized(rasterize)

    ax[0].add_collection(collection0)
    ax[0].set_xlim(xlim[0], xlim[1])
    ax[0].set_ylim(ylim[0], ylim[1])
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[0].set_title(subtitle1)
    ax[0].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[0].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[0].set_xticklabels([''] * len(ax[0].get_xticks()))
    ax[0].set_yticklabels([''] * len(ax[0].get_yticks()))

    collection1 = PolyCollection(polygons1, array=Q1, cmap=cmap, norm=norm, edgecolor=edgecolors[1])
    print('finished PolyCollection')

    collection1.set_rasterized(rasterize)

    ax[1].add_collection(collection1)
    ax[1].set_xlim(xlim[0], xlim[1])
    ax[1].set_ylim(ylim[0], ylim[1])
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(ylabel)
    ax[1].set_title(subtitle2)
    ax[1].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[1].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[1].set_xticklabels([''] * len(ax[1].get_xticks()))
    ax[1].set_yticklabels([''] * len(ax[1].get_yticks()))

    cbar = fig.colorbar(collection0, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label(cbar_label)

    # optional save plot
    if save:
        plt.savefig('figures/' + save_name + '.pdf', dpi = dpi, bbox_inches='tight')

    plt.show()


# function to do a 2D plot of the mesh with the colormap according to Q
def plot_2Dx3b(s0, s1, s2, q0, q1, q2, polygons0, polygons1, polygons2, Q0, Q1, Q2, cmap='viridis', vmin=0, vmax=1, 
              edgecolor1='face', edgecolor2='face', edgecolor3='face', 
              cbar_label='Q_value', subtitle1='', subtitle2='', subtitle3='', 
              xlabel="", ylabel="", save=True, save_name='image2D', figsize=(18, 10), 
              logscale=False, logmin=1e-18, xlim=(0, 1), ylim=(0, 1), rasterize = True, edge_widths = [1, 1, 1], quiver_scale = 20, dpi = 1000):

    # optionally prepare Q for logscale
    if logscale:
        Q0 = np.where(np.isinf(Q0), -np.inf, np.log10(np.maximum(Q0, np.zeros(len(Q0)) + logmin)))
        Q1 = np.where(np.isinf(Q1), -np.inf, np.log10(np.maximum(Q1, np.zeros(len(Q1)) + logmin)))
        Q2 = np.where(np.isinf(Q2), -np.inf, np.log10(np.maximum(Q2, np.zeros(len(Q2)) + logmin)))

    # define plot, norm, poly collection
    fig, ax = plt.subplots(1, 3, figsize=figsize)
    plt.subplots_adjust(wspace=0.08)
    norm1 = mcolors.Normalize(vmin=vmin[0], vmax=vmax[0])
    norm2 = mcolors.Normalize(vmin=vmin[1], vmax=vmax[1])
    norm3 = mcolors.Normalize(vmin=vmin[2], vmax=vmax[2])

    collection0 = PolyCollection(polygons0, array=Q0, cmap=cmap, norm=norm1, edgecolor=edgecolor1, linewidth = edge_widths[0])
    print('finished PolyCollection')

    collection0.set_rasterized(rasterize)
    
    ax[0].add_collection(collection0)
    ax[0].set_xlim(xlim[0], xlim[1])
    ax[0].set_ylim(ylim[0], ylim[1])
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[0].set_title(subtitle1)
    ax[0].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[0].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[0].set_xticklabels([''] * len(ax[0].get_xticks()))
    ax[0].set_yticklabels([''] * len(ax[0].get_yticks()))
    
    step = len(s0) // 2000
    x_quiver0 = s0[::step, 0]
    y_quiver0 = s0[::step, 1]
    u_quiver0 = q0[::step, 2]
    v_quiver0 = q0[::step, 3]

    quiver0 = ax[0].quiver(x_quiver0, y_quiver0, u_quiver0, v_quiver0, angles='xy', scale_units='xy', scale=quiver_scale, color='gray')
    quiver0.set_rasterized(rasterize)
    cbar0 = fig.colorbar(collection0, ax=ax[0], orientation='vertical', pad=0.02)
    #cbar0.set_label(cbar_label)

    collection1 = PolyCollection(polygons1, array=Q1, cmap=cmap, norm=norm2, edgecolor=edgecolor2, linewidth = edge_widths[1])
    print('finished PolyCollection')

    collection1.set_rasterized(rasterize)
    
    ax[1].add_collection(collection1)
    ax[1].set_xlim(xlim[0], xlim[1])
    ax[1].set_ylim(ylim[0], ylim[1])
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(ylabel)
    ax[1].set_title(subtitle2)
    ax[1].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[1].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[1].set_xticklabels([''] * len(ax[1].get_xticks()))
    ax[1].set_yticklabels([''] * len(ax[1].get_yticks()))

    step = len(s1) // 2000
    x_quiver1 = s1[::step, 0]
    y_quiver1 = s1[::step, 1]
    u_quiver1 = q1[::step, 2]
    v_quiver1 = q1[::step, 3]

    quiver1 = ax[1].quiver(x_quiver1, y_quiver1, u_quiver1, v_quiver1, angles='xy', scale_units='xy', scale=quiver_scale, color='gray')
    quiver1.set_rasterized(rasterize)
    cbar1 = fig.colorbar(collection1, ax=ax[1], orientation='vertical', pad=0.02)
    #cbar1.set_label(cbar_label)

    collection2 = PolyCollection(polygons2, array=Q2, cmap=cmap, norm=norm3, edgecolor=edgecolor3, linewidth = edge_widths[2])
    print('finished PolyCollection')

    collection2.set_rasterized(rasterize)

    ax[2].add_collection(collection2)
    ax[2].set_xlim(xlim[0], xlim[1])
    ax[2].set_ylim(ylim[0], ylim[1])
    ax[2].set_xlabel(xlabel)
    ax[2].set_ylabel(ylabel)
    ax[2].set_title(subtitle3)
    ax[2].set_xticks([0.2, 0.4, 0.6, 0.8])
    ax[2].set_yticks([0.2, 0.4, 0.6, 0.8])
    ax[2].set_xticklabels([''] * len(ax[2].get_xticks()))
    ax[2].set_yticklabels([''] * len(ax[2].get_yticks()))

    step = len(s2) // 2000
    x_quiver2 = s2[::step, 0]
    y_quiver2 = s2[::step, 1]
    u_quiver2 = q2[::step, 2]
    v_quiver2 = q2[::step, 3]

    quiver2 = ax[2].quiver(x_quiver2, y_quiver2, u_quiver2, v_quiver2, angles='xy', scale_units='xy', scale=quiver_scale, color='gray')
    quiver2.set_rasterized(rasterize)
    cbar2 = fig.colorbar(collection2, ax=ax[2], orientation='vertical', pad=0.02)
    #cbar2 = fig.colorbar(collection2, ax=ax, orientation='vertical', pad=0.02)
    cbar2.set_label(cbar_label)

    # optional save plot
    if save:
        plt.savefig('figures/' + save_name + '.pdf', dpi = dpi, bbox_inches='tight')
        #plt.savefig('figures/' + save_name + '.png', dpi=500)

    plt.show()


# function to do a animation of the mesh evoulution in 2D
def animation2D(file_name, frames, fps=30, animation_name='animation2D', cbar_label='Q-value', cmap='viridis', edgecolor='face', vmin=0, vmax=1, title='', xlabel = "", ylabel = "", lim=(0, 1), figsize=(12, 10), logscale=False, logmin=1e-18, quantity_index = 1, plot_seeds = False, do_quiver = False):
    
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
    ax.set_xticklabels([''] * len(ax.get_xticks()))
    ax.set_yticklabels([''] * len(ax.get_yticks()))

    # optional set title
    if title != '':
        ax.set_title(title)

    # set colorbar
    cbar = fig.colorbar(collection, ax=ax)
    cbar.set_label(cbar_label)

    # add empty scatter plot to eventually later on show seeds
    scatter_plot = ax.scatter([], [], color = 'tab:red', s = 10)

    quiver_plot = None

    # animation update function
    def update(frame):
        nonlocal quiver_plot

        # open file for this frame
        file_path = "files/" + file_name + str(frame) + ".csv"
        seeds, polygons, Q = process_file(file_path)
        s = seeds
        p = polygons
        q = Q
        
        # get Q and polygons for collection
        Q = Q[:, quantity_index]
        if logscale:
            Q = np.where(np.isinf(Q), -np.inf, np.log10(np.maximum(Q, np.zeros(len(Q)) + logmin)))
        collection.set_paths(polygons)
        collection.set_array(Q)

        # optional plot seeds
        if plot_seeds:
            scatter_plot.set_offsets(seeds)

        if do_quiver:
            if quiver_plot:
                quiver_plot.remove()
            step = len(p) // 3000
            x_quiver = s[::step, 0]
            y_quiver = s[::step, 1]
            u_quiver = q[::step, 2]
            v_quiver = q[::step, 3]

            quiver_plot = ax.quiver(x_quiver, y_quiver, u_quiver, v_quiver, angles='xy', scale_units='xy', scale=10, color='gray')
                
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


### 1D FV PLOTTING OPTIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# function to do 1D-plot of vmesh and cartesian Q values
def plot_1D(filenames, labels, sort_option, title = '', x_label = '', y_label = 'Q-Value', xlim = (0,1), ylim = (0, 1), save_name = '', bin_size = 0, logscale = False, quantity_index = 1, analytic_solution = "",
            x0 = 0.5, gamma = 5.0/3.0, left_state = (1,1,0), right_state = (0.1, 0.125, 0), velocity = 0.5, a = 0, b = 0.1, 
                hl = 2, hr= 1, g = 1, rasterized = False):
    
    # go through all filenames
    for i in range(len(filenames)):

        # open file
        seed, polygons, Q = process_file("files/" + filenames[i] + ".csv", sort_option)

        # plot Q[:, quantity_index] sorted by x or y
        if sort_option[i] == 'y':
            plt.scatter(seed[:, 1], Q[:, quantity_index[i]], marker = '.', label= labels[i], rasterized = rasterized)
        elif sort_option[i] == 'x':
            #plt.plot(seed[:, 0], Q[:, quantity_index[i]], marker = '+', label= labels[i], linewidth = 1, rasterized = rasterized)
            plt.scatter(seed[:, 0], Q[:, quantity_index[i]], label= labels[i], s = 1, rasterized = rasterized)
        elif sort_option[i] == 'diagonal':
            # plot a projection onto diagonal here
            plt.scatter(((1/np.sqrt(2)) * (seed[:, 0] + seed[:, 1])), Q[:, quantity_index[i]], label= labels[i], s=0.8, rasterized = rasterized)
        elif sort_option[i] == 'diagonal_filter':
            # calculate distance of points to diagonal
            distance_to_diagonal = np.abs(seed[:, 0] - seed[:, 1]) / np.sqrt(2)
            
            # filter indices by distance limit
            filtered_indices = distance_to_diagonal < 0.05
            filtered_seed = seed[filtered_indices]
            filtered_Q = Q[filtered_indices]

            # plot filtered seeds projected onto diagonal
            proj_seed = (1/np.sqrt(2)) * (filtered_seed[:, 0] + filtered_seed[:, 1])
            plt.scatter(proj_seed, filtered_Q[:, quantity_index[i]], label=labels[i], s=0.8, rasterized = rasterized)

        # optional binned plot in x direction
        if bin_size != 0:
            avg_seed = [np.mean(seed[i:i+bin_size, 0]) for i in range(0, int(len(seed[:, 1])), bin_size)]
            avg_Q = [np.mean(Q[i:i+bin_size, quantity_index]) for i in range(0, int(len(Q[:, quantity_index])), bin_size)]
            plt.plot(avg_seed, avg_Q, label = labels[i] + 'avg', color = 'black')

        if analytic_solution == "shock_tube" and i == range(len(filenames))[-1]:
                
            ### For the analytical solution to the shock tube we use a 
            ### package written by Jerko Škifić. Feel free to check it out: https://github.com/ibackus/sod-shocktube
            import sodshock  # Copyright (c) 2015 Jerko Škifić

            positions, regions, values = sodshock.solve(left_state, right_state, (xlim[0], xlim[1], x0), Q[0, 0], gamma, 1000)

            val_names = ['time - cannot plot analytical solution for time', 'rho', 'u', 'v  - cannot plot analytical solution for y-velocity', 'e  - cannot plot analytical solution for energy density', 'p'] # only rho, u, p will work though (only ones that we can directly compare)

            plt.plot(values['x'], values[val_names[quantity_index[0]]], color = 'red', label = 'analytic', linewidth = 1, rasterized = rasterized)

        # optional add analytic solution for advecting step
        if analytic_solution == "adv_step" and i == range(len(filenames))[-1]:
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
            plt.plot(lsp, h_analytic, label = 'analytic', color = 'grey', linewidth = 1, alpha = 0.5, rasterized = rasterized)

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
        plt.savefig("figures/" + save_name + ".pdf")

    #plt.show()


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


### FV L1 ERROR OPTIONS ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
    #plt.show() 


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
    plt.xlabel("N in [a.u.]")
    plt.ylabel("L1 error in [a.u.]")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig("figures/L1_error_over_N.pdf")
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
    plt.savefig("figures/"+ save_name +".pdf")
    plt.show()    


### CREATE DIFFERENT STRUCTFILES ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# function to create structure file of an airfoil that can be used in a voronoi mesh simulation
def struct_airfoil(angle_of_attack = 0, n_pts = 300, rd_perturb_size = 0.0001, pt0 = (0.3, 0.5), save_name = 'struct', plot = False):

    # import airfoil
    from aeropy.airfoils.shapes import naca
    x, y = naca("2412", n = n_pts)

    # shrink to right size
    x = 0.5*np.array(x[::-1][:-1])
    y = 0.5*np.array(y[::-1][:-1])

    # add small random pertubation to avoid degeneracies in voronoi mesh generation
    x = np.array([x_i + 0.0001 * random.random() for x_i in x])
    y = np.array([y_i + 0.0001 * random.random() for y_i in y])

    # turn wing to correct angle of attack
    # rotation matrix
    angle = - np.deg2rad(angle_of_attack)
    rot = np.array([[np.cos(angle), -np.sin(angle)],
                     [np.sin(angle), np.cos(angle)]])
    x_r, y_r = np.dot(rot, np.vstack((x,y)))
    x_r += pt0[0]
    y_r += pt0[1]
    points = list(zip(x_r,y_r))

    # save points into csv file
    with open('files/' + save_name + '.csv', 'w', newline='') as file:
        f = csv.writer(file)
        for pt in points:
            f.writerow(pt)
    print("struct file for airfoil succesfully generated")

    # optional plot the airfoil shape
    if plot:
        plt.plot(np.append(x_r, x_r[0]), np.append(y_r, y_r[0]))
        plt.axis('equal')
        plt.title(rf"NACA 2412 Airfoil, attack angle $\alpha = {angle_of_attack}$")
        plt.show()


# function to create structure file of a circle that can be used in a voronoi mesh simulation
def struct_circle(radius = 0.05, n_pts = 200, p0 = (0.3, 0.5), rd_perturb_size = 0.0005, save_name = 'struct'):
    
    # generate points
    points = [(radius * np.cos((-2 * np.pi * i)/n_pts) + p0[0] + random.random()*rd_perturb_size,
               radius * np.sin((-2 * np.pi * i)/n_pts) + p0[1] + random.random()*rd_perturb_size) for i in range(n_pts)]
    
    # store in csv
    with open('files/' + save_name + '.csv', 'w', newline='') as file:
        f = csv.writer(file)
        for pt in points:
            f.writerow(pt)
    print("struct file for circle succesfully generated")


# function to create structure file for a given .shp file (e.g. Africa)
def struct_shp(filepath = 'shapefiles_coastlines/Africa.shp', shrink_factor = 110, p0 = (0.45, 0.55), sample_step = 40, crop = (0, 124), plot = False, save_name = 'struct'):
    
    # import geopandas and load coastline data
    import geopandas as gpd
    coastline = gpd.read_file(filepath)
    coords = coastline.get_coordinates()

    # scale coordinates accordingly
    x = coords["x"].to_numpy()/shrink_factor + p0[0]
    y = coords["y"].to_numpy()/shrink_factor + p0[1]

    x_r = []
    y_r = []

    # sample coordinates
    for i in range(1, len(x), sample_step):
        x_r.append(x[i])
        y_r.append(y[i])
    
    print('original size:', len(x), 'size now:', len(x_r))

    # crop the coast
    x_r = x_r[crop[0]: crop[1]]
    y_r = y_r[crop[0]: crop[1]]
    print('after cropping:', len(x_r))

    # make sure that everything is inside of boundaries of mesh generator (0,1)
    x_s = []#[0.995, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25]
    y_s = []#[0.005001, 0.005003, 0.0050002, 0.00500023, 0.005000014, 0.00500042, 0.00500045, 0.0050000123 ,0.005000119 ,0.00500001234 ,0.0050001001 ,0.005000102 ,0.0050000405 ,0.005000003]
    for i in range(len(x_r)):
        if not (x_r[i] <= 0.005 or x_r[i] >= 0.995 or y_r[i] <= 0.005 or y_r[i] >= 0.995):
            x_s.append(x_r[i])
            y_s.append(y_r[i])
    print('after excluding outside:', len(x_s))

    #x_s.append(0.9949)
    #y_s.append(0.37)
    #x_s.append(0.99495)
    #y_s.append(0.15)

    # optional plot the shape
    if plot:
        #plt.figure(figsize=(6.1,6))
        plt.scatter(x, y, s = 1, label = 'original data', color = "tab:blue")
        plt.plot(x_s, y_s, color = "tab:orange", label = 'cropped and sampled')
        plt.legend()
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title("loaded shapefile")
        #plt.show()

    # store points in csv file
    points = [(x_s[i], y_s[i]) for i in range(len(x_s))]
    with open('files/' + save_name + '.csv', 'w', newline='') as file:
        f = csv.writer(file)
        for pt in points:
            f.writerow(pt)
    print('struct file for shape succesfully generated')


### 1D DG PLOT OPTIONS ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# basis functions on 1D reference element, length of Q vector determines polynomial order
def monomial_basis_func(x, Q):
    y  = 0

    # Q = [time, Q_0, Q_1, Q_2] ... -> Q[1] = Q_0
    for i in range(1, len(Q)):
        y += Q[i] * x**(i-1)

    return y

# function to plot 1D Discontinous Galerkin functions for monomial basis functions
def DG_plot_1D_monomial(filenames, colors, linestyles, labels, title = "1D DG scalar upwind advection", xlabel = "x in [a.u.]", ylabel = "Q-Value in [a.u.]", save_name = "DG_plot_1D", element_resolution = 100, xlim = (0, 1), ylim = (-0.1, 1.1), figsize = None, legend_loc = 1):

    # set plot options
    plt.figure(figsize=figsize)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])

    # linspace on reference element
    lsp_reference = np.linspace(-1, 1, element_resolution)

    # plot all files
    for i in range(len(filenames)):
        
        s, p, q = process_file("files/" + filenames[i] + ".csv")

        for j in range(len(s)):
            dx = p[i][2][0]-p[i][1][0]
            lsp_local = np.linspace(j*dx, (j+1)*dx, element_resolution)
            y = [monomial_basis_func(x, q[j]) for x in lsp_reference]

            # plots first element plus label and then all the other elements without label in the same color
            if j == 0:
                plt.plot(lsp_local, y, color = colors[i], linestyle = linestyles[i], label = labels[i])
            else:
                plt.plot(lsp_local, y, color = colors[i], linestyle = linestyles[i])
        print("plotted " + filenames[i])

    # show legend and save plot
    plt.legend(loc = legend_loc)
    plt.savefig('figures/' + save_name + ".pdf")
    #plt.show()


### 1D DG L1 ERROR ----------------------
# get L1 error of a single file for given initial conditions 
def get_L1_DG_bar(filename, integration_resolution = 100, h1 = 1, h2 = 0, x0 = 0.25):

    # load file
    s, p, q = process_file("files/" + filename + ".csv")
    t = q[0][0]
    N_elements = np.shape(q)[0]

    # function to get L1 error of single element
    def get_element_error(index):

        # local and reference spaces
        local_lsp = np.linspace(p[index][1][0]+1e-16, p[index][2][0]-1e-16, integration_resolution)
        reference_lsp = np.linspace(-1, 1, integration_resolution)

        # for linspaces get simulated and analytical values
        y_theo = np.array([analytic_Q(x, t, h1, h2, x0) for x in local_lsp])
        y_sim = np.array([monomial_basis_func(x, q[index]) for x in reference_lsp])

        # sum up errors and divide by resolution to get element L1 error
        L_elt = (np.sum(np.abs(y_theo - y_sim)))/integration_resolution
        return L_elt

    L1_sum = 0

    # sum up element errors
    for i in range(N_elements):
        L1_sum += get_element_error(i)
    
    # divide by element number
    return t, L1_sum/N_elements


# function to plot L1 error over time
def DG_1D_L1_over_time(filenames, filerange, labels, integral_resolution = 100, h1 = 1, h2 = 0, x0 = 0.25):

    # loop through all filenames
    for i in range(len(filenames)):

        times = []
        L1s = []

        # loop through all times
        for j in filerange:
            t, L1 = get_L1_DG_bar(filenames[i] + str(j), integral_resolution, h1, h2, x0)
            times.append(t)
            L1s.append(L1)

        # plot L1 over time
        plt.plot(times, L1s, label = labels[i])
    
    # plot style
    plt.legend()
    plt.title("DG 1D L1 error over time")
    plt.ylabel("L1 error")
    plt.xlabel("time t in [a.u.]")
    plt.savefig("figures/L1_DG_1D_over_time.pdf")
    plt.show()


# plot L1 error over N
def DG_1D_L1_over_N(filenames, N_xs, dataname = "", integral_resolution = 100, h1 = 1, h2 = 0, x0 = 0.25):

    L1s = []

    # go through given filenames
    for i in range(len(filenames)):
        t, L1 = get_L1_DG_bar(filenames[i], integral_resolution, h1, h2, x0)
        L1s.append(L1)

    # fit 1/x function onto L1 over N data
    x = N_xs
    y = L1s
    lsp = np.linspace(min(np.log(x)), max(np.log(x)), 10)
    popt, pcov = curve_fit(fit_function, np.log(x), np.log(y))
    print(popt)

    # scatterplot and plot style
    plt.scatter(N_xs, L1s, label = dataname + f" order: {np.abs(popt[0]):.2f}", marker = 'o')
    plt.plot(np.exp(lsp), np.exp(fit_function(lsp, *popt)), color = 'grey')
    plt.title("DG 1D L1 error over resolution")
    plt.xlabel("N_x")
    plt.ylabel("L1 error")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig("figures/DG_L1_error_over_N.pdf")



### 2D DG PLOT OPTIONS ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# basis function on 2D cartesian
def legendre_basis_func(x, y, Q):
    length = len(Q)
    if length == 2:
        return Q[1]
    elif length == 4:
        return Q[1] + x*Q[2] + y*Q[3]
    elif length == 7:
        return Q[1] + x*Q[2] + y*Q[3] + x*y*Q[4] + (1/2)*(3*x**2 - 1) * Q[5] +  (1/2)*(3*y**2 - 1) * Q[6]
    else:
        print("Error: len(Q)-1 != 1, 3 or 6. That is not implemented!")
        return - np.inf


# calculates 2D image and returns it as a np.array
def DG_calc_2D_image(filename, element_resolution = 10):
    
    # load file
    s, p, q = process_file("files/" + filename + ".csv")

    # image resolution
    N_row = int(np.sqrt(len(s)))
    pxl = element_resolution * N_row

    image = np.zeros((pxl, pxl))

    # load image
    for i in range(0, pxl, element_resolution):
        for j in range(0, pxl, element_resolution):
            for k1 in range(element_resolution):
                for k2 in range(element_resolution):
                    Q = q[int((N_row*i + j)/element_resolution)]
                    x = 2*(k2/element_resolution) - 1
                    y = 2*(k1/element_resolution) - 1
                    image[i+k1, j+k2] = legendre_basis_func(x, y, Q)
    return image


# plots 2D DG cartesian
def DG_plot_2D_legendre(filename, cmap = "inferno", vmin = -0.2, vmax = 1.2, save_name = "DG_2D_image", element_resolution = 10, figsize = (5, 4), xticks = [0.2, 0.4, 0.6, 0.8], yticks = [0.2, 0.4, 0.6, 0.8]):

    # load image
    image = DG_calc_2D_image(filename, element_resolution)

    # plot settings
    fig, ax = plt.subplots(figsize = figsize)
    ax.set_xticks(np.array(xticks) * len(image))
    ax.set_yticks(np.array(yticks) * len(image))
    ax.set_xticklabels([''] * len(ax.get_xticks()))
    ax.set_yticklabels([''] * len(ax.get_yticks()))

    # plot image
    implot = plt.imshow(image, origin="lower", cmap = cmap, vmin = vmin, vmax = vmax,  interpolation="none")
    implot.set_rasterized(True)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("figures/" + save_name + ".pdf", dpi = 1000)
    plt.show()

                    
### 2D DG L1 ERROR -----------------------
# analytical solution to advected square in 2D
def step_func_2D(x, y, t, v = [0.5, 0.5], p0 = [0, 0], a = 0.3, b = 0.3):
    
    if (x - v[0]*t >= p0[0] and x - v[0]*t < p0[0] + a and y - v[1]*t >= p0[1] and y - v[1]*t < p0[1] + b):
        return 1
    else:
        return 0
    
# analytical solution to advected gaussian in 2D
def gaussian_2D(x, y, t, v = [0.5, 0.5], p0 = [0.5, 0.5], A = 1, sigma = 0.25):
    if x < v[0]*t or y < v[1]*t:
        return 0
    else:
        return A * np.exp(-((((x - v[0]*t) - p0[0])*((x - v[0]*t) - p0[0]) + ((y - v[1]*t) - p0[1])*((y - v[1]*t) - p0[1]))/(2*sigma*sigma)))


# calculates L1 error for 2D advected step function
def get_L1_DG_2D(filename, analytic_sol = 'step', integration_resolution = 100, v = [0.5, 0.5], p0 = [0, 0], a = 0.3, b = 0.3, A = 1, sigma = 0.25):
    
    s, p, q = process_file("files/" + filename + ".csv")
    t = q[0][0]
    N_elements = np.shape(q)[0]

    # function to get L1 error of single element
    def get_2D_element_error(index):

        # local and reference spaces
        local_lsp_x = np.linspace(p[index][1][0]+1e-16, p[index][2][0]-1e-16, integration_resolution)
        local_lsp_y = np.linspace(p[index][0][1]+1e-16, p[index][1][1]-1e-16, integration_resolution)
        reference_lsp = np.linspace(-1, 1, integration_resolution)

        # for linspaces get simulated and analytical values
        y_theo = []
        if analytic_sol == "step":
            y_theo = np.array([[step_func_2D(x, y, t, v, p0, a, b) for y in local_lsp_y] for x in local_lsp_x])
        elif analytic_sol == "gaussian":
            y_theo = np.array([[gaussian_2D(x, y, t, v, p0, A, sigma) for y in local_lsp_y] for x in local_lsp_x])
        else:
            print("choose step or gaussian as analytical solution")
        y_sim = np.array([[legendre_basis_func(x, y, q[index]) for y in reference_lsp] for x in reference_lsp])

        # sum up erros and divide by resolution to get element L1 error
        L_elt = (np.mean(np.abs(y_theo - y_sim)))
        return L_elt

    L1_sum = 0

    # sum up element errors
    for i in range(N_elements):
        L1_sum += get_2D_element_error(i)

    # divide by element number
    return t, L1_sum/N_elements


# function to plot L1 error over time 2D
def DG_2D_L1_over_time(filenames, filerange, labels, analytic_sol = 'step', integral_resolution = 100, v = [0.5, 0.5], p0 = [0, 0], a = 0.3, b = 0.3, A = 1, sigma = 0.25):

    # loop through all filenames
    for i in range(len(filenames)):
        
        times = []
        L1s = []

        # loop through all times
        for j in filerange:
            t, L1 = get_L1_DG_2D(filenames[i] + str(j), analytic_sol, integral_resolution, v, p0, a, b, A, sigma)

            times.append(t)
            L1s.append(L1)

        # plot L1 over time
        plt.plot(times, L1s, label = labels[i])

    # plot style
    plt.legend()
    plt.title("DG 2 L1 error over time")
    plt.ylabel("L1 error")
    plt.xlabel("time t in [a.u.]")
    plt.savefig("figures/L1_DG_2D_over_time.pdf")
    plt.show()


# plot 2D L1 error over N
def DG_2D_L1_over_N(filenames, N_xs, analytic_sol = 'step', dataname = "", integral_resolution = 100, v = [0.5, 0.5], p0 = [0, 0], a = 0.3, b = 0.3, A= 1, sigma = 0.25, color = None):

    ts = []
    L1s = []

    # go through given filenames
    for i in range(len(filenames)):
        t, L1 = get_L1_DG_2D(filenames[i], analytic_sol, integral_resolution, v, p0,a, b, A, sigma)
        ts.append(t)
        L1s.append(L1)

    print(t)

    # fit 1/x function onto L1 over N data
    x = N_xs
    y = L1s
    lsp = np.linspace(min(np.log(x)), max(np.log(x)), 10)
    popt, pcov = curve_fit(fit_function, np.log(x), np.log(y))
    print(popt)


    # scatterplot and plot style
    plt.scatter(N_xs, L1s, label = dataname + f" order: {np.abs(popt[0]):.2f}", marker = 'o', color = color)
    plt.plot(np.exp(lsp), np.exp(fit_function(lsp, *popt)), color = 'grey')
    plt.title("DG 2D L1 error over resolution")
    plt.xlabel("N_x")
    plt.ylabel("L1 error")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig("figures/DG_2D_L1_error_over_N.pdf")