#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module for meshing 2D images, used at CINEMAX summer schoool.
Created on Mon Jul  1 23:31:10 2019

@author: vand@dtu.dk 2019
"""

#%% LOADING PACKAGES

import numpy as np
import scipy.interpolate
import triangle  # https://rufat.be/triangle/installing.html
import matplotlib.path

#%% LOADING FUNCTIONS

def distribute_points(curve,ael):
    '''
    Distributes points of a curve acording to the average edge length. This is
    a helping function used by resample_contours function.
    
    Parameters
    ----------
    curve : A (N,2) numpy arrays containing 2D coordinates of one contour. 
        Closed contour has the same coordinate in the first and the last element.
    ael : Desired average edge length in pixels.

    Returns
    -------
    curve : A (n,2) numpy arrays containing 2D coordinates of the contour
        resampled according to ael.
    '''
    
    # point distances
    d = np.sqrt(np.sum(np.square(np.diff(curve, axis=0)), axis=1)) 
    # length parametrization of old points, ensuring that it starts with 0
    p = np.cumsum(np.hstack((0,d)))
    # number of points in the new contour, one larger than nr_segments
    num_pts = int(p[-1]/ael) + 1 
    f = scipy.interpolate.interp1d(p, curve, axis=0)
    curve = f(np.linspace(0, p[-1], num_pts)) 
    return curve


def resample_contours(contours,ael):
    '''
    Resamples contours according to the average edge length ael.
    
    Parameters
    ----------
    contours : A list of (N,2) numpy arrays. Each array contains 2D 
        coordinates of one contour. Closed contours have the same coordinate 
        in the first and the last element. This is the data structure returned
        by skimage.measure.find_contours.
    ael : Desired average edge length in pixels.

    Returns
    -------
    contours : Same dat structure as input contours but resampled according to 
        ael.
    '''
    contours_resampled = [] 
    for c in contours:
        c_new = distribute_points(c,ael)
        s = c_new.shape[0]
        i = np.array_equal(c_new[0],c_new[-1]) # a flag indicating connected curve
        if s>2+i: # discarding small components
            contours_resampled.append(c_new)
    return contours_resampled


def contours_to_segments(contours):
    '''
    Collects all line segments from the list of contours.
    
    Parameters
    ----------
    contours : A list of (N,2) numpy arrays with  coordinates of the contours.

    Returns
    -------
    vertices : A (V,2) numpy array with coordinates of the vertices.
    edges : A (E,2) numpy array with edges represented as a pair of vertex 
        indices. Returned edges are oriented.
    '''
    
    vertices = np.empty((0,2), dtype = float)
    edges = np.empty((0,2), dtype = int)
    for c in contours:
        s = c.shape[0] 
        i = np.array_equal(c[0],c[-1]) # a flag indicating connected curve
        a = np.arange(s-1)
        e = np.array([a,np.mod(a+1,s-i)]).T + vertices.shape[0]
        edges = np.vstack((edges,e))
        vertices = np.vstack((vertices,c[:s-i]))
    return vertices, edges


def remove_small_contours(contours, min_nr):
    '''
    Removes small contours from the contour list according to minimal number of 
    points.

    Parameters
    ----------
    contours : A list of (N,2) numpy arrays with coordinates of the contours.
    min_nr : Minimal number oof contour points. Smaller contours will be removed.

    Returns
    -------
    contours_new : A list of (n,2) numpy arrays with coordinates of the contours.

    '''
    contours_new = [c for c in contours if c.shape[0] >= min_nr]
    return contours_new


def triangulate_domain(vertices, segments, max_area, domain_shape):
    '''
    Triangulates an image domain using conforming constrained Delaunay 
    triangulation, with segments giving enforced connections. An input segment
    may be broken into smaller parts.
        
    Parameters
    ----------
    vertices : A (v,2) numpy array with coordinates of the vertices.
    segments : A (E,2) numpy array with edges represented as a pair of vertex 
        indices. Segments do not need to be oriented.
    max_area : Maximal triangle area for resulting triangluation.
    domainn_shape: Shape of the image domain to be triangulated. 

    Returns
    -------
    vertices : A (V,2) numpy array with coordinates of the vertices. Input 
        vertices constitute the firs v vertices of output set.
    triangles : A (T,3) numpy array with triangles represented as triplets of
        vertex indices.
    edges : A (E,2) numpy array with edges represented as pairs of vertex 
        indices. Returned edges are NOT oriented, even if input segments are.
        Some of the input segments may be broken into multiple edges by adding
        new vertices.
    '''
    
    corners = (np.array(domain_shape)-1) * np.array([[0,0], [0,1], [1,1], [1,0]])
    vertices = np.vstack((vertices, corners))
    B = triangle.triangulate(dict(vertices=vertices, segments=segments), 'pcqa' + 
                             str(max_area))
    vertices = B['vertices']
    triangles = B['triangles']
    segments = B['segments']
    return vertices, triangles, segments


def label_triangles_from_intensity(vertices, triangles, image, threshold):
    '''
    Lables triangles in a mesh according to the image intensity. 
    
    Parameters
    ----------
    vertices : A (V,2) numpy array with coordinates of the vertices.
    triangles : A (T,3) numpy array with triangles represented as triplets of
        vertex indices.
    image : A (r,c) numpy array with pixel intenisties.
    threshold : Intensity threshold determining whether label is 0 or 1.

    Returns
    -------
    labels : A (T,) numpy array with boolean triangle labels.
    '''
    
    nr,nc = image.shape
    xgrid, ygrid = np.mgrid[:nr, :nc]
    xy = np.vstack((xgrid.ravel(), ygrid.ravel())).T
    image = image.ravel()
    labels = np.array([np.mean(image[matplotlib.path.Path(vertices[t]).
                    contains_points(xy)])>threshold for t in triangles])    
    return labels


def label_triangles_from_contours(vertices, triangles, edges):

    '''
    Lables triangles in a mesh according to the edges of the contours. 
    Triangle adjacent to (0,0) will be labeled 0. Other triangles will be labeled
    with either 0 or 1, depending on whether the interface given by edges has
    been crossed.
    
    Parameters
    ----------
    vertices : A (V,2) numpy array with coordinates of the vertices.
    triangles : A (T,3) numpy array with triangles represented as triplets of
        vertex indices.
     edges : A (E,2) numpy array with edges represented as pairs of vertex 
        indices. Edges do not need to be oriented.

    Returns
    -------
    labels : A (T,) numpy array with boolean triangle labels.
    '''
    
    # Breadth-first search on triangles in a mesh, with labels switching 
    # when interface is crossed.
    corner = np.flatnonzero(np.all(vertices==0,axis=1))
    corner_triangles = np.flatnonzero(np.any(triangles==corner,axis=1))
    # the first triangle to process has vertex in (0,0)
    ct = corner_triangles[0] 
    se = np.sort(edges, axis=1) # preparation for easier interface-check
    ti = np.hstack((triangles, np.arange(triangles.shape[0]).reshape((-1,1))))
    # preparation for easier adjacency check
    ei = np.vstack((ti[:,[1,0,3]], ti[:,[2,1,3]], ti[:,[0,2,3]])) 
    labels = -np.ones(triangles.shape[0], dtype=int) # label -1 is for unlabeled
    queu = []
    labels[ct] = 0 # first triangle labeled with 0, others will be either 0 or 1
    queu.append(ct) 
    while len(queu)>0:
        # an array containing vertices of a triangle 
        tt = triangles[queu[0]]  
        # a list contining edges of a triangle
        es = [tt[[0,1]], tt[[1,2]], tt[[2,0]]] 
        
        # a list containing indices to adjacent triangles  
        ats = [ei[np.all(e == ei[:,0:2], axis=1), 2] for e in es]  
        # a list containing interface information
        ais = [np.any(np.all(np.sort(e)==se, axis=1)) for e in es] 
        for t, i in zip(ats, ais):
            # checking whether t exist and is unlabeled
            if (t.size) and (labels[t]<0): 
                labels[t] = np.logical_xor(labels[queu[0]], i) # label t
                queu.append(t.item()) 
        queu.pop(0)
    return labels.astype(bool)


def label_triangles_from_curves(vertices, triangles, edges_oriented, edges):
    '''
    Lables triangles in a mesh according to the oriented edges of the contours.
    Triangle adjacent to edges_oriented will be labeled witht 0 on the one side
    of the edges and 1 on the other side. This labeling will be propagated but
    not crossing the interfacee given by edges. Triangles left unllabelled 
    have label -1. 
    
    Parameters
    ----------
    vertices : A (V,2) numpy array with coordinates of the vertices.
    triangles : A (T,3) numpy array with triangles represented as triplets of
        vertex indices.
    edges_oriented : A (e,2) numpy array with oriented edges represented as 
        pairs of vertex indices.
    edges : A (E,2) numpy array with edges represented as pairs of vertex 
        indices. Edges do not need to be oriented. Usually edges_oriented is
        a subset of edges.

    Returns
    -------
    labels : A (T,) numpy array with integer triangle labels.
    '''
    
    # Labeling triangles adjacent to interface edges.
    ti = np.hstack((triangles, np.arange(triangles.shape[0]).reshape((-1,1))))
     # preparation for easier adjacency check
    ei = np.vstack((ti[:,[1,0,3]], ti[:,[2,1,3]], ti[:,[0,2,3]]))
    se = np.sort(edges, axis=1) # preparation for easier interface-check
    
    labels = -1*np.ones(triangles.shape[0], dtype=int)
    queu = []
    
    for e in edges_oriented:
        # an array containing index to adjacent triangle, or empty array
        t_next = ei[np.all(e == ei[:, 0:2], axis=1), 2] 
        # an index to opposite triangle
        t_opposite = ei[np.all(e[::-1] == ei[:, 0:2], axis=1), 2] 
        if t_next.size: # not empty    
            labels[t_next] = 1
            queu.append(t_next.item())
        if t_opposite.size:
            labels[t_opposite] = 0     
            queu.append(t_opposite.item())
            
    # Propagating labels without crossing interface   
    while len(queu)>0:
        # an array containing vertices of a triangle  
        tt = triangles[queu[0]] 
        # a list contining edges of a triangle
        es = [tt[[0,1]], tt[[1,2]], tt[[2,0]]]  
        
        # a list containing indices to adjacent triangles  
        ats = [ei[np.all(e == ei[:,0:2], axis=1), 2] for e in es]  
        # a list containing interface information
        ais = [np.any(np.all(np.sort(e)==se, axis=1)) for e in es] 
        for t, i in zip(ats, ais):
            # checking whether t exist, is unlabeled, and not over interface
            if (t.size) and (labels[t]<0) and (not i): 
                labels[t] = labels[queu[0]]
                queu.append(t.item())
        queu.pop(0)
    return labels

    
def save_mesh(filename, vertices, triangles, labels):
    '''
    Saves a file containing vertices, triangles and labels of a mesh.
    
    Parameters
    ----------
    filename : A filenime of the text file to be created.
    vertices : A (V,2) numpy array with coordinates of the vertices.
    triangles : A (T,3) numpy array with triangles represented as triplets of
        vertex indices.
    labels : A (T,) numpy array with integer (or boolean) triangle labels.
    '''
    
    file = open(filename,'w')
    for v in vertices:
        file.write('v {0:.2f} {1:.2f}\n'.format(v[0], v[1]))
    for t in triangles:
        file.write('f {0:d} {1:d} {2:d}\n'.format(t[0], t[1], t[2]))
    for l in labels:
        file.write('l {0:d}\n'.format(l))
    file.close() 
    
    
def load_mesh(filename):
    '''
    Loads a file containing vertices, triangles and labels of a mesh.
    
    Parameters
    ----------
    filename : A filenime of the text file containing mesh information.
        
    Returns
    -------
    vertices : A (V,2) numpy array with coordinates of the vertices.
    triangles : A (T,3) numpy array with triangles represented as triplets of
        vertex indices.
    labels : A (T,) numpy array with integer (or boolean) triangle labels.
    '''
    
    vertices = []
    triangles = []
    labels =[]
    file = open(filename, "r")
    for line in file:
        if line[0]=='v':
            vertices.append(np.array([float(n) for n in line[2:].split()]))
        if line[0]=='f':
            triangles.append(np.array([int(n) for n in line[2:].split()]))
        if line[0]=='l':
            labels.append(np.array([int(n) for n in line[2:].split()]))
    file.close()
    return np.array(vertices), np.array(triangles), np.array(labels)
 


