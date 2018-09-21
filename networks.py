# -*- coding: utf-8 -*-


"""
Created on Tue Jun 12 2018
Last edit: Tue Jun 12 2018
@author: kfinc

"""

import numpy as np
import pandas as pd

def calculate_lsn_edges (A, labels):
    """Function calculates number of edges between and within predefined large-scale networks (LSNs).
    The function takes binary symetrical adjacency matrix, module assignment of each ROI and calculate number of edges between and within each
    large-scale network.

    Parameters
    ------------
    array: N x N binary ajdacency matrix
    array: N-length vector with module assignment for each node

    Returns
    ------------
    array: M x M matrix with number of edges between each module

    """
    columns = np.unique(labels)
    lsn_matrix = np.zeros((len(labels), len(columns)))
    lsn_edges = np.zeros((len(columns), len(columns)))

    for col in range(len(columns)):
        module = columns[col, ]
        for row in range(len(labels)):
            if (labels[row, ] == module):
                lsn_matrix[row, col] = 1
            else:
                lsn_matrix[row, col] = 0

    lsn_edges = lsn_matrix.T @ A @ lsn_matrix
    return lsn_edges


def allegiance_matrix(M):
    """Calculates M x M allegiance matrix from M x N matrix, where M represents nodes, and N represents time windows. 
    Each value on allegiance matrix represents the probability that node i and node j have been assigned to the same 
    functional community (Bassett et al., 2014).
    
    Parameters
    ------------
    array: M x N module assignment matrix 

    Returns
    ------------
    array: M x M allegiance matrix
        
    """
    
    roi_n = len(M[:,0])
    A = np.zeros((roi_n, roi_n))
    
    for i in range(roi_n):
        for j in range(i):
            if i == j:
                continue
            else:
                vector = M[i, :] == M[j, :]
                p = vector.mean()
                A[i, j] = p
    A = A + A.T
    return A


def allegiance_matrices_4d(M):
    """Calculates 4D array composed of allegiance matrices for each subject and each condition/session.
    
    Parameters
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x M(window) 

    Returns
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x N(node), where N x N is allegiance matrix 
        
    """
    
    AM = np.zeros((46, 4, 264, 264))
    sub_n = len(M[:,0,0,0])
    ses_n = len(M[0,:,0,0])
    
    for sub in range(sub_n):
        for ses in range(ses_n):
            AM[sub, ses, :, :] = allegiance_matrix(M[sub, ses, :, :])
    return AM    
    

def sort_matrices_4d(M, idx):
    """Sorts matrices according to predefinded index.
    
    Parameters
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x N(node) (unsorted)
    array: N-length vector with index to sort matrix

    Returns
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x N(node) (sorted)
    
    """
    M1 = M[:,:,:,idx]
    M2 = M1[:,:,idx,:]
    return M2


def dumming(labels):
    """Gennerate vectors with dummy variables from single vector with categorical variable
    
    Parameters
    ------------
    array: np.array, N-length vector with categorical variable

    Returns
    ------------
    array: 3D array with M rows and N columns with binary variables
    
    
    """

    columns = np.unique(labels)
    dummies = np.zeros((len(labels), len(columns)))

    for col in range(len(columns)):
        module = columns[col]
        for row in range(len(labels)):
            if (labels[row] == module):
                dummies[row, col] = 1
            else:
                dummies[row, col] = 0
    return dummies


def dumming_pd(labels):
    """Gennerate vectors with dummy variables from single vector with categorical variable
    
    Parameters
    ------------
    array: np.array, N-length vector with categorical variable

    Returns
    ------------
    array: pd.DataFrame with M rows and N columns with binary variables
    
    """
    columns = np.unique(labels)
    return pd.DataFrame(dumming(labels), columns=columns)

def upper_tri_masking(A):
    """Getting values of upper triangle of matrix without diagonal"""
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:,None] < r
    return A[mask]


def fc_cartography(M, modules):
    """Function which calculates mean integration and recruitment values from sorted allegiance matrices"""
    dummy_networks = dumming(sorted(modules))
    roi_n = len(dummy_networks)
    net = np.size(dummy_networks,1)
    diagnostics = np.zeros((net, net))

    for i in range(net):
        for j in range(net):
            vec1 = dummy_networks[:,i].astype('bool')
            vec2 = dummy_networks[:,j].astype('bool')

            L = M[vec1,:]
            P = L[:,vec2]

            if i == j:
                m = upper_tri_masking(P).mean()
            else:
                m = P.mean()

            diagnostics[i, j] = m
    return diagnostics


def fc_cartography_4d(M, modules):
    """Function which calculates mean integration and recruitment values from sorted allegiance matrices. 4D version"""
    sub_n = len(M[:,0,0,0])
    ses_n = len(M[0,:,0,0])
    mod_n = len(np.unique(modules))

    fc_cart = np.zeros((sub_n, ses_n, mod_n, mod_n))

    for i in range(sub_n):
        for j in range(ses_n):
            X = M[i,j,:,:]
            fc_cart[i, j, :, :] = fc_cartography(X, modules)
    return fc_cart

