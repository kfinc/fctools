import numpy as np

def m_wei(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov_wei(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m_wei(x, w)) * (y - m_wei(y, w))) / np.sum(w)

def corr_wei(x, y, w):
    """Weighted Correlation"""
    return cov_wei(x, y, w) / np.sqrt(cov_wei(x, x, w) * cov_wei(y, y, w))
