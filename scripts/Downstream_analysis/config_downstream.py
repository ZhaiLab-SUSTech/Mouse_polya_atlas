# -*- coding: utf-8 -*-
"""
Shared configuration variables for the downstream analysis pipeline.
"""

# The primary list of tissues, defining the order for most operations.
TISSUE_LIST = [
    'brain', 'thyroid', 'thymus', 'heart', 'lung', 'liver', 'spleen', 
    'pancreas', 'stomach', 'small', 'large', 'adrenal', 
    'kidney', 'muscle', 'adipose', 'bone', 'testis', 'sperm'
]


# Dictionary to map short tissue names to full names for plotting.
SHORT_TO_LONG_NAME = {
    'brain': 'Brain',
    'thyroid': 'Thyroid',
    'thymus': 'Thymus',
    'heart': 'Heart',
    'lung': 'Lung',
    'liver': 'Liver',
    'spleen': 'Spleen',
    'pancreas': 'Pancreas',
    'stomach': 'Stomach',
    'small': 'Small intestine',
    'large': 'Large intestine',
    'adrenal': 'Adrenal gland',
    'kidney': 'Kidney',
    'muscle': 'Muscle',
    'adipose': 'Adipose',
    'bone': 'Bone marrow',
    'testis': 'Testis',
    'sperm': 'Sperm'
}