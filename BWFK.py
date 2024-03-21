# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 15:29:53 2024

@author: Abdulkerim Capar
"""
import os
from PIL import Image
import numpy as np
from scipy.ndimage import distance_transform_edt
import statsmodels.stats.inter_rater as ir 
import Utilities

def create_boundary_image(binary_image):
    rows, cols = binary_image.shape
    # boundary_image = np.zeros_like(binary_image)
    boundary_image = np.ones_like(binary_image) * 1

    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current pixel has opposite value neighbors
            if binary_image[i, j] == 1:
                if binary_image[i, j] != binary_image[i-1, j] or \
                   binary_image[i, j] != binary_image[i+1, j] or \
                   binary_image[i, j] != binary_image[i, j-1] or \
                   binary_image[i, j] != binary_image[i, j+1]:
                       boundary_image[i, j] = 0

    return boundary_image

def fleissKappaWeighted(rate, n, weights):
    N = len(rate)
    k = len(rate[0])
    # print("#raters = ", n, ", #subjects = ", N, ", #categories = ", k)
    
    weights = N * weights / sum(weights)

    # Calculate PA
    sum_pa = 0
    for row, w in zip(rate, weights):
        sum_row = 0
        for i in row:
            sum_row += i**2
        pi = (sum_row - n) / (n * (n - 1))
        sum_pa += w * pi
    PA = sum_pa / N
    # print("PA = ", PA)
    
    # Calculate PE
    sum_pe = 0
    for i in range(k):
        sum_pe_j = 0
        for row, w in zip(rate, weights):
            pj = row[i] / (N * n)
            sum_pe_j += w * pj
        sum_pe += sum_pe_j**2
    PE = sum_pe
    # print("PE =", PE)


    kappa = -float("inf")
    try:
        kappa = (PA - PE) / (1 - PE)
        kappa = float("{:.3f}".format(kappa))
    except ZeroDivisionError:
        print("Expected agreement = 1")

    # print("Fleiss' Kappa =", kappa)

    return kappa


def RunBWFK(image_paths, resizeW, resizeH, DL):
    pix_matrix = []
    distance_transforms = []
    
    for i, path in enumerate(image_paths):

        image = Image.open(path)
        image2 = image.resize((resizeW,resizeH), Image.NEAREST)        
        
        image_data = np.array(image2)  
        image_data[image_data == 0] = 0 
        image_data[image_data == 1] = 0 
        image_data[image_data == Utilities.STROMA_PIX_VALUE] = 1 #stroma
        
        boundary_data = create_boundary_image(image_data)

        distance_transform = distance_transform_edt(boundary_data)
        
        distance_transforms.append(distance_transform) 
        
        pixel_vector = image_data.flatten()
        pix_matrix.append(pixel_vector)

    
    
    average_distance_transform = np.mean(distance_transforms, axis=0)
    
    average_distance_transform[average_distance_transform > DL] = DL
    
    pix_matrix =  np.column_stack(pix_matrix)
    rater_data = ir.aggregate_raters(pix_matrix)    

    weights = average_distance_transform.flatten()

    boundary_weighted_fleiss_kappa = fleissKappaWeighted(rater_data[0], len(pix_matrix[0]), weights)  
    
    return boundary_weighted_fleiss_kappa

        

experiment_dir = r"c:\Users\DELL\anaconda\inter_observer\github\sample_data\experiment1"
xml_file_name = "lymphocyte.xml"
mask_file_name = "mask.png"
observer_list = ["observer1", "observer2", "observer3", "observer4" ]

files_mask = Utilities.ReadFilePathsFromExperimentFolder(experiment_dir, mask_file_name, observer_list)
print(files_mask)

BWFK_results = []
for index, image_masks in enumerate(files_mask):
    boundary_weighted_fleiss_kappa = RunBWFK(image_masks[1], 300, 300, 5)
    BWFK_results.append(boundary_weighted_fleiss_kappa)
    print(boundary_weighted_fleiss_kappa)

mean_BWFK = np.mean(BWFK_results)
print(mean_BWFK)


