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

def calcAgreement(P, obIndex, centers_all, cell_diameter):
    agreement = 1
    for i, centers in enumerate(centers_all):
        if i==obIndex:
            continue
        minDistance = 9999
        for center in centers:
            distance = Utilities.euclidean_distance(P[0], P[1], center[0],center[1])
            if minDistance > distance:
                minDistance = distance
        if minDistance < cell_diameter:
            agreement+=1           
    
    return agreement


def ReadAllCenters(image_xmls, image_masks, stroma_pix_value):
    centers_all = []     
    
    for index, (obs_xml_file, obs_mask_file) in enumerate(zip(image_xmls, image_masks)):
                
        centers = Utilities.readCentersFromXmlWithMaskImage(obs_xml_file, obs_mask_file, stroma_pix_value)
        
        centers_all.append(centers)        
    
    return  centers_all

def RunDBCAA(centers_all, cell_diameter):
    predict_num = 0
    predict_num = 0
    for centers in centers_all:
        predict_num += len(centers)    
    
    obs_num = len(centers_all)
    agreements_cell = []
    
    
    for i, centers in enumerate(centers_all):        
        for center in centers:            
            agreement = calcAgreement(center, i, centers_all, cell_diameter)
            agreements_cell.append(agreement)
            
           
    distance_based_cell_agreement = np.mean(agreements_cell) / obs_num    
    
    return distance_based_cell_agreement




experiment_dir = r"c:\Users\DELL\anaconda\inter_observer\github\sample_data\experiment1"
xml_file_name = "lymphocyte.xml"
mask_file_name = "mask.png"
observer_list = ["observer1", "observer2", "observer3", "observer4" ]

experiments_coor = Utilities.ReadFilePathsFromExperimentFolder(experiment_dir, xml_file_name, observer_list)
experiments_mask = Utilities.ReadFilePathsFromExperimentFolder(experiment_dir, mask_file_name, observer_list)
print(experiments_coor)

DBCAA_results = []
for index, (image_coors, image_masks) in enumerate(zip(experiments_coor, experiments_mask)):
    centers_all = ReadAllCenters(image_coors[1], image_masks[1], Utilities.STROMA_PIX_VALUE)
    distance_based_cell_agreement = RunDBCAA(centers_all, 20);
    DBCAA_results.append(distance_based_cell_agreement)
    print(distance_based_cell_agreement)
    

DBCAA_results = np.mean(DBCAA_results)
print(DBCAA_results)



