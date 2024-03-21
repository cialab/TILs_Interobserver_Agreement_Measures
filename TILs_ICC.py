# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 15:29:53 2024

@author: Abdulkerim Capar
"""
import os
import math
from PIL import Image
import numpy as np
import xml.etree.ElementTree as ET
import icc
import Utilities


def CircleArea(radius):
    area = math.pi * radius**2
    return area

def CountPix(gray_image, pix_val):
    image_array = np.array(gray_image)
    count = np.sum(image_array == pix_val)
    return count

def CalculateTILsScore(xml_file_path, mask_file_path):
    try:
        mask_image = Image.open(mask_file_path)
        width, height = mask_image.size
        stromal_area = CountPix(mask_image, Utilities.STROMA_PIX_VALUE)               
        
        cell_area = CircleArea(Utilities.CELL_RADIUS)
        tot_cell_area = 0
        
        bounding_boxes = Utilities.readCentersFromXml(xml_file_path, width, height)
        
        for box in bounding_boxes:
            center_x, center_y = box
            
            #check the cell inside stroma
            if (center_x < width) and (center_y < height) and (mask_image.getpixel((center_x, center_y)) == Utilities.STROMA_PIX_VALUE):
                tot_cell_area = tot_cell_area + cell_area 
        
        
        tils_score = 0
        if stromal_area > 0:
           tils_score = 100 * tot_cell_area /  stromal_area
        
        return tils_score

    except Exception:
        pass

def RunScoreAgreementAll(experimentsXml, experimentsMask):
        
    result_list = []
    result_matrix = []
    image_list = []
    
    
    for index1, (image_xmls, image_masks) in enumerate(zip(experimentsXml, experimentsMask)):
        image_tils_scores = []
        image_list.append(image_xmls[0])
        for index2, (obs_xml, obs_mask) in enumerate(zip(image_xmls[1], image_masks[1])): 
            
            tils_score = CalculateTILsScore(obs_xml, obs_mask)
            image_tils_scores.append(tils_score)
            result_list.append((index1, index2, tils_score))
        
        result_matrix.append(image_tils_scores)
    
    
    return np.array(result_matrix), np.array(result_list), image_list 
    
    

experiment_dir = r"c:\Users\DELL\anaconda\inter_observer\github\sample_data\experiment1"
xml_file_name = "lymphocyte.xml"
mask_file_name = "mask.png"
observer_list = ["observer1", "observer2", "observer3", "observer4" ]


files_xml = Utilities.ReadFilePathsFromExperimentFolder(experiment_dir, xml_file_name, observer_list)
files_mask = Utilities.ReadFilePathsFromExperimentFolder(experiment_dir, mask_file_name, observer_list)
# print(files_xml)


result_matrix, result_list, image_list = RunScoreAgreementAll(files_xml, files_mask)


icc_val, Fvalue, lbound, ubound = icc.icc_func(result_matrix) 

print(result_list)
print(icc_val)



