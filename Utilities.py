# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 15:05:07 2024

@author: Abdulkerim Capar
"""
import os
from PIL import Image
import xml.etree.ElementTree as ET
import math

# Constants
STROMA_PIX_VALUE = 2
CELL_RADIUS = 12

def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

# Read box annotation coordinates form ASAP formatted xml file, without mask control
def readCentersFromXml(xml_file_path, imgW, imgH):
    
    try:
        centers = []
        
        tree = ET.parse(xml_file_path)
        root = tree.getroot()    
        for annotation in root.iter('Annotation'):
            xmin = ymin = float('inf')
            xmax = ymax = float('-inf')
            for coordinate in annotation.iter('Coordinate'):
                x, y = int(float(coordinate.get('X'))), int(float(coordinate.get('Y')))
                xmin = min(xmin, x)
                ymin = min(ymin, y)
                xmax = max(xmax, x)
                ymax = max(ymax, y)
                
            center_x = int((xmin + xmax) / 2)
            center_y = int((ymin + ymax) / 2)
            
            if (center_x < imgW) and (center_y < imgH):                
                centers.append((center_x, center_y))
        return centers

    except Exception:
        pass

# Read box annotation coordinates form ASAP formatted xml file checking that the box centers are on foreground
def readCentersFromXmlWithMaskImage(xml_file_path, mask_file_path, foreground_pix_val):
    
    try:
        centers = []
        mask_image = Image.open(mask_file_path)
        width, height = mask_image.size
        
        tree = ET.parse(xml_file_path)
        root = tree.getroot()    
        for annotation in root.iter('Annotation'):
            xmin = ymin = float('inf')
            xmax = ymax = float('-inf')
            for coordinate in annotation.iter('Coordinate'):
                x, y = int(float(coordinate.get('X'))), int(float(coordinate.get('Y')))
                xmin = min(xmin, x)
                ymin = min(ymin, y)
                xmax = max(xmax, x)
                ymax = max(ymax, y)
                
            center_x = int((xmin + xmax) / 2)
            center_y = int((ymin + ymax) / 2)
            
            if (center_x < width) and (center_y < height) and (mask_image.getpixel((center_x, center_y)) == foreground_pix_val):            
                centers.append((center_x, center_y))
        return centers

    except Exception:
        pass

# Read <file_name> named images in <image_level_folder_path> folder for each <observer_list> folder 
def ReadFilePathsFromImageFolder(image_level_folder_path, file_name, observer_list):
    files_all = []
    for index, observer_level_folder in enumerate(observer_list):
        observer_level_folder_path = os.path.join(image_level_folder_path, observer_level_folder)

        if os.path.isdir(observer_level_folder_path):
            file_path = os.path.join(observer_level_folder_path, file_name)       
            files_all.append(file_path) 
    
    return files_all

# Read <file_name> named images in <experiment_level_folder_path> folder for each image folder and for each <observer_list> folder 
def ReadFilePathsFromExperimentFolder(experiment_level_folder_path, file_name, observer_list):
    files_all = []
    for index, image_level_folder in enumerate(os.listdir(experiment_level_folder_path)):
        image_level_folder_path = os.path.join(experiment_level_folder_path, image_level_folder)
        
        if not os.path.isdir(image_level_folder_path):
            continue
        
        files = ReadFilePathsFromImageFolder(image_level_folder_path, file_name, observer_list)
        files_all.append((image_level_folder, files))
    
    return files_all



# experiment_dir = r"c:\Users\DELL\anaconda\inter_observer\github\sample_data\experiment1"
# xml_file_name = "lymphocyte.xml"
# mask_file_name = "mask.png"
# observer_list = ["observer1", "observer2", "observer3", "observer4" ]

# files_xml = ReadFilePathsFromExperimentFolder(experiment_dir, xml_file_name, observer_list)
# print(files_xml)
