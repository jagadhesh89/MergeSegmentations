import numpy as np
from annoy import AnnoyIndex
import numpy as np
from matplotlib import pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')
import joblib
from tiatoolbox.utils.visualization import  overlay_prediction_contours
import matplotlib as mpl
from tiatoolbox.utils.visualization import overlay_prediction_contours
from tiatoolbox.wsicore.wsireader import WSIReader
import argparse
import pickle

def get_centroid_types(pred):
    centroids = []
    cell_types = []
    box = []
    contour = []
    prob = []
    pred_id = []
    for eachpred in pred:
        centroids.append(pred[eachpred]["centroid"])
        cell_types.append(pred[eachpred]["type"])
        box.append(pred[eachpred]["box"])
        contour.append(pred[eachpred]["contour"])
        prob.append(pred[eachpred]["prob"])
        pred_id.append(eachpred)
    return centroids, cell_types, box, contour, prob, pred_id


def change_mon_types(mon_types):
    updated_mon_types = []
    for eachtype in mon_types:
        if eachtype == 0:
            updated_mon_types.append(0)
        #mon is lymph and pan is inflam
        if eachtype == 1:
            updated_mon_types.append(1)
        #mon is macro and pan is inflam
        if eachtype == 2:
            updated_mon_types.append(1)
        #mon is macro and pan is inflam
        if eachtype == 3:
            updated_mon_types.append(1)
    return updated_mon_types

def check_type_for_merge(pan_type, mon_type):
    matchflag = 0
    if pan_type == mon_type:
        matchflag = 1
    return matchflag


def merge_coordinates(coords_method1, coords_method2, radius_threshold, pan_types, mon_types, pan_prob, mon_prob, updated_mon_types):
    merged_coords = []
    dim = len(coords_method1[0])

    # Build Annoy index with the coordinates from the second method
    t = AnnoyIndex(dim, 'euclidean')
    for i, coord in enumerate(coords_method2):
        t.add_item(i, coord)
    t.build(10)


    merged_indices_method1 = []
    merged_indices_method2 = []
    merged_types = []
    merged_source = []

    # For each coordinate in method1, find coordinates in method2 within the radius threshold
    for i, coord1 in enumerate(coords_method1):
        #print(len(coords_method1))
        indices = t.get_nns_by_vector(coord1, n=50, search_k=-1, include_distances=True)
        #print(len(indices))
        for j, distance in zip(*indices):
            #print(j, distance)
            if distance <= radius_threshold:
                pan_type = pan_types[i]
                mon_type = mon_types[j]
                merged_indices_method1.append(i)
                merged_indices_method2.append(j)
                coord2 = coords_method2[j]
                merged_coords.append(((coord1[0] + coord2[0]) / 2, (coord1[1] + coord2[1]) / 2))
                #print(updated_mon_types[j])
                #break

                matchflag = check_type_for_merge(pan_type, updated_mon_types[j])
                if matchflag == 0:
                    if pan_prob[i] > mon_prob[j]:
                        merged_types.append(pan_type)
                        merged_source.append(0)
                        continue
                    else:
                        merged_types.append(int(mon_type)+5)
                        merged_source.append(1)
                        continue
                else:
                    merged_types.append(int(mon_type)+5)
                    merged_source.append(1)



    # non-merged coordinates from both methods
    non_merged_coords1 = [coord for i, coord in enumerate(coords_method1) if i not in merged_indices_method1]
    non_merged_coords2 = [coord for i, coord in enumerate(coords_method2) if i not in merged_indices_method2]
    #non_merged_coords1 = []
    #non_merged_coords2 = []


    return merged_coords, non_merged_coords1, non_merged_coords2, merged_indices_method1, merged_indices_method2, merged_types, merged_source

def main():
    parser = argparse.ArgumentParser(description='merge_preds')

    parser.add_argument("-p", "--model1pred",
                        dest="model1",
                        help="The path to the dat file of model1 (pannuke)",
                        required=True)
    parser.add_argument("-m", "--model2pred",
                        dest="model2",
                        help="The path to the dat file of model2 (monusac)",
                        required=True)
    parser.add_argument("-o", "--outputpath",
                        dest="outdir",
                        help="The path to output directory",
                        required=True)

    inpArgs = parser.parse_args()
    pan_preds = os.path.abspath(inpArgs.model1)
    mon_preds = os.path.abspath(inpArgs.model2)
    
    slidename=inpArgs.pannuke.split("/")[-2].replace(".dat","")

    wsi_pred_mon = joblib.load(mon_preds)
    wsi_pred_pan = joblib.load(pan_preds)

    mon_centroids, mon_types, mon_box, mon_countour, mon_prob, mon_id = get_centroid_types(wsi_pred_mon)
    pan_centroids, pan_types, pan_box, pan_countour, pan_prob, pan_id = get_centroid_types(wsi_pred_pan)

    updated_mon_types = change_mon_types(mon_types)

    mon_pan_mapping = {0:0, 0:4, 1:6, 2:7, 3:8}

    mon_map = {0:"epi", 1:"lymphocyte", 2:"macrophage", 3:"neutrophil"}
    pan_map = {0:"neo-epi", 1:"inflam", 2:"connective", 3:"dead", 4:"non-neo-epi"}

    radius_threshold = 10

    merged_coords, non_merged_coords1, non_merged_coords2, merged_indices_method1,merged_indices_method2, merged_types, merged_source = merge_coordinates(pan_centroids, mon_centroids, radius_threshold, pan_types, mon_types, pan_prob, mon_prob, updated_mon_types)

    combined_dict = {}
    combined_pan_ids = []
    combined_mon_ids = []
    merged_indices_method1 = list(merged_indices_method1)
    merged_indices_method2 = list(merged_indices_method2)
    for i in range(0,len(merged_coords)):
        index_m1 = merged_indices_method1[i]
        index_m2 = merged_indices_method2[i]
        pan_id_req = pan_id[index_m1]
        mon_id_req = mon_id[index_m2]
        combined_pan_ids.append(pan_id_req)
        combined_mon_ids.append(mon_id_req)
        box_req = wsi_pred_pan[pan_id_req]['box']
        centroid_req = list(merged_coords[i])
        contour = wsi_pred_pan[pan_id_req]['contour']
        source = merged_source[i]
        if source == "0":
            prob_req = wsi_pred_pan[pan_id_req]['prob']
            type_req = wsi_pred_pan[pan_id_req]['type']
        else:
            prob_req = wsi_pred_mon[mon_id_req]['prob']
            type_req = int(wsi_pred_mon[mon_id_req]['type'])
        combined_dict[pan_id_req] = {'box':box_req, 'centroid':centroid_req,'contour':contour, 'prob':prob_req, 'type':type_req}

    
    typestest = {}
    for eachpred in wsi_pred_mon:
        if eachpred not in combined_mon_ids:
            '''combined_dict[eachpred] = wsi_pred_mon[eachpred]
            typestest[combined_dict[eachpred]['type']] = 1'''
            box_req = wsi_pred_mon[eachpred]['box']
            centroid_req = wsi_pred_mon[eachpred]['centroid']
            contour = wsi_pred_mon[eachpred]['contour']
            prob_req = wsi_pred_mon[eachpred]['prob']
            type_req = int(wsi_pred_mon[eachpred]['type']) + 5
            combined_dict[eachpred] = {'box':box_req, 'centroid':centroid_req,'contour':contour, 'prob':prob_req, 'type':type_req}
    #print(typestest)

    for eachpred in wsi_pred_pan:
        if eachpred not in combined_pan_ids:
            box_req = wsi_pred_pan[eachpred]['box']
            centroid_req =  wsi_pred_pan[eachpred]['centroid']
            contour = wsi_pred_pan[eachpred]['contour']
            prob_req = wsi_pred_pan[eachpred]['prob']
            type_req = int(wsi_pred_pan[eachpred]['type'])
            combined_dict[eachpred] = {'box':box_req, 'centroid':centroid_req,'contour':contour, 'prob':prob_req, 'type':type_req}
            
    with open(outdir + "/" + slidename+ ".dat", 'wb') as file:
        pickle.dump(combined_dict, file)

if __name__ =="__main__":
    main()   