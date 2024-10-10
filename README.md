# MergeSegmentations

Usage: 


python merge_preds.py -m <method1_predictions> -p <method2_predictions>

The input to the scripts are prediction files from the monusac and pannuke models. This can be tweaked to take inputs from any similar models as long as the outputs follow a datastructure that is similar to the format described as follows:

The file is a pickle dump with underlying datastructure as a dictionary. Where each key corresponds to a unique nuclei id. 

For example:
{
  'id1':{'box':[],'centroid':[],'contour':[], 'prob':[], 'type':[]}
  'id2':{'box':[],'centroid':[],'contour':[], 'prob':[], 'type':[]}  
}

The box describes the bounding box of the nuclei
The centroid describes the centroid x,y coordinates of the nuclei
The countour describes the boundary of the nuclei in coordinates
The prob details the probability of the prediction of nuclei
The type details the type of the nuclei, i.e epithelial etc. 

The primary function of this code is the "merge_coordinates" function in the script that uses [Annoy](https://github.com/spotify/annoy) to merge the predictions.

The output is a ".dat" file which is a pickle file, which has all the combined/integrated/merged predictions in a datastructure described in the example above. 
