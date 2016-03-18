north_pole = (90,0)
weight_limit = 1100
sleigh_weight = 10

import pandas as pd
import numpy as np
from haversine import haversine

def weighted_trip_length(stops): 
    tuples = [tuple(x) for x in stops[['Latitude','Longitude']].values]
    weights = stops.Weight.tolist();
    # adding the last trip back to north pole, with just the sleigh weight
    tuples.append(north_pole)
    weights.append(sleigh_weight)
    
    dist = 0.0
    prev_stop = north_pole
    prev_weight = sum(weights)
    for i, tup in enumerate(tuples):        
        dist += haversine(tup, prev_stop) * prev_weight
        prev_stop = tup
        prev_weight -= weights[i]
    return dist

def weighted_reindeer_weariness(all_trips):
    uniq_trips = all_trips.TripId.unique()
    
    if any(all_trips.groupby('TripId').Weight.sum() > weight_limit - sleigh_weight):
        print("One of the sleighs over weight limit!")
 
    dist = 0
    for t in uniq_trips:
        this_trip = all_trips[all_trips.TripId==t]
        dist += weighted_trip_length(this_trip);
    
    return dist    


if __name__ == '__main__':
    datadir = sys.argv[1];
    gifts = pd.read_csv(datadir+'gifts.csv')
    sample_sub = pd.read_csv(datadir+'sample_submission.csv')

    all_trips = sample_sub.merge(gifts, on='GiftId')

    print(weighted_reindeer_weariness(all_trips))
