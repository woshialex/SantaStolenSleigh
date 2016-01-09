#Follow these steps:

##put gifts.csv file under the following directory and all the results will be saved to this directory too.
datadir='/home/qiliu/Documents/Data/santa/';

#1. generate initial clusters: generate the initial clusters which is based on the longitude value of the gifts, all the trips basically go straight from south to north, with maximum capacity of 980. The lat>-60 and lat<-60 are treated as two separate region.

##uncomment this line to run step 1!!!
#python initial.py $datadir

#2. compile the code

##uncomment this line to run step 2!!!
#make

#3. run the parallel C++ code to do global optimization over all trips
# cool_speed=1.0002 with 4 cores takes about 4 days to get score 1.23890e13
BETA=4;
cool_speed=1.0002;
Nthread=4;

##uncomment this line to run step 3!!!!
#./santa.x $datadir goodinitial_980.csv $BETA $cool_speed 0.01 0.02 500 result.csv $Nthread

#to get our leaderboard score 1.238507e13, run with cool_speed=1.00002 on 4 cores for 10 days,
#cool_speed=1.00002;
#./santa.x $datadir goodinitial_980.csv $BETA $cool_speed 0.01 0.02 500 result.csv $Nthread
#then continue the solution that is saved with score 1.238667e13 at beta=36.5 during the run, use speed=1.00005 to run for about another 5 days to converge.
#./santa.x $datadir 1238667.csv 36.5 $cool_speed 0.01 0.02 500 final_result.csv $Nthread

## -- 4. make the submission file

##uncomment this line to run step 4!!!
#cut -d, -f1,5 $datadir/result.csv > $datadir/submission.csv
