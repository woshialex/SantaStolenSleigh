/*
 * Santa Stolen Sleigh competition - greedy + Simulated Annealing algorithm
 * Qi Liu, Mirsad Buljubasic
 */


// ADD LICENCE PART



/********************************************* SHORT ALGORITHM DESCRIPTION *****************************************************

-INPUT: Starting solution to the problem
-Local Search (Simulated Annealing) algorithm is performed starting from the input solution. No problem
decomposition is made.
-Only feasible moves considered during entire procedure
-Only moves between near trips (routes) are performed: defined by "lat_width" and "lon_width" parameters
-Moves inside a single route and moves between the routes

Moves inside a single route:
1) SWAP(i,i+1) - swap positions of 2 subsequent gifts i and i + 1
2) Permutation(i, N) - evaluate all possible permutations of gifts between positions i and i+N (N = 5 is used)
   and take the best one
3) Shift(i,j) - shift a gift at position i to position j

Moves between the routes:
1) Shift - shift gift from one route to another: Gift to shift and route to shift to are chosen randomly
   (respecting constraints and range) and the best position to shift to is chosen (in terms of objective function)
2) Swap - swap two gifts between the routes. Two routes and first gift are chosen randomly, while the second gift
   is chosen to minimize the objective.
2) ShiftBlock (or ShiftN) - shift a block of gifts (gifts between positions i and i + N) from one route to another: routes and
   beginning position i are chosen randomly as before.
3) SwapBest - choose two gifts g1 and g2 from routes r1 and r2 respectively, drop them from the current routes and then
   insert g1 to best position in r2 and g2 to best position in r1 (can be seen as generalization of Swap move)
4) CROSS (or SwapTails) - exchange tail blocks between two trips


The last move is
CUT - cutting a single trip into two. Cutting is simply done by breaking the trip at the best position.

For  Permutation(i, N) and  Shift(i,j) only the moves improving the objective are accepted,
while all other moves are accepted according to the current state (temperature).

Initial temperature and cooling speed are chosen experimentally.

Moves evaluation is done in parallel, with "NMB_THREADS" being number of threads utilized

*****************************************************************************************************************/



#include "gift.h"
#include <ctime>


// list of all trips contained in solution
// filled initially by reading from the input solution file
vector<Trip> trips;

// data structure that stores all trip pairs
// local search move between two trips can be performed only if pair belongs to the list
// construction is based on mean longitude of the trips
vector<pair<int,int> > trip_pairs;

// local moves between two trips are allowed only if difference between mean latitudes of trips does not exceed given value (lat_width)
// this limitation has been made in order to speed up the search
// expressed in radians
double latitude_difference_threshold;
double longitude_difference_threshold;

//TODO
// indicate which moves will be performed
// will be set in main
vector<double> prob(4,1.0);



void parallel_job(void *sucp, int id)
{

	// between routes optimization
	vector<int> *suc = (vector<int>*) sucp;
	for(int j = 0; j < 100; j++) //don't use too many threads
	{
		vector<int> sucOne = local_search_between_routes(trips, trip_pairs, latitude_difference_threshold,prob,id);
		mtx[MAXTRIP-1].lock();
		for(int k = 0; k < suc->size(); k++) (*suc)[k] += sucOne[k];
		mtx[MAXTRIP-1].unlock();
	}

	//inside route optimization
	int sucopt = selfupdate(trips, id);
	int exactN = 5;
	double ssuc = selfupdateExact(trips, exactN, id);
	double s1move = selfupdate1move(trips,latitude_difference_threshold, id);

	//swapping tails
	double tsscore = swaptails(trips, trip_pairs, id);
}



int main(int argc, char *argv[])
{

	srand(1);
	cout.precision(10);

	//-------------------- program parameters ---------------------------
	// 10 parameters required
	if(argc != 10)
	{
		cerr << "need 10 parameters" << endl;
		cerr << "santa.x datadir start_file_name BETA cool_speed lon_width lat_width stopping_cond save_file_name N_THREADS" << endl;
		return 0;
	}

	string datadir = string(argv[1]); // data directory

	// read input solution
	trips = load_data(datadir + string(argv[2]));
	for(auto &t:trips) t.update_info();
	cout << "total trips " << trips.size() << endl;
	double ss = totalscore(trips);
	cout << "initial score: " << ss << endl;
	int N = trips.size();

	// initial temperature and cooling rate
	BETA = atof(argv[3]);                  //annealing starting beta = 1/T (temperature)
	double cooling_speed = atof(argv[4]);  //annealing cooling speed

	// local search range limitation
	longitude_difference_threshold = atof(argv[5]);
	latitude_difference_threshold = atof(argv[6]);

	// stopping condition: stop when number of performed moves in the last iteration is less than a given value (stopping_suc)
	int stopping_suc = atoi(argv[7]);
	bool sep = false;

	//--------------------------------------------------------------------

	//initialize the random generator
	for(int i = 0; i < 32; i++) rs[i] = rand();
	
	prob[1] = prob[2] = prob[3] = 1;
	NMB_THREADS = atoi(argv[9]);
	assert(NMB_THREADS <= 16);    //don't use too many cores or else wasted


	double goodresult = min(1.240e13, ss) - 1e7;
	time_t start = time(NULL);
	vector<int> suc(prob.size(), 0);

	// Simulated Annealing procedure
	int nmbIters = 5000000;
	for(int i = 1; i < nmbIters; i++)
	{
		assert(trips.size() < MAXTRIP);          //mtx array size limit
		for(auto &t:trips) t.update_range();
		int Npair = set_working_pairs(trip_pairs, trips, sep, longitude_difference_threshold);

		// do it in parallel
		vector<thread> jobs(NMB_THREADS);
		for(int l = 0; l < NMB_THREADS; l++) jobs[l] = thread(parallel_job, &suc, l);
		for(int l = 0; l < NMB_THREADS; l++) jobs[l].join();

		if(i % 10 == 0)
		{
			for(auto &t:trips){ t.update_info(); }
			//recalculate auxiliary variables to avoid accumulated errors
			//remove empty trips
			int b = trips.size();
			trips.erase(remove_if(trips.begin(),trips.end(), [](Trip &t){return t.size()<=2;}), trips.end());
			int e = trips.size();
			if(b != e) cout << " empty trips removed : " << b - e << endl;

			double cutscore = selfupdateCut(trips, 0);

			time_t now = time(NULL);
			cout << " - Estimated time remaining: " << log(2000/BETA) / log(cooling_speed) * difftime(now,start)/i/3600.0 << " hours" << endl;
		}
		
		if(i % 2 == 0)
		{
			ss = totalscore(trips);
			cout << i << " :";
			for(auto x:suc) cout << x << ' ';
			cout << shift11_count << ' ';
			shift11_count = 0;
			cout << " score: " << ss << endl;
			if(ss < goodresult)
			{
				cout << "save at beta = " << BETA << endl;
				save_result(trips,datadir + to_string(int(ss/1e7))+".csv");
				goodresult = ss - 1e7;
			}

			if(suc[0] < stopping_suc)break;
			for(int k = 0; k < suc.size(); k++) suc[k] = 0;
		}

		BETA *= cooling_speed;

	}


	//end

	time_t end = time(NULL);
	cout << "time: " << difftime(end,start) << " seconds" << endl;

	double ss2 = totalscore(trips,true);
	cout << "check true total score " << ss2 << ' ' << ss2 - ss << endl;
	int count = 0;
	for(auto &trip:trips) if(trip.size() > 3) count++;
	cout << "total trips remained" << count << endl;

	save_result(trips,datadir+string(argv[8]));
	return 0;
}
