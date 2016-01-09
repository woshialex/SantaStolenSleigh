#include<iostream>
#include<cmath>
#include<ctime>
#include<vector>
#include<cstring>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<cassert>
#include<utility>
#include<algorithm>
#include<iomanip>
#include<thread>
#include<mutex>
using namespace std;

const double PI = 3.14159265358979323846;
const double DEG_TO_RAD = PI/180.0;
const double EARTH_RADIUS = 6371000.0;
const double BASE = 10.0;
const double MAXW = 1000.0;
const double SOUTH = -60*DEG_TO_RAD;
const double SPOLE = -85*DEG_TO_RAD;

static int shift11_count = 0;
const int MAXTRIP = 1700;

mutex mtx[MAXTRIP];//max number of trips
double BETA;
unsigned int rs[32];

int NMB_THREADS;//should be less than 16

class Trip;
inline bool overlap(const Trip & A, const Trip & B, const double delta);

class Gift{
public:
	int id;
	double lat,lon,weight,cumdis,cumweight;
	double x,y,z;
	Gift(int gid, double lat, double lon, double weight){
		this->id = gid;
		this->lat = lat * DEG_TO_RAD;
		this->lon = (lon+180) * DEG_TO_RAD;
		this->weight = weight;
		this->cumdis = 0.0;
		this->cumweight = weight;
		x = 0.5*cos(this->lat)*sin(this->lon);
		y = 0.5*cos(this->lat)*cos(this->lon);
		z = 0.5*sin(this->lat);
	}
	bool operator<(const Gift& other) const{
		return id<other.id;
	}

	double haversine(const Gift & other) const{
		double dx = x - other.x;
		double dy = y - other.y;
		double dz = z - other.z;
		double r = sqrt(dx*dx+dy*dy+dz*dz);
		return 2.0*asin(r);
	}
	void show(){
		cout<<id<<' '<<lat<<' '<<lon<<' '<<weight<<' '<<cumdis<<endl;
	}
};

const Gift NorthPole(0, 90, 0, BASE);

class Trip{
public:
	vector<Gift> gifts;
	long double ws;
	double min_lat;

	double min_lon;
	double max_lon;

	double i_de, i_l1, i_l2; //try_insert auxilary information
	double r_de, r_dis;//try_remove auxilary information

	Trip(){
		ws = 0.0;
		min_lat = 0.0;
	}

	Trip(vector<Gift>::iterator first, vector<Gift>::iterator last){
		gifts.push_back(NorthPole);
		gifts.insert(gifts.end(),first,last);
		gifts.push_back(NorthPole);
		update_info();
	}
	Trip(vector<Gift>::reverse_iterator first, vector<Gift>::reverse_iterator last){
		gifts.push_back(NorthPole);
		gifts.insert(gifts.end(),first,last);
		gifts.push_back(NorthPole);
		update_info();
	}

	void update_info(){
		ws = 0.0;
		for(unsigned int i=1;i<gifts.size();i++){
			gifts[i].cumweight = gifts[i-1].cumweight + gifts[i].weight;
			gifts[i].cumdis = gifts[i-1].cumdis + gifts[i].haversine(gifts[i-1]);
			ws += gifts[i].weight * gifts[i].cumdis;
		}
		update_range();

		if(gifts.back().cumweight > MAXW + 2*BASE){
			cout<<"EXCEED WEIGHT LIMIT!! WARN!! "<<gifts.back().cumweight - 2*BASE<<endl;
		}
	}

	void update_range(){
		if(gifts.size()<3)return;
		min_lon = gifts[1].lon;
		max_lon = gifts[1].lon;
		min_lat = gifts[1].lat;
		for(unsigned int i=2;i<gifts.size()-1;i++){
			if(gifts[i].lat<min_lat)min_lat = gifts[i].lat;
			if(gifts[i].lat<SPOLE)continue;
			double l = gifts[i].lon;
			if(max_lon<PI/2 && l>PI*3/2)l-=2*PI;
			if(min_lon>PI*3.0/2 && l<PI/2)l+=2*PI;
			if(max_lon<l)max_lon=l;
			if(min_lon>l)min_lon=l;
		}
	}

	Gift& operator[](int i){
		return gifts[i];
	}

	int size(){
		return gifts.size();
	}

	Gift& back(){
		return gifts.back();
	}

	void append(Gift x){
		if (!gifts.empty()){
			Gift &last = gifts.back();
			x.cumdis = last.cumdis + x.haversine(last);
			x.cumweight = last.cumweight + x.weight;
			ws += x.cumdis * x.weight;
		}
		if(x.lat<min_lat)min_lat=x.lat;
		gifts.push_back(x);
	}

	void show() const{
		for(auto &g:gifts){
			cout<<g.id<<','<<g.lat<<','<<g.lon<<endl;
		}
	}

	int cuttrip(int ri){
		//some trip may be better to cut to two trips
		//PI/2-lat; is the haversine to the northpole
		if(gifts.size()<=10)return 0;
		int pos = 0;
		double cost = 1.0e9;
		const double totalW = gifts.back().cumweight;
		for(unsigned int i=2;i<gifts.size()-3;i++){
			double w = totalW - gifts[i].cumweight;
			double de = (PI/2-gifts[i].lat) * BASE + (PI/2-gifts[i+1].lat)*(BASE + w) -
				w * gifts[i+1].cumdis;
			if(de<cost){
				cost = de;
				pos = i;
			}
		}
		if(pos==0)return 0;
		if(exp(-BETA*cost)>(double)rand_r(&rs[ri])/RAND_MAX)return pos;
		else return 0;
	}


	int selfoptimize2opt(int ri)
	{

		if(gifts.size() <= 5) return 0;
		int count = 0;
		for(int iter = 0; iter < 5; iter++)
		{

			int i = rand_r(&rs[ri]) % (gifts.size() - 3);
			int j = rand_r(&rs[ri]) % (gifts.size() - 3);
			while(fabs(j - i) < 1.001) j = rand_r(&rs[ri]) % (gifts.size() - 3);
			if(j < i) { int t = i; i = j; j = t; }

			double de = 0.0;

			double CW = gifts.back().cumweight - gifts[i].cumweight;
			double hav1 = gifts[i].haversine(gifts[j]);
			de += CW * hav1;

			CW -= gifts[j].weight;
			for(int l = j; l >= i + 2; l--)
			{
				de += CW * (gifts[l].cumdis - gifts[l - 1].cumdis);
				CW -= gifts[l - 1].weight;
			}

			double hav2 = gifts[i + 1].haversine(gifts[j + 1]);
			de += CW * hav2;

			for(int l = i; l <= j ; l++)
				de -= (gifts.back().cumweight - gifts[l].cumweight) * (gifts[l + 1].cumdis - gifts[l].cumdis);

			if (de < 0 ||  exp(-BETA*de) > (double)rand_r(&rs[ri])/RAND_MAX)
			{

				//cout << " 2opt " << score() << " " << strictscore() << endl;

				// vector of new cum. distances and weights for j, j - 1,..., i + 1
				vector<double> newCumDis, newCumW;

				double hav3 = gifts[i + 1].cumdis - gifts[i].cumdis;
				double hav4 = gifts[j + 1].cumdis - gifts[j].cumdis;

				newCumDis.push_back(gifts[i].cumdis + hav1);
			    newCumW.push_back(gifts[i].cumweight + gifts[j].weight);
				for(int l = j - 1; l >= i + 1; l--)
				{
					newCumDis.push_back(newCumDis.back() + (gifts[l + 1].cumdis - gifts[l].cumdis));
					  newCumW.push_back(newCumW.back() + gifts[l].weight);
				}

				reverse(gifts.begin() + i + 1, gifts.begin() + j + 1);

				for(int l = i + 1; l <= j; l++)
				{
					gifts[l].cumdis = newCumDis[l - (i + 1)];
					gifts[l].cumweight = newCumW[l - (i + 1)];
				}

				for(unsigned int l = j + 1; l < gifts.size(); l++) gifts[l].cumdis += (hav1 + hav2 - hav3 - hav4);

				count++;
				ws += de;
			}

		}
		return count;
	}

	double subroute_cost(const vector<Gift> &subroute, double wleft){
		double cost = 0.0;
		double cumdis = 0.0;
		for(int i=1;i<subroute.size();i++){
			cumdis += subroute[i].haversine(subroute[i-1]);
			cost += subroute[i].weight *cumdis;
		}
		cost += cumdis * wleft;
		return cost;
	}
	double enumerate_exact(vector<Gift>& subroute, double wleft){//begin and end point are fixed
		//return the best solution and cost change
		const double begincost = subroute_cost(subroute,wleft);
		double cost = begincost;
		vector<Gift> tmproute(subroute);
		sort(tmproute.begin()+1,tmproute.end()-1);
		do{
			double c = subroute_cost(tmproute,wleft);
			if(c<cost){//update route
				cost = c;
				subroute.assign(tmproute.begin(),tmproute.end());
			}
		}while(next_permutation(tmproute.begin()+1,tmproute.end()-1));
		return cost-begincost;
	}

	double selfoptimize1move(double delta,int ri){//randomly pick one and move it to the optimum place
		if(gifts.size()<5)return false;
		double totalcost = 0.0;
		int N = gifts.size()-2;
		for(int i=0;i<N/2;i++){//do N/2 times
			int ii = rand_r(&rs[ri])%N+1;
			//maybe slow version but working hopefully
			Trip A(*this);
			A.try_remove(ii);
			A.remove(ii);
			int pos = A.try_insert(gifts[ii],delta);
			if(pos<0)continue;
			double de = A.ws + A.i_de - ws;
			if(de<0){//do the move
				A.insert(gifts[ii],pos);
				totalcost += de;
				*this = A;
			}

		}
		return totalcost;
	}

	double selfoptimizeNexact(int N, int ri){//return improvement amount, N>=3
		if(gifts.size()<N+2)return 0.0;
		int skip = N/2;
		double imp = 0.0;
		for(int i=0;i<gifts.size()-N-2;i+=skip){
			//in order, optimize N point, if N is large, i can skip
			//int start = i;
			//random order
			int start = rand_r(&rs[ri])%(gifts.size()-N-2);
			int end = start+N+2;
			//enumrate all possibilities of the N gifts in between! can be improved with some exact solution method like bound and branch so N can be larger
			vector<Gift> subroute(gifts.begin()+start, gifts.begin()+end);
			double wleft = gifts.back().cumweight - subroute.back().cumweight;
			double de = enumerate_exact(subroute, wleft);
			if(de<0){//update the route information
				imp += de;
				for(int j=1;j<subroute.size()-1;j++){
					gifts[start+j] = subroute[j];
				}
			}
		}
		//update cumdis and cumweight information
		for(int i=1;i<gifts.size();i++){
			gifts[i].cumdis = gifts[i-1].cumdis + gifts[i-1].haversine(gifts[i]);
			gifts[i].cumweight = gifts[i-1].cumweight + gifts[i].weight;
		}
		ws += imp;
		return imp;
	}

	int selfoptimize(int ri){
		if(gifts.size()<=3)return 0;
		int count=0;
		for(int i=0;i<5;i++){
			int pos = rand_r(&rs[ri])%(gifts.size()-3);
			//swap pos+1 and pos + 2 ?
			//0-1-2-3-... -> 0-2-1-3-..
			double d01 = gifts[pos+1].cumdis - gifts[pos].cumdis;
			double d23 = gifts[pos+3].cumdis - gifts[pos+2].cumdis;
			double d12 = gifts[pos+2].cumdis - gifts[pos+1].cumdis;
			double d02 = gifts[pos+2].haversine(gifts[pos]);
			double d13 = gifts[pos+3].haversine(gifts[pos+1]);
			
			double de = (gifts.back().cumweight - gifts[pos+2].cumweight) * (d13+d02-d23-d01);
			de += gifts[pos+1].weight*(d12+d02-d01);
			de -= gifts[pos+2].weight*(d01+d12-d02);

			if (de<0 ||  exp(-BETA*de) > (double)rand_r(&rs[ri])/RAND_MAX){ //accept
				Gift tmp = gifts[pos+1];
				gifts[pos+1] = gifts[pos+2];
				gifts[pos+1].cumweight -= tmp.weight;
				gifts[pos+1].cumdis -= d01+d12-d02;
				tmp.cumweight += gifts[pos+1].weight;
				tmp.cumdis += d02+d12-d01;
				gifts[pos+2] = tmp;
				for(int j=pos+3;j<gifts.size();j++){
					gifts[j].cumdis += d02+d13-d01-d23;
				}
				count++;
				ws += de;
			}
		}
		return count;
	}

	//shiftN0
	bool shiftN0(Trip &other, double delta, double maxdis, int N, int ri){
		//do this only for N>=2, for N<2 there are other things need to be take care of
		//maxdis: dis between 1 and N, only do so if they are very close to each other to reduce calculation
		if(gifts.size()<2+N || other.size()<3)return false;
		double dis_1_N = 0.0;
		int ii = 0;
		for(int loop=0;loop<N*2;loop++){//try many times to find a small group
			ii = rand_r(&rs[ri])%(gifts.size()-1-N)+1;
			dis_1_N = gifts[ii+N-1].cumdis - gifts[ii].cumdis;
			if(dis_1_N<maxdis)break;
		}
		if(dis_1_N>maxdis)return false;

		double weight_1_N = gifts[ii+N-1].cumweight - gifts[ii-1].cumweight;

		if(weight_1_N + other.back().cumweight>MAXW+2*BASE)return false;
		
		//calculate cost reduced if removed
		double ddis = gifts[ii+N].cumdis - gifts[ii-1].cumdis - gifts[ii-1].haversine(gifts[ii+N]);
		double de1 = weight_1_N * gifts[ii].cumdis + (gifts.back().cumweight-gifts[ii+N-1].cumweight)*ddis;
		//calculate cost added by adding to other
		int pos = other.insert_posN(gifts[ii],gifts[ii+N-1],dis_1_N,weight_1_N,delta);

		double l1 = gifts[ii].haversine(other[pos]);
		double ddiso = dis_1_N + l1 + gifts[ii+N-1].haversine(other[pos+1]) - (other[pos+1].cumdis - other[pos].cumdis);
		double de2 = (other.back().cumweight - other[pos].cumweight)*ddiso + weight_1_N*(other[pos].cumdis+l1);

		double de = de2-de1;

		if (de<0 ||  exp(-BETA*de)> (double)rand_r(&rs[ri])/RAND_MAX){ //accept

			//add to the other
			vector<Gift> & ogifts = other.gifts;

			for(int i=pos+1;i<ogifts.size();i++){
				ogifts[i].cumweight += weight_1_N;
				ogifts[i].cumdis += ddiso;
			}
			ogifts.insert(other.gifts.begin()+pos+1,gifts.begin()+ii,gifts.begin()+ii+N);
			
			double changedis = 0.0 - ogifts[pos+1].cumdis + ogifts[pos].cumdis + l1;
			for(int i=pos+1;i<=pos+N;i++){
				ogifts[i].cumweight = ogifts[i-1].cumweight + ogifts[i].weight;
				ogifts[i].cumdis += changedis;
			}
			other.ws += de2;

			//remove from self
			for(int i=ii+N;i<gifts.size();i++){
				gifts[i].cumweight -= weight_1_N;
				gifts[i].cumdis -= ddis;
			}
			gifts.erase(gifts.begin()+ii,gifts.begin()+ii+N);
			ws -= de1;//did not count self energy, but as long as the total is correct, it is fine. It is not difficult to calculate package self energy.!!!

			return true;
		}
		else
			return false;
	}

	//shift10
	bool shift10(Trip &other, double delta, int ri){
		if(gifts.size()<3 || other.size()<3)return false;
		int ii = rand_r(&rs[ri])%(gifts.size()-2)+1;

		Gift &g = gifts[ii];
		if(g.weight + other.back().cumweight>MAXW+2*BASE){
			//return false;
			return shift11(other,delta,ii, ri);
		}
		//calculate cost -reduced if removed
		try_remove(ii);
		//calculate cost added by adding to other
		int pos = other.try_insert(g, delta);
		if(pos<0)return false;
		double de = other.i_de + r_de;

		//too full then reject it out, try may not work
		//double remain = MAXW+2*BASE-gifts.back().cumweight;
		//bool tmp = gifts.back().cumweight - g.weight - other.back().cumweight>0;
		//double factor = (5.0-remain)/2.0*(tmp?1:0.0);
		//if(remain<0)factor=100;
		//if(factor<1)factor=1;
		double factor = 1.0;
		//double small = (300-gifts.back().cumweight>0 && tmp)?0.2:1.0;//try! don't know whether it works
		double small = 1.0;

		if (de<0 ||  exp(-BETA*de/factor*small)> (double)rand_r(&rs[ri])/RAND_MAX){ //accept

			//add to the other
			other.insert(g,pos);
			//remove from self
			this->remove(ii);

			return true;
		}
		else
			return false;

	}

	void try_remove(int pos){
		r_dis = -1.0*(gifts[pos+1].cumdis - gifts[pos-1].cumdis - gifts[pos-1].haversine(gifts[pos+1]));
		r_de = -1.0*gifts[pos].weight * gifts[pos].cumdis + (gifts.back().cumweight-gifts[pos].cumweight)*r_dis;
	}

	void remove(int pos){
		for(int i=pos+1;i<gifts.size();i++){
			gifts[i].cumweight -= gifts[pos].weight;
			gifts[i].cumdis += r_dis;
		}
		gifts.erase(gifts.begin()+pos);
		ws += r_de;
	}

	int insert_posN(const Gift & gb, const Gift &ge, double selfdis, double selfweight, double delta=0.08){
		int N = gifts.size();
		double mincost = 1000.0;
		int pos = 0;//will insert after pos
		const double latplus=gb.lat+delta;
		const double latminus = gb.lat -delta;
		for(int i=0;i<N-1;i++){
			double lat = gifts[i].lat;
			double nlat = gifts[i+1].lat;
			if((lat < latplus && lat>latminus) ||
					(lat>latplus && nlat <latplus) ||
					(lat<latminus && nlat>latminus)){
				double l1 = gb.haversine(gifts[i]);
				double ddiso = selfdis + l1 + ge.haversine(gifts[i+1]) - (gifts[i+1].cumdis - gifts[i].cumdis);
				double c = (gifts.back().cumweight - gifts[i].cumweight)*ddiso + selfweight*(gifts[i].cumdis+l1);
				if(c<mincost){
					mincost = c;
					pos = i;
				}
			}
		}
		return pos;
	}

	void insert(Gift &g, int pos, bool calculated = true){
		//i_l1,i_l2,i_de shared variable, set when do try_insert

		if(!calculated){
			i_l1 = g.haversine(gifts[pos]);
			i_l2 = g.haversine(gifts[pos+1]);
			double ddis = i_l1 + i_l2 - (gifts[pos+1].cumdis - gifts[pos].cumdis);
			i_de = (gifts.back().cumweight - gifts[pos].cumweight)*ddis + g.weight*(gifts[pos].cumdis+i_l1);
		}

		g.cumdis = gifts[pos].cumdis + i_l1;
		g.cumweight = gifts[pos].cumweight + g.weight;
		double ddis = i_l1+i_l2-(gifts[pos+1].cumdis-gifts[pos].cumdis);
		for(int i=pos+1;i<gifts.size();i++){
			gifts[i].cumweight += g.weight;
			gifts[i].cumdis += ddis;
		}
		gifts.insert(gifts.begin()+pos+1,g);
		ws += i_de;
	}


	int try_insert(const Gift & g, double delta){
		int N = gifts.size();
		double mincost = 1.0e9;
		int pos = -1;//will insert after pos
		const double latplus=g.lat+delta;
		const double latminus = g.lat -delta;
		for(int i=0;i<N-1;i++){
			double lat = gifts[i].lat;
			double nlat = gifts[i+1].lat;
			if((lat < latplus && lat>latminus) ||
					(lat>latplus && nlat <latplus) ||
					(lat<latminus && nlat>latminus)){

				double l1 = g.haversine(gifts[i]);
				double l2 = g.haversine(gifts[i+1]);
				double ddis = l1 + l2 - (gifts[i+1].cumdis - gifts[i].cumdis);
				double de = (gifts.back().cumweight - gifts[i].cumweight)*ddis + g.weight*(gifts[i].cumdis+l1);
				if(de<mincost){//save information for insert!!
					mincost = de;
					pos = i;
					this->i_de = de;
					this->i_l1 = l1;
					this->i_l2 = l2;
				}
			}
		}
		return pos;
	}

	double score(){
		return ws;
	}

	double strictscore(){
		double ss = 0.0;
		double w = gifts.back().weight;
		for(unsigned i = gifts.size()-1; i-- >0;){
			ss += w * gifts[i].haversine(gifts[i+1]);
			w += gifts[i].weight;
		}
		return ss;
	}

	//---------- shift11 -----------------
	bool shift11(Trip &other, double delta, int ii, int ri){
		if(gifts.size()<3 || other.size()<3)return false;

		int k=0, jj=0;
		for(k=0;k<5;k++){//try at most 5 times
			jj = rand_r(&rs[ri]) % (other.size() -2) +1;
			if((other[jj].weight - gifts[ii].weight + gifts.back().cumweight <= MAXW + 2 * BASE)&&(gifts[ii].weight - other[jj].weight + other.back().cumweight <= MAXW + 2 * BASE)) break;
		}
		if(k==5)return false;

		Trip A(*this);
		Trip B(other);

		A.try_remove(ii);
		A.remove(ii);
		B.try_remove(jj);
		B.remove(jj);
		int jj_pos = A.try_insert(other[jj],delta);
		if(jj_pos<0)return false;
		int ii_pos = B.try_insert(gifts[ii],delta);
		if(ii_pos<0)return false;
		//need to insert first then do remove, or else the variables are ineffective

		double de = A.r_de + A.i_de + B.r_de + B.i_de;
		if (de<0 ||  exp(-BETA*de)> (double)rand_r(&rs[ri])/RAND_MAX){ //accept
			A.insert(other[jj], jj_pos);
			B.insert(gifts[ii], ii_pos);

			*this = A;
			other = B;
			//assert(ws == A.ws && gifts[0].id == A[0].id);
			shift11_count++;
			return true;
		}
		return false;
	}

	//-------------------------- SWAP -----------------------------------
	bool updateSWAP(Trip &other,  double delta, int ri)
	{

		if(gifts.size() < 3 || other.size() < 3) return false;
		int ii = rand_r(&rs[ri]) % (gifts.size() - 2) + 1;

		// evaluate swap move (cost and feasibility)
		// find a position on other route to swap with g
		// both, cost and feasibility have to be checked
		double de;
		int pos = other.swap_pos(*this, ii, de, delta);

		if(pos < 0) return false;

		if (de < 0 ||  exp(-BETA * de) > (double) rand_r(&rs[ri]) / RAND_MAX)
		{
			Gift tmpG = gifts[ii];
			performSWAP(ii, other[pos]);
			other.performSWAP(pos, tmpG);
			return true;
		}
		else return false;

	}

	// swap this[ii] with other[jj]
	void performSWAP(int pos, Gift g)
	{
		double l1 = g.haversine(gifts[pos - 1]);
		double l2 = g.haversine(gifts[pos + 1]);
		double l3 = gifts[pos].cumdis - gifts[pos - 1].cumdis;
		double l4 = gifts[pos + 1].cumdis - gifts[pos].cumdis;
		double ddiso = l1 + l2 - (l3 + l4);
		double de = (gifts.back().cumweight - gifts[pos - 1].cumweight) * ddiso +
				    (g.weight - gifts[pos].weight) * (gifts[pos - 1].cumdis + l1) -
				    gifts[pos].weight * (l2 - l4);

		for(int i = pos + 1; i < gifts.size(); i++)
		{
			gifts[i].cumweight += (g.weight - gifts[pos].weight);
			gifts[i].cumdis += ddiso;
		}

		gifts[pos] = g;
		ws += de;

		gifts[pos].cumdis = gifts[pos - 1].cumdis + l1;
		gifts[pos].cumweight = gifts[pos - 1].cumweight + g.weight;
	}


	// cost difference with swapping gift at position "pos" with gift g
	double swap_cost(const Gift &g, int pos)
	{
		double l1 = g.haversine(gifts[pos - 1]);
		double l2 = g.haversine(gifts[pos + 1]);
		double l3 = gifts[pos].cumdis - gifts[pos - 1].cumdis;
		double l4 = gifts[pos + 1].cumdis - gifts[pos].cumdis;
		double ddiso = l1 + l2 - (l3 + l4);
		double de = (gifts.back().cumweight - gifts[pos].cumweight)*ddiso + 
			g.weight * (gifts[pos-1].cumdis + l1) -
			gifts[pos].weight * gifts[pos].cumdis;
		return de;
	}


	int swap_pos(Trip& other, int ii, double &diff, double delta = 0.08)
	{
		int N = gifts.size();
		double mincost = 1.0e12;
		int pos = 0;//will insert after pos
		const double latplus=other[ii].lat+delta;
		const double latminus = other[ii].lat - delta;
		for(int i=1;i<N-1;i++){
			double lat = gifts[i].lat;
			if((lat < latplus && lat>latminus)){
				if(other[ii].weight - gifts[i].weight + gifts.back().cumweight > MAXW + 2 * BASE) continue;
				if(gifts[i].weight - other[ii].weight + other.back().cumweight > MAXW + 2 * BASE) continue;
				double c = swap_cost(other[ii], i) + other.swap_cost(gifts[i], ii);
				if(c<mincost){
					mincost = c;
					pos = i;
				}
			}
		}
		diff = mincost;
		if(diff > 999) return -1;
		return pos;
	}
	//------------------------------------------------------------------------------------
	

	double swap_tail(Trip &other, int ri){
		if(gifts.size()<3 || other.size()<3)return 0.0;
		if(!overlap(*this,other,0.0))return false;
		if(((min_lat<SOUTH) != (other.min_lat<SOUTH)) && (rand_r(&rs[ri])%5>0))return 0.0;
		//when they are separated at -60, try less frequently

		double cost = 1.0e13;
		int ii=-1,jj=0;
		//swap the tail part starting at i+1 & j+1
		for(int i=0;i<gifts.size()-1;i++){
			double w1 = gifts.back().cumweight - gifts[i].cumweight;
			//double l1 = gifts[i+1].cumdis - gifts[i].cumdis;
			for(int j=0;j<other.size()-1;j++){
				if((i==0 && j==0) || (i==gifts.size()-2 && j==other.size()-2))continue;
				double w2 = other.back().cumweight - other[j].cumweight;
				//need to satisfy weight constraint
				if(gifts[i].cumweight + w2 >MAXW + 2*BASE || other[j].cumweight + w1>MAXW+2*BASE)continue;

				double de = w1 * (gifts[i+1].haversine(other[j]) + other[j].cumdis - gifts[i+1].cumdis) +
					w2 * (gifts[i].haversine(other[j+1]) + gifts[i].cumdis - other[j+1].cumdis);
				if(de<cost){
					cost = de;
					ii = i;
					jj = j;
				}
			}
		}
		if(ii<0)return 0.0;

		if(exp(-cost*BETA)>(double)rand_r(&rs[ri])/RAND_MAX){//do the tail-swap
			vector<Gift> selftail(gifts.begin()+ii+1,gifts.end());
			gifts.erase(gifts.begin()+ii+1,gifts.end());
			vector<Gift> othertail(other.gifts.begin()+jj+1,other.gifts.end());
			other.gifts.erase(other.gifts.begin()+jj+1,other.gifts.end());
			gifts.insert(gifts.end(),othertail.begin(),othertail.end());
			other.gifts.insert(other.gifts.end(),selftail.begin(),selftail.end());

			this->update_info();
			other.update_info();
			return cost;
		}
		else return 0.0;
	}
};


inline bool overlap(const Trip & A, const Trip & B, const double delta){
	if(( A.min_lon <= B.max_lon + delta && B.min_lon<= A.max_lon +delta) || ( A.min_lat <SPOLE && B.min_lat <SPOLE ))return true;
	if(B.min_lon>PI*1.5 && A.max_lon<PI*0.5)return (A.min_lon+2*PI<=B.max_lon+delta && B.min_lon<=A.max_lon+2*PI+delta);
	if(A.min_lon>PI*1.5 && B.max_lon<PI*0.5)return (B.min_lon+2*PI<=A.max_lon+delta && A.min_lon<=B.max_lon+2*PI+delta);
	return false;
}



vector<Trip> load_data(string filename){
	ifstream infile;
	infile.open(filename.c_str());
	string line;
	int gid,tid;
	double lat,lon,weight;
	char c;
	vector<Trip> data;
	int prev_id = -1234;
	Gift start(0, 90, 0, BASE);
	if(infile.is_open()){
		getline(infile,line);//skip header
		int count = 0;
		while(getline(infile,line)){
			count++;
			stringstream x(line);
			x>>gid>>c>>lat>>c>>lon>>c>>weight>>c>>tid;
			if(false){
				cout<<gid<<','<<lat<<','<<lon<<','<<weight<<','<<tid<<endl;
				if(count>10)break;
			}
			if ( tid != prev_id){
				if ( !data.empty() ){
					data.back().append(start);//end of the trip
				}
				prev_id = tid;
				data.push_back(Trip());
				data.back().append(start); //start of the trip
			}
			Trip& last = data.back();
			Gift g(gid,lat,lon,weight);
			last.append(g);
		}
		data.back().append(start);//end of the trip
		infile.close();
	}
	else cerr<<"unable to open input file "<<filename<<endl;

	return data;
}

double totalscore(vector<Trip> trips, bool strict=false){
	double ss = 0.0;
	if(strict){
		for(auto &i:trips){
			ss += i.strictscore();
		}

	}else{
		for(auto &i:trips){
			ss += i.score();
		}
	}
	return ss * EARTH_RADIUS;
}

void save_result(vector<Trip> trips, string filename){
	ofstream file;
	file.open(filename.c_str());
	file<<"GiftId,Latitude,Longitude,Weight,TripId"<<endl;
	int c=0;
	for(auto &trip:trips){
		c += 1;
		for(int i=1;i<trip.size()-1;i++){
			Gift &x = trip[i];
			file<<x.id<<','<<x.lat/DEG_TO_RAD<<','<<x.lon/DEG_TO_RAD-180<<','<<x.weight<<','<<c<<endl;
		}
	}
	file.close();
}

double selfupdateExact(vector<Trip> &trips, int Nex, int ri){
	clock_t start = clock();
	if(Nex<2){
		cerr<<"exact solution looking for at least 2 gifts"<<endl;
		return 0.0;
	}
	if(Nex>8){
		cout<<"WARNING, N is too large N! will be infinite"<<endl;
	}
	double imp = 0.0;
	int N = trips.size();
	int s = int(N/NMB_THREADS)+1;
	int e = min(ri*s+s,N);
	for(int i=ri*s;i<e;i++){
		mtx[i].lock();
		imp += trips[i].selfoptimizeNexact(Nex, ri);
		mtx[i].unlock();
	}
	imp *= EARTH_RADIUS;
	clock_t end = clock();
	//cout<<"selfupdate Exact "<<Nex<<" point improved :"<<imp<<", takes "<<setw(4)<<(end-start)/(double)(CLOCKS_PER_SEC)<<" second"<<endl;
	return imp;
}

double selfupdateCut(vector<Trip> &trips, int ri){
	clock_t starttime = clock();
	int N = trips.size();
	Gift start(0, 90, 0, BASE);
	double de = 0.0;
	//return de;
	for(int i=0;i<N;i++){
		int pos = trips[i].cuttrip(ri);
		if(pos==0)continue;
		else cout<<"cut! "<<i<< "at "<<pos<<endl;
		double orig_score = trips[i].ws;
		//cut the trip after pos
		trips.push_back(Trip());
		vector<Gift> &nt = trips.back().gifts;
		nt.push_back(start);
		vector<Gift> &t = trips[i].gifts;
		trips.back().gifts.insert(trips.back().gifts.end(),t.begin()+pos+1,t.end());
		t.erase(t.begin()+pos+1,t.end()-1);

		//update information for both trips
		trips.back().update_info();
		trips[i].update_info();

		de += trips[i].ws + trips.back().ws - orig_score;
	}
	de *= EARTH_RADIUS;
	//cout<<"Try cutting trips"<<endl;
	clock_t end = clock();
	//cout<<"selfupdate Cut trips : "<<de<<", takes "<<setw(4)<<(end-starttime)/(double)(CLOCKS_PER_SEC)<<" second"<<endl;
	cout<<"New # of trips "<<trips.size()<<endl;
	return de;
}

double swaptails(vector<Trip> &trips, const vector<pair<int,int> > &trip_pair, int ri){
	clock_t start = clock();

	double cost = 0.0;
	int N = trip_pair.size();
	for(int k=0;k<(N/NMB_THREADS);k++){
		int x = rand_r(&rs[ri])%N;//random pairs
		int i = trip_pair[x].first;
		int j = trip_pair[x].second;
		if(i<j){
			if(!mtx[i].try_lock())continue;
			if(!mtx[j].try_lock()){mtx[i].unlock();continue;}
		}else{
			if(!mtx[j].try_lock())continue;
			if(!mtx[i].try_lock()){mtx[j].unlock();continue;}
		}
		cost += trips[i].swap_tail(trips[j], ri);
		mtx[j].unlock();
		mtx[i].unlock();
	}
	cost *= EARTH_RADIUS;
	clock_t end = clock();
	//cout<<"swaptails improved : "<<cost<<", takes "<<setw(4)<<(end-start)/(double)(CLOCKS_PER_SEC)<<" seconds"<<endl;
	return cost;
}

int selfupdate(vector<Trip> &trips, int ri){
	clock_t start = clock();
	int suc = 0;
	int N = trips.size();
	int s = int(N/NMB_THREADS)+1;
	int e = min(ri*s+s,N);
	for(int i=ri*s;i<e;i++){
		mtx[i].lock();
		suc += trips[i].selfoptimize(ri);
		mtx[i].unlock();
//		suc += trips[i].selfoptimize2opt(ri); //it has a bug need to fix
	}
	clock_t end = clock();
	//cout<<"self update 2opt times : "<< suc <<", takes"<<setw(4)<<(end-start)/(double)(CLOCKS_PER_SEC)<<" seconds"<<endl;
	return suc;
}

double selfupdate1move(vector<Trip> &trips, double delta, int ri){
	clock_t start = clock();
	double cost = 0.0;
	int N = trips.size();
	int s = int(N/NMB_THREADS)+1;
	int e = min(ri*s+s,N);
	for(int i=ri*s;i<e;i++){
		mtx[i].lock();
		cost += trips[i].selfoptimize1move(delta, ri);
		mtx[i].unlock();
	}
	cost *= EARTH_RADIUS;
	clock_t end = clock();
	//cout<<"selfupdate 1move improved :"<<cost<<", takes "<<setw(4)<<(end-start)/(double)(CLOCKS_PER_SEC)<<" seconds"<<endl;
	return cost;
}


vector<int> local_search_between_routes(vector<Trip> & trips, const vector<pair<int,int> > &trip_pair, double lat_range, const vector<double> &prob, int ri)
{
	vector<int> suc(prob.size(),0);
	int N = trip_pair.size();
	for(int k=0;k<(N/NMB_THREADS);k++){
		int x = rand_r(&rs[ri])%N;//random pairs
		int i = trip_pair[x].first;
		int j = trip_pair[x].second;
		if(i<j){
			if(!mtx[i].try_lock())continue;
			if(!mtx[j].try_lock()){mtx[i].unlock();continue;}
		}else{
			if(!mtx[j].try_lock())continue;
			if(!mtx[i].try_lock()){mtx[j].unlock();continue;}
		}
		if(trips[i].shift10(trips[j], lat_range,ri))suc[0]++;
		for(int loopN = 1; loopN<prob.size()-1; loopN++){//loopN -> loopN-1 shifts
			if((double)rand_r(&rs[ri])/RAND_MAX<prob[loopN])
				if(trips[i].shiftN0(trips[j], lat_range, 0.02*sqrt(loopN), loopN,ri))suc[loopN]++;
		}
		if((double)rand_r(&rs[ri])/RAND_MAX<prob.back())if(trips[i].updateSWAP(trips[j], lat_range, ri))suc[prob.size()-1]++;
		mtx[j].unlock();
		mtx[i].unlock();
	}
	return suc;
}


int set_working_pairs(vector<pair<int,int> > & trip_pair, vector<Trip> &trips, const bool sep, const double lon_range){
	trip_pair.clear();
	int N = trips.size();
	if(!sep){
		for(int i=0;i<N;i++){
			if(trips[i].size()<3)continue;
			for(int j=0;j<N;j++){
				if(j==i)continue;
				if(trips[j].size()<3)continue;

				if(overlap(trips[i], trips[j],lon_range)) {
					trip_pair.push_back(make_pair(i,j));
				}
			}
		}
	}else{
		for(int i=0;i<N;i++){
			if(trips[i].size()<3)continue;
			bool mlat = (trips[i].min_lat>SOUTH);
			for(int j=0;j<N;j++){
				if(j==i)continue;
				if(trips[j].size()<3)continue;
				//seprate to two regions to optimize
				bool jlat = (trips[j].min_lat>SOUTH);
				if(mlat!=jlat)continue;

				if(overlap(trips[i],trips[j],lon_range)) {
					trip_pair.push_back(make_pair(i,j));
				}
			}
		}

	}
	return trip_pair.size();
}
