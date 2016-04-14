#include<iostream>
#include<cmath>
#include <map>
#include <utility>
#include<cstdlib>
#include<tclap/CmdLine.h>
#include<boost/random/uniform_int.hpp>
#include<boost/random/uniform_real.hpp>
#include<boost/random/poisson_distribution.hpp>
#include<boost/random/exponential_distribution.hpp>
#include<boost/random/variate_generator.hpp>
#include<boost/random/mersenne_twister.hpp>
#include<ctime>
#include<string>
#include<algorithm>

using namespace std;

typedef boost::mt19937 prgType;

#define BIGINT 9999999

//This structure contains all the information about a branch in the coalescent tree
struct Branch_t {
    int num;
    int current;
    int external;
    int left;
    int upper_left;
    int upper_right;
    int right;
    int population;
    double total_time_to_event;
    double length;
    vector<double> variant_locs;
    vector<double> seg_starts;
    vector<double> seg_ends;
    vector<double> migration_times;
};

//This is the function which "shifts" time for exponential growth
double exponentialGrowth(double beta,double total_time,double basic_coal_time) {
    if (1+beta*basic_coal_time*exp(-beta*total_time) <=0) {
        return -99;
    }
    else {
    double etimeshift = (1/beta)*log(1+beta*basic_coal_time*exp(-beta*total_time));
    return etimeshift;
    }
}

//This is the function which "shifts" time for faster than exponential growth
double fteGrowth(double bigN,double alpha,double ttot,double basic_coal_time,double fte_exp) {
    double ftetimeshift = (pow((basic_coal_time*pow(2*bigN,(fte_exp-1))*alpha*fte_exp)+pow((1+pow(2*bigN,(fte_exp-1))*(fte_exp-1)*alpha*ttot),(fte_exp/(fte_exp-1))),(fte_exp-1)/fte_exp)-(1+pow(2*bigN,(fte_exp-1))*(fte_exp-1)*alpha*ttot))/(pow(2*bigN,(fte_exp-1))*alpha*(fte_exp-1));
    return ftetimeshift;
}

//This function generates the exponential time to the next coalescent event and then does whatever "shifting" is necessary for the current growth state
double timeToCoal(int num_rem,int grow,double beta,double bigN,prgType& rng,double ttot,double new_beta,double fte_exp, double lambda) {
        double standard_rate = (1/lambda)*(num_rem*(num_rem-1.0))/2;
        if (standard_rate == 0) {
                return 99999;
        }
        boost::exponential_distribution<> exp_dist(standard_rate);
        boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
        double x = getexp(); 
        if (grow == -1) {
            return x;
        }
        if (grow == 0) {
            double new_x = new_beta*x;
            return new_x;
        }
        if (grow == 1) { 
            double value = exponentialGrowth(beta,ttot,x);
            return value;
        }
        if (grow == 2) {
            double value = fteGrowth(bigN,beta,ttot,x,fte_exp);
            return value;
        }
        cout << "Growth model not recognized, please enter new -g value\n";
        exit(219); 
}

//Generating the exponential time to the next migration event
double timeToMig(double pop_size,prgType& rng,int pop,double M) {
        double standard_rate = (M*pop_size)/2;
        if (standard_rate == 0) {
                return BIGINT;
        }
        boost::exponential_distribution<> exp_dist(standard_rate);
        boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
        double x = getexp();
        return x;
}

//Generating the exponential time to next recombination event     
double timeToRecomb(double total_seglength,double rec_rate,prgType& rng,int num_rem) {
    double recomb_rate = (rec_rate/2)*total_seglength;
     //   cout << "Recombination Rate error, rec_rate = " << rec_rate << " total_seglength = " << total_seglength << " Num_Rem = " << num_rem <<endl;
        
    if (recomb_rate>0) {
        boost::exponential_distribution<> exp_dist(recomb_rate);
        boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
        double x = getexp();
        return x;
    }
    if (recomb_rate==0) {
        return BIGINT;
    }
} 

//After growth is complete or when an instantaneous growth event occurs this function changes the population size scaling parameter to reflect current populations size
double getNewPopSize(int grow,double time_grow,double beta,double bigN,double old_beta, double fte_exp, double pop_size) {
    if (grow==-1) {
        return pop_size;
    }
    if (grow==0) {
        return old_beta;
    }
    if (grow==1) {
       double new_beta = exp(-beta*time_grow);
       return new_beta;
   }
   if (grow==2) { 
        double new_beta = 1/(pow(1+pow(2*bigN,(fte_exp-1))*(fte_exp-1)*beta*time_grow,(1/(fte_exp-1))));
        return new_beta;
    }
    if (grow==9) {
        double new_beta = beta;
        return new_beta;
    }
}

//After coalescence (or recombination although it should be redundant there) this recalculates the remaining total sequence length in the current active branches
double getLength(vector<Branch_t> & Tree,int line) {
    if (Tree[line].left == -1 && Tree[line].right == -1) {
        return Tree[line].length;
    }
    if (Tree[line].left != -1 ) {
        double more_length = getLength(Tree,Tree[line].left);
        return Tree[line].length + more_length;
    }
    else {
        double more_length = getLength(Tree,Tree[line].right);
        return Tree[line].length + more_length;
    }
}

//This function looks at the branches above the current branch and grabs the mutations that carryover onto the current branch
void propogateMutation(vector<Branch_t> & Current_Tree,int line) {
    int sizzler = Current_Tree.size();
    int s = Current_Tree[line].seg_starts.size();
    int up_left = Current_Tree[line].upper_left;
    int up_right = Current_Tree[line].upper_right;
    if(up_left != -1) {
        int up_vars = Current_Tree[up_left].variant_locs.size();
        int exist_vars = Current_Tree[line].variant_locs.size();
        for (int j=0;j<up_vars;++j) {
            int seen = 0;
            for (int k=0;k<exist_vars;++k) {
                if (Current_Tree[line].variant_locs[k]==Current_Tree[up_left].variant_locs[j]) {
                    seen++;
                    break;
                }
            }
            if (seen==0) {
                for (int h=0;h<s;++h) {
                    if (Current_Tree[line].seg_starts[h] <= Current_Tree[up_left].variant_locs[j] && Current_Tree[line].seg_ends[h] >= Current_Tree[up_left].variant_locs[j]) {
                        Current_Tree[line].variant_locs.push_back(Current_Tree[up_left].variant_locs[j]);
                        break;
                    }
                }
            }
        }
    }
    if(up_right != -1) {
        int up_vars = Current_Tree[up_right].variant_locs.size();
        int exist_vars = Current_Tree[line].variant_locs.size();
        for (int j=0;j<up_vars;++j) {
            int seen = 0;
            for (int k=0;k<exist_vars;++k) {
                if (Current_Tree[line].variant_locs[k]==Current_Tree[up_right].variant_locs[j]) {
                    seen++;
                    break;
                }
            }
            if (seen==0) {
                for (int h=0;h<s;++h) {
                    if (Current_Tree[line].seg_starts[h] <= Current_Tree[up_right].variant_locs[j] && Current_Tree[line].seg_ends[h] >= Current_Tree[up_right].variant_locs[j]) {
                        Current_Tree[line].variant_locs.push_back(Current_Tree[up_right].variant_locs[j]);
                        break;
                    }
                }
            }
        }
    }
}

//Going branch by branch from the final branch and generating the number of mutations on each branch, their location, and then calling the propogate function to grab mutations higher on the tree
void populateMutations(vector<Branch_t> & Current_Tree,prgType& rng,vector<double> & Variant_Pos,double theta) {
    double n = Current_Tree.size();
    for (int i=Current_Tree.size()-1;i>-1;--i) {
        double length = Current_Tree[i].length;        
        double mut_rate = (theta/2)*length;
        int num_muts=0;
        boost::uniform_real<> cl(0,1); 
        boost::variate_generator<prgType&, boost::uniform_real<> > getunif(rng,cl);
        if (mut_rate !=0) {
            boost::poisson_distribution<> pois_dist(mut_rate);
            boost::variate_generator<prgType&, boost::poisson_distribution<> > getpois(rng,pois_dist);
            num_muts = getpois();
      //      cout << "Branch: " << i;
        }
        for (int j=0;j<num_muts;++j) {
            double where = getunif();
            int covered = 0;
            int s = Current_Tree[i].seg_starts.size();
            for (int h=0;h<s;++h) {
                if (where >= Current_Tree[i].seg_starts[h] && where <= Current_Tree[i].seg_ends[h]) {
                    ++covered;
                }
            }
            if (covered==1) {
            
                Current_Tree[i].variant_locs.push_back(where);
                Variant_Pos.push_back(where);
            }
        }
        propogateMutation(Current_Tree,i);
    }
} 

//This function carries out a recombination event once it has been decided that a recombination is the next event
int recombinationEvent(vector<Branch_t> & Current_Tree,double which,double & total_seglength,double tevent,prgType& rng,int* count,double ttot,vector<double> & Recomb_Pos) {
    double tracker=0;
    int recomb_line = -1;
    double recomb_line_length = 0;
    int n = Current_Tree.size();
    for (int i=0;i<n;++i) {
        double min_prob = -1;
        double prob_length = 0;
        double max_prob = -1;
        //Choosing which branch has the recombination
        if (Current_Tree[i].current==1) {
            int h = Current_Tree[i].seg_starts.size();
            for (int j=0;j<h;++j) {
                prob_length = prob_length + Current_Tree[i].seg_ends[j]-Current_Tree[i].seg_starts[j];
            }
            
            double prob_length_div = prob_length/total_seglength;
            min_prob = tracker;
            max_prob = tracker + prob_length_div;
            if (which > min_prob && which <= max_prob) {
                recomb_line = i;
                recomb_line_length = prob_length;
                break;
   
            }
            else {
                tracker = max_prob;
            }
        }
    }
    
    //Making the branch with the recombination not current and selecting recombination location
    Branch_t & recomb_branch = Current_Tree[recomb_line];
    int current_pop = recomb_branch.population;

    recomb_branch.length = tevent + ttot - getLength(Current_Tree,recomb_line);
    recomb_branch.total_time_to_event = tevent + ttot;
    recomb_branch.current = 0;
    boost::uniform_real<> cl(0,recomb_line_length);
    boost::variate_generator<prgType&, boost::uniform_real<> > getunif2(rng,cl);
    double where = getunif2();
    int s = recomb_branch.seg_starts.size();
    
    double running_tot = 0;
    double rec_loc = -1;
    int recomb_seg = -1;
    for (int i=0;i<s;++i) {
        double seg_length = recomb_branch.seg_ends[i]-recomb_branch.seg_starts[i];
        if (where < (running_tot + seg_length) && where > running_tot)  {
            rec_loc = where - running_tot;
            recomb_seg = i;
            break;
        }
        else {
            running_tot = running_tot + seg_length;
        }
     }
     double new_break = recomb_branch.seg_starts[recomb_seg] + rec_loc;
     Recomb_Pos.push_back(new_break);
     //Creating first new branch from recombination  
     Branch_t new_branch1;
     new_branch1.num = *count;
     ++*count;
     new_branch1.current = 1;
     new_branch1.length = 0;
     new_branch1.population = current_pop;
     new_branch1.total_time_to_event = 0;
     new_branch1.external = 0;
     new_branch1.left = -1;
     new_branch1.upper_left = -1;
     new_branch1.upper_right = -1;
     new_branch1.right = recomb_branch.num;  

     if (recomb_seg>0) {
        for (int i=0;i<recomb_seg;++i) {
            new_branch1.seg_starts.push_back(recomb_branch.seg_starts[i]);
            new_branch1.seg_ends.push_back(recomb_branch.seg_ends[i]);
        }
        new_branch1.seg_starts.push_back(recomb_branch.seg_starts[recomb_seg]);
        new_branch1.seg_ends.push_back(new_break);
     }

     if (recomb_seg == 0 ) {
        new_branch1.seg_starts.push_back(recomb_branch.seg_starts[recomb_seg]);
        new_branch1.seg_ends.push_back(new_break);
     }
     
     //Making recombining branch point to new branch1
     recomb_branch.upper_left = new_branch1.num;
     
     //Creating Second new branch from recombination
     Branch_t new_branch2;
     new_branch2.num = *count;
     ++*count;
     new_branch2.current = 1;
     new_branch2.length = 0;
     new_branch2.population = current_pop;
     new_branch2.total_time_to_event = 0;
     new_branch2.external = 0;
     new_branch2.upper_left = -1;
     new_branch2.upper_right = -1;
     new_branch2.left = recomb_branch.num;
     new_branch2.right = -1;  
     if (recomb_seg < s-1) {
         new_branch2.seg_starts.push_back(new_break);
         new_branch2.seg_ends.push_back(recomb_branch.seg_ends[recomb_seg]);
         for (int i=recomb_seg+1;i<s;++i) {
             new_branch2.seg_starts.push_back(recomb_branch.seg_starts[i]);
             new_branch2.seg_ends.push_back(recomb_branch.seg_ends[i]);
         }
     }
     
     if (recomb_seg==s-1) {
         new_branch2.seg_starts.push_back(new_break);
        //double final_end = recomb_branch.seg_ends[0];
        double check_end = recomb_branch.seg_ends[recomb_seg];
         new_branch2.seg_ends.push_back(check_end);
     }
     
    //Making recombining branch point to newbranch2
    recomb_branch.upper_right = new_branch2.num; 
     
    Current_Tree.push_back(new_branch1);
    Current_Tree.push_back(new_branch2);
    
    int new_n = Current_Tree.size();
    total_seglength=0;
    for (int j=0;j<new_n;++j) {
        if (Current_Tree[j].current==1) {
            int blocks = Current_Tree[j].seg_starts.size();
            for (int h=0;h<blocks;++h) {
                total_seglength = total_seglength + (Current_Tree[j].seg_ends[h]-Current_Tree[j].seg_starts[h]);
            }
        }
    }
    return current_pop;
}

//Carrying out a Coalescent Event    
void performCoal(vector<Branch_t> & Current_Tree,int line1,int line2,int pop,int* count,double tevent,double ttot,double & total_seglength) { 

    //Creating the new current branch            
     Branch_t newcoalbranch;
    newcoalbranch.num = *count;
    ++*count;
    newcoalbranch.current =1;
    newcoalbranch.external=0;
    newcoalbranch.length=0;
    newcoalbranch.upper_left = -1;
    newcoalbranch.upper_right = -1;
    newcoalbranch.population = pop;
    newcoalbranch.total_time_to_event =0;
    
    Current_Tree[line1].upper_left = newcoalbranch.num;
    Current_Tree[line2].upper_right = newcoalbranch.num;
    Current_Tree[line1].current=0;
    Current_Tree[line2].current=0;
     
    newcoalbranch.left = Current_Tree[line1].num;
    newcoalbranch.right = Current_Tree[line2].num;

    Current_Tree[line1].length = ttot + tevent - getLength(Current_Tree,line1);
    Current_Tree[line1].total_time_to_event = ttot + tevent;
    Current_Tree[line2].length = ttot + tevent - getLength(Current_Tree,line2);
    Current_Tree[line2].total_time_to_event = ttot + tevent;

    //Merging the different segment portions in the new branch
    //First create a pair of vectors with all the segment starting and ending positions
    int s1 = Current_Tree[line1].seg_starts.size();
    int s2 = Current_Tree[line2].seg_starts.size();
    
    vector<double> all_seg_starts;
    vector<double> all_seg_ends;
    
    for (int i=0;i<s1;++i) {
        all_seg_starts.push_back(Current_Tree[line1].seg_starts[i]);
        all_seg_ends.push_back(Current_Tree[line1].seg_ends[i]);
    }
    for (int i=0;i<s2;++i) {
        all_seg_starts.push_back(Current_Tree[line2].seg_starts[i]);
        all_seg_ends.push_back(Current_Tree[line2].seg_ends[i]);
    }

    int all_size = all_seg_starts.size();
 //   for (int i=0;i<all_size;++i) {
 //       if (all_seg_ends[i] <= all_seg_starts[i]) {
      //      cout << "Start: " << all_seg_starts[i] << " End: " << all_seg_ends[i] << endl;
 //       }
//    }
        
    //Next create a multimap structure with each segment start and end as a key paired with the same integer value
    multimap<double, int> merge_map;
    
    for(int i=0;i<all_size;++i) {
        merge_map.insert(pair<double, int>(all_seg_starts[i],i));
        merge_map.insert(pair<double, int>(all_seg_ends[i],i));
    }

    int temp_count = -1;
    double first = -1;
    double last = -1;
    vector<int> values;
    int seen=0;
    
    //Loop through the multimap, keeping track of overlapping segments in order to correctly merge them
    for (multimap<double, int>::iterator it = merge_map.begin();it != merge_map.end();++it)  {
        if (temp_count==-1) {
            first = (*it).first;
            values.push_back((*it).second);
            seen++;
            temp_count++;
        }
        if (seen==0) {
            if ((*it).first==last) {
                values.push_back((*it).second);
                seen++;
                temp_count=0;
            }
            if ((*it).first!=last) {
                newcoalbranch.seg_starts.push_back(first);
                newcoalbranch.seg_ends.push_back(last);            
                values.clear();
                first=(*it).first;
                last=-1;
                values.push_back((*it).second);
                seen++;
                temp_count=0;
            }
        }
        if (temp_count>0) {                
            int in_size = values.size();
            int in_seen = 0;
            for (int j=0;j<in_size;++j) {
                if ((*it).second == values[j]) {
                    in_seen++;
                    break;
                }
            }
            if (in_seen>0) {seen--;}
            if (in_seen==0) {
                values.push_back((*it).second);
                seen++;
            }
            if (seen==0) {
                last = (*it).first;
            }
        }     
        temp_count++;
    }
    newcoalbranch.seg_starts.push_back(first);
    newcoalbranch.seg_ends.push_back(last);    
        
    int end_size = newcoalbranch.seg_starts.size();    
    values.clear();
    merge_map.clear();    
     
    Current_Tree.push_back(newcoalbranch);
     
     int new_n = Current_Tree.size();
     total_seglength=0;
     for (int j=0;j<new_n;++j) {
        if (Current_Tree[j].current==1) {
            int blocks = Current_Tree[j].seg_starts.size();
            for (int h=0;h<blocks;++h) {
                total_seglength = total_seglength + (Current_Tree[j].seg_ends[h]-Current_Tree[j].seg_starts[h]);
            }
        }
    }
     
}      

//Initial stage of a coalescent event once it has been decided the next event is a coalescent event         
void coalEvent(vector<Branch_t> & Current_Tree,double which1,double which2,double tevent,int* count,double ttot,double & total_seglength,bool subdivision,int pop) {
    double tracker1 = 0;
    double tracker2 = 0;
    int n = Current_Tree.size();
    vector<int> active;
    
    //Choosing which two lines coalesce
    if (!subdivision) {
        for (int i=0;i<n;++i) {
            if (Current_Tree[i].current==1) {
                active.push_back(Current_Tree[i].num);
            }
        }
    }
    if (subdivision) {
        for (int i=0;i<n;++i) {
            if (Current_Tree[i].current==1 && Current_Tree[i].population == pop) {
                active.push_back(Current_Tree[i].num);
            }
        }
    }
    double num_active = active.size();
    double prob_choosefirst = (1/num_active);
    double prob_choosesecond = (1/(num_active-1));
    int first_line = -1;
    int second_line = -1;
    for (int i=0;i<num_active;++i) {
        if (which1 > tracker1 && (which1 < (tracker1 + prob_choosefirst))) {
            first_line = active[i];
            active.erase (active.begin()+i);
            break;
        }
        else {
            tracker1 = tracker1 + prob_choosefirst;
        }
    }
    for (int i=0;i<(num_active-1);++i) {
        if (which2 > tracker2 && which2 < (tracker2 + prob_choosesecond)) {
            second_line = active[i];
            break;
        }
        else {
            tracker2 = tracker2 + prob_choosesecond;
        }
    }        
        
    //Call the function which finishes the coalescent event
    performCoal(Current_Tree,first_line,second_line,pop,count,tevent,ttot,total_seglength);   
}

//Carry out a migration event when it has been decided the next event is a migration
void migrationEvent(vector<Branch_t> & Tree,double which,double tevent, double ttot,int origin_pop,int dest_pop) {
    //select which branch to migrate
    int n = Tree.size();
    double tracker = 0;
    vector<int> active;
    for (int i=0;i<n;++i) {
        if (Tree[i].current==1 && Tree[i].population == origin_pop) {
            active.push_back(Tree[i].num);
        }
    }
    int num_active = active.size();
    double prob_choose = 1.0/num_active;
    int mig_line = -1;
    for (int i=0;i<num_active;++i) {
        if (which > tracker && which <= (tracker + prob_choose)) {
            mig_line = active[i];
            break;
        }
        else {
            tracker = tracker + prob_choose;
        }
    }   
    
    //perform migration
   
    Tree[mig_line].population = dest_pop;
    Tree[mig_line].migration_times.push_back(ttot + tevent);
}

//If via recombination a mutation has been propogated to every haplotype in the sample this function removes these locations from the variable positions in the final haplotypes
void removeFullMutations(vector<Branch_t> & Tree,vector<double> & Pos) {
    vector<int> external;
    int num_branches = Tree.size();
    for (int i=0;i<num_branches;++i) {
        if (Tree[i].external==1) {
                external.push_back(Tree[i].num);
        }
    }
    int num_ext = external.size();
    int num_muts = Pos.size();
    vector<int> times_seen;
    for (int i=0;i<num_muts;++i) {
        times_seen.push_back(0);
    }
    for (int i=(num_ext-1);i>-1;--i) {
        int personal_muts = Tree[i].variant_locs.size();
        for (int j=0;j<personal_muts;++j) {
            for(int h=0;h<num_muts;++h) {
                if (Tree[i].variant_locs[j]==Pos[h]) {
                    times_seen[h]++;
                    break;
                }
            }
        }
    }
    for (int i=(num_muts-1);i>-1;--i) {
        if (times_seen[i]==num_ext) {
            Pos.erase (Pos.begin()+i);
        }
    }
} 

//This function checks to see if all segments have individually reached their MRCA
int allMRCA(vector<Branch_t> & Tree,double & total_seglength,bool subdivision,int & pop1_size,int & pop2_size,int & num_rem) {
    
    multimap<int, int> branch_count_match;
    multimap<double, int> all_segs_per_branch;
    vector<int> current;
    int num_branches = Tree.size();
    
    for (int i=0;i<num_branches;++i) {
        if (Tree[i].current==1) {
            current.push_back(Tree[i].num);
        }
    }
    
    int num_current = current.size();
    int tracker = 0;
    int tracker2 = 0;
    for (int i=0;i<num_current;++i) {
        int s_length = Tree[current[i]].seg_starts.size();
        branch_count_match.insert(pair<int, int>(current[i], tracker2));
        for(int j=0;j<s_length;++j) {
            all_segs_per_branch.insert(pair<double, int>(Tree[current[i]].seg_starts[j],tracker2));
            all_segs_per_branch.insert(pair<double, int>(Tree[current[i]].seg_ends[j],tracker2));
            tracker++;
        }
        tracker2++;
    }
    
    vector<int> score;
    vector<int> current_open;
    vector<int> future_open;
    for (int i=0;i<num_current;++i) {
        score.push_back(0);
        current_open.push_back(0);
        future_open.push_back(0);
    } 
    double current_break = -1;
    int temp_count  =-1 ;
    int current_score = 0;
    for (multimap<double, int>::iterator it = all_segs_per_branch.begin();it != all_segs_per_branch.end();++it)  {
        if (temp_count==-1) {
            current_break = (*it).first;
            temp_count++;
        }
        if ((*it).first==current_break) {
            if (current_open[(*it).second]==0) {
                future_open[(*it).second] = 1;
            }
            if (current_open[(*it).second]==1) {
                future_open[(*it).second] = 0;
                if (current_score > 1) {
                    for (int h=0;h<num_current;++h) {
                        if (current_open[h]==1) {
                            score[h]++;
                        }
                    }
                    if (total_seglength>2) {
                        return 0;
                    }
                }
            }
        } else {
            current_open = future_open;
            current_break = (*it).first;
            current_score = 0;
            for (int j=0;j<num_current;++j) {
                current_score = current_score + current_open[j];
            }
            if(current_open[(*it).second]==0) {
                future_open[(*it).second] = 1;
            }
            if(current_open[(*it).second]==1) {
                future_open[(*it).second] = 0;
                if (current_score > 1) {
                    for (int h=0;h<num_current;++h) {
                        if (current_open[h]==1) {
                            score[h]++;
                        }
                    }
                    if (total_seglength>2) {
                        return 0;
                    }
                }
            }
        }
    }
    
    int mrca = 0;
    for (int i=0;i<num_current;++i) {
        if (score[i]>0) {
            mrca++;
        }
        if (score[i]==0) {
            for (multimap<int, int>::iterator it = branch_count_match.begin();it != branch_count_match.end();++it)  {
                if ((*it).second==i) {
                    if (subdivision) {
                        if(Tree[(*it).first].population==1) {
                            pop1_size--;
                        }
                        if(Tree[(*it).first].population==2) {
                            pop2_size--;
                        }
                        Tree[(*it).first].current=0;
                        break;
                    } else {
                        num_rem--;
                        Tree[(*it).first].current=0;
                        break;
                    }
                
                }
            }
        }
    }
    
    vector<int> active;
    for (int i=0;i<num_branches;++i) {
           if (Tree[i].current==1) {
               active.push_back(Tree[i].num);
           }
    }
    int num_active = active.size();
    if (num_active < 1) {
        int we;
    }
    
    
    total_seglength = 0;
    for (int j=0;j<num_branches;++j) {
        if (Tree[j].current==1) {
            int blocks = Tree[j].seg_starts.size();
            for (int h=0;h<blocks;++h) {
                total_seglength = total_seglength + (Tree[j].seg_ends[h]-Tree[j].seg_starts[h]);
            }
        }
    } 
    if (mrca>0) {
        return 0;
    } else {
        return 1;
    }
}     

//This function goes through the tree and generates a vector with the external branch numbers reachable from the current haplotype segment branch      
void everyNode(vector<Branch_t> & Tree,double point,int current,vector<int> & found) {
    int s1 = Tree[current].seg_starts.size();
    int f1 = found.size();
    int left = Tree[current].left;
    int right = Tree[current].right;
    if(Tree[current].external==1) {
        int seen = 0;
        for(int i=0;i<s1;++i) {
            if (point >= Tree[current].seg_starts[i] && point <= Tree[current].seg_ends[i]) {
                ++seen;
                break;
            }
        }
        if (seen>0) {
            int in_seen = 0;
            for(int i=0;i<f1;++i) {
                if (current==found[i]) {
                    ++in_seen;
                    break;
                }
            }
            if (in_seen==0) {
                found.push_back(current);
                int number = found.size();
            }
        }
    } else {
        if(left!=-1) {
            int l1 = Tree[left].seg_starts.size();
            int l_seen = 0;
            for (int i=0;i<l1;++i) {
                if (point >= Tree[left].seg_starts[i] && point <= Tree[left].seg_ends[i]) {
                    ++l_seen;
                    break;
                }
            }
            if (l_seen>0) {
                everyNode(Tree,point,left,found);
            }
        }
        if(right!=-1) {
            int r1 = Tree[right].seg_starts.size();
            int r_seen = 0;
            for (int i=0;i<r1;++i) {
                if (point >= Tree[right].seg_starts[i] && point <= Tree[right].seg_ends[i]) {
                    ++r_seen;
                    break;
                }
            }
            if (r_seen>0) {
                everyNode(Tree,point,right,found);
            }
        }
    }
}
    
//This function will return the MRCA for a given segment of haplotype
double getMRCA(vector<Branch_t> & Tree,int current) {
    double MRCA = 0;
    int keep_going = 0;
    while(keep_going==0) {
        if (Tree[current].left!=-1) {
            current = Tree[current].left;
            MRCA = MRCA + Tree[current].length;
        } else {
            if (Tree[current].right!=-1) {
                current = Tree[current].right;
                MRCA = MRCA + Tree[current].length;
            } else {
                keep_going=1;
            }
        }
    }
    return MRCA;
}

//This function finds the MRCA of each haplotype segment and uses this to calculate the total Tree length and the TMRCA
void getLengths(vector<Branch_t> & Tree,vector<double> & Recomb_Pos,vector<double> & lengths,int num_ext) {
    std::sort(Recomb_Pos.begin(),Recomb_Pos.end());    
    int size = Recomb_Pos.size();
    int tree_size = Tree.size();
    double TMRCA = 0;
    double total_length = 0;
    for (int i=0;i<(size-1);++i) {
        int j = i + 1;
        double weight = Recomb_Pos[j] - Recomb_Pos[i];
        double midpoint = (Recomb_Pos[j] + Recomb_Pos[i])/2;
        double current_length = 0;
        for (int h=0;h<tree_size;++h) {
            int full = -1;
            vector<int> found;
            everyNode(Tree,midpoint,h,found);
            int found_size = found.size();
            if (found_size==num_ext) {full=1;}
            if (found_size<num_ext) {full=0;}
            if (found_size>num_ext) {
                int we_have_problem=1;
            }
            
            //use found vector here
//            cout << "Current Branch: " << h << " Found vector: ";
//            for (int k=0;k<found_size;++k) {
//               cout << found[k] << " ";
//            }
//            cout << endl;
            
            if (full==0) {
                int s1 = Tree[h].seg_starts.size();
                for (int k=0;k<s1;++k) {
                    if (midpoint >= Tree[h].seg_starts[k] && midpoint <= Tree[h].seg_ends[k]) {
                        current_length = current_length + Tree[h].length;
                        break;
                    }
                }
            }
            if (full==1) {
                total_length = total_length + (current_length*weight);
                double current_MRCA = getMRCA(Tree,h);
                TMRCA = TMRCA + (current_MRCA*weight);
                break;
            }
        }
    }
    lengths.push_back(TMRCA);
    lengths.push_back(total_length);
}

int infTest(double time_grow,double alpha,double bigN,double fte_exp) {
    double infvalue = 1/(-alpha*(fte_exp-1)*pow(bigN,fte_exp-1));
    if (infvalue < time_grow) {
        return 1;
    } else {
        return 0;
    }
}

int tinyTest(int gro,double time_grow,double alpha,double bigN,double fte_exp) {
    double size = 2;
    if (gro==1) {
        size = bigN*exp(-alpha*time_grow);
    }
    if (gro==2) {
        double dN = pow(bigN,fte_exp-1);
        size = pow((dN/(1+(alpha*(fte_exp-1)*dN*time_grow))),(1/(fte_exp-1)));
    }
    if (size<1) {
        return 1;
    } else {
        return 0;
    }
}


////////////////////////////////////Branch Descendant Crap///////////////////////////////////



int numberDescendants(int current,vector<Branch_t> & Tree) {
    int descendants = 0;
    if (Tree[current].left!=-1) {
        if (Tree[Tree[current].left].external==1) {
            ++descendants;
        } else {
            int leftdescendants = numberDescendants(Tree[current].left,Tree);
            descendants = descendants + leftdescendants;
        }
    }
    if (Tree[current].right!=-1) {
        if (Tree[Tree[current].right].external==1) {
            ++descendants;
        } else {
            int leftdescendants = numberDescendants(Tree[current].right,Tree);
            descendants = descendants + leftdescendants;
        }
    }
    return descendants;
}
        
int linesWithDescendants(vector<Branch_t> & Tree, int interest) {
    int num_branches = Tree.size();
    int number = 0;
    for (int i=0;i<num_branches;++i) {
        if (Tree[i].current==1) {
            int descendants = numberDescendants(i,Tree);
            if (descendants==interest) {
                ++number;
            }
        }
    }    
    return number;
}

void split(const string& str, const string& delimiters, vector<string>& tokens) {
        string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        string::size_type pos = str.find_first_of(delimiters, lastPos);
        while (string::npos != pos || string::npos != lastPos) {
            tokens.push_back(str.substr(lastPos,pos-lastPos));
            lastPos=str.find_first_not_of(delimiters,pos);
            pos = str.find_first_of(delimiters,lastPos);
        }
}


void newGetLengths(vector<Branch_t> & Tree,vector<double> & Recomb_Pos,int num_ext, double bigN,std::vector<int>& sizes) {
    std::sort(Recomb_Pos.begin(),Recomb_Pos.end());    
    int size = Recomb_Pos.size();
    int tree_size = Tree.size();
    double total_length = 0;
    
    for (int h=0;h<tree_size;++h) {
        
        map<std::string, double> desc_map;

        for (int i=0;i<(size-1);++i) {
            int j = i + 1;
            double weight = Recomb_Pos[j] - Recomb_Pos[i];
            double midpoint = (Recomb_Pos[j] + Recomb_Pos[i])/2;
            int run = 0;
            int s1 = Tree[h].seg_starts.size();
            for (int k=0;k<s1;++k) {
                if (midpoint >= Tree[h].seg_starts[k] && midpoint <= Tree[h].seg_ends[k]) {
                    ++run;
                    break;
                }
            }
            if (run>0) {
                vector<int> found;
                everyNode(Tree,midpoint,h,found);
                int desc = found.size();
                sort(found.begin(),found.end());
                
                std::stringstream ss;
                
                ss << desc << ":";
                for (int k=0;k<(desc-1);++k) {
                    ss << found[k] << ",";
                }
                ss << found[(desc-1)];
                
                string value = ss.str();
                double width = weight;
                int seen = 0;
                for (map<std::string, double>::iterator it=desc_map.begin();it!=desc_map.end();++it) {
                    string current = (*it).first;
                    if (current==value) {
                        width = width + (*it).second;
                        ++seen;
                    }
                }
                if (seen>0) {
                    desc_map.erase(value);
                }
                desc_map.insert(pair<std::string, double>(value,width));
            }
        }
        
        for (map<std::string, double>::iterator it=desc_map.begin();it!=desc_map.end();++it) {
            vector<string> kinfo;
            string delim = ":";
            string tinfo = (*it).first;
            split(tinfo,delim,kinfo);

            
            double second = atof(kinfo[0].c_str());
            double b_length = 2*bigN*exp(log((*it).second) + log(Tree[h].length));
            total_length+=b_length;
	    for (int ii=0;ii<sizes.size();++ii) {
	      if (second==sizes[ii]) {
		cout << h << "," << kinfo[0] << "," << (*it).second << "," << Tree[h].length << "," << b_length << endl;    
		break;
	      }
	    }
            kinfo.clear();            
        }
       
        
        desc_map.clear();   
    }
     cout << "Total Length: " << total_length << endl;
}
                       
////////////////////////////
///////////////////////////
            
            
                                           
int main(int argc, char** argv) {
    try {
        TCLAP::CmdLine cmd("For more detailed parameter descriptions please see the users guide", ' ', "1.39");
        TCLAP::ValueArg<double> recombRate("r","recombination_rate","Recombination rate r, this is the per generation probability a recombination occurs in the region be simulated, along with value entered for -N FTEC calculates the value rho as 4Nr",false,0,"double");
        cmd.add(recombRate);
        TCLAP::SwitchArg subDiv("s","subdivision","Is the population subdivided?",false);
        cmd.add(subDiv);
        TCLAP::ValueArg<double> mutRate("m", "mutation_rate", "REQUIRED: Mutation rate m, this is the per base mutation rate multiplied by the number of bases being simulated, along with value entered for -N FTEC calculates the value theta as 4Nm",true,-1,"double");
        cmd.add(mutRate);
        TCLAP::ValueArg<int> numExt("n", "number_samples", "REQUIRED: Total number of samples to simulate",true,-1,"double");
        cmd.add(numExt);
        TCLAP::ValueArg<int> pop2Size("","p2n", "Number of simulated samples drawn from population 2",false,-1,"double");
        cmd.add(pop2Size);
        TCLAP::ValueArg<double> pop2Prop("","p2prop", "Size of population 2 relative to population 1",false,0,"double");
        cmd.add(pop2Prop);
        TCLAP::ValueArg<double> popMerge("","merge_time","Time until populations merge into a single ancestral population",false,-1,"double");
        cmd.add(popMerge);
        TCLAP::MultiArg<double> migRate("","migration","Migration rate m, proportion of chromosomes leaving a population each generatino, if a value is entered for --migration2 this value becomes the migration rate from population 1 to population 2",false,"double");
        cmd.add(migRate);
        TCLAP::MultiArg<double> migRate2("","migration2","If differential migration rates are desired, this value is the migration rate from population 2 to population 1",false,"double");
        cmd.add(migRate2);  
        TCLAP::MultiArg<double> migTime("","migration_time","Time(s) for migration rate changes",false,"double ");
        cmd.add(migTime);
        TCLAP::MultiArg<double> migTime2("","migration_time2","If differential migration rates are used, this parameter gives the times at which the rate from population 2 to population 1 changes",false,"double");
        cmd.add(migTime2); 
        TCLAP::ValueArg<int> typeGrowth("g","growth_type","Type of continuous size change",false,-1,"int");
        cmd.add(typeGrowth);
        TCLAP::ValueArg<int> repetitions("x","repetitions","Number of simulation repetitions",false,1,"int");
        cmd.add(repetitions);
        TCLAP::ValueArg<double> timeGrow("t","time_to_grow","Time, in coalescent units, backwards from the present during which the population size is continuously changing",false,-1,"double");
        cmd.add(timeGrow);
        TCLAP::ValueArg<double> Beta("a","growth_constant_alpha","Constant alpha",false,0,"double");
        cmd.add(Beta);
        TCLAP::MultiArg<double> iG("b","intant_prop","Proportion of population remaining after instantaneous growth/contraction event",false,"double ");
        cmd.add(iG);
        TCLAP::MultiArg<double> tInstant("","t_instant","Time in the past, in coalescent units, at which an instantaneous growth/contraction event occured",false,"double ");
        cmd.add(tInstant);
        TCLAP::MultiArg<double> iG2("","b2","Proportion of population 2  remaining after instantaneous growth/contraction event",false,"double ");
        cmd.add(iG2);
        TCLAP::MultiArg<double> tInstant2("","t_instant2","Time in the past, in coalescent units, at which an instant growth/contraction event occured in population 2",false,"double ");
        cmd.add(tInstant2);
        TCLAP::ValueArg<double> B2("","a2","Growth constant alpha for population 2",false,0,"double");
        cmd.add(B2);
        TCLAP::ValueArg<double> timeGrow2("","t2","Time in the past, in coalescent units, during which population 2 was continuously changing",false,-1,"double");
        cmd.add(timeGrow2);
        TCLAP::ValueArg<int> typeGrowth2("","g2","Type of growth in population 2",false,-1,"int");
        cmd.add(typeGrowth2);
        TCLAP::ValueArg<double> bN("N","pop_size","REQUIRED: Current population size",true,0,"double");
        cmd.add(bN);
        TCLAP::SwitchArg timed("","time","Should program print running times?",false);
        cmd.add(timed);
        TCLAP::ValueArg<double> fte_constant("c","fte_exp","Faster than exponential growth exponent",false,-99,"double");
        cmd.add(fte_constant);
        TCLAP::ValueArg<double> fte_constant2("","c2","Faster than exponential growth exponent in population 2",false,-99,"double");
        cmd.add(fte_constant2);
        TCLAP::SwitchArg lengthed("L","length","Should program output total tree length and TMRCA (in 2N)?",false);
        cmd.add(lengthed);
        TCLAP::ValueArg<double> mSize("","merge_size","Size of merged population as proportion of population 1 starting size",false,-9,"double");
        cmd.add(mSize);
        TCLAP::ValueArg<unsigned int> rSeed("","seed","Random number generator seed",false,1,"unsigned int");
        cmd.add(rSeed);
	TCLAP::ValueArg<std::string> oSizes("","osizes","What branch sizes should be output, comma separated list",true,"2","string");
        cmd.add(oSizes);
	
         
        cmd.parse(argc,argv);      
        
        cout << "FTEC ";
        int num_ext = numExt.getValue();
        cout << num_ext << " ";
        int reps = repetitions.getValue();
        cout << reps << " ";
        double mutation = mutRate.getValue();
        cout << "-m " << mutation << " ";
        double rec = recombRate.getValue();
        cout << "-r " << rec << " ";
        double bigN = bN.getValue();
        cout << "-N " << bigN << " ";
        double time_grow = timeGrow.getValue();
        if (time_grow!=-1) {cout << "-t " << time_grow << " ";}
        int gro = typeGrowth.getValue();
        if (gro!=-1) {cout << "-g " << gro << " ";}
        double beta = Beta.getValue();
        if (beta!=0) {cout << "-a " << beta << " ";}
        double fte_exp = fte_constant.getValue();
        if (fte_exp!=-99) {cout << "-c " << fte_exp << " ";}
        vector<double> instgrow = iG.getValue();
	vector<std::string> o_sizes;
	std::string presizes = oSizes.getValue();


	std::string d1 = ",";
	split(presizes,d1,o_sizes);
	std::vector<int> a_sizes;
	for (int ii=0;ii<o_sizes.size();++ii) {
	  a_sizes.push_back(atoi(o_sizes[ii].c_str()));
	}
	  

       int inst_empty = 0;
        if (instgrow.empty()) {
            ++inst_empty;
            instgrow.push_back(-1);
        }
        if (inst_empty==0) {
            cout << "-b ";
            for (int i=0;i<instgrow.size();++i) {
                cout << instgrow[i] << " ";
            }
        } 
        vector<double> instgrow2 = iG2.getValue();
        int inst_empty2 = 0;
        if (instgrow2.empty()) {
            ++inst_empty2;
            instgrow2.push_back(-1);
        }
        if (inst_empty2==0) {
            cout << "--b2 ";
            for (int i=0;i<instgrow2.size();++i) {
                cout << instgrow2[i] << " ";
            }
        }       
        vector<double> t_instant = tInstant.getValue();
        int t_inst_empty = 0;
        if (t_instant.empty()) {
            ++t_inst_empty;
            t_instant.push_back(-1);
        }
        if (t_inst_empty==0) {
            cout << "--t_instant ";
            for (int i=0;i<t_instant.size();++i) {
                cout << t_instant[i] << " ";
            }
        }
        vector<double> t_instant2 = tInstant2.getValue();
        int t_inst_empty2 = 0;
        if (t_instant2.empty()) {
            ++t_inst_empty2;
            t_instant2.push_back(-1);
        }
        if (t_inst_empty2==0) {
            cout << "--t_instant2 ";
            for (int i=0;i<t_instant2.size();++i) {
                cout << t_instant2[i] << " ";
            }
        }
        bool outersubdivision = subDiv.getValue();
        if (outersubdivision) {cout << "-s ";}
        vector<double> mig = migRate.getValue();
        int was_empty1 = 0;
        if (mig.empty()) {
            ++was_empty1;
            mig.push_back(0);
        }
        if (was_empty1==0) {
            int mig_size = mig.size();
            cout << "--migration ";
            for (int i=0;i<mig_size;++i) {
                cout << mig[i] << " ";
            }
        }
        vector<double> mig2 = migRate2.getValue();
        int was_empty2 = 0;
        if (mig2.empty()) {
            ++was_empty2;
            mig2.push_back(0);
        }
        if (was_empty2==0) {
            int mig2_size = mig2.size();
            cout << "--migration2 ";
            for (int i=0;i<mig2_size;++i) {
                cout << mig2[i] << " ";
            }
        }
        vector<double> mtimes = migTime.getValue();
        int was_empty3 = 0;
        if (mtimes.empty()) {
            ++was_empty3;
            mtimes.push_back(-1);
        }
        if (was_empty3==0) {
            int mtimesize = mtimes.size();
            cout << "--migration_times ";
            for (int i=0;i<mtimesize;++i) {
                cout << mtimes[i] << " ";
            }
        }
        vector<double> mtimes2 = migTime2.getValue();
        int was_empty4 = 0;
        if (mtimes2.empty()) {
            ++was_empty4;
            mtimes2.push_back(-1);
        }
        if (was_empty4==0) {
            int mtimesize2 = mtimes2.size();
            cout << "--migration_times2 ";
            for (int i=0;i<mtimesize2;++i) {
                cout << mtimes2[i] << " ";
            }
        }
        double tao = popMerge.getValue();
        if (tao!=-1) {cout << "--merge_time " << tao << " ";}
        int outer_pop2_size = pop2Size.getValue();
        if (outer_pop2_size!=-1) {cout << "--p2n " << outer_pop2_size << " ";}
        double pop2_prop = pop2Prop.getValue(); 
        if (pop2_prop!=0) {cout << "--p2prop " << pop2_prop << " ";}
        double time_grow2 = timeGrow2.getValue();
        if (time_grow2!=-1) {cout << "--t2 " << time_grow2 << " ";}
        int gro2 = typeGrowth2.getValue();
        if (gro2!=-1) {cout << "--g2 " << gro2 << " ";}
        double beta2 = B2.getValue();
        if (beta2!=0) {cout << "--a2 " << beta2 << " ";}
        double fte_exp2 = fte_constant2.getValue();
        if (fte_exp2!=-99) {cout << "--c2 " << fte_exp2 << " ";}
        double mergesize = mSize.getValue();
        if (mergesize!=-9) {cout << "--merge_size " << mergesize << " ";}
        bool report_lengths = lengthed.getValue();
        if (report_lengths) {cout << "-L ";}
        bool timer = timed.getValue();
        if (timer) {cout << "--time ";}
        cout << endl;
        int mig_size = mig.size();
        int mtimes_size = mtimes.size();
        int mig2_size = mig2.size();
        int mtimes2_size = mtimes2.size();
        unsigned int randseed = rSeed.getValue();
        
        clock_t start = clock();

        //Uniform Number Generator setup
        prgType rng;
        unsigned int seed = time(0);
        if (randseed != 1) {
            seed = randseed;
        }
        rng.seed(seed);
        boost::uniform_real<> cl(0,1);
        boost::variate_generator<prgType&, boost::uniform_real<> > getunif(rng,cl);
    
//    cout << time(0) << endl<<endl;
        
        //Checking the parameters entered on the command line for unsupported or illegitimate values
        if ((!outersubdivision) && ((pop2_prop!=0) || (outer_pop2_size!=-1))) {
            cout << "Populations subdivision requires subdivision flag -s\n";
            exit(274);
        }
        if (gro!=-1 && t_inst_empty==0) {
            for (int i=0;i<t_instant.size();++i) {
                if (t_instant[i] < time_grow) {
                    cout << "Program doesn't support instant size changes while population is growing, please enter new -t or --t_instant values\n";
                    exit(273);
                }
            }
        }
        if (gro2!=-1 && t_inst_empty2==0) {
            for (int i=0;i<t_instant2.size();++i) {
                if (t_instant2[i] < time_grow2) {
                    cout << "Program doesn't support instant size changes while population is growing, please enter new --t2 or --t_instant2 values\n";
                    exit(273);
                }
            }
        }
        if (inst_empty==0) {
            if (instgrow.size()!=t_instant.size()) {
                cout << "Number of instant size change times and new populations size proportions need to match, please enter an equal number of --t_instant and -b values\n";
                exit(272);
            }
        }
        if (inst_empty2==0) {
            if (instgrow2.size()!=t_instant2.size()) {
                cout << "Number of instant size change times and new populations size proportions need to match for population 2, please enter an equal number of --t_instant2 and -b values\n";
                exit(272);
            }
        }
        if (t_inst_empty2==0 && tao!=-1) {
            for (int i=0;i<t_instant2.size();++i) {
                if (t_instant2[i] > tao) {
                    cout << "Cannot have instant size change in population 2 after populations have merged, please enter new --t_instant2 or --merge_time value\n";
                    exit(271);
                }
            }
        }            
        if (t_inst_empty==0 && (t_instant.size() > 1)) {
            for (int i=0;i<t_instant.size();++i) {
                if (t_instant[i] <= 0) {
                    cout << "Demographic events can only occur in the past, please enter positive --t_instant value(s)\n";
                    exit(248);
                }
                for (int j=(i+1);j<t_instant.size();++j) {
                    if (t_instant[j] <= t_instant[i]) {
                        cout << "Please enter instant population size change times in chronological order\n";
                        exit(270);
                    }
                }
            }
        }
        if (t_inst_empty2==0 && (t_instant2.size() > 1)) {
            for (int i=0;i<t_instant2.size();++i) {
                if (t_instant[i] <= 0) {
                    cout << "Demographic events can only occur in the past, please enter positive --t_instant2 value(s)\n";
                    exit(248);
                }
                for (int j=(i+1);j<t_instant2.size();++j) {
                    if (t_instant2[j] <= t_instant2[i]) {
                        cout << "Please enter instant population 2 size change times in chronological order\n";
                        exit(270);
                    }
                }
            }
        }    
        if ((inst_empty==0 && t_inst_empty!=0) || (inst_empty!=0 && t_inst_empty==0)) {
            cout << "Both a time and new population size are required for instant population size change, please enter both valid -b and --t_instant values\n";
            exit(269);
        }
        if ((inst_empty2==0 && t_inst_empty2!=0) || (inst_empty2!=0 && t_inst_empty2==0)) {
            cout << "Both a time and new population 2 size are required for instant population size change, please enter both valid --b2 and --t_instant2 values\n";
            exit(269);
        }
        if (inst_empty==0) {
            for (int i=0;i<instgrow.size();++i) {
                if (instgrow[i] <= 0) {
                    cout << "Population must have size greater than 0, please enter valid -b value\n";
                    exit(268);
                }
            }
        }
        if (inst_empty2==0) {
            for (int i=0;i<instgrow2.size();++i) {
                if (instgrow2[i] <= 0) {
                    cout << "Population 2 must have size greater than 0, please enter valid --b2 value\n";
                    exit(268);
                }
            }
        }
        if ((outersubdivision) && pop2_prop<=0) {
            cout << "\nPopulation 2 must have postive relative size, pleae enter valid --p2prop value\n";
            exit(267);
        }
        if ((outersubdivision) && outer_pop2_size<0) {
            cout << "\nPopulation 2 needs non-negative size, please enter valid --p2n value\n";
            exit(266);
        }
        if (was_empty3==0) {
            for (int i=0;i<mtimes_size;++i) {
                if (mtimes[i]<0) {
                    cout << "\nMigration rate change times need to be positive, please enter new --migration_time values\n";
                    exit(265);
                }
            }
        }
        if (was_empty4==0) {
            for (int i=0;i<mtimes2_size;++i) {
                if (mtimes2[i]<0) {
                    cout << "\nMigration2 rate change times need to be positive, please enter new --migration_time2 values\n";
                    exit (265);
                }
            }
        }
        if (was_empty3==0) {
            int temp_j = 0;
            while (temp_j<mtimes_size) {
                for (int i=1;i<mtimes_size;++i) {
                    if (mtimes[i]<=mtimes[temp_j]) {
                        cout << "\nMigration rate change times in incorrect order, please enter change times in chronological order\n";
                        exit(264);
                    }
                }
                ++temp_j;
            }
        }
        if (was_empty4==0) {
            int temp_j = 0;
            while (temp_j<mtimes2_size) {
                for (int i=1;i<mtimes2_size;++i) {
                    if (mtimes2[i]<=mtimes2[temp_j]) {
                        cout << "\nMigration2 rate change times in incorrect order, please enter change times in chronological order\n";
                        exit(264);
                    }
                }
                ++temp_j;
            }
        }
        for (int i=0;i<mig_size;++i) {
            if (mig[i]<0) {
                cout << "\nMigration rates need to be positive, please enter new --migration values\n";
                exit(263);
            }
        }
        for (int i=0;i<mig2_size;++i) {
            if (mig2[i]<0) {
                cout << "\nMigration2 rates need to be positive, please enter new --migration values\n";
                exit(263);
            }
        }
        if (mig2_size!=(mtimes_size+1) && mtimes2[0]!=-1) {
            cout << "\nIncorrect Number of migration parameters and times to change migration rate for population 2, please enter new --migration2 and --migration_time2 value\n";
            exit(262);
        }
        if (mig_size!=(mtimes_size+1) && mtimes[0]!=-1) {
            cout << "\nIncorrect Number of migration parameters and times to change migration rate, please enter new --migration and --migration_time values\n";
            exit(262);
        }  
        if (tao!=-1) {
            for (int i=0;i<mtimes_size;++i) {
                if (mtimes[i] >= tao) {
                    cout << "\nInvalid time for migration rate change, occurs after populations merge.  Please enter new --migration_time or --merge_time values\n";
                    exit(260);
                }
            }
            for (int i=0;i<mtimes2_size;++i) {
                if (mtimes2[i] >= tao) {
                    cout << "\nInvalid time for migration rate change, occurs after populations merge.  Please enter new --migration_time2 or --merge_time values\n";
                    exit(260);
                }
            }
        }
        if ((gro == 2) && (fte_exp ==1)) {
            cout << "\nFor exponential growth please enter a 1 at the -g command\n";
            exit(259);
        }
        
        if ((gro == 2) && (fte_exp <= 0)) {
            cout << "\nFTEC currently requires an exponent value greater than 0\n";
            exit(258);
        }
        if ((gro2 == 2) && (fte_exp2 <= 0)) {
            cout << "\nFTEC currently requires an exponent value greater than 0, please check value for population 2\n";
            exit(258);
        }
        if (reps <= 0) {
            cout << "\nNumber of repetitions less than or equal to 0, performing no simulations\n";
            exit(254);
        }
        if ((gro != -1) && (gro != 1 ) && (gro != 2)) {
            cout << "\nGrowth type not currently supported, please enter either nothing, 1, or 2 at the -g prompt\n";
            exit(253);
        }
        
        if (gro2 != -1 && gro2 != 1 && gro2 != 2) {
            cout << "\nGrowth type in second population not currently supported, please enter either nothing, 1, or 2 at the -g prompt\n";
            exit(252);
        }
        if ((gro != -1) && (beta==0)) {
            cout << "\nMissing growth constant, population size remaining constant\n";
            gro = -1;
        }
        if ((gro2 != -1) && (beta2==0)) {
            cout << "\nMissing growth constant for second population, population size remaining constant\n";
            gro2 = -1;
        }
        if (mergesize <= 0 && mergesize != -9 ) {
            cout << "\nMerged population must have size greater than 0\n";
            exit(247);
        }
        if ((outersubdivision) && tao==-1 && (was_empty1>0) && (was_empty2>0)) {
            cout << "\nWith no merge time and no migration simulation will never reach MRCA\n";
            exit(245);
        }
        if ((time_grow > tao && tao > 0) || (time_grow2 > tao && tao > 0)) {
	        cout << "\nProgram does not support growth after populations merge, please set new --merge_time or -t and/or --t2 values\n";
	        exit(261);
	    }
        for (int re=0;re<reps;++re) {
            
        if (re==0) {
            cout.precision(0);
       //     cout << "fte " << num_ext << " " << reps << endl;
            cout << fixed << seed << endl << endl;
            cout.precision(6);
        }    
        
        int lines_remaining = 0;    
        bool subdivision = outersubdivision;
        //if (subdivision) {cout << "true" << endl;}
        double total_seglength = num_ext;
        int num_rem = num_ext;
        bool ever_subdivision = subdivision;
        vector<Branch_t> Tree;
        vector<double> Variant_Pos;
        vector<double> Recomb_Pos;
        Recomb_Pos.push_back(0);
        Recomb_Pos.push_back(1);
        int pop2_size = outer_pop2_size;
        int pop1_size = 0;
        double ttot = 0;
        double theta = 4*bigN*mutation;
        double recomb = 4*bigN*rec;
        double M1 = 4*bigN*mig[0];
        double M2 = M1;
        if (was_empty2==0) {
            M2 = 4*bigN*mig2[0];
        }
        double tadd_subdiv = -1;
        double tadd_growth = -1;
        double tadd_growth2 = -1;
        double tadd_instgrow = -1;
        double tadd_instgrow2 = -1;
        double new_beta = -1;
        double bigN2 = pop2_prop*bigN;
        int growth = gro;
        int growth2 = gro2;
        double new_beta2 = -1;
        int rep_counter = 0;
        int mig_rate_count = 0;
        double t_add_mig_change1 = 0;
        int mig2_rate_count = 0;
        double t_add_mig_change2 = 0;
        int inst_add_count = 0;
        int inst_add_count2 = 0;
        
        if (gro==2 && beta < 0) {
            int inf_size = infTest(time_grow,beta,bigN,fte_exp);
            if (inf_size==1) {
                cout << "\nWith current parameters population 1 would have been infinite in size in the past, please enter new parameters.\n";
                exit(275);
            }
        }
        if (gro2==2 && beta2 < 0) {
            int inf_size = infTest(time_grow2,beta2,bigN2,fte_exp2);
            if (inf_size==1) {
                cout << "\nWith current parameters population 2 would have been infinite in size in the past, please enter new parameters.\n";
                exit(275);
            }
        }
        
        if (gro!= -1 && beta > 0) {
            int tiny_size = tinyTest(gro,time_grow,beta,bigN,fte_exp);
            if (tiny_size==1) {
                cout << "\nWith current parameters population 1 shrinks too far and would have size < 1 individual in the past, please enter new parameters.\n";
                exit(276);
            }
        }
        if (gro2!= -1 && beta2 > 0) {
            int tiny_size = tinyTest(gro,time_grow2,beta2,bigN2,fte_exp2);
            if (tiny_size==1) {
                cout << "\nWith current parameters population 2 shrinks too far and would have size < 1 individual in the past, please enter new parameters.\n";
                exit(276);
            }
        }
        
        if (time_grow==-1 && growth>=0) {
            cout << "\nPlease enter the amount of time during which the population was growing\n";
            exit(255);
        }
        if (theta <= 0) {
            cout << "\nThe mutation rate must be greater than 0\n";
            exit(256);
        }
        if (recomb < 0) {
            cout << "\nThe recombination rate must be greater than or equal to 0\n";
            exit(257);
        }
        

	    int count = 0;
	    if(subdivision) {
	        pop1_size = num_ext - pop2_size;
	        if (pop1_size <=0) {
	            cout << "\nPopulation 2 too large for total population size\n";
	            exit(262);
	        }
	        
	        //Create external branches for each population
	        for (int i=0;i<pop1_size;++i) {
	            Branch_t extBranch;
	            extBranch.num = i;
	            extBranch.current=1;
	            extBranch.length=0;
	            extBranch.population = 1;
	            extBranch.external=1;
	            extBranch.left= -1;
	            extBranch.upper_left=-1;
	            extBranch.upper_right=-1;
	            extBranch.right= -1;
	            extBranch.total_time_to_event=0;
	            extBranch.seg_starts.push_back(0);
	            extBranch.seg_ends.push_back(1);
	            Tree.push_back(extBranch);
	            ++count;
	        }  
	    
	        for (int i=0;i<pop2_size;++i) {
	            Branch_t extBranch;
	            extBranch.num = count;
	            extBranch.current=1;
	            extBranch.length=0;
	            extBranch.population = 2;
	            extBranch.external=1;
	            extBranch.left= -1;
	            extBranch.right= -1;
	            extBranch.upper_left=-1;
	            extBranch.upper_right=-1;
	            extBranch.total_time_to_event=0;
	            extBranch.seg_starts.push_back(0);
	            extBranch.seg_ends.push_back(1);
	            Tree.push_back(extBranch);
	            ++count;
	        }  
	        
	        //while the two populations are seperate generate times until next event
	        while (subdivision) {
 
               double tevent = 0;
	           int event = -1;
	           vector<double> times;
	           double t_coal_pop1 = timeToCoal(pop1_size,growth,beta,bigN,rng,ttot,new_beta,fte_exp,1);
	           times.push_back(t_coal_pop1);
               double t_coal_pop2 = timeToCoal(pop2_size,growth2,beta2,bigN,rng,ttot,new_beta2,fte_exp2,pop2_prop);
	           times.push_back(t_coal_pop2);
               double t_mig_pop1 = timeToMig(pop1_size,rng,1,M1);
	           times.push_back(t_mig_pop1);
               double t_mig_pop2 = timeToMig(pop2_size,rng,2,M2);
	           times.push_back(t_mig_pop2); 
               double t_rec_event = timeToRecomb(total_seglength,recomb,rng,num_rem);
	           times.push_back(t_rec_event);
	           
               double min_time = BIGINT;
	        //select earliest of possible events
	            for (int i=0;i<5;++i) {
	                if (times[i] < min_time) {
	                    tevent = times[i];
	                    event = i;
	                    min_time = tevent;
	                }
	           }

	           //if time to next event is longer than time until populations merge => end subdivision, end growth, calculate new population size, move to non-subdivided portion of program below
	           if (ttot + tevent > tao && tao > 0) {
	               tadd_subdiv = tao - ttot;
    //cout << "Pops Merge\n";
    //cout << "Pop 1 Size : " << pop1_size << " Pop 2 Size : " << pop2_size << endl;
	               num_rem = pop1_size + pop2_size; 
                   if (mergesize==-9) {
                       double old_beta = new_beta;
                       double old_beta2 = new_beta2;
                       double temp_beta1 = getNewPopSize(growth,ttot,beta,bigN,old_beta,fte_exp,1);
                       double temp_beta2 = getNewPopSize(growth2,ttot,beta2,bigN2,old_beta2,fte_exp2,pop2_prop);
                       new_beta = temp_beta1 + temp_beta2;
                   } else {
                       new_beta = mergesize;
                   }
	               growth = 0;
	               growth2 = 0;
	               subdivision= 0;
	               break;
	           }
               if (mig_rate_count<mtimes_size) {
                if(ttot + tevent > mtimes[mig_rate_count] && mtimes[mig_rate_count] > 0) {
                    t_add_mig_change1 = mtimes[mig_rate_count] - ttot;
                    ++mig_rate_count;
                    M1 = 4*bigN*mig[mig_rate_count];
                    if (was_empty2!=0) {
                       M2=M1;
                    }
                    continue;
                }
               }
               if (mig2_rate_count <mtimes2_size) {
                if(ttot + tevent > mtimes2[mig2_rate_count] && mtimes2[mig2_rate_count] > 0) {
                    t_add_mig_change2 = mtimes2[mig2_rate_count] - ttot;
                    ++mig2_rate_count;
                    M2 = 4*bigN*mig2[mig2_rate_count];
                    continue;
                }
               }
               if (t_add_mig_change1 > 0) {
                   tevent = tevent + t_add_mig_change1;
                   t_add_mig_change1 = 0;
               }
               if (t_add_mig_change2 > 0) {
                   tevent = tevent + t_add_mig_change2;
                   t_add_mig_change2 = 0;
               }
	           //if instant growth has occured add appropriate extra time to the time to next event
	          if (tadd_instgrow > 0) {
	              tevent = tevent + tadd_instgrow;
	              tadd_instgrow = 0;
	          }
              //same as above, except for population 2
              if (tadd_instgrow2 > 0) {
                  tevent = tevent + tadd_instgrow2;
                  tadd_instgrow2 = 0;
              }     
	          // if population growth in population 1 has ended add appropriate extra time to the time to next event
	          if (tadd_growth > 0) {
	              tevent = tevent + tadd_growth;
	              tadd_growth = 0;
	          }
	          // if population growth in population 2 has ended add appropriate extra time to the the time to next event
	          if (tadd_growth2 > 0) {
	              tevent = tevent + tadd_growth2;
	              tadd_growth2 = 0;
	          }
	          //if time to next event exceeds the time for population growth for the first time => end growth, calculate new population size, additional time to next event
	          if (ttot + tevent > time_grow && time_grow > 0 && tadd_growth < 0) {
	              tadd_growth = time_grow - ttot;
	              double old_beta = new_beta;
	              new_beta = getNewPopSize(growth,time_grow,beta,bigN,old_beta,fte_exp,1);
                  growth = 0;  
	              continue;
	          }
	          //same as above except for population 2
	          if (ttot + tevent > time_grow2 && time_grow2 > 0 && tadd_growth2 < 0) {
	              tadd_growth2 = time_grow2 - ttot;
	              double old_beta2 = new_beta2;
	              new_beta2 = getNewPopSize(growth2,time_grow2,beta2,bigN2,old_beta2,fte_exp2,pop2_prop);
	              growth2 = 0;
	              continue;
	          }
	          //if the time to next event exceeds time until population size change for the first time => end growth, calculate new pop size, additional time to next event
	          if (inst_add_count < t_instant.size()) {
                if (ttot + tevent > t_instant[inst_add_count] && t_instant[inst_add_count] > 0) {
	                tadd_instgrow = t_instant[inst_add_count] - ttot;
	                double old_beta = new_beta;
	                new_beta = getNewPopSize(9,t_instant[inst_add_count],instgrow[inst_add_count],bigN,old_beta,fte_exp,1);
	                growth = 0;
                    inst_add_count++;
	                continue;
	            }
              }
             
              if (inst_add_count2 < t_instant2.size()) {
                if (ttot + tevent > t_instant2[inst_add_count2] && t_instant2[inst_add_count2] > 0) {
	                tadd_instgrow2 = t_instant2[inst_add_count2] - ttot;
	                double old_beta2 = new_beta2;
	                new_beta2 = getNewPopSize(9,t_instant2[inst_add_count2],instgrow2[inst_add_count2],bigN2,old_beta2,fte_exp2,pop2_prop);
	                growth2 = 0;
                    inst_add_count2++;
	                continue;
	            }
              }              
              
	          //if event = 0 then we have a coalescent event in pop1
	           if (event==0) {
	              double which_first = getunif();
	              double which_second = getunif();
	              coalEvent(Tree,which_first,which_second,tevent,&count,ttot,total_seglength,subdivision,1);
  // cout << "Coal Pop 1\n";
                  if (rep_counter >= num_ext && recomb>0) {
                        int mrca = allMRCA(Tree,total_seglength,subdivision,pop1_size,pop2_size,num_rem);
	                    if (mrca==1) {
                            num_rem=1;
                            subdivision=0;
                            break;
	                    }
                    }                  
	              pop1_size = pop1_size -1;
	          }
              
	          //if event = 1 coalescent event in pop2
	          if (event==1) {
	            double which_first = getunif();
	            double which_second = getunif();
	            coalEvent(Tree,which_first,which_second,tevent,&count,ttot,total_seglength,subdivision,2);
 //  cout << "Coal Pop 2\n";
                if (rep_counter >= num_ext && recomb>0) {
                        int mrca = allMRCA(Tree,total_seglength,subdivision,pop1_size,pop2_size,num_rem);
	                    if (mrca==1) {
                            subdivision=0;
                            num_rem=0;
                            break;
	                    }
                }
                pop2_size = pop2_size -1;
	          }
	          //if event = 2 migration pop1 => pop2
	          if (event==2) {
	              double which = getunif();
	              migrationEvent(Tree,which,tevent,ttot,1,2);
	              --pop1_size;
 //   cout << "Mig Pop1 -> Pop 2\n";
	              ++pop2_size;
	          }
	          //if event = 3 migration pop2 => pop1
	          if (event==3) {
	              double which = getunif();
	              migrationEvent(Tree,which,tevent,ttot,2,1);
	              ++pop1_size;
//   cout << "Mig Pop2 -> Pop1 \n";               
	              --pop2_size;
	          }
	          //if event = 4 recombination event
	          if (event==4) {
	              double which = getunif();
                  
 //  cout << "Recombination\n";
	              int capture = recombinationEvent(Tree,which,total_seglength,tevent,rng,&count,ttot,Recomb_Pos);
	              if (capture==1) {pop1_size = pop1_size+1;}
	              if (capture==2) {pop2_size = pop2_size+1;}
	          }
              ++rep_counter;
	          ttot = ttot + tevent;
              int still_remaining = pop1_size + pop2_size;
              if (still_remaining < 2) {
                  num_rem=1;
                  subdivision = 0;
                  break;
              }            
	      }
	  }
	                     
	    //if there was not ever subdivision this initializes the external branches on the coalescent tree                
	    if(!ever_subdivision) {
	        for (int i=0;i<num_rem;++i) {
	            Branch_t extBranch;
	            extBranch.num = i;
	            extBranch.current=1;
	            extBranch.length=0;
	            extBranch.external=1;
	            extBranch.population = 1;
	            extBranch.left= -1;
	            extBranch.right= -1;
	            extBranch.upper_left=-1;
	            extBranch.upper_right=-1;
	            extBranch.total_time_to_event=0;
	            extBranch.seg_starts.push_back(0);
	            extBranch.seg_ends.push_back(1);
	            Tree.push_back(extBranch);
	            ++count;
	        }
	     }
	     
	     if(!subdivision) {

	    //get time to next event, coalescent or recombination and select the earlier of the two
	     while (num_rem>1) { 
 
	          double tevent = 0;
	          int event = -1;
	          double t_coal_event = timeToCoal(num_rem,growth,beta,bigN,rng,ttot,new_beta,fte_exp,1);
	          if (t_coal_event==-99) {
	              cout << "Population size too small under current settings, rerun with new growth time or growth rate\n"; 
	              abort();
	          }
	          double t_rec_event = timeToRecomb(total_seglength,recomb,rng,num_rem);
              
	          if (t_coal_event <= t_rec_event) {
	              event = 1;
	              tevent = t_coal_event;
	          }
	          else {
	              event = 2;
	              tevent = t_rec_event;
	          }
              //checking to see if population has shrunk too far
              double current_size = getNewPopSize(growth,ttot,beta,bigN,new_beta,fte_exp,1);
              if (current_size*bigN <= 1) {
                  cout << "Population size too small under current settings, rerun with new growth time or growth rate\n";
                  exit(265);
              }
              
	          //These next conditions are the same as above
	          if (tadd_instgrow > 0) {
	              tevent = tevent + tadd_instgrow;
	              tadd_instgrow = 0;
	          }     
	          if (tadd_growth > 0) {
	              tevent = tevent + tadd_growth;
	              tadd_growth = 0;
	          }
	          if (tadd_subdiv > 0) {
	              tevent = tevent + tadd_subdiv;
	              tadd_subdiv = 0;
	          }
	          if (ttot + tevent > time_grow && time_grow > 0 && tadd_growth < 0) {
	              tadd_growth = time_grow - ttot;
	              double old_beta = new_beta;
	              new_beta = getNewPopSize(growth,time_grow,beta,bigN,old_beta,fte_exp,1);
                  lines_remaining = num_rem;
	              growth = 0;
//                  cout << endl << "Number of Lines Remaining at end of growth with" << endl;
//                  for (int mac=0;mac<num_ext;++mac) {
//                    int var_w_mac = linesWithDescendants(Tree,mac);
//                    if (var_w_mac!=0) {
//                        cout << mac << " Descendants: " << var_w_mac << endl;
//                    }
//                  }
//                  cout << endl;
                  continue;
	          } 
              if (inst_add_count < t_instant.size()) {
	            if (ttot + tevent > t_instant[inst_add_count] && t_instant[inst_add_count] > 0) {
	                tadd_instgrow = t_instant[inst_add_count] - ttot;
	                double old_beta = new_beta;
	                new_beta = getNewPopSize(9,t_instant[inst_add_count],instgrow[inst_add_count],bigN,old_beta,fte_exp,1);
                    lines_remaining = num_rem;
                    growth = 0;
                    inst_add_count++;
	                continue;
	            }
              }
	          //event 1 = coalescence
	          if (event ==1) {
	              double which_first = getunif();
	              double which_second = getunif();

	              coalEvent(Tree,which_first,which_second,tevent,&count,ttot,total_seglength,subdivision,1);
//                cout << "Tree Size: " << Tree.size() << " Total Seglength: " << total_seglength << endl;
                  if (rep_counter >= num_ext && recomb>0) {
                        int mrca = allMRCA(Tree,total_seglength,subdivision,pop1_size,pop2_size,num_rem);
	                    if (mrca==1) {
                            num_rem=0;
                            break;
	                    }
                    } 
	              num_rem = num_rem -1;
	          }
	          //event = 2 recombination
	          if (event == 2) {
	              double which = getunif();
	             
	              int capture = recombinationEvent(Tree,which,total_seglength,tevent,rng,&count,ttot,Recomb_Pos);
	              num_rem = num_rem + 1;
	          } 
              ++rep_counter;    
	          ttot = ttot + tevent;   
	      }   
	    

        //This calls the function which adds mutations to the final tree
//	     cout << "Adding Mutations\n";         
//        populateMutations(Tree,rng,Variant_Pos,theta);
       
        
         //This calls the function which removes any mutations which, through recombination, have populated onto every haplotype in the sample
//	     cout << "Removing Mutations\n";
//	    removeFullMutations(Tree,Variant_Pos);
	                     
        
//	    	     Next is all final output stuff
//	     cout << "Formating and pring final output\n";
	     int final_size = Tree.size();
//         std::sort(Variant_Pos.begin(),Variant_Pos.end());
//	     int num_sites = Variant_Pos.size();
//	     cout << "\n//\nsegsites: " << num_sites << endl << "positions: ";

//        cout << "//" << endl;
                         
//	     for (int i=0;i<num_sites;++i) {    
//	         cout << Variant_Pos[i] << " ";
//        }
//	     cout << endl;

//        multimap<int, double> desc_map; 
//
//        for (int i=0;i<final_size;++i) {
//            if (Tree[i].external==1) {
//                desc_map.insert(pair<int, double>(1,Tree[i].length));
//                
//            } else {
//                int descendants = numberDescendants(i,Tree);
//                if (descendants < 6) {
//                    
//                    double l1 = Tree[i].length;
//                    double l2 = 0;
//                    for (int no=0;no<Tree[i].seg_starts.size();++no) {
//                        l2 = l2 + (Tree[i].seg_ends[no] - Tree[i].seg_starts[no]);
//                    }
//                    
//                    double fl = l1 * l2;
//                                        
//                    desc_map.insert(pair<int, double>(descendants,fl));
//                }
//            }
//        }
//        
//        vector<int> seen_before;
//            
//        for (multimap<int, double>::iterator it1 = desc_map.begin();it1!=desc_map.end();++it1) {            
//            int money = (*it1).first;
//            int should_run = 0;
//            for (int mo=0;mo<seen_before.size();++mo) {
//                int prev = seen_before[mo];
//                if (money==prev) {
//                    ++should_run;
//                }
//            }
//            if (should_run==0) {
//                seen_before.push_back(money);
//            
//                pair<multimap<int, double>::iterator, multimap<int, double>::iterator> d_num_pair;
//                d_num_pair = desc_map.equal_range(money);
//            
////                int track = 0;
////                double sum1 = 0;
////                double square1 = 0;        
//            
//                
//                cout << money << ":";
//                for (multimap<int, double>::iterator it2 = d_num_pair.first; it2 != d_num_pair.second;++it2) {
//                    double current1 = (*it2).second;
//                    cout << " " << current1;
////                    double currentsquare = current1*current1;
////                    sum1 = sum1 + current1;
////                    square1 = square1 + currentsquare;
////                    track++;  
//                }
//                cout << endl;
//            
//                //double avg = sum1/track;
//                //double var = (square1 - (track*avg*avg))/(track - 1);
//            
// //               cout << money << " " << track << " " << sum1 << " " << square1 << endl;
//            }
//            
////            cout << endl; 
//        }
          
          if (re==0) {
              cout << "Branch,Descendants,Prop_of_branch,Branch_length,Prop_x_length\n";
          } 

          newGetLengths(Tree,Recomb_Pos,num_ext,bigN,a_sizes);  
          



	   
	     
//	     cout << Tree[Tree[count-1].left].total_time_to_event << endl;
	     
//             for (int i=0;i<final_size;++i) {
////            cout << "\nNext Branch\n";
//           cout << "Branch Number " << Tree[i].num << endl;
//           cout << "Branch Length " << Tree[i].length << endl;
//             cout << "Total Time to Event " << Tree[i].total_time_to_event << endl;
//           cout << "Branch Left: " << Tree[i].left << " Branch Right: " << Tree[i].right << endl;
//         cout << "Branch Left: " << Tree[i].left << " Pop " << Tree[Tree[i].left].population << " Branch Right: " << Tree[i].right << " Pop " << Tree[Tree[i].right].population << endl;
////
////           cout << "Branch Upper Left " << Tree[i].upper_left << " Branch Upper Right " << Tree[i].upper_right << endl;    
////             cout << "Branch Variants ";
////             int s = Tree[i].variant_locs.size();
////             for (int j=0;j<s;++j) {
////                 cout << Tree[i].variant_locs[j] << " ";
////             }
//           cout << "\nBranch Segments ";
//           int s2 = Tree[i].seg_starts.size();
//         for (int j=0;j<s2;++j) {
//              cout << "[" << Tree[i].seg_starts[j] << "," << Tree[i].seg_ends[j] << "] "; 
//           }
//           cout << endl;
//       }
//         if (num_sites>0) {
//            for(int i=0;i<final_size;++i) {       
//                if (Tree[i].external == 1) {
////               cout << "Haplotype" << Tree[i].num << " ";
//                    int num_varis = Tree[i].variant_locs.size();
//                    for (int j=0;j<num_sites;++j) {
//                        int seen = 0;
//                        for (int h=0;h<num_varis;++h) {
//                            if (Tree[i].variant_locs[h]==Variant_Pos[j]) {
//                              ++seen;
//                            }
//                        }
//                        if (seen > 0) {cout << 1;}
//                        else {cout << seen;}
//                    }
//                    cout << endl;
//                }
//            }
//        }

   //    cout << "Lines Remaining at End of Growth: " << lines_remaining << endl;
         if (report_lengths) {
            vector<double> length_vector;
            getLengths(Tree,Recomb_Pos,length_vector,num_ext);
            cout << "\nTotal Tree Length: " << length_vector[1] <<  " TMRCA: " << length_vector[0] << endl;
        }
       
       }
       Variant_Pos.clear();
	   Recomb_Pos.clear();
	   Tree.clear();
	   }
        
        if (timer) {
	        clock_t end = clock();
	        double duration = (double)(end - start)/CLOCKS_PER_SEC;
	        cout << "Elapsed time was " << duration << " sec\n";
        }
    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    return 0;
}
    



