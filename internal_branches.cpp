#include<iostream>
#include<cmath>
#include <map>
#include <utility>
#include<cstdlib>
#include<tclap/CmdLine.h>
#include<eigen/Eigen/Core>
#include<eigen/Eigen/Dense>
#include<boost/random/uniform_int.hpp>
#include<boost/random/uniform_real.hpp>
#include<boost/random/poisson_distribution.hpp>
#include<boost/random/exponential_distribution.hpp>
#include<boost/random/variate_generator.hpp>
#include<boost/random/mersenne_twister.hpp>
#include<ctime>
#include<string>
#include<sstream>
#include<fstream>

using namespace std;
using Eigen::MatrixXd;

typedef boost::mt19937 prgType;

void split(const string& str, const string& delimiters, vector<string>& tokens) {

  if (!str.empty()) {

    size_t d_w = delimiters.length();
    size_t s_w = str.length();

    //string::size_type lastPos = str.find(delimiters);
    string::size_type lastPos = 0;  
    string::size_type pos = str.find(delimiters);
    if (pos==string::npos) {
      tokens.push_back(str);
    } else {
      while (string::npos != pos) {
	tokens.push_back(str.substr(lastPos,pos-lastPos));
	lastPos=pos+d_w;
	pos = str.find(delimiters,lastPos);
      }
      if (lastPos < s_w) {
	tokens.push_back(str.substr(lastPos,s_w));
      }
    }
  }

}

void fillQs(MatrixXd& cQs,int& nn) {
  
  for (int i=0;i<cQs.rows();++i) {
    double current_val = 0;
    cQs(i,i) = 1;
    for (int j=i+1;j<cQs.cols();++j) {
      current_val+= log(nn - j - 1) - log(nn - j + 1);
      cQs(i,j) = exp(current_val);
    }
  }

}
    
MatrixXd getStartTimes(int& size,unsigned int& nn,bool& emit,std::string& sfile,std::string& output) {

  vector<MatrixXd> start_probs;
  vector<MatrixXd> Qs;
  
  MatrixXd preiM(1,1);
  preiM.setZero();
  start_probs.push_back(preiM);

  MatrixXd initialMatrix(1,1);
  initialMatrix.setZero();
  initialMatrix(0,0) = 1;
  start_probs.push_back(initialMatrix);

  MatrixXd Qs1(1,1);
  Qs1(0,0) = 1; 
  MatrixXd Qs2(1,1);
  Qs2(0,0) = 1;
  Qs.push_back(Qs1);Qs.push_back(Qs2);

  if (sfile.compare("[none]")==0) {
  
    MatrixXd secondMatrix(2,2);
    secondMatrix.setZero();
    secondMatrix(0,0) = 1;
    secondMatrix(1,1) = 1;
    start_probs.push_back(secondMatrix);
  
    int current_n = 3;
    
    while (current_n <= nn) {
      
      MatrixXd cQs(current_n-1,current_n-1);
      cQs.setZero();
      fillQs(cQs,current_n);
      Qs.push_back(cQs);
         
      MatrixXd c_probs(size,current_n);
      if (current_n <= size) {
	c_probs.resize(current_n,Eigen::NoChange);
      }
      c_probs.setZero();
      c_probs(0,0) = 1;
      if (current_n <= size) {
	c_probs(current_n-1,current_n-1) = 1;
      }
      
      double c_n = (double) current_n;
      
      for (int csize=2;csize<=c_probs.rows();++csize) {
	for (int cbegin=(csize-1);cbegin<c_probs.cols()-1;++cbegin) {
	  if (csize==2 && cbegin==1) {
	    c_probs((csize-1),cbegin) = 1;
	  } else {
	    double cb_n = (double) cbegin;
	    double coal_prob = 2.0/((c_n - cb_n)*(c_n - (cb_n-1)));
	    double val = 0;
	    for (int mm=1;2*mm<=csize;++mm) {
	      
	      
	      if (csize==2) {
		int yy=0;int zz=0;
		MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
		MatrixXd inQ2 = Qs[(current_n-2)].row(zz).col(cbegin-1);
		val+= ((c_n*(c_n-1.0))/2.0)*c_probs((mm-1),yy)*c_probs((csize-mm)-1,zz)*inQ1(0,0)*inQ2(0,0);
	      } else {
		if (mm==1) {
		  int yy=0;
		  for (int zz=((csize-mm)-1);zz<cbegin;++zz) {
		    MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
		    MatrixXd inQ2 = Qs[(current_n-2)].row(zz).col(cbegin-1);
		    MatrixXd inStart1 = start_probs[(current_n-1)].row((csize-mm)-1).col(zz);
		    val+= c_n*c_probs(mm-1,yy)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
		  }
		} else {
		  if (mm==(csize-mm)) {
		    for (int yy=mm-1;yy<cbegin;++yy) {
		      for (int zz=(yy+1);zz<cbegin;++zz) {
			MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
			MatrixXd inQ2 = Qs[(current_n-2)].row(zz).col(cbegin-1);
			MatrixXd inStart1 = start_probs[(current_n - mm)].row((csize-mm)-1).col((zz - mm) + 1);
			val += c_probs((mm-1),yy)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
		      }
		    }
		  } else {
		    for (int yy=mm-1;yy<cbegin;++yy) {
		      for (int zz=((csize-mm)-1);zz<cbegin;++zz) {
			if (yy!=zz) {
			  if (yy > zz) {
			    MatrixXd inQ1 = Qs[(current_n-1)].row(zz).col(cbegin-1);
			    MatrixXd inQ2 = Qs[(current_n-2)].row(yy).col(cbegin-1);
			    MatrixXd inStart1 = start_probs[(current_n - (csize-mm))].row(mm-1).col(yy - (csize-mm) + 1);
			    val += c_probs((csize-mm)-1,zz)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
			  } else {
			    MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
			    MatrixXd inQ2 = Qs[(current_n - 2)].row(zz).col(cbegin-1);
			    MatrixXd inStart1 = start_probs[(current_n - mm)].row((csize-mm)-1).col(zz - mm + 1);
			    val+= c_probs((mm-1),yy)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	    c_probs(csize-1,cbegin) = coal_prob*val;
	  }
	}
      }

      // for (int i=0;i<c_probs.rows();++i) { 
      //   for (int j=0;j<c_probs.cols();++j) {
      // 	std::cout << c_probs(i,j) << " ";
      //   }
      //   std::cout << std::endl;
      // }
      // std::cout << std::endl;
    
      start_probs.push_back(c_probs);
  
      //for (int i=0;i<cQs.rows();++i) { 
      //   for (int j=0;j<cQs.cols();++j) {
      // 	std::cout << cQs(i,j) << " ";
      //   }
      //   std::cout << std::endl;
      // }
      // std::cout << std::endl;

      ++current_n;
    
    }
  } else {
       
    for (int cn=2;cn<=nn;++cn) {
      MatrixXd junk(1,1);
      start_probs.push_back(junk);
    }
        
    const char* fn = sfile.c_str();
    std::ifstream myfile (fn);
    std::string line;
    std::vector<int> done;
    bool exact_done = false;
          
    if (myfile.is_open()) {
      
      while(getline(myfile,line)) {

	std::vector<std::string> pre_in_line;
	std::string predelim = "=";
	std::string delim1 = " ";
	if (!line.empty()) {

	  split(line,predelim,pre_in_line);
	  if (pre_in_line[0].compare("n")==0) {
	  
	    int current_n = atoi(pre_in_line[1].c_str());
	    if (current_n==nn) {
	      exact_done = true;
	    }
	    done.push_back(current_n);

	    MatrixXd c_probs(size,current_n);
	    if (current_n <= size) {
	      c_probs.resize(current_n,Eigen::NoChange);
	    }
      
	    for (int csize=0;csize<c_probs.rows();++csize) {
	      std::string in_line;
	      std::vector<std::string> break_line;
	    
	      getline(myfile,in_line);
	      if (!in_line.empty()) {
		split(in_line,delim1,break_line);
		if ((break_line.size()==c_probs.cols()) || (break_line.size()==c_probs.cols()+1)) {
		  int end_val = break_line.size();
		  if (break_line.size()==c_probs.cols() + 1) {
		    if (break_line[(break_line.size()-1)].empty()) {
		      end_val--;
		    } else {
		      std::cerr << "Incorrect number of columns for entry " << line << std::endl;
		    }
		  }
		  for (int cbegin = 0;cbegin<end_val;++cbegin) {
		    c_probs(csize,cbegin) = atof(break_line[cbegin].c_str());
		  }   
		} else {
		  std::cerr << "Incorrect number of rows for entry " << line << std::endl;
		  exit(1);
		}
	      } else {
		std::cerr << "Incorrect number of rows for entry " << line << std::endl;
		exit(1);
	      }
	    }
	  
	    start_probs[current_n] = c_probs;
	  }
	}
      }
    } else {
      std::cerr << "Unable to open start file " << sfile << std::endl;
      exit(1);
    }

    if ((exact_done==false) || (emit==true)) {

      std::sort (done.begin(), done.end());
      
      int current_n = 2;
      int done_index = 0;

      while (current_n <= nn) {
	
	if (current_n > 2) {
	  MatrixXd cQs(current_n-1,current_n-1);
	  cQs.setZero();
	  fillQs(cQs,current_n);
	  Qs.push_back(cQs);
	}
	
	if (done[done_index]==current_n) {
	  ++done_index;
	  ++current_n;
	} else {
    
	  MatrixXd c_probs(size,current_n);
	  if (current_n <= size) {
	    c_probs.resize(current_n,Eigen::NoChange);
	  }
	  c_probs.setZero();
	  c_probs(0,0) = 1;
	  if (current_n <= size) {
	    c_probs(current_n-1,current_n-1) = 1;
	  }

	  double c_n = (double) current_n;
      
	  for (int csize=2;csize<=c_probs.rows();++csize) {
	    for (int cbegin=(csize-1);cbegin<c_probs.cols()-1;++cbegin) {
	      if (csize==2 && cbegin==1) {
		c_probs((csize-1),cbegin) = 1;
	      } else {
		double cb_n = (double) cbegin;
		double coal_prob = 2.0/((c_n - cb_n)*(c_n - (cb_n-1)));
		double val = 0;
		for (int mm=1;2*mm<=csize;++mm) {
	      
	      
		  if (csize==2) {
		    int yy=0;int zz=0;
		    MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
		    MatrixXd inQ2 = Qs[(current_n-2)].row(zz).col(cbegin-1);
		    val+= ((c_n*(c_n-1.0))/2.0)*c_probs((mm-1),yy)*c_probs((csize-mm)-1,zz)*inQ1(0,0)*inQ2(0,0);
		  } else {
		    if (mm==1) {
		      int yy=0;
		      for (int zz=((csize-mm)-1);zz<cbegin;++zz) {
			MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
			MatrixXd inQ2 = Qs[(current_n-2)].row(zz).col(cbegin-1);
			MatrixXd inStart1 = start_probs[(current_n-1)].row((csize-mm)-1).col(zz);
			val+= c_n*c_probs(mm-1,yy)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
		      }
		    } else {
		      if (mm==(csize-mm)) {
			for (int yy=mm-1;yy<cbegin;++yy) {
			  for (int zz=(yy+1);zz<cbegin;++zz) {
			    MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
			    MatrixXd inQ2 = Qs[(current_n-2)].row(zz).col(cbegin-1);
			    MatrixXd inStart1 = start_probs[(current_n - mm)].row((csize-mm)-1).col((zz - mm) + 1);
			    val += c_probs((mm-1),yy)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
			  }
			}
		      } else {
			for (int yy=mm-1;yy<cbegin;++yy) {
			  for (int zz=((csize-mm)-1);zz<cbegin;++zz) {
			    if (yy!=zz) {
			      if (yy > zz) {
				MatrixXd inQ1 = Qs[(current_n-1)].row(zz).col(cbegin-1);
				MatrixXd inQ2 = Qs[(current_n-2)].row(yy).col(cbegin-1);
				MatrixXd inStart1 = start_probs[(current_n - (csize-mm))].row(mm-1).col(yy - (csize-mm) + 1);
				val += c_probs((csize-mm)-1,zz)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
			      } else {
				MatrixXd inQ1 = Qs[(current_n-1)].row(yy).col(cbegin-1);
				MatrixXd inQ2 = Qs[(current_n - 2)].row(zz).col(cbegin-1);
				MatrixXd inStart1 = start_probs[(current_n - mm)].row((csize-mm)-1).col(zz - mm + 1);
				val+= c_probs((mm-1),yy)*inStart1(0,0)*inQ1(0,0)*inQ2(0,0);
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
		c_probs(csize-1,cbegin) = coal_prob*val;
	      }
	    }
	  }
      
	  // for (int i=0;i<c_probs.rows();++i) { 
	  //   for (int j=0;j<c_probs.cols();++j) {
	  // 	std::cout << c_probs(i,j) << " ";
	  //   }
	  //   std::cout << std::endl;
	  // }
	  // std::cout << std::endl;
    
	  start_probs[current_n] = c_probs;
  
	  //for (int i=0;i<cQs.rows();++i) { 
	  //   for (int j=0;j<cQs.cols();++j) {
	  // 	std::cout << cQs(i,j) << " ";
	  //   }
	  //   std::cout << std::endl;
	  // }
	  // std::cout << std::endl;

	  ++current_n;
    
	}
      }
    }
  }

  if (output.compare("start_probabilities")==0) {
    if (!emit) {
      std::cout << "n=" << nn << std::endl;
      MatrixXd finalm = start_probs[nn];
      for (int ii=0;ii<finalm.rows();++ii) {
	std::cout << finalm(ii,0);
	for (int jj=1;jj<finalm.cols();++jj) {
	  std::cout << " " << finalm(ii,jj);
	}
	std::cout << std::endl;
      }
    } else {
      for (int ln=2;ln<=nn;++ln) {
	std::cout << "n=" << ln << std::endl;
	MatrixXd finalm = start_probs[ln];
	for (int ii=0;ii<finalm.rows();++ii) {
	  std::cout << finalm(ii,0);
	  for (int jj=1;jj<finalm.cols();++jj) {
	    std::cout << " " << finalm(ii,jj);
	  }
	  std::cout << std::endl;
	}  
	std::cout << std::endl;
      }
    }
  }
  return start_probs[nn];
}

MatrixXd getEndTimes(unsigned int& nn,std::string& file,bool& emit,std::string& output) {
  	
  if (file.compare("[none]")==0) {

    if (!emit) {
    
      MatrixXd end_probs(nn-1,nn);
      end_probs.setZero();

      for (int start=0;start<(nn-1);++start) {
	double val = 1;
	for (int end = start+1;end<nn;++end) {
	  end_probs(start,end) = val* (2.0/(nn - (end-1)));
	  val = val * (nn - (end+1))/(nn - (end-1));
	}
      }

      if (output.compare("end_probabilities")==0) {
	std::cout << "n=" << nn << std::endl;
	for (int i=0;i<end_probs.rows();++i) {
	  std::cout << end_probs(i,0);
	  for (int j=1;j<end_probs.cols();++j) {
	    std::cout << " " << end_probs(i,j);
	  }
	  std::cout << endl;
	}
      }
      return end_probs;
    } else {
      
      for (int current_n=2;current_n<=nn;++current_n) {

	MatrixXd end_probs(current_n-1,current_n);
	end_probs.setZero();

	for (int start=0;start<(current_n-1);++start) {
	  double val = 1;
	  for (int end = start+1;end<current_n;++end) {
	    end_probs(start,end) = val* (2.0/(current_n - (end-1)));
	    val = val * (current_n - (end+1))/(current_n - (end-1));
	  }
	}

	if (output.compare("end_probabilities")==0) {
	  std::cout << "n=" << current_n << std::endl;
	  for (int i=0;i<end_probs.rows();++i) {
	    std::cout << end_probs(i,0);
	    for (int j=1;j<end_probs.cols();++j) {
	      std::cout << " " << end_probs(i,j);
	    }
	    std::cout << endl;
	  }
	  std::cout << endl;
	}
	if (current_n==nn) {
	  return end_probs;
	}
      }
    }
  } else {

    if (!emit) {

      const char* fn = file.c_str();
      std::ifstream myfile (fn);
      std::string line;
      std::vector<int> done;
      
      bool found = false;
      MatrixXd end_probs(nn-1,nn);
      
      if (myfile.is_open()) {
      
	while(getline(myfile,line)) {
	  
	  std::vector<std::string> pre_in_line;
	  std::string predelim = "=";
	  std::string delim1 = " ";
	  
	  if (!line.empty()) {
	    split(line,predelim,pre_in_line);
	    if (pre_in_line[0].compare("n")==0) {
	      if (atoi(pre_in_line[1].c_str())==nn) {
		found = true;

		for (int csize=0;csize<nn-1;++csize) {
		  std::string in_line;
		  std::vector<std::string> break_line;
	    
		  getline(myfile,in_line);
		  if (!in_line.empty()) {
		    split(in_line,delim1,break_line);
		    if ((break_line.size()==nn) || (break_line.size()==(nn+1))) {
		      int end_val = break_line.size();
		      if (break_line.size()==(nn + 1)) {
			if (break_line[(break_line.size()-1)].empty()) {
			  end_val--;
			} else {
			  std::cerr << "Incorrect number of columns for entry " << line << std::endl;
			}
		      }
		      if (output.compare("end_probabilities")==0) {
			std::cout << atof(break_line[0].c_str());
		      }
		      for (int cbegin = 1;cbegin<end_val;++cbegin) {
			if (output.compare("end_probabilities")==0) {
			  std::cout << " " << atof(break_line[cbegin].c_str());
			}
			end_probs(csize,cbegin) = atof(break_line[cbegin].c_str());
		      } 
		      if (output.compare("end_probabilities")==0) {
			std::cout << std::endl;
		      } 
		    }else {
		      std::cerr << "Incorrect number of rows for entry " << line << std::endl;
		      exit(1);
		    }
		  } else {
		    std::cerr << "Incorrect number of rows for entry " << line << std::endl;
		    exit(1);
		  }
		}
		break;
	      }
	    }
	  }
	}
	myfile.close();
      } else {
	std::cerr << "Unable to open start file " << file << std::endl;
	exit(1);
      }

      if (found==false) {
	std::cerr << "Value not found in provided end file, generating from scratch\n";
	end_probs.setZero();

	for (int start=0;start<(nn-1);++start) {
	  double val = 1;
	  for (int end = start+1;end<nn;++end) {
	    end_probs(start,end) = val* (2.0/(nn - (end-1)));
	    val = val * (nn - (end+1))/(nn - (end-1));
	  }
	}

	if (output.compare("end_probabilities")==0) {
	  std::cout << "n=" << nn << std::endl;
	  for (int i=0;i<end_probs.rows();++i) {
	    std::cout << end_probs(i,0);
	    for (int j=1;j<end_probs.cols();++j) {
	      std::cout << " " << end_probs(i,j);
	    }
	    std::cout << endl;
	  }
	}
      }
      return end_probs;
    } else {

      std::vector<MatrixXd> end_probs;
      std::vector<int> done;
    
      for (int cn=2;cn<=nn;++cn) {
	MatrixXd junk1(1,1);
	end_probs.push_back(junk1);
      }

      const char* fn = file.c_str();
      std::ifstream myfile (fn);
      std::string line;

      
      if (myfile.is_open()) {
      
	while(getline(myfile,line)) {
	  
	  std::vector<std::string> pre_in_line;
	  std::string predelim = "=";
	  std::string delim1 = " ";
	  
	  if (!line.empty()) {
	    split(line,predelim,pre_in_line);
	    if (pre_in_line[0].compare("n")==0) {

	      int current_n = atoi(pre_in_line[1].c_str());
	      MatrixXd c_probs(current_n-1,current_n);
	      
	      for (int csize=0;csize<current_n-1;++csize) {
		std::string in_line;
		std::vector<std::string> break_line;
	    
		getline(myfile,in_line);
		if (!in_line.empty()) {
		  split(in_line,delim1,break_line);
		  if ((break_line.size()==current_n) || (break_line.size()==(current_n+1))) {
		    int end_val = break_line.size();
		    if (break_line.size()==(current_n + 1)) {
		      if (break_line[(break_line.size()-1)].empty()) {
			end_val--;
		      } else {
			std::cerr << "Incorrect number of columns for entry " << line << std::endl;
		      }
		    }
		    
		    for (int cbegin = 0;cbegin<end_val;++cbegin) {
		      c_probs(csize,cbegin) = atof(break_line[cbegin].c_str());
		    }   
		    
		  } else {
		    std::cerr << "Incorrect number of rows for entry " << current_n << std::endl << line << std::endl;
		    exit(1);
		  }
		} else {
		  std::cerr << "Incorrect number of rows for entry " << current_n << std::endl << line << std::endl;
		  exit(1);
		}
	      }
	      end_probs[current_n] = c_probs;
	      done.push_back(current_n);
	    }
	  }
	}
	myfile.close();
      } else {
	std::cerr << "Unable to open start file " << file << std::endl;
	exit(1);
      }

      std::sort (done.begin(), done.end());
    
      int current_n = 2;
      int done_index = 0;
      
      while (current_n <= nn) {

	if (output.compare("end_probabilities")==0) {
	  std::cout << "n=" << current_n << std::endl; 
	}
	if (done[done_index]==current_n) {

	  if (output.compare("end_probabilities")==0) {
	
	    MatrixXd inprob = end_probs[current_n];
	    for (int i=0;i<inprob.rows();++i) {
	      std::cout << inprob(i,0);
	      for (int j=1;j<inprob.cols();++j) {
		std::cout << " " << inprob(i,j);
	      }
	      std::cout << std::endl;
	    }
	    std::cout << std::endl;
	  }
	  ++done_index;
	  ++current_n;
	
	} else {
	  
	  MatrixXd end_probs(current_n-1,current_n);
	  end_probs.setZero();

	  for (int start=0;start<(current_n-1);++start) {
	    double val = 1;
	    for (int end = start+1;end<current_n;++end) {
	      end_probs(start,end) = val* (2.0/(current_n - (end-1)));
	      val = val * (current_n - (end+1))/(current_n - (end-1));
	    }
	  }

	  if (output.compare("end_probabilities")==0) {
	    std::cout << "n=" << current_n << std::endl;
	    for (int i=0;i<end_probs.rows();++i) {
	      std::cout << end_probs(i,0);
	      for (int j=1;j<end_probs.cols();++j) {
		std::cout << " " << end_probs(i,j);
	      }
	      std::cout << endl;
	    }
	  }
	  ++current_n;
	}
      }
      return end_probs[nn];
    }
  }
}

//The probability of a branch with n descendants arising from a coalescence between branches with i, n-i descendants 
double splitProb(int n,int i) {
    int num_val = 0;
    if ((n % 2)==0) {
        if (2*i==n) {
            num_val = 1;
        }
    }
    double val = (double) (2-num_val)/(n-1);
    return val;
} 


//Iterates over every possible split of m descendants into n, m-n and calculates new matrix of probabilities for Tree entry m 
void fillMatrix(vector<double>& probs, int m, vector<vector<double> > Tree, int a) {
  for (int n=m-1;2*n>=m;--n) {
    int i = m-n;
    double p_split = splitProb(m,n);
    for (int col=0;col<m/a+1;++col) {
            
      double pval = 0;
      for (int ind1 = col;ind1>=0;--ind1) {
	int ind2 = col-ind1;
	if (ind1<(n/a)+1 && ind2<(i/a)+1) {
	  pval+= Tree[n][ind1] * Tree[i][ind2];                                          
	}
      }
      probs[col] = probs[col] + p_split*pval;
    }
  }
}
    

vector<double> getBnums(int& size,unsigned int& nn,std::string& file,bool& emit,std::string& output) {

   
  vector<vector<double> > Tree;

  for (int i=0;i<size;++i) {
    vector<double> firstBranch;
    firstBranch.push_back(1);
    for (int j=0;j<(i/size)+1;++j) {
      firstBranch.push_back(0);
    }
    Tree.push_back(firstBranch);
  }

  vector<double> firstBranch;
  firstBranch.push_back(0);
  firstBranch.push_back(1);
  Tree.push_back(firstBranch);

  if (file.compare("NA")==0) {
    for (int m=size+1;m<=nn;++m) {
      
      vector<double> c_prob;
      for (int i=0;i<(m/size)+1;++i) {
	c_prob.push_back(0);
      }            
      
      fillMatrix(c_prob,m,Tree,size);            
      Tree.push_back(c_prob);
    }
    
  } else {
    
    for (int m=size+1;m<=nn;++m) {
      
      vector<double> c_prob;
      for (int i=0;i<(m/size)+1;++i) {
	c_prob.push_back(0);
      }            
      
      Tree.push_back(c_prob);
    }
    
    std::vector<int> done;
    const char* fn = file.c_str();
    std::ifstream myfile (fn);
    std::string line;
      
    if (myfile.is_open()) {
      
	while(getline(myfile,line)) {
	  
	  std::vector<std::string> pre_pre_in_line;
	  std::vector<std::string> pre_in_line;
	  std::vector<std::string> pre_in_line2;
	  std::string predelim = "=";
	  std::string predelim2 = ",";
	  std::string delim1 = " ";
	  
	  if (!line.empty()) {
	    split(line,predelim2,pre_pre_in_line);
	    split(pre_pre_in_line[0],predelim,pre_in_line);
	    if (pre_in_line[0].compare("n")==0) {
	      int current_n = atoi(pre_in_line[1].c_str());
	      if (current_n > size) {
		split(pre_pre_in_line[1],predelim,pre_in_line2);
		int current_size = atoi(pre_in_line2[1].c_str());
		if (current_size==size) {
		  std::string in_line;
		  std::vector<std::string> in_in_line;
		  getline(myfile,in_line);
		  split(in_line,delim1,in_in_line);
		  if (in_in_line.size()==Tree[current_n].size()) {
		    for (int ii=0;ii<in_in_line.size();++ii) {
		      Tree[current_n][ii] = atof(in_in_line[ii].c_str());
		    }
		    done.push_back(current_n);
		  } else {
		    std::cerr << "Number of probabilties in file for size " << size << " Current_N " << current_n << " does not match expectations: " << in_in_line.size() << " in file " << Tree[current_n].size() << " expected\n";
		    exit(1);
		  }
		}
	      } else {
		std::string in_line;
		getline(myfile,in_line);
	      }
	    }
	  }
	}
	myfile.close();
    } else {
      std::cerr << "Unable to open branch length file " << file << std::endl;
      exit(1);
    }

    std::sort (done.begin(), done.end());
    
    int done_ind = 0;
    
    for (int m=size+1;m<=nn;++m) {
      
      if (m==done[done_ind]) {
	++done_ind;
	std::cout << "Using file of " << m << std::endl;
      } else {

	vector<double> c_prob;
	for (int i=0;i<(m/size)+1;++i) {
	  c_prob.push_back(0);
	}            
      
	fillMatrix(c_prob,m,Tree,size);            
	Tree[m] = c_prob;
	std::cout << "Calculating for " << m << std::endl;
      }
    }
  }

  if (output.compare("branch_numbers")==0) {    
    if (!emit) {
      vector<double> finalm = Tree[nn];
      for (int ii=0;ii<finalm.size();++ii) {
	std::cout << nn << " " << size << " " << ii << " " << finalm[ii] << std::endl;
      } 
    } else {
      for (int ln=2;ln<=nn;++ln) {
	for (int ii=0;ii<Tree[ln].size();++ii) {
	  std::cout << ln << " " << size << " " << ii << " " << Tree[ln][ii] << std::endl;
	}
      }
    }
  }
  return Tree[nn];
}

void getLengthDistribution(unsigned int& samplesize,int& csize,MatrixXd& eprobs,MatrixXd& sprobs,std::string& prebreaks,prgType& rng,bool expected,unsigned int& gnum) {


  boost::uniform_real<> cl(0,1);
  boost::variate_generator<prgType&, boost::uniform_real<> > getunif(rng,cl);
  
  std::vector<double> breaks;
  std::vector<std::string> bvals;
  std::string delim = ",";

  split(prebreaks,delim,bvals);
  if (bvals.size() > 1) {
    for (int i=0;i<bvals.size();++i) {
      breaks.push_back(atof(bvals[i].c_str()));
    }
  } else {
    bvals.clear();
    if (atof(prebreaks.c_str())!=0) {
      breaks.push_back(atof(prebreaks.c_str()));
    } else {
      const char* fn = prebreaks.c_str();
      std::ifstream myfile (fn);
      std::string line;
      if (myfile.is_open()) {
	getline(myfile,line);
	split(line,delim,bvals);
	for (int i=0;i<bvals.size();++i) {
	  breaks.push_back(atof(bvals[i].c_str()));
	}
      } else {
	std::cerr << "Unable to open interpret breaks statement as comma seperated value: " << prebreaks << std::endl;
	exit(1);
      }
      myfile.close();
    }
  }
  std::sort(breaks.begin(),breaks.end());
 
  if (breaks[0]!=0.0) {
    std::vector<double>::iterator it;
    it= breaks.begin();
    breaks.insert(it,0);
  }

  // std::vector<double> c_probs;  
  std::vector<double> c_counts;
  for (int i=0;i<breaks.size();++i) {
    c_counts.push_back(0);
  }

  std::vector<double> rstarts;
  double start_sum = 0;
  for (int i=0;i<sprobs.row(csize-1).size();++i) {
    start_sum+=sprobs(csize-1,i);
  }
  for (int i=0;i<sprobs.row(csize-1).size();++i) {
    rstarts.push_back(sprobs(csize-1,i)/start_sum);
  }

  std::vector<double> expected_lengths;
  if (expected==true) {
    double base = (double) samplesize;
    while (base>=2) {
      double ex_val = 2.0/(base*(base-1.0));
      expected_lengths.push_back(ex_val);
      --base;
    }
  }

  unsigned int gene = 0;
  while (gene < gnum) {

    std::vector<double> val_vec;

    std::vector<double> current_lengths;
    if (expected==true) {
      current_lengths = expected_lengths;
    } else {
      double base = (double) samplesize;
      while (base>=2) {
	double rate = (base*(base-1.0))/2.0;
	boost::exponential_distribution<> exp_dist(rate);
	boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
	current_lengths.push_back(getexp());
	--base;
      }
    }
    int evals = 0;
  
    while ((evals < samplesize) && (evals < 2000)) {
      double stp = getunif();
      double enp = getunif();
      int stint = -1;
      int stent = -1;
      double start = 0;
      double val = 0;
      
      for (int ii=0;ii<rstarts.size();++ii) {
	if ((stp > start) && (stp <= (start + rstarts[ii]))) {
	  stint = ii;
	  break;
	} else {
	  start+=rstarts[ii];
	}
      }
      start = 0;
      for (int ii=stint;ii<eprobs.row(stint).size();++ii) {
	if ((enp > start) && (enp <= (start + eprobs(stint,ii)))) {
	  stent = ii;
	  break;
	} else {
	  start+=eprobs(stint,ii);
	}
      }
      if ((stint==-1) || (stent==-1)) {
	std::cerr << "Error choosing branch start/end\n";
	exit(1);
      }
      
      for (int ii=stint;ii<stent;++ii) {
	val+=current_lengths[ii];	
      }

      val_vec.push_back(val);
      ++evals;
    }
    gene+=evals;

    std::sort(val_vec.begin(),val_vec.end());

    unsigned int cc=val_vec.size();
    unsigned int ind = (val_vec.size()-1);
    for (int ii=(breaks.size()-1);ii>=0;--ii) {
      if (val_vec[ind] < breaks[ii]) {
	c_counts[ii]+=cc;
      } else {
	while ((val_vec[ind] > breaks[ii]) && (ind > 0)) {
	  cc--;
	  ind--;
	  if (ind==0) {
	    break;
	  }
	}
	if (ind>0) {
	  c_counts[ii]+=cc;
	} else {
	  if (val_vec[ind] < breaks[ii]) {
	    c_counts[ii]+=cc;
	  }
	}
      }
    }
    val_vec.clear();
  
  }

  for (int i=0;i<breaks.size();++i) {
    std::cout << csize << " " << samplesize << " " << breaks[i] << " ";
    printf("%5.4f",c_counts[i]/((double) gene));
    std::cout << std::endl;
  }
}

void getR2probs(int& csize,unsigned int& samplesize,bool expected,MatrixXd& eprobs,MatrixXd& sprobs,vector<double>& bnums,prgType& rng,unsigned int& gnum) {

  boost::uniform_real<> cl(0,1);
  boost::variate_generator<prgType&, boost::uniform_real<> > getunif(rng,cl);

  double final_val = 0;
  double total_evals = 0;

  std::vector<double> rstarts;
  double start_sum = 0;
  for (int i=0;i<sprobs.row(csize-1).size();++i) {
    start_sum+=sprobs(csize-1,i);
  }
  for (int i=0;i<sprobs.row(csize-1).size();++i) {
    rstarts.push_back(sprobs(csize-1,i)/start_sum);
  }

  std::vector<double> expected_lengths;
  if (expected==true) {
    double base = (double) samplesize;
    while (base>=2) {
      double ex_val = 2.0/(base*(base-1.0));
      expected_lengths.push_back(ex_val);
      --base;
    }
  }


  std::vector<double> current_lengths;
  if (expected==true) {
    current_lengths = expected_lengths;
  }  else {
     double base = (double) samplesize;
     while (base>=2) {
       double rate = (base*(base-1.0))/2.0;
       boost::exponential_distribution<> exp_dist(rate);
       boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
       current_lengths.push_back(getexp());
       --base;
     }
   }

  double evals = 0;
  
  for (unsigned int gene=0;gene<gnum;++gene) {
    
    if (expected!=true) {
      if ((evals > samplesize) || (evals > 2000.0)) {
	current_lengths.clear();
	double base = (double) samplesize;
	while (base>=2) {
	  double rate = (base*(base-1.0))/2.0;
	  boost::exponential_distribution<> exp_dist(rate);
	  boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
	  current_lengths.push_back(getexp());
	  --base;
	}
	evals = 0;
      }
    }
    
    std::vector<double> branches;
    double btp = getunif();
    double branch_start = 0;
    int32_t br_num = -1;
    for (int i=0;i<bnums.size();++i) {
      if ((btp > branch_start) && (btp <= (branch_start + bnums[i]))) {
	br_num = i;
	break;
      } else {
	branch_start+=bnums[i];
      }
    }
    if (br_num==-1) {
      std::cerr << "Error selecting branch number" << std::endl;
      exit(1);
    }
    evals+=(double) br_num;
    if (br_num > 0) {
      for (int ii=0;ii<br_num;++ii) {
	double stp = getunif();
	double enp = getunif();
	int stint = -1;
	int stent = -1;
	double start = 0;
	double val = 0;
	
	for (int ii=0;ii<rstarts.size();++ii) {
	  if ((stp > start) && (stp <= (start + rstarts[ii]))) {
	    stint = ii;
	    break;
	  } else {
	    start+=rstarts[ii];
	  }
	}
	start = 0;
	for (int ii=stint;ii<eprobs.row(stint).size();++ii) {
	  if ((enp > start) && (enp <= (start + eprobs(stint,ii)))) {
	    stent = ii;
	    break;
	  } else {
	    start+=eprobs(stint,ii);
	  }
	}
	if ((stint==-1) || (stent==-1)) {
	  std::cerr << "Error choosing branch start/end (2)\n";
	  exit(1);
	}
      
	for (int ii=stint;ii<stent;++ii) {
	  val+=current_lengths[ii];	
	}

	branches.push_back(val);
      }    
      double numerator = 0;
      double denominator = 0;
      for (int ii=0;ii<branches.size();++ii) {
	numerator+=branches[ii]*branches[ii];
	denominator+=branches[ii];
      }
      denominator = denominator*denominator;
      final_val+=(numerator/denominator);     
    }
  }
  final_val = final_val/((double) gnum);

  std::cout << csize << " " << samplesize << " ";
  printf("%5.4f",final_val);
  std::cout << std::endl;
  
}  

void getGenePortions(int& csize,unsigned int& samplesize,bool expected,MatrixXd& eprobs,MatrixXd& sprobs,vector<double>& bnums,prgType& rng,unsigned int& gnum) {

  boost::uniform_real<> cl(0,1);
  boost::variate_generator<prgType&, boost::uniform_real<> > getunif(rng,cl);

  std::vector<double> rstarts;
  double start_sum = 0;
  for (int i=0;i<sprobs.row(csize-1).size();++i) {
    start_sum+=sprobs(csize-1,i);
  }
  for (int i=0;i<sprobs.row(csize-1).size();++i) {
    rstarts.push_back(sprobs(csize-1,i)/start_sum);
  }

  std::vector<double> expected_lengths;
  if (expected==true) {
    double base = (double) samplesize;
    while (base>=2) {
      double ex_val = 2.0/(base*(base-1.0));
      expected_lengths.push_back(ex_val);
      --base;
    }
  }

  double evals = 0;    
  std::vector<double> current_lengths;
  if (expected==true) {
    current_lengths = expected_lengths;
  } else {
    double base = (double) samplesize;
    while (base>=2) {
      double rate = (base*(base-1.0))/2.0;
      boost::exponential_distribution<> exp_dist(rate);
      boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
      current_lengths.push_back(getexp());
      --base;
    }
  }

  for (unsigned int gene=0;gene<gnum;++gene) {

    std::vector<double> branches;
    double btp = getunif();
    double branch_start = 0;
    int32_t br_num = -1;
    for (int i=0;i<bnums.size();++i) {
      if ((btp > branch_start) && (btp <= (branch_start + bnums[i]))) {
	br_num = i;
	break;
      } else {
	branch_start+=bnums[i];
      }
    }
    if (br_num==-1) {
      std::cerr << "Error selecting branch number" << std::endl;
      exit(1);
    }

    if (expected==false) {
      if ((br_num + evals > samplesize) || (br_num + evals > 1000)) {
	current_lengths.clear(); 
	double base = (double) samplesize;
	 while (base>=2) {
	   double rate = (base*(base-1.0))/2.0;
	   boost::exponential_distribution<> exp_dist(rate);
	   boost::variate_generator<prgType&, boost::exponential_distribution<> > getexp(rng,exp_dist);
	   current_lengths.push_back(getexp());
	   --base;
	 }
      }
    }

    for (int ii=0;ii<br_num;++ii) {
      double stp = getunif();
      double enp = getunif();
      int stint = -1;
      int stent = -1;
      double start = 0;
      double val = 0;
	
      for (int ii=0;ii<rstarts.size();++ii) {
	if ((stp > start) && (stp <= (start + rstarts[ii]))) {
	  stint = ii;
	  break;
	} else {
	  start+=rstarts[ii];
	}
      }
      start = 0;
      for (int ii=stint;ii<eprobs.row(stint).size();++ii) {
	if ((enp > start) && (enp <= (start + eprobs(stint,ii)))) {
	  stent = ii;
	  break;
	} else {
	  start+=eprobs(stint,ii);
	}
      }
      if ((stint==-1) || (stent==-1)) {
	std::cerr << "Error choosing branch start/end (2)\n";
	exit(1);
      }
      
      for (int ii=stint;ii<stent;++ii) {
	val+=current_lengths[ii];	
      }
      
      branches.push_back(val);
    }    
    
    evals+=br_num;
    
    double sum_length = 0;
    double squared = 0;
    double variance = 0;
    for (int i=0;i<branches.size();++i) {
      sum_length+=branches[i];
      squared+=branches[i]*branches[i];
    }
    double mean = 0;
    if (branches.size()>0) {
      double bsizer = (double) branches.size();
      mean = sum_length/bsizer;
      if (branches.size()>1) {
	variance = (squared - bsizer*mean*mean)/(bsizer - 1.0);
      }
    }
    
    std::cout << csize << " " << samplesize << " ";
    printf("%lu %4.3f %4.3f",branches.size(),sum_length,variance);
    for (int i=0;i<branches.size();++i) {
      std::cout << " " << branches[i];
    }
    std::cout << std::endl;
    
  }

}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Fun times", ' ', "0.00");
    TCLAP::ValueArg<std::string> startProb("","starts","file with starting probabilities, in format output by start_probabilities command",false,"[none]","string");
    cmd.add(startProb);
    TCLAP::ValueArg<std::string> endProb("","ends","Comma separated list of files with ending probabilities, same order as \"s\" command. If entering only partial list of files, replace missing entries with \"NA\".",false,"[none]","string");
    cmd.add(endProb);
    TCLAP::ValueArg<std::string> bNum("","bnum","Comma separated list of files with number of branches in genealogies, same order as \"s\" command. If entering only partial list of files, replace missing entries with \"NA\".",false,"[none]","string");
    cmd.add(bNum);
    TCLAP::ValueArg<std::string> sampleSize("n","sample_size","Comma separated list of total sample sizes to output values for",true,"[none]","string");
    cmd.add(sampleSize);
    TCLAP::ValueArg<std::string> siZes("s","sizes","Comma separated list of branch sizes to generate values for. With \"start_probabilities\" command the maximum value from this list is used",true,"[none]","string");
    cmd.add(siZes);
    TCLAP::ValueArg<std::string> outPut("o","output","Desired output. Options: \"start_probabilities\" \"end_probabilities\" \"branch_numbers\" \"length_distribution\" \"branch_lengths\" \"r2_prob\" \"genealogy_portions\"",true,"[none]","string");
    cmd.add(outPut);
    TCLAP::ValueArg<int> branchToGen("","bn","For each size given, this number of random branch lengths will be generated, output is two columns, first with size, second with length",false,1,"int");
    cmd.add(branchToGen);
    TCLAP::ValueArg<double> rSeed("","seed","Random number generator seed",false,-9,"double");
    cmd.add(rSeed);
    TCLAP::SwitchArg expectedValues("","expected","Use expected values instead of randomly generated coalescent times",false);
    cmd.add(expectedValues);
    TCLAP::ValueArg<std::string> distBreaks("","breaks","Either a comma separated list of values to output cdf probabilities for, or a filename that contains such a list",false,"[none]","string");
    cmd.add(distBreaks);
    TCLAP::SwitchArg emitAll("","emit_all","For \"start_probabilities\" and \"end_probabilities\" should the program emit all sample sizes, or only final sample size",false);
    cmd.add(emitAll);
    TCLAP::ValueArg<unsigned int> GeneNum("x","reps","For \"length_distribution\" this is the number of individual branches used to estimate the distribution of lengths. For \"r2_prob\" and \"genealogy_portions\" this is the number of genealogies realized [default=1]",false,1,"unsigned int");
    cmd.add(GeneNum);

    cmd.parse(argc,argv);

    std::string sfiles = startProb.getValue();
    std::string efiles = endProb.getValue();
    std::string prebfiles = bNum.getValue();
    std::string presizes = siZes.getValue();
    double randseed = rSeed.getValue();
    std::string output = outPut.getValue();
    std::string presamplesize = sampleSize.getValue();
    bool expected = expectedValues.getValue();
    std::string prebreaks = distBreaks.getValue();
    bool emit = emitAll.getValue();
    int br_n = branchToGen.getValue();
    unsigned int gnum = GeneNum.getValue();

    std::string delim1 = ",";
    std::vector<int> sizes;
    std::vector<std::string> pre_s_vec;
    split(presizes,delim1,pre_s_vec);
    for (int ii=0;ii<pre_s_vec.size();++ii) {
      sizes.push_back(atoi(pre_s_vec[ii].c_str()));
    }

    //std::vector<std::string> sfiles;
    //std::vector<std::string> efiles;
    std::vector<std::string> bfiles;

    unsigned int max_ssize = 1;
    std::vector<unsigned int> sample_sizes;
    std::vector<std::string> pre_ss;
    split(presamplesize,delim1,pre_ss);
    for (int ii=0;ii<pre_ss.size();++ii) {
      sample_sizes.push_back(atoi(pre_ss[ii].c_str()));
      if (atoi(pre_ss[ii].c_str()) > max_ssize) {
	max_ssize = atoi(pre_ss[ii].c_str());
      }
    }

    if (prebfiles.compare("[none]")!=0) {
      split(prebfiles,delim1,bfiles);
      if (bfiles.size() != sizes.size()) {
	std::cerr << "Incorrect number of files for number of branches given number of sizes entered. Enter NAs in string for missing files\n";
     	exit(1);
      }
    } else {
      for (int ii=0;ii<sizes.size();++ii) {
     	bfiles.push_back("NA");
      }
    }

    prgType rng;
    double seed = time(0);
    if (randseed != -9) {
      seed = randseed;
    }
    rng.seed(seed);
  
    bool did_something = false;

    int max_size = -1;
    for (int i=0;i<sizes.size();++i) {
      if (sizes[i] > max_size) {
	max_size = sizes[i];
      }
    }
        
    bool header = true;

    if (output.compare("end_probabilities")==0) {
      if (emit==true) {
	MatrixXd eprobs = getEndTimes(max_ssize,efiles,emit,output);    
      } else {
	for (int ss=0;ss<sample_sizes.size();++ss) {
	  unsigned int samplesize = sample_sizes[ss];
	  MatrixXd eprobs = getEndTimes(samplesize,efiles,emit,output);    
	}
      }
      did_something=true;
    } else {
      if (output.compare("start_probabilities")==0) {  
	if (emit==true) {
	  MatrixXd sprobs = getStartTimes(max_size,max_ssize,emit,sfiles,output);
	} else {
	  for (int ss=0;ss<sample_sizes.size();++ss) {
	    unsigned int samplesize = sample_sizes[ss];
	    MatrixXd sprobs = getStartTimes(max_size,samplesize,emit,sfiles,output);
	  }
	}
	did_something = true;
      } else {
	for (int ss=0;ss<sample_sizes.size();++ss) {
	  unsigned int samplesize = sample_sizes[ss];
	  if (output.compare("branch_numbers")==0) {
	    for (int re=0;re<sizes.size();++re) {
	      int csize = sizes[re];
	      vector<double> bnums = getBnums(csize,samplesize,bfiles[re],emit,output);
	      did_something=true;
	    }
	  } else {
	    for (int re=0;re<sizes.size();++re) {
	      int csize = sizes[re];
	      MatrixXd eprobs = getEndTimes(samplesize,efiles,emit,output);    
	      MatrixXd sprobs = getStartTimes(csize,samplesize,emit,sfiles,output);
	      if (output.compare("length_distribution")==0) {
		if (header==true) {
		  std::cout << "BranchSize SampleSize Break CDF\n";
		  header=false;
		}
		getLengthDistribution(samplesize,csize,eprobs,sprobs,prebreaks,rng,expected,gnum);
		did_something=true;
	      } else {
		if (output.compare("r2_prob")==0) {
		  if (header==true) {
		    std::cout << "BranchSize SampleSize ProbR2" << std::endl;
		    header = false;
		  }
		  vector<double> bnums = getBnums(csize,samplesize,bfiles[re],emit,output);
		  getR2probs(csize,samplesize,expected,eprobs,sprobs,bnums,rng,gnum);
		  did_something=true;
		} else {
		  if (output.compare("genealogy_portions")==0) {
		    if (header==true) {
		      std::cout << "BranchSize SampleSize NumBranches SumLength Variance IndvBranches\n";
		      header=false;
		    }
		    vector<double> bnums = getBnums(csize,samplesize,bfiles[re],emit,output);
		    getGenePortions(csize,samplesize,expected,eprobs,sprobs,bnums,rng,gnum);
		    did_something=true;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  
    if (did_something==false) {
      std::cerr << "Invalid output option " << output << "\nSee help for valid options\n";
      exit(1);
    }
  

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
  return 0;
}
