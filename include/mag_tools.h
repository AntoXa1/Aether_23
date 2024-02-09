
#define sign(x) ((x) > 0) ? 1 : ((x) < 0) ? -1 : 0

#define GET_NAME(v) #v
// returns the name of the variable
#define SHOW(x) std::cout << #x"=" << x << std::endl;
// prints "varible name" = variable

#define PLOT(v) PrintToFileAndRunCompanionPlotScript(v, #v);

#define SPLOT(x,y,n) PrintToFileAndRunCompanionScatterPlotScript(x,y,n);

#include <vector>

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{
  std::vector<double> linspaced;
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) {
      linspaced.push_back(start);
      return linspaced;
  }
  double delta = (end - start) / (num - 1);
  for(int i=0; i < num-1; ++i){
      linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end); 
  // make start and end exactly the same as the input
  return linspaced;
}


template <typename T>
void PrintToFileAndRunCompanionPlotScript(T x, const std::string& name) { 
  // std::string cwd = "/Users/dorodnitsyn/WORK/PROJECTS/Aether_pub/Aether/run.test_mgrid";

  char cwd_[100];
  std::string  cwd(getcwd(cwd_, sizeof(cwd_)));       
  auto fileName = cwd + "/" + "_" + name+"_dbg.dat";   

  SHOW(fileName)

  int n_rows=0, n_cols=0, n_depth=0;

  if constexpr (std::is_same<T,arma::Cube<float>>::value) {      
    SHOW("arma_cube")
    n_rows=x.n_rows;
    n_cols=x.n_cols;
    n_depth=x.n_slices;
  }else if constexpr (std::is_same<T,arma::Mat<float>>::value) {
    // subview of the cube should know about cols and rows
    std::cout << typeid(x).name() << std::endl;       
    SHOW("arma_mat")
    n_rows=x.n_rows;
    n_cols=x.n_cols;
  }else{
    SHOW("not a cube or mat:")
    std::cout << typeid(x).name() << std::endl;
    n_rows=x.n_rows;
    n_cols=x.n_cols;  
  }
    SHOW(n_rows) 
    SHOW(n_cols)
    SHOW(n_depth)     

    // try {
    std::ofstream output_file(fileName);
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j+=1) {          
          output_file << x(i,j) << "\t";          
       }
    output_file << std::endl;
    }
    output_file.close();

  const std::string pyFileToExecute = cwd + "/" + "./_myplotfrom.py";
  SHOW(pyFileToExecute)
  
  system(pyFileToExecute.c_str());

}

// char cwd[100];
// std::string fdir(getcwd(cwd, sizeof(cwd)));        

//   std::string  fdir(getcwd(cwd, sizeof(cwd)));
//   auto fileName = fdir + "/" + "_data.dat";   
//   SHOW(fileName)
//   exit(10);
//     SHOW(len)

//     std::ofstream output_file(fileName);

//     for (int i = 0; i < len; i++) {        
//           output_file << x[i] <<"\t"<< y[i]<< std::endl;           
//     }
//     output_file.close();

//   const std::string pyFileToExecute = fdir + "/" + "./_myscat.py";
//   SHOW(pyFileToExecute)

template <typename T>
void PrintToFileAndRunCompanionScatterPlotScript(T x, T y, int len) { 
  //print two columns of the length, len
  char cwd1[100];
  std::string fdir(getcwd(cwd1, sizeof(cwd1)));
  auto fileName = fdir + "/" + "_data.dat";   
  SHOW(fileName) 
  SHOW(len)
    std::ofstream output_file(fileName);
    for (int i = 0; i < len; i++) {        
          output_file << x[i] <<"\t"<< y[i]<< std::endl;           
    }
    output_file.close();
    // exit(10);

  const std::string pyFileToExecute = fdir + "/" + "./_myscat.py";
  // const std::string pyFileToExecute = fdir + "/" + "./plot_from_cpp.jl";
  SHOW(pyFileToExecute)
  
  system(pyFileToExecute.c_str());
}

// class DipoleLine{
// public:
//   int numElem;
  
//   DipoleLine(int numElemIn,int tPow); 

// // friend class Grid;

//   std::vector<double> xx;
//   std::vector<double> zz;
//   std::vector<double> qq;
//   std::vector<double> rr;
//   std::vector<double> tt;
//   double tPower;

// };

