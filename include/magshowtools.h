
#define GET_NAME(v) #v
// returns the name of the variable

#define SHOW(x) std::cout << #x"=" << x << std::endl;
// prints "varible name" = variable

#define PLOT(v) PrintToFileAndExecCompanionPlotScript(v, #v);

template <typename T>

void PrintToFileAndExecCompanionPlotScript(T x, const std::string& name) { 
  std::string cwd = "/Users/dorodnitsyn/WORK/PROJECTS/Aether_pub/Aether/run.test_mgrid";
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
